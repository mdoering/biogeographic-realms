import math
from pathlib import Path

import geopandas as gpd
import numpy as np
import polars as pl
import pyproj
import pytest
import shapely
from shapely import Polygon
from test_grid import POLE_CUT, POLE_WHOLE, A, B, C, S, antimeridian_cells, north_cells, write_grid

import realms.boundaries
from realms.boundaries import (
    Rotation,
    boundary,
    boundary_km_over_land,
    edge_adjacency,
    edge_lengths,
    edge_share,
    gate,
    near_cells,
    rotations,
)
from realms.grid import read_centres
from realms.hexmaps import _densified

ROOT = Path(__file__).parent.parent
RES5 = ROOT / "data/reference/grids/ISEA3H-5.gpkg"
GEOD = pyproj.Geod(ellps="WGS84")
SPHERE = pyproj.Geod(a=6371008.8, f=0)
GLOBE = shapely.box(-180, -90, 180, 90)


def vertices(shape) -> list[tuple[float, float]]:
    return [xy for part in shapely.get_parts(shape) for xy in part.exterior.coords[:-1]]


def same_point(p, q) -> bool:
    return abs(p[1] - q[1]) < 1e-9 and abs((p[0] - q[0] + 180) % 360 - 180) < 1e-9


def shared_vertices(a, b) -> list[tuple[float, float]]:
    found = []
    for p in vertices(a):
        if any(same_point(p, q) for q in vertices(b)) and not any(same_point(p, f) for f in found):
            found.append(p)
    return found


def km_between(p, q, geod=GEOD) -> float:
    return geod.inv(p[0], p[1], q[0], q[1])[2] / 1e3


def test_shared_edge_length_of_two_hexagons(tmp_path):
    """A and C are west of the antimeridian, B east of it, and S is cut in two by it, as the grid files hold such a cell."""
    cells = antimeridian_cells()
    edges = edge_lengths(write_grid(tmp_path / "grid.gpkg", cells), cache_dir=tmp_path / "cache")
    assert set(edges) == {tuple(sorted(p)) for p in [(A, S), (A, C), (C, S), (B, S)]}
    shapes = {cell: shape for cell, _, _, shape in cells}
    for (a, b), km in edges.items():
        assert a < b
        p, q = shared_vertices(shapes[a], shapes[b])
        # the band of SNAP degrees that finds the shared edge adds a few centimetres at its ends
        assert km == pytest.approx(km_between(p, q), rel=1e-5)
    # about 111 km for a degree of meridian, S to B across the line as much as A to S
    assert edges[tuple(sorted((B, S)))] == pytest.approx(edges[tuple(sorted((A, S)))], rel=1e-9)
    assert 100 < edges[tuple(sorted((B, S)))] < 120


@pytest.mark.parametrize("noise", [0.0, 1.6e-6, -1.6e-6])
def test_an_edge_over_the_pole_is_the_geodesic_over_it(tmp_path, noise):
    """The cells either side of the north pole at resolution 5 share the edge from (-168.75, 88.8) to (11.25, 88.8), which runs over the pole."""
    north_a, north_b = north_cells(noise)
    gpkg = write_grid(tmp_path / "grid.gpkg", [("1", -78.75, 87.655, north_a), ("2", 101.25, 87.655, north_b)])
    edges = edge_lengths(gpkg, cache_dir=tmp_path / "cache")
    over_the_pole = km_between((-168.75, 88.7955), (11.25, 88.7955))
    assert 2 * (90 - 88.7955) * 111 * 0.99 < over_the_pole < 2 * (90 - 88.7955) * 112
    assert edges == {("1", "2"): pytest.approx(over_the_pole, rel=1e-4)}


def test_an_edge_over_the_pole_cut_by_the_antimeridian_is_whole(tmp_path):
    gpkg = write_grid(tmp_path / "grid.gpkg", [("1", 11.25, 88.795, POLE_WHOLE), ("2", -168.75, 88.795, POLE_CUT)])
    edges = edge_lengths(gpkg, cache_dir=tmp_path / "cache")
    assert edges == {("1", "2"): pytest.approx(km_between((-78.75, 89.217516), (101.25, 89.217516)), rel=1e-4)}


def test_edges_are_cached_while_the_grid_file_is_unchanged(tmp_path, monkeypatch):
    gpkg = write_grid(tmp_path / "grid.gpkg", antimeridian_cells())
    calls = []
    measure = realms.boundaries._measure_edges
    monkeypatch.setattr(realms.boundaries, "_measure_edges", lambda *args: calls.append(1) or measure(*args))
    first = edge_lengths(gpkg, cache_dir=tmp_path / "cache")
    assert (tmp_path / "cache" / "edges-grid.parquet").exists()
    assert edge_lengths(gpkg, cache_dir=tmp_path / "cache") == first
    assert len(calls) == 1
    # a different grid under the same name is measured again
    write_grid(gpkg, antimeridian_cells()[:3])
    assert set(edge_lengths(gpkg, cache_dir=tmp_path / "cache")) == {tuple(sorted((A, C)))}
    assert len(calls) == 2


@pytest.mark.skipif(not RES5.exists(), reason=f"the resolution 5 grid is not at {RES5} (scripts/fetch_reference.py)")
def test_the_edges_of_every_cell_add_up_to_its_perimeter(tmp_path):
    """Every cell of resolution 5, those across the antimeridian and over the poles included, and a closed tiling has 3N - 6 edges."""
    edges = edge_lengths(RES5, cache_dir=tmp_path)
    cells, shapes, _ = realms.boundaries._cells(RES5)
    assert len(edges) == 3 * len(cells) - 6
    around = dict.fromkeys(cells, 0.0)
    for (a, b), km in edges.items():
        around[a] += km
        around[b] += km
    perimeter = {cell: GEOD.geometry_area_perimeter(shape)[1] / 1e3 for cell, shape in zip(cells, shapes, strict=True)}
    worst = max(abs(around[c] / perimeter[c] - 1) for c in cells)
    assert worst < 1e-3


def test_parts_of_one_cell_are_never_a_boundary():
    edges = {("1", "2"): 10.0, ("2", "3"): 20.0, ("1", "3"): 5.0}
    labels = {"1:land": 1, "1:sea": 2, "2:land": 1, "2:sea": 1, "3:sea": 2}
    assert boundary(labels, "land", edges).height == 0
    sea = boundary(labels, "sea", edges)
    assert sea.schema == {"cell_a": pl.String, "cell_b": pl.String, "km": pl.Float64}
    assert sea.rows() == [("1", "2", 10.0), ("2", "3", 20.0)]
    with pytest.raises(ValueError, match="kind"):
        boundary(labels, "coast", edges)


def squares(columns: int, rows: int) -> list[tuple[str, float, float, Polygon]]:
    """One-degree squares from 0 E and 0 N, ids `<row><column>` starting at 10."""
    return [
        (f"{10 + r}{10 + c}", c + 0.5, r + 0.5, shapely.box(c, r, c + 1, r + 1))
        for r in range(rows)
        for c in range(columns)
    ]


def regions_of(boxes) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame({"code": [str(i) for i in range(len(boxes))]}, geometry=list(boxes), crs="EPSG:4326")


def test_a_boundary_along_a_reference_edge_has_share_one(tmp_path):
    gpkg = write_grid(tmp_path / "grid.gpkg", squares(8, 4))
    edges = edge_lengths(gpkg, cache_dir=tmp_path / "cache")
    adjacency = edge_adjacency(edges)
    regions = regions_of([shapely.box(0, 0, 3, 4), shapely.box(3, 0, 8, 4)])
    near = near_cells(regions, gpkg, adjacency)
    # columns 2 and 3 meet both regions, and 1 and 4 are their neighbours
    assert near == {f"{10 + r}{10 + c}" for r in range(4) for c in (1, 2, 3, 4)}

    def split(*cuts):
        return {f"{cell}:land": sum(lon > cut for cut in cuts) for cell, lon, _, _ in squares(8, 4)}

    along = boundary(split(3), "land", edges)
    assert along.height == 4
    assert edge_share(along, near) == 1.0
    assert edge_share(boundary(split(7), "land", edges), near) == 0.0
    both = boundary(split(3, 7), "land", edges)
    assert edge_share(both, near) == pytest.approx(0.5, rel=1e-3)
    assert math.isnan(edge_share(boundary({}, "land", edges), near))


def test_near_cells_reach_across_the_antimeridian(tmp_path):
    """S is drawn from 179 to 181; a region border along 180.5, held as pieces at -180 to -179.5, crosses it."""
    gpkg = write_grid(tmp_path / "grid.gpkg", antimeridian_cells())
    adjacency = edge_adjacency(edge_lengths(gpkg, cache_dir=tmp_path / "cache"))
    west = shapely.MultiPolygon([shapely.box(170, 0, 180, 20), shapely.box(-180, 0, -179.5, 20)])
    regions = regions_of([west, shapely.box(-179.5, 0, -170, 20)])
    assert near_cells(regions, gpkg, adjacency) == {A, B, C, S}
    # one region, even in two pieces, is no border
    assert near_cells(regions_of([shapely.union_all([west, shapely.box(-179.5, 0, -170, 20)])]), gpkg, adjacency) == set()


def sphere_km2(shape) -> float:
    def ring_km2(ring) -> float:
        return abs(SPHERE.geometry_area_perimeter(Polygon(ring))[0]) / 1e6

    return sum(ring_km2(p.exterior) - sum(ring_km2(h) for h in p.interiors) for p in shapely.get_parts(shape))


def assert_drawn_on_the_map(shape):
    assert shape.is_valid and not shape.is_empty
    west, south, east, north = shape.bounds
    assert -180 <= west and east <= 180 and -90 <= south and north <= 90


def test_a_rotated_region_over_a_pole_or_the_antimeridian_stays_valid_and_keeps_its_area():
    cap = shapely.box(-180, -90, 180, -70)  # as the region file draws Antarctica, down to the pole
    across = shapely.MultiPolygon([shapely.box(170, 0, 180, 20), shapely.box(-180, 0, -170, 20)])
    lake = Polygon(shapely.box(0, 60, 30, 80).exterior.coords, [shapely.box(10, 65, 20, 70).exterior.coords])
    rotations_ = [Rotation.of(0, 0, 0), Rotation.of(40, 10, 20), Rotation.of(200, 300, 29), Rotation.of(95, 45, 7)]
    # Edges are straight in longitude and latitude, so the areas are measured with them cut short. pyproj cannot measure the cap as
    # drawn, whose edge along 70 S starts and ends at one point.
    area = {id(cap): 2 * math.pi * 6371.0088**2 * (1 - math.sin(math.radians(70))), id(across): sphere_km2(_densified(across, 0.05)), id(lake): sphere_km2(_densified(lake, 0.05))}
    for rotation in rotations_:
        for shape, inside in ((cap, (30, -80)), (across, (175, 10)), (across, (-175, 10)), (lake, (5, 75))):
            rotated = rotation.geometry(shape)
            assert_drawn_on_the_map(rotated)
            assert sphere_km2(rotated) == pytest.approx(area[id(shape)], rel=2e-3)
            lon, lat = rotation.points([inside[0]], [inside[1]])
            assert rotated.intersects(shapely.Point(lon[0], lat[0]))
        hole_lon, hole_lat = rotation.points([15], [67.5])
        assert not rotation.geometry(lake).intersects(shapely.Point(hole_lon[0], hole_lat[0]))

    # a region that the rotation moves over the north pole covers it from every side
    box = shapely.box(0, 65, 20, 85)
    for axis in (100.0, 280.0):
        rotation = Rotation.of(0, axis, 15)
        if rotation.points([10], [75])[1][0] > 89.9:
            break
    rotated = rotation.geometry(box)
    assert_drawn_on_the_map(rotated)
    assert all(rotated.intersects(shapely.Point(lon, 89.9)) for lon in (-179, -90, 0, 90, 179))
    assert sphere_km2(rotated) == pytest.approx(sphere_km2(_densified(box, 0.05)), rel=2e-3)

    regions = regions_of([cap, across])
    moved = rotations_[1].apply(regions)
    assert list(moved["code"]) == ["0", "1"] and moved.crs == regions.crs


def exact_km_over_land(regions: gpd.GeoDataFrame, land) -> float:
    """The union of the regions' boundaries, without the cut at the antimeridian, intersected with the land: geodesic km."""
    # each region is rotated on its own, so the edge two regions share differs in the last digits: the snap makes it one line
    lines = shapely.intersection(shapely.union_all(shapely.boundary(regions.geometry.values), grid_size=1e-6), land)
    km = 0.0
    for part in shapely.get_parts(lines):
        if part.geom_type != "LineString":
            continue
        xy = np.asarray(part.coords)
        a, b = xy[:-1], xy[1:]
        keep = ~((np.abs(a[:, 0]) > 180 - 1e-6) & (np.abs(b[:, 0]) > 180 - 1e-6))
        km += GEOD.inv(a[keep, 0], a[keep, 1], b[keep, 0], b[keep, 1])[2].sum() / 1e3
    return km


def test_rotations_keep_the_land_length_within_tolerance():
    regions = regions_of([shapely.box(0, 0, 20, 20), shapely.box(20, 0, 40, 20), shapely.box(0, 20, 20, 40), shapely.box(20, 20, 40, 40)])
    land = shapely.box(-30, -10, 70, 60)
    unrotated = boundary_km_over_land(regions, land=land)
    assert unrotated == pytest.approx(exact_km_over_land(regions, land), rel=1e-3)

    kept = rotations(regions, 5, seed=3, land_tolerance=0.2, land=land)
    assert len(kept) == 5
    assert kept == rotations(regions, 5, seed=3, land_tolerance=0.2, land=land)
    for rotation in kept:
        assert 0 <= rotation.shift_deg < 360 and 0 <= rotation.tilt_deg <= 30
        assert abs(rotation.land_km - unrotated) <= 0.2 * unrotated
        assert rotation.land_km == pytest.approx(boundary_km_over_land(regions, rotation, land=land))
        # the rotated polygons agree with the rotated boundary that was measured
        assert exact_km_over_land(rotation.apply(regions), land) == pytest.approx(rotation.land_km, rel=0.01)

    # the kept draws are those of all draws that met the tolerance, in order, and some did not
    every = rotations(regions, 250, seed=3, land_tolerance=math.inf, land=land)
    within = [r for r in every if abs(r.land_km - unrotated) <= 0.2 * unrotated]
    assert [r.matrix for r in kept] == [r.matrix for r in within[:5]]
    assert len(within) < len(every)

    with pytest.raises(ValueError, match="only"):
        rotations(regions, 5, seed=3, land_tolerance=0.2, land=land, max_attempts=4)


# A world for the gate: the grid of resolution 5, every cell on land, and reference regions in sectors of uneven width, each cut at
# uneven latitudes, so that no rotation maps them onto themselves.
SECTORS = [-180.0, -110.0, -35.0, 25.0, 95.0, 180.0]
BANDS = [[-15.0, 40.0], [5.0], [-40.0, 25.0], [10.0, 60.0], [-25.0]]


def reference_regions() -> gpd.GeoDataFrame:
    boxes = []
    for k, (west, east) in enumerate(zip(SECTORS[:-1], SECTORS[1:], strict=True)):
        cuts = [-90.0, *BANDS[k], 90.0]
        boxes += [shapely.box(west, south, east, north) for south, north in zip(cuts[:-1], cuts[1:], strict=True)]
    return regions_of(boxes)


def reference_label(lon: float, lat: float) -> int:
    k = min(int(np.searchsorted(SECTORS, lon, side="right")) - 1, len(BANDS) - 1)
    return 10 * k + int(np.searchsorted(BANDS[k], lat, side="right"))


def grid_label(lon: float, lat: float, meridians: list[float], parallels: list[float]) -> int:
    sector = int(np.searchsorted(meridians, lon, side="right")) % len(meridians)
    return 10 * sector + int(np.searchsorted(parallels, lat, side="right"))


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    if not RES5.exists():
        pytest.skip(f"the resolution 5 grid is not at {RES5} (scripts/fetch_reference.py)")
    edges = edge_lengths(RES5, cache_dir=tmp_path_factory.mktemp("edges"))
    centres = read_centres(RES5).rows()
    regions = reference_regions()

    def labels(rule, kind="land"):
        return {f"{cell}:{kind}": rule(lon, lat) for cell, lat, lon in centres}

    maps = {
        "off": labels(lambda lon, lat: grid_label(lon, lat, [-150, -60, 10, 80, 140], [15])),
        "copies": labels(reference_label),
        "independent": labels(lambda lon, lat: grid_label(lon, lat, [-130, -50, 50, 120], [-35, 50])),
    }
    # sea parts that copy the regions, which the gate must not see
    maps["independent"] |= labels(reference_label, kind="sea")
    drawn = rotations(regions, 39, seed=11, land=GLOBE)
    result = gate(maps, "off", regions, drawn, edges=edges, gpkg=RES5, adjacency=edge_adjacency(edges))
    return {row["mode"]: row for row in result.iter_rows(named=True)}, result


def test_a_map_that_copies_the_reference_regions_fails_the_gate(world):
    rows, result = world
    assert result.columns == ["mode", "share", "raise", "null_q95", "passes"]
    assert result["mode"].to_list() == ["copies", "independent"]
    copies = rows["copies"]
    assert copies["share"] == 1.0
    assert copies["raise"] > copies["null_q95"] + 0.2
    assert not copies["passes"]


def test_a_map_independent_of_the_reference_passes_the_gate(world):
    rows, _ = world
    independent = rows["independent"]
    assert independent["share"] < 0.6
    assert independent["raise"] <= independent["null_q95"]
    assert independent["passes"]


def test_the_gate_needs_its_base_and_a_rotation(tmp_path):
    gpkg = write_grid(tmp_path / "grid.gpkg", squares(2, 1))
    edges = edge_lengths(gpkg, cache_dir=tmp_path / "cache")
    regions = regions_of([shapely.box(0, 0, 1, 1), shapely.box(1, 0, 2, 1)])
    maps = {"native_only": {"1010:land": 1, "1011:land": 2}}
    with pytest.raises(ValueError, match="base"):
        gate(maps, "off", regions, [Rotation.of(10, 0, 5)], edges=edges, gpkg=gpkg, adjacency=edge_adjacency(edges))
    with pytest.raises(ValueError, match="rotation"):
        gate(maps | {"off": maps["native_only"]}, "off", regions, [], edges=edges, gpkg=gpkg, adjacency=edge_adjacency(edges))


def test_without_a_land_tolerance_every_draw_is_kept_and_no_land_is_read(tmp_path):
    regions = regions_of([shapely.box(0, 0, 20, 20), shapely.box(20, 0, 40, 20)])
    land = shapely.box(-30, -10, 70, 60)
    every = rotations(regions, 7, seed=5, land_tolerance=math.inf, land=land)
    unfiltered = rotations(regions, 7, seed=5, land_tolerance=None, land=tmp_path / "no such land.zip")
    assert [r.matrix for r in unfiltered] == [r.matrix for r in every]
    assert all(math.isnan(r.land_km) for r in unfiltered)


def test_the_starting_land_tolerance_is_one_half():
    import inspect

    assert inspect.signature(rotations).parameters["land_tolerance"].default == 0.5
    assert inspect.signature(rotations).parameters["max_tilt_deg"].default == 30.0
