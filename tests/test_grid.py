import math
import zipfile

import geopandas as gpd
import numpy as np
import polars as pl
import pyproj
import pytest
import shapely
from shapely import Polygon

import realms.grid
from realms.grid import Grid, load_grid, part_adjacency, part_id, parts

GEOD = pyproj.Geod(ellps="WGS84")
R = 1.0
STEP = math.sqrt(3) * R  # between the centres of two hexagons side by side


def km2(geometry) -> float:
    return sum(abs(GEOD.geometry_area_perimeter(p)[0]) for p in shapely.get_parts(geometry)) / 1e6


def hexagon(lon: float, lat: float) -> Polygon:
    """A pointy-top hexagon of circumradius 1 degree, so that its neighbours sit at (+-STEP, 0) and (+-STEP/2, +-1.5)."""
    return Polygon([(lon + R * math.cos(math.radians(a)), lat + R * math.sin(math.radians(a))) for a in (90, 150, 210, 270, 330, 30)])


def split_at_antimeridian(shape: Polygon):
    """The hexagon, with longitudes up to 181, cut at 180 into two pieces the way the grid files hold it."""
    east = shapely.clip_by_rect(shape, -999, -90, 180, 90)
    west = shapely.clip_by_rect(shape, 180, -90, 999, 90)
    return shapely.MultiPolygon([east, shapely.transform(west, lambda c: c - [360.0, 0.0])])


# Four ids that no float can tell apart
A, B, C, S = "658282526011250001", "-658282526011250002", "658282526011250003", "658282526011250004"
LAT = 10.0
WEST_OF_S = hexagon(180 - STEP, LAT)  # A
EAST_OF_S = hexagon(-180 + STEP, LAT)  # B, seen from the other side of the line
ABOVE = hexagon(180 - STEP / 2, LAT + 1.5 * R)  # C touches A and S, and has an edge on the line itself
SPLIT = hexagon(180.0, LAT)  # S


def write_grid(path, cells, ids_as_integers=False):
    """`cells`: (id, centre lon, centre lat, geometry). Columns as in the GBIF files."""
    ids = [int(i) if ids_as_integers else i for i, *_ in cells]
    frame = gpd.GeoDataFrame(
        {
            "seqnum": [float(n) for n in range(1, len(cells) + 1)],
            "lat_c": [lat for _, _, lat, _ in cells],
            "lon_c": [lon for _, lon, _, _ in cells],
            "res": [6.0] * len(cells),
            "geogrid_id": ids,
        },
        geometry=[g for *_, g in cells],
        crs="EPSG:4326",
    )
    frame.to_file(path, driver="GPKG")
    return path


def write_land(path, shapes):
    """A zipped shapefile, as Natural Earth's land comes."""
    folder = path.with_suffix("")
    folder.mkdir()
    gpd.GeoDataFrame({"featurecla": ["Land"] * len(shapes)}, geometry=list(shapes), crs="EPSG:4326").to_file(folder / "land.shp", driver="ESRI Shapefile")
    with zipfile.ZipFile(path, "w") as archive:
        for file in sorted(folder.iterdir()):
            archive.write(file, file.name)
    return path


def antimeridian_cells(centre_lon=180.0):
    return [
        (A, 180 - STEP, LAT, WEST_OF_S),
        (B, -180 + STEP, LAT, EAST_OF_S),
        (C, 180 - STEP / 2, LAT + 1.5 * R, ABOVE),
        (S, centre_lon, LAT, split_at_antimeridian(SPLIT)),
    ]


@pytest.fixture
def no_land(tmp_path):
    return write_land(tmp_path / "land.zip", [Polygon([(0, 40), (1, 40), (1, 41), (0, 41)])])


@pytest.fixture
def grid(tmp_path, no_land) -> Grid:
    return load_grid(write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), no_land)


def test_ids_stay_text(tmp_path, no_land):
    for as_integers in (False, True):
        grid = load_grid(write_grid(tmp_path / f"{as_integers}.gpkg", antimeridian_cells(), as_integers), no_land)
        assert grid.cells["cell"].to_list() == [A, B, C, S]
        assert set(grid.adjacency) == {A, B, C, S}


def test_floating_ids_are_refused(tmp_path, no_land):
    path = write_grid(tmp_path / "grid.gpkg", antimeridian_cells())
    frame = gpd.read_file(path)
    frame["geogrid_id"] = frame["geogrid_id"].astype(float)
    frame.to_file(path, driver="GPKG")
    with pytest.raises(ValueError, match="geogrid_id"):
        load_grid(path, no_land)


def test_the_centres_are_read_without_the_geometry_as_the_grid_has_them(tmp_path, no_land, grid):
    path = write_grid(tmp_path / "grid.gpkg", antimeridian_cells())
    centres = realms.grid.read_centres(path)
    assert centres.schema == {"cell": pl.String, "lat": pl.Float64, "lon": pl.Float64}
    assert centres.equals(grid.cells.select("cell", "lat", "lon"))
    for as_integers in (False, True):  # ids that no float can tell apart stay text
        assert realms.grid.read_centres(write_grid(tmp_path / f"{as_integers}.gpkg", antimeridian_cells(), as_integers))["cell"].to_list() == [A, B, C, S]


def test_floating_ids_are_refused_by_the_centres_too(tmp_path):
    path = write_grid(tmp_path / "grid.gpkg", antimeridian_cells())
    frame = gpd.read_file(path)
    frame["geogrid_id"] = frame["geogrid_id"].astype(float)
    frame.to_file(path, driver="GPKG")
    with pytest.raises(ValueError, match="geogrid_id"):
        realms.grid.read_centres(path)


def test_cells(grid):
    assert grid.resolution == 6
    assert grid.cells.columns == ["cell", "lat", "lon", "area_km2", "land_km2", "sea_km2"]
    assert grid.cells.schema["cell"] == pl.String
    row = grid.cells.filter(pl.col("cell") == A).row(0, named=True)
    assert (row["lat"], row["lon"]) == (LAT, 180 - STEP)


def test_adjacency_is_symmetric_and_has_the_right_neighbours(grid):
    adjacency = grid.adjacency
    assert all(cell not in near for cell, near in adjacency.items())
    assert all(cell in adjacency[other] for cell, near in adjacency.items() for other in near)
    # B touches S only, across the line; A and C touch each other and S; A and B are two cells apart
    assert adjacency == {A: {C, S}, B: {S}, C: {A, S}, S: {A, B, C}}


@pytest.mark.parametrize("centre_lon", [180.0, -180.0])
def test_a_split_cell_finds_neighbours_on_both_sides(tmp_path, no_land, centre_lon):
    grid = load_grid(write_grid(tmp_path / "grid.gpkg", antimeridian_cells(centre_lon)), no_land)
    # A and C are drawn at lon 178 to 180 and B at -180 to -178, with S in two pieces between them
    assert grid.adjacency[S] == {A, B, C}


def test_cells_meeting_in_a_point_are_not_neighbours(tmp_path, no_land):
    """An L of three squares: 1 and 2 only share the corner (1, 1). Square 3 is 4e-9 degrees off square 1, as neighbours are in the files."""
    cells = [
        ("1", 0.5, 0.5, Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])),
        ("2", 1.5, 1.5, Polygon([(1, 1), (2, 1), (2, 2), (1, 2)])),
        ("3", 1.5, 0.5, Polygon([(1 + 4e-9, 0), (2, 0), (2, 1), (1 + 4e-9, 1)])),
    ]
    assert load_grid(write_grid(tmp_path / "grid.gpkg", cells), no_land).adjacency == {"1": {"3"}, "2": {"3"}, "3": {"1", "2"}}


def test_areas_match_pyproj(grid):
    expected = {A: km2(WEST_OF_S), B: km2(EAST_OF_S), C: km2(ABOVE), S: km2(SPLIT)}
    for row in grid.cells.iter_rows(named=True):
        assert row["area_km2"] == pytest.approx(expected[row["cell"]], rel=1e-3)
        assert row["land_km2"] + row["sea_km2"] == pytest.approx(row["area_km2"])


# The two cells either side of the north pole at resolution 5, as the file draws them: both have an edge along 88.8 N from lon -168.75
# to 11.25 that really runs over the pole, and the second one is drawn the long way round, so it crosses itself.
def north_cells(noise: float):
    """With `noise` degrees on the longitude of the first vertex, which both cells share. Polar longitudes in the files are off by up to 2e-6."""
    v0 = (-168.75 - noise, 88.7955)
    a = Polygon([v0, (-124.608, 86.6306), (-93.283, 85.1627), (-64.217, 85.1627), (-32.892, 86.6306), (11.25, 88.7955)])
    b = Polygon([v0, (11.25, 88.7955), (55.392, 86.6306), (86.717, 85.1627), (115.783, 85.1627), (147.108, 86.6306)])
    return a, b


@pytest.mark.parametrize("noise", [0.0, 1.6e-6, -1.6e-6])
def test_two_cells_meeting_across_a_pole_are_neighbours(tmp_path, no_land, noise):
    north_a, north_b = north_cells(noise)
    assert not north_b.is_valid
    grid = load_grid(write_grid(tmp_path / "grid.gpkg", [("1", -78.75, 87.655, north_a), ("2", 101.25, 87.655, north_b)]), no_land)
    assert grid.adjacency == {"1": {"2"}, "2": {"1"}}
    assert grid.cells["area_km2"].to_list() == pytest.approx([km2(north_a), km2(north_b)], rel=1e-3)


# The same at resolution 6, where the cell on the far side is split at the antimeridian: its pole edge is cut in two at 180.
POLE_WHOLE = Polygon([(-41.083541, 88.0238), (-6.737647, 87.464256), (29.237647, 87.464256), (63.583541, 88.0238), (101.25, 89.217516), (-78.75, 89.217516)])
POLE_CUT = shapely.MultiPolygon(
    [
        Polygon([(173.262353, 87.464256), (138.916459, 88.0238), (101.25, 89.217516), (180.0, 89.217516), (180.0, 87.464256)]),
        Polygon([(-78.75, 89.217516), (-116.416459, 88.0238), (-150.762353, 87.464256), (-180.0, 87.464256), (-180.0, 89.217516)]),
    ]
)


def test_a_pole_edge_cut_by_the_antimeridian_is_whole_again(tmp_path, no_land):
    cells = [("1", 11.25, 88.795, POLE_WHOLE), ("2", -168.75, 88.795, POLE_CUT)]
    grid = load_grid(write_grid(tmp_path / "grid.gpkg", cells), no_land)
    assert grid.adjacency == {"1": {"2"}, "2": {"1"}}
    # the cells are the same size, which drawn in pieces they are not: the pieces leave out the cap over the pole
    assert km2(POLE_CUT) < 0.95 * km2(POLE_WHOLE)
    assert grid.cells["area_km2"].to_list() == pytest.approx([km2(POLE_WHOLE)] * 2, rel=1e-3)


def test_land_is_the_geodesic_area_of_the_overlap(tmp_path):
    # B lies wholly on land. S has its west piece wholly on land, and a box inside its east piece. A and C have none, but C's edge
    # on the antimeridian touches the land that starts at -180.
    covers_b_and_west_of_s = Polygon([(-180, 8), (-176, 8), (-176, 13), (-180, 13)])
    inside_east_of_s = Polygon([(179.5, 9.6), (180, 9.6), (180, 10.4), (179.5, 10.4)])
    land = write_land(tmp_path / "land.zip", [covers_b_and_west_of_s, inside_east_of_s])
    grid = load_grid(write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), land)
    cells = {row["cell"]: row for row in grid.cells.iter_rows(named=True)}

    assert cells[A]["land_km2"] == 0
    assert cells[A]["sea_km2"] == pytest.approx(cells[A]["area_km2"])
    assert cells[C]["land_km2"] == pytest.approx(0, abs=1e-6)
    assert cells[B]["land_km2"] == pytest.approx(cells[B]["area_km2"], rel=1e-6)
    assert cells[B]["sea_km2"] == pytest.approx(0, abs=1e-6)
    west_piece = shapely.clip_by_rect(SPLIT, 180, -90, 999, 90)
    assert cells[S]["land_km2"] == pytest.approx(km2(west_piece) + km2(inside_east_of_s), rel=1e-3)
    assert cells[S]["sea_km2"] == pytest.approx(cells[S]["area_km2"] - cells[S]["land_km2"])


def test_land_never_exceeds_the_cell(tmp_path):
    everything = Polygon([(-180, -90), (180, -90), (180, 90), (-180, 90)])
    grid = load_grid(write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), write_land(tmp_path / "land.zip", [everything]))
    assert (grid.cells["land_km2"] <= grid.cells["area_km2"]).all()
    assert (grid.cells["sea_km2"] >= 0).all()
    assert grid.cells["land_km2"].to_list() == pytest.approx(grid.cells["area_km2"].to_list(), rel=1e-6)


def test_land_is_kept_within_the_cell(tmp_path, no_land, monkeypatch):
    """Overlap areas and cell areas come from different polygons, so the overlap can come out a little over or under."""
    monkeypatch.setattr(realms.grid, "_land_km2", lambda shapes, land: np.array([1e12, -5.0, 0.0, 10.0]))
    grid = load_grid(write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), no_land)
    cells = {row["cell"]: row for row in grid.cells.iter_rows(named=True)}
    assert cells[A]["land_km2"] == cells[A]["area_km2"]
    assert cells[A]["sea_km2"] == 0
    assert cells[B]["land_km2"] == 0
    assert cells[B]["sea_km2"] == cells[B]["area_km2"]
    assert cells[S]["land_km2"] == 10


def test_part_id():
    assert part_id("-658282526011250001", "land") == "-658282526011250001:land"
    assert part_id("7", "sea") == "7:sea"
    with pytest.raises(ValueError, match="coast"):
        part_id("7", "coast")


@pytest.fixture
def coastal(tmp_path):
    """A strip of land along 10 N through A, S and B, so that each of them has a land part and a sea part. C has no land."""
    strip = [Polygon([(177, 9.8), (180, 9.8), (180, 10.2), (177, 10.2)]), Polygon([(-180, 9.8), (-177, 9.8), (-177, 10.2), (-180, 10.2)])]
    land = write_land(tmp_path / "land.zip", strip)
    return load_grid(write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), land)


def test_parts(coastal):
    present = [part_id(A, "land"), part_id(A, "sea"), part_id(B, "sea"), part_id(A, "sea")]
    table = parts(coastal, present)
    assert table.columns == ["part", "cell", "kind", "area_km2"]
    assert table["part"].to_list() == sorted({part_id(A, "land"), part_id(A, "sea"), part_id(B, "sea")})
    cell_a = coastal.cells.filter(pl.col("cell") == A).row(0, named=True)
    area = dict(zip(table["part"], table["area_km2"], strict=True))
    assert area[part_id(A, "land")] == pytest.approx(cell_a["land_km2"])
    assert area[part_id(A, "sea")] == pytest.approx(cell_a["sea_km2"])
    assert table.filter(pl.col("part") == part_id(B, "sea"))["kind"].to_list() == ["sea"]
    assert parts(coastal, []).columns == ["part", "cell", "kind", "area_km2"]


def test_parts_refuse_what_is_not_a_part(coastal):
    with pytest.raises(ValueError, match="nowhere"):
        parts(coastal, ["nowhere:land"])
    with pytest.raises(ValueError, match="coast"):
        parts(coastal, [f"{A}:coast"])
    with pytest.raises(ValueError, match="not a part"):
        parts(coastal, [A])


def test_part_adjacency(coastal):
    present = [part_id(c, k) for c in (A, S, C) for k in ("land", "sea")] + [part_id(B, "sea")]
    near = part_adjacency(coastal, present)
    assert set(near) == set(present)
    # the same kind in a neighbouring cell, and the other part of the same cell
    assert near[part_id(A, "land")] == {part_id(S, "land"), part_id(C, "land"), part_id(A, "sea")}
    assert part_id(S, "sea") not in near[part_id(A, "land")]
    assert part_id(C, "sea") not in near[part_id(A, "land")]
    # B has only a sea part here, so nothing of kind land is next to it, and it meets S's sea
    assert near[part_id(B, "sea")] == {part_id(S, "sea")}
    assert all(part in near[other] for part, others in near.items() for other in others)


def test_part_adjacency_only_joins_what_is_present(coastal):
    near = part_adjacency(coastal, [part_id(A, "land"), part_id(S, "sea"), part_id(C, "land")])
    assert near == {part_id(A, "land"): {part_id(C, "land")}, part_id(C, "land"): {part_id(A, "land")}, part_id(S, "sea"): set()}
    assert part_adjacency(coastal, []) == {}
