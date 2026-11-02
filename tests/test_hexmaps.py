from pathlib import Path

import geopandas as gpd
import pytest
import shapely
from shapely import Polygon
from test_grid import POLE_CUT, POLE_WHOLE, A, B, C, S, antimeridian_cells, km2, write_grid, write_land

import realms.grid
from realms.grid import part_id
from realms.hexmaps import EQUAL_EARTH, part_geometries, plot_parts

ROOT = Path(__file__).parent.parent
RES5 = ROOT / "data/reference/grids/ISEA3H-5.gpkg"
LAND = ROOT / "data/reference/natural_earth/ne_10m_land.zip"
# Longitude covered by the widest piece we accept: a band around the globe spans 360, a real cell at resolution 5 under 20
MAX_DEGREES = 30


@pytest.fixture
def coastal(tmp_path):
    """A strip of land along 10 N through A, S and B, so each has a land and a sea part. S straddles the antimeridian, C has no land."""
    strip = [Polygon([(177, 9.8), (180, 9.8), (180, 10.2), (177, 10.2)]), Polygon([(-180, 9.8), (-177, 9.8), (-177, 10.2), (-180, 10.2)])]
    return write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), write_land(tmp_path / "land.zip", strip)


def every_part(cells):
    return [part_id(c, kind) for c in cells for kind in ("land", "sea")]


def flat_km2(geometry) -> float:
    """Area on the equal-area map, with the edges cut short enough to follow the curves that straight lines in longitude and latitude become."""
    return gpd.GeoSeries([shapely.segmentize(geometry, 0.01)], crs="EPSG:4326").to_crs(EQUAL_EARTH).area.iloc[0] / 1e6


def assert_no_band(geoms, over_a_pole=False):
    """No piece runs round the globe. A cell over a pole does span half the longitudes, along the pole itself, so it is only held to the map."""
    for part, shape in zip(geoms["part"], geoms.geometry, strict=True):
        for piece in shapely.get_parts(shape):
            west, _, east, _ = piece.bounds
            assert -180 <= west and east <= 180
            assert over_a_pole or east - west < MAX_DEGREES, f"{part} spans {east - west} degrees of longitude"
    if over_a_pole:
        return
    projected = geoms.to_crs(EQUAL_EARTH)
    # 30 degrees of longitude at the equator, with a margin, where a band would be the width of the map (about 34,000 km)
    for part, shape in zip(projected["part"], projected.geometry, strict=True):
        for piece in shapely.get_parts(shape):
            west, _, east, _ = piece.bounds
            assert east - west < 4e6, f"{part} is {(east - west) / 1e3:.0f} km wide on the map"


def test_a_cell_across_the_antimeridian_draws_no_band_around_the_globe(coastal):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, every_part([A, B, C, S]))
    assert geoms.crs == "EPSG:4326"
    assert list(geoms.columns) == ["part", "geometry"]
    assert_no_band(geoms)

    by_part = dict(zip(geoms["part"], geoms.geometry, strict=True))
    sea_of_s = by_part[part_id(S, "sea")]
    # S is split into a piece either side of the line, and nothing is lost by cutting it
    assert {round(p.bounds[0]) for p in shapely.get_parts(sea_of_s)} == {-180, 179}
    cell_s = shapely.union_all([by_part[part_id(S, "land")], sea_of_s])
    assert km2(cell_s) == pytest.approx(km2(antimeridian_cells()[3][3]), rel=1e-6)


def test_a_cell_over_a_pole_draws_a_cap_not_a_band(tmp_path, no_land_file):
    cells = [("1", 11.25, 88.795, POLE_WHOLE), ("2", -168.75, 88.795, POLE_CUT)]
    geoms = part_geometries(write_grid(tmp_path / "pole.gpkg", cells), no_land_file, every_part(["1", "2"]))
    geoms = geoms[~geoms.geometry.is_empty]
    assert sorted(geoms["part"]) == ["1:sea", "2:sea"]
    assert_no_band(geoms, over_a_pole=True)
    def centred(cell, shape, lon_c):
        """The cell as `grid.py` draws it, moved to centre longitude 0 so that it is all within -180 and 180 and can be projected whole."""
        return shapely.transform(realms.grid._planar(cell, shape, lon_c), lambda xy: xy - [lon_c, 0.0])

    drawn = {"1:sea": centred("1", POLE_WHOLE, 11.25), "2:sea": centred("2", POLE_CUT, -168.75)}
    for part, shape in zip(geoms["part"], geoms.geometry, strict=True):
        assert shape.bounds[3] == pytest.approx(90)  # reaches the pole
        # the cell is whole, cap included, which cutting it at the antimeridian must not lose. Its edges are straight in longitude and
        # latitude here, so the area is not the geodesic one, but it is the same for the cell drawn whole and the cell drawn in pieces
        assert flat_km2(shape) == pytest.approx(flat_km2(drawn[part]), rel=1e-6), part
        assert flat_km2(shape) > 1.05 * flat_km2(POLE_WHOLE)  # without the cap over the pole it would be this
    assert flat_km2(geoms.geometry.iloc[0]) == pytest.approx(flat_km2(geoms.geometry.iloc[1]), rel=1e-3)


@pytest.fixture
def no_land_file(tmp_path):
    return write_land(tmp_path / "land.zip", [Polygon([(0, 40), (1, 40), (1, 41), (0, 41)])])


def test_land_and_sea_parts_tile_their_cell(coastal):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, every_part([A, B, C, S]))
    by_part = dict(zip(geoms["part"], geoms.geometry, strict=True))
    cells = {c: shape for c, _, _, shape in antimeridian_cells()}
    land_seen = 0
    for cell in (A, B, C, S):
        land_part, sea_part = by_part[part_id(cell, "land")], by_part[part_id(cell, "sea")]
        assert km2(shapely.intersection(land_part, sea_part)) == pytest.approx(0, abs=1e-6)
        assert km2(land_part) + km2(sea_part) == pytest.approx(km2(cells[cell]), rel=1e-6)
        land_seen += km2(land_part) > 0
    assert land_seen == 3 and km2(by_part[part_id(C, "land")]) == 0
    # land is what grid.py counts as land
    grid = realms.grid.load_grid(grid_file, land)
    for row in grid.cells.iter_rows(named=True):
        assert km2(by_part[part_id(row["cell"], "land")]) == pytest.approx(row["land_km2"], abs=1e-3)
        assert km2(by_part[part_id(row["cell"], "sea")]) == pytest.approx(row["sea_km2"], abs=1e-3)


@pytest.mark.skipif(not (RES5.exists() and LAND.exists()), reason="needs data/reference/grids/ISEA3H-5.gpkg and the Natural Earth land")
def test_a_coastal_cell_of_the_real_grid_is_tiled_by_its_parts():
    centres = realms.grid.read_centres(RES5)
    # near the coast of Lower Saxony, the sea and the land both within one cell
    nearest = centres.with_columns(((centres["lon"] - 8.0) ** 2 + (centres["lat"] - 53.8) ** 2).alias("d")).sort("d").row(0, named=True)["cell"]
    geoms = part_geometries(RES5, LAND, every_part([nearest]))
    land_part, sea_part = geoms.geometry
    assert km2(land_part) > 0 and km2(sea_part) > 0
    assert km2(shapely.intersection(land_part, sea_part)) == pytest.approx(0, abs=1e-3)
    frame = gpd.read_file(RES5)
    row = frame[frame["geogrid_id"].astype(str) == nearest].iloc[0]
    cell = realms.grid._planar(nearest, row.geometry, row["lon_c"])
    # on the map, where cutting an edge in two does not change the area as it does for the geodesic one
    assert flat_km2(land_part) + flat_km2(sea_part) == pytest.approx(flat_km2(cell), rel=1e-5)
    assert_no_band(geoms)


def test_parts_that_are_not_in_the_grid_are_refused(coastal):
    grid_file, land = coastal
    with pytest.raises(ValueError, match="not in the grid"):
        part_geometries(grid_file, land, ["12345:land"])
    with pytest.raises(ValueError, match="not a part id"):
        part_geometries(grid_file, land, [A])


# Realms for the four cells of `coastal`: 1 is the land of A, S and B, 2 their sea, 3 the cell above
LABELS = {part_id(c, "land"): 1 for c in (A, S, B)} | {part_id(c, "sea"): 2 for c in (A, S, B)} | {part_id(C, "sea"): 3}


def test_the_dissolved_layer_has_one_feature_per_realm(coastal, tmp_path):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, list(LABELS))
    out = tmp_path / "maps" / "realms"
    plot_parts(geoms, LABELS, out, "Realms", sub_realms={part_id(A, "sea"): 1, part_id(B, "sea"): 2}, imputed={part_id(C, "sea")})

    assert (tmp_path / "maps" / "realms.png").stat().st_size > 5000
    parts = gpd.read_file(tmp_path / "maps" / "realms.gpkg", layer="parts")
    dissolved = gpd.read_file(tmp_path / "maps" / "realms.gpkg", layer="realms")
    assert parts.crs == "EPSG:4326" == dissolved.crs
    assert sorted(parts["part"]) == sorted(LABELS)
    assert {"part", "realm", "sub_realm", "imputed"} <= set(parts.columns)
    by_part = parts.set_index("part")
    assert by_part.loc[part_id(A, "sea"), "sub_realm"] == 1 and by_part.loc[part_id(B, "sea"), "sub_realm"] == 2
    assert by_part["sub_realm"].isna().sum() == len(parts) - 2
    assert by_part["imputed"].to_dict() == {p: p == part_id(C, "sea") for p in LABELS}
    assert by_part["realm"].to_dict() == LABELS

    assert sorted(dissolved["realm"]) == [1, 2, 3]
    assert dissolved["dissolved"].tolist() == [3, 3, 1]
    for shape in dissolved.geometry:
        assert shape.geom_type in ("Polygon", "MultiPolygon") and shape.is_valid
    # nothing is lost or counted twice by dissolving
    assert km2(shapely.union_all(dissolved.geometry)) == pytest.approx(km2(shapely.union_all(parts.geometry)), rel=1e-6)
    for realm, shape in zip(dissolved["realm"], dissolved.geometry, strict=True):
        assert km2(shape) == pytest.approx(km2(shapely.union_all(parts.geometry[parts["realm"] == realm])), rel=1e-6)


def test_sub_realms_and_imputed_default_to_none_and_the_picture_is_the_same_every_time(coastal, tmp_path):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, list(LABELS))
    for name in ("a", "b"):
        plot_parts(geoms, LABELS, tmp_path / name, "Realms")
    parts = gpd.read_file(tmp_path / "a.gpkg", layer="parts")
    assert parts["sub_realm"].isna().all() and not parts["imputed"].any()
    assert (tmp_path / "a.png").read_bytes() == (tmp_path / "b.png").read_bytes()


def test_a_part_with_no_label_is_refused(coastal, tmp_path):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, list(LABELS))
    with pytest.raises(ValueError, match="no realm"):
        plot_parts(geoms, {p: 1 for p in list(LABELS)[1:]}, tmp_path / "x", "Realms")


def test_realms_are_coloured_apart_from_their_neighbours_by_the_adjacency_given(coastal, tmp_path, monkeypatch):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, list(LABELS))
    seen = {}
    import realms.hexmaps as hexmaps

    real = hexmaps._separating_colours

    def spy(realms_, neighbours):
        seen["n"] = {k: set(v) for k, v in neighbours.items() if v}
        return real(realms_, neighbours)

    monkeypatch.setattr(hexmaps, "_separating_colours", spy)
    adjacency = {p: set() for p in LABELS}
    adjacency[part_id(A, "land")] = {part_id(A, "sea")}
    adjacency[part_id(A, "sea")] = {part_id(A, "land")}
    plot_parts(geoms, LABELS, tmp_path / "m", "Realms", adjacency=adjacency)
    assert seen["n"] == {1: {2}, 2: {1}}


def test_a_part_that_vanishes_in_projection_is_refused(coastal, tmp_path, monkeypatch):
    """A part with area that projects to nothing would be a silent hole in the map: plot_parts raises instead."""
    import realms.hexmaps as hexmaps

    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, list(LABELS))
    target = dict(zip(geoms["part"], geoms.geometry, strict=True))[part_id(A, "sea")]
    real = hexmaps._densified

    def losing(shape, step=hexmaps.SEGMENT):
        return shapely.MultiPolygon() if shape.equals(target) else real(shape, step)

    monkeypatch.setattr(hexmaps, "_densified", losing)
    with pytest.raises(ValueError, match="vanish") as error:
        plot_parts(geoms, LABELS, tmp_path / "m", "Realms")
    assert part_id(A, "sea") in str(error.value)
    assert not (tmp_path / "m.png").exists()


# A coastal sea part at resolution 5 that shapely.segmentize (GEOS 3.14) turned into an empty polygon before `_densified` replaced it
VANISHED = "4919636839102606602:sea"


@pytest.mark.skipif(not (RES5.exists() and LAND.exists()), reason="needs data/reference/grids/ISEA3H-5.gpkg and the Natural Earth land")
def test_the_part_that_segmentize_emptied_is_drawn_and_would_be_refused(tmp_path, monkeypatch):
    import realms.hexmaps as hexmaps

    geoms = part_geometries(RES5, LAND, [VANISHED])
    assert not geoms.geometry.iloc[0].is_empty
    plot_parts(geoms, {VANISHED: 1}, tmp_path / "drawn", "One part")  # the densifying in use keeps it
    assert (tmp_path / "drawn.png").exists()
    monkeypatch.setattr(hexmaps, "_densified", lambda shape, step=hexmaps.SEGMENT: shapely.segmentize(shape, step))
    with pytest.raises(ValueError, match=f"vanish.*{VANISHED}"):
        plot_parts(geoms, {VANISHED: 1}, tmp_path / "lost", "One part")
