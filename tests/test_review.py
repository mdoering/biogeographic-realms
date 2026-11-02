import geopandas as gpd
import numpy as np
import polars as pl
import pytest
import shapely
from helpers import grid, matrix, presence
from shapely import Polygon
from test_grid import A, B, C, S, antimeridian_cells, write_grid, write_land

from realms.boundaries import _line_km, edge_lengths
from realms.candidates import candidate_maps, merge_sequence
from realms.cluster import Partition
from realms.endemism import realm_endemism
from realms.grid import part_id
from realms.hexmaps import part_geometries
from realms.maps import HUES
from realms.matrix import cell_taxa
from realms.review import (
    distinct_colours,
    edge_changes,
    edge_lines,
    endemic_species,
    markdown,
    match_colours,
    overlap,
    part_areas,
    part_neighbours,
    plot_boundary_change,
    plot_realm,
    realm_neighbours,
    realm_table,
    sub_realm_table,
    unfloored_maps,
    view,
)

# ------------------------------------------------------------------------------------------------ colours

AREAS = {p: 1.0 for p in "abcdefgh"}


def test_a_realm_takes_the_colour_of_its_largest_overlap_parent_or_child():
    k7 = dict(zip("abcdefgh", [0, 0, 0, 1, 1, 2, 2, 2]))
    colours = {0: "red", 1: "blue", 2: "green"}
    # finer: realm 0 splits into 0 (two parts) and 9 (one part); 9 is new and takes an unused colour
    k8 = dict(zip("abcdefgh", [0, 0, 9, 1, 1, 2, 2, 2]))
    found = match_colours(k7, colours, k8, AREAS, palette=["red", "blue", "green", "pink"])
    assert found == {0: "red", 1: "blue", 2: "green", 9: "pink"}
    # coarser: 1 and 2 merge into 1; it holds more of 2 than of 1, so it takes 2's colour, one to one
    k6 = dict(zip("abcdefgh", [0, 0, 0, 1, 1, 1, 1, 1]))
    assert match_colours(k7, colours, k6, AREAS) == {0: "red", 1: "green"}


def test_the_split_off_realm_avoids_its_neighbours_colours_when_the_palette_runs_out():
    k7 = dict(zip("abcd", [0, 0, 1, 1]))
    k8 = dict(zip("abcd", [0, 5, 1, 1]))
    found = match_colours(k7, {0: "red", 1: "blue"}, k8, AREAS, neighbours={5: {0}}, palette=["red", "blue"])
    assert found == {0: "red", 1: "blue", 5: "blue"}
    assert match_colours(k7, {0: "red", 1: "blue"}, k8, AREAS, neighbours={5: {0, 1}}, palette=["red", "blue"])[5] == "red"


def test_distinct_colours_go_to_the_largest_realms_first():
    labels = dict(zip("abcd", [3, 1, 1, 2]))
    areas = {"a": 5.0, "b": 1.0, "c": 1.0, "d": 2.0}
    assert distinct_colours(labels, areas, ["x", "y", "z"]) == {3: "x", 1: "y", 2: "z"}
    with pytest.raises(ValueError, match="only 2 colours"):
        distinct_colours(labels, areas, ["x", "y"])


def test_overlap_counts_area_over_the_parts_both_maps_hold():
    assert overlap({"a": 0, "b": 1, "x": 3}, {"a": 5, "b": 5, "y": 4}, {"a": 2.0, "b": 3.0}) == {(0, 5): 2.0, (1, 5): 3.0}


def test_part_neighbours_are_the_same_kind_next_door_and_the_other_kind_in_the_cell():
    near = part_neighbours(["1:land", "1:sea", "2:land", "3:sea"], {"1": {"2", "3"}, "2": {"1"}, "3": {"1"}})
    assert near == {"1:land": {"1:sea", "2:land"}, "1:sea": {"1:land", "3:sea"}, "2:land": {"1:land"}, "3:sea": {"1:sea"}}
    assert realm_neighbours({"1:land": 0, "1:sea": 1, "2:land": 0, "3:sea": 1}, near) == {0: {1}, 1: {0}}


# ------------------------------------------------------------------------------------------------ tables


def test_endemic_species_follow_realm_share_on_the_species_columns():
    rows = [("a", "s1"), ("b", "s2"), ("b", "w"), ("a", "w2"), ("c", "w"), ("d", "s3")]
    labels = {"a": 0, "c": 0, "b": 1, "d": 1}
    found = endemic_species(labels, matrix(rows))
    # realm 0 holds s1, w2, w: s1 and w2 endemic; realm 1 holds s2, w, s3: s2 and s3
    assert found == {0: (2, 3), 1: (2, 3)}
    shares = realm_endemism(labels, matrix(rows))
    assert {r: e / p for r, (e, p) in found.items()} == pytest.approx(shares)
    # a species with 19 of its 20 units in a realm is endemic at 95%, with 18 it is not
    many = [(f"u{i}", "x") for i in range(20)] + [(f"u{i}", "y") for i in range(20)]
    at = {f"u{i}": 0 if i < 19 else 1 for i in range(20)}
    assert endemic_species(at, matrix(many)) == {0: (2, 2), 1: (0, 2)}
    at = {f"u{i}": 0 if i < 18 else 1 for i in range(20)}
    assert endemic_species(at, matrix(many))[0] == (0, 2)


def test_endemic_species_count_only_species_columns():
    taxa = pl.DataFrame({"species_id": ["s1", "s2"], "genus_id": ["G", "G"], "family_id": ["F", "F"]})
    ct = cell_taxa(presence([("a", "s1"), ("b", "s2")]), taxa)
    assert endemic_species({"a": 0, "b": 1}, ct) == {0: (1, 1), 1: (1, 1)}


def _map():
    """Two cells, each with land and sea; realm 0 holds both land parts (units u1, u2, with u2 a sub-realm), realm 1 the sea, one part imputed."""
    parts = pl.DataFrame(
        {
            "part": ["1:land", "2:land", "1:sea", "2:sea"],
            "k": 7,
            "realm": [0, 0, 1, 1],
            "unit": ["u1:land", "u2:land", "u3:sea", None],
            "imputed": [False, False, False, True],
            "selected": True,
        }
    )
    candidates = pl.DataFrame(
        {"unit": ["u1:land", "u2:land", "u3:sea"], "k": 7, "realm": [0, 0, 1], "sub_realm": [None, 40, None], "sub_realm_of": [None, 0, None]},
        schema={"unit": pl.String, "k": pl.Int64, "realm": pl.Int64, "sub_realm": pl.Int64, "sub_realm_of": pl.Int64},
    )
    cells = pl.DataFrame({"cell": ["1", "2"], "lat": [0.0, 1.0], "lon": [0.0, 1.0], "area_km2": [10.0, 20.0], "land_km2": [4.0, 15.0], "sea_km2": [6.0, 5.0]})
    return parts, candidates, cells


def test_the_realm_table_adds_up_parts_units_endemics_and_sub_realms():
    parts, candidates, cells = _map()
    areas = part_areas(cells)
    assert areas == {"1:land": 4.0, "1:sea": 6.0, "2:land": 15.0, "2:sea": 5.0}
    table = realm_table(parts, candidates, areas, {"u1:land": 4.0, "u2:land": 15.0, "u3:sea": 6.0}, {0: (3, 10), 1: (1, 2)})
    rows = {r["realm"]: r for r in table.iter_rows(named=True)}
    assert table["realm"].to_list() == [0, 1]  # largest first
    assert rows[0] | {} == {
        "realm": 0, "area_km2": 19.0, "land_km2": 19.0, "sea_km2": 0.0, "area_share": pytest.approx(19 / 30), "units": 2, "land_units": 2,
        "sea_units": 0, "parts": 2, "imputed_parts": 0, "species": 10, "endemic_species": 3, "sub_realms": 1, "sub_realm_units": 1, "sub_realm_km2": 15.0,
    }
    assert rows[1]["sea_km2"] == 11.0 and rows[1]["imputed_parts"] == 1 and rows[1]["units"] == 1 and rows[1]["sub_realms"] == 0
    assert rows[1]["sub_realm_km2"] == 0.0
    assert sub_realm_table(candidates, {"u1:land": 4.0, "u2:land": 15.0, "u3:sea": 6.0}).to_dicts() == [
        {"realm": 0, "sub_realm": 40, "units": 1, "area_km2": 15.0}
    ]


def test_the_realm_table_refuses_parts_and_candidates_that_disagree():
    parts, candidates, cells = _map()
    wrong = candidates.with_columns(pl.when(pl.col("unit") == "u1:land").then(1).otherwise(pl.col("realm")).alias("realm"))
    with pytest.raises(ValueError, match="not those of the candidates"):
        realm_table(parts, wrong, part_areas(cells), {}, {})
    with pytest.raises(ValueError, match="one k"):
        realm_table(pl.concat([parts, parts.with_columns(pl.col("k") + 1)]), candidates, part_areas(cells), {}, {})


def test_markdown_formats_columns_and_leaves_nulls_empty():
    table = pl.DataFrame({"k": [6, 7], "x": [0.12345, None], "note": ["a", "b"]})
    assert markdown(table, {"x": ".2f"}) == "| k | x | note |\n|---|---|---|\n| 6 | 0.12 | a |\n| 7 |  | b |\n"


# ------------------------------------------------------------------------------------------------ without the floor


def test_unfloored_maps_merge_the_modules_without_a_floor():
    cells, adjacency = grid(4, 2)
    # four modules of two cells each, the columns of the grid; module 3's cells are tiny, so a floor would absorb it
    module = {c: int(c.split(",")[0]) for c in cells}
    rows = [(c, f"m{module[c]}") for c in cells] + [("0,0", "x01"), ("1,0", "x01"), ("2,0", "x23"), ("3,0", "x23")]
    ct = matrix(rows)
    partition = Partition(paths={c: (module[c], 0) for c in cells}, unassigned=[])
    units, maps = unfloored_maps(partition, ct, adjacency, [2, 4, 9])
    assert units == partition.leaf_labels()  # the modules are contiguous, so resolving pieces moves nothing
    assert sorted(maps) == [2, 4]  # no map has more realms than there are modules
    assert maps[4] == units
    assert maps[2] == merge_sequence(units, ct, k_max=4, k_min=2)[2]
    # with the floor, module 3 is a sub-realm and there is no map of 4 realms; without it there is
    areas = {c: (0.1 if module[c] == 3 else 10.0) for c in cells}
    floored = candidate_maps(partition, ct, areas, adjacency, k_max=4, k_min=2, floor_km2=1.0)
    assert sorted(floored.maps) == [2, 3] and set(floored.sub_realms) == {units["3,0"]}
    assert unfloored_maps(partition, ct, adjacency, [5])[1] == {}


# ------------------------------------------------------------------------------------------------ boundaries and figures


def test_edge_changes_sort_edges_into_kept_added_and_removed():
    base = pl.DataFrame({"cell_a": ["1", "1", "2"], "cell_b": ["2", "3", "3"], "km": [1.0, 2.0, 3.0]})
    other = pl.DataFrame({"cell_a": ["1", "3"], "cell_b": ["2", "4"], "km": [1.0, 4.0]})
    found = edge_changes(base, other)
    assert found.rows() == [("1", "2", 1.0, "kept"), ("1", "3", 2.0, "removed"), ("2", "3", 3.0, "removed"), ("3", "4", 4.0, "added")]


@pytest.fixture
def coastal(tmp_path):
    strip = [Polygon([(177, 9.8), (180, 9.8), (180, 10.2), (177, 10.2)]), Polygon([(-180, 9.8), (-177, 9.8), (-177, 10.2), (-180, 10.2)])]
    return write_grid(tmp_path / "grid.gpkg", antimeridian_cells()), write_land(tmp_path / "land.zip", strip)


def test_edge_lines_are_the_shared_edges_as_edge_lengths_measures_them(coastal, tmp_path):
    grid_file, _ = coastal
    edges = edge_lengths(grid_file, cache_dir=tmp_path / "cache")
    pairs = sorted(edges)
    lines = edge_lines(grid_file, pairs)
    km = _line_km(np.array(lines, dtype=object))
    for pair, length in zip(pairs, km.tolist(), strict=True):
        assert length == pytest.approx(edges[pair], rel=1e-6), pair
    for line in lines:
        west, _, east, _ = line.bounds
        assert -181 < west and east < 181 and east - west < 30


LABELS = {part_id(c, "land"): 1 for c in (A, S, B)} | {part_id(c, "sea"): 2 for c in (A, S, B)} | {part_id(C, "sea"): 3}


def test_the_realm_and_boundary_figures_are_drawn(coastal, tmp_path):
    grid_file, land = coastal
    geoms = part_geometries(grid_file, land, list(LABELS))
    parts = pl.DataFrame(
        {"part": list(LABELS), "k": 7, "realm": list(LABELS.values()), "unit": [p if p != part_id(C, "sea") else None for p in LABELS],
         "imputed": [p == part_id(C, "sea") for p in LABELS], "selected": True}
    )
    units = parts.filter(pl.col("unit").is_not_null())
    candidates = units.select("unit", "k", "realm").with_columns(
        pl.when(pl.col("unit") == part_id(A, "sea")).then(77).otherwise(None).cast(pl.Int64).alias("sub_realm"), pl.lit(None, dtype=pl.Int64).alias("sub_realm_of")
    )
    drawn = plot_realm(geoms, parts, candidates, 2, HUES[0], tmp_path / "realm", "Realm 2", floored={u: r for u, r in units.select("unit", "realm").iter_rows()})
    assert drawn == {"realm": 2, "parts": 3, "sub_realms": 1, "imputed": 0, "units_above_floor": [2], "outside_view": 0}
    assert (tmp_path / "realm.png").stat().st_size > 5000
    with pytest.raises(ValueError, match="no parts"):
        plot_realm(geoms, parts, candidates, 5, HUES[0], tmp_path / "none", "Realm 5")

    edges = edge_lengths(grid_file, cache_dir=tmp_path / "cache")
    pairs = sorted(edges)
    changes = pl.DataFrame({"cell_a": [a for a, _ in pairs], "cell_b": [b for _, b in pairs], "km": [edges[p] for p in pairs],
                            "change": ["kept", "added", "removed", "added", "kept", "removed"][: len(pairs)]})
    regions = gpd.GeoDataFrame({"code": ["X"]}, geometry=[shapely.box(170, 0, 179, 20)], crs="EPSG:4326")
    plot_boundary_change(geoms, LABELS, {1: HUES[0], 2: HUES[1], 3: HUES[2]}, changes, edge_lines(grid_file, pairs), regions,
                         "+proj=laea +lat_0=10 +lon_0=180", (170.0, 0.0, 180.0, 25.0), tmp_path / "change", "Change", focus={part_id(A, "land")})
    assert (tmp_path / "change.png").stat().st_size > 5000


def test_the_view_holds_most_of_the_area_and_leaves_a_far_speck_out():
    squares = [shapely.box(x, y, x + 1, y + 1) for x in range(10) for y in range(10)] + [shapely.box(500, 500, 500.1, 500.1)]
    flat = gpd.GeoDataFrame({"part": [str(i) for i in range(len(squares))]}, geometry=squares, crs="EPSG:8857")
    assert view(flat) == (0.0, 0.0, 10.0, 10.0)
    assert view(flat, share=1.0) == (0.0, 0.0, 500.1, 500.1)
