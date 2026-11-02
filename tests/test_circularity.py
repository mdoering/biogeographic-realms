import geopandas as gpd
import polars as pl
import pytest
import shapely

from realms.boundaries import EDGE_SCHEMA, boundary
from realms.circularity import NONE, centre_regions, common_k, compare, headline, labels_at, moved, moves_by_region

# five cells in a row, a - b - c - d - e, each with a land and a sea part; every shared edge is 10 km, except d-e at 5 km
EDGES = {("a", "b"): 10.0, ("b", "c"): 10.0, ("c", "d"): 10.0, ("d", "e"): 5.0}


def parts(rows):
    return pl.DataFrame(rows, schema={"part": pl.String, "k": pl.Int64, "realm": pl.Int64, "selected": pl.Boolean}, orient="row")


def test_labels_at_takes_one_k_of_parts():
    table = parts([("a:land", 2, 0, False), ("b:land", 2, 1, False), ("a:land", 3, 5, True), ("b:land", 3, 6, True)])
    assert labels_at(table, 3) == {"a:land": 5, "b:land": 6}
    with pytest.raises(ValueError, match="k=4"):
        labels_at(table, 4)


def test_compare_counts_shared_and_lone_parts_and_scores_the_shared_ones():
    a = {"a:land": 0, "b:land": 0, "c:land": 1, "d:land": 1, "x:sea": 3}
    b = {"a:land": 7, "b:land": 7, "c:land": 9, "d:land": 9, "y:sea": 1, "z:sea": 1}
    found = compare(a, b)
    assert found == {"v_measure": pytest.approx(1.0), "shared": 4, "only_a": 1, "only_b": 2}
    split = compare(a, {"a:land": 0, "b:land": 1, "c:land": 0, "d:land": 1})
    assert split["v_measure"] == pytest.approx(0.0) and split["shared"] == 4


def test_common_k_is_the_shared_k_most_runs_selected():
    candidates = {"off": [5, 6, 7, 8], "introduced": [5, 6, 7, 8, 9], "native_only": range(5, 11), "no_layer2": range(5, 13), "no_griis": range(5, 10)}
    selected = {"off": 7, "introduced": 7, "native_only": 7, "no_layer2": 12, "no_griis": 9}
    assert common_k(candidates, selected) == 7
    # a tie goes to the coarser k, and a selected k that not every run has does not count
    assert common_k({"a": [5, 6], "b": [5, 6, 9]}, {"a": 6, "b": 9}) == 6
    assert common_k({"a": [5, 6], "b": [5, 6]}, {"a": 6, "b": 5}) == 5
    assert common_k({"a": [5, 6], "b": [5, 6, 9]}, {"a": 9, "b": 9}) == 5
    with pytest.raises(ValueError, match="no k"):
        common_k({"a": [5], "b": [6]}, {"a": 5, "b": 6})


def test_moved_is_the_symmetric_difference_of_the_boundary_edges():
    one = boundary({"a:land": 0, "b:land": 0, "c:land": 1, "d:land": 1, "e:land": 1, "a:sea": 0, "b:sea": 1}, "land", EDGES)
    # the b-c boundary stays (realms renamed), c-d is new, and d-e is new because e changes realm
    two = boundary({"a:land": 4, "b:land": 4, "c:land": 5, "d:land": 6, "e:land": 4}, "land", EDGES)
    found = moved(one, two)
    assert found.schema == pl.Schema(EDGE_SCHEMA)
    assert found.rows() == [("c", "d", 10.0), ("d", "e", 5.0)]
    assert moved(two, one).rows() == found.rows()
    assert moved(one, one).height == 0
    # the land and sea parts of a cell are never one boundary, so a map's sea edges are not among its land edges
    assert moved(boundary({"a:sea": 0, "b:sea": 1}, "sea", EDGES), boundary({"a:sea": 0, "b:sea": 0}, "sea", EDGES)).rows() == [("a", "b", 10.0)]


def test_moves_by_region_takes_the_first_cell_then_the_second_then_none():
    edges = pl.DataFrame([("a", "b", 10.0), ("b", "c", 4.0), ("c", "d", 3.0), ("d", "e", 1.0)], schema=EDGE_SCHEMA, orient="row")
    region_of = {"a": "Europe", "b": None, "c": "Africa", "e": None}  # d has no entry at all
    found = moves_by_region(edges, region_of)
    assert found.columns == ["region", "km", "edges"]
    assert found.rows() == [("Europe", 10.0, 1), ("Africa", 7.0, 2), (NONE, 1.0, 1)]
    assert moves_by_region(edges.head(0), region_of).height == 0


def test_centre_regions_holds_only_the_centres_inside_a_region():
    regions = gpd.GeoDataFrame(
        {"LEVEL3_COD": ["BBB", "AAA"], "LEVEL1_COD": [2, 1]},
        geometry=[shapely.box(10, 0, 20, 10), shapely.box(0, 0, 10, 10)],
        crs="EPSG:4326",
    )
    centres = pl.DataFrame({"cell": ["in1", "in2", "edge", "out", "far"], "lon": [5.0, 15.0, 10.0, 25.0, 15.0 - 360.0], "lat": [5.0, 5.0, 5.0, 5.0, 5.0]})
    found = centre_regions(centres, regions)
    # a centre on the edge between two goes to the lower level 3 code (AAA, Europe); a turn of longitude is the same place
    assert found == {"in1": "Europe", "in2": "Africa", "edge": "Europe", "out": None, "far": "Africa"}


def test_headline_is_the_strongest_wcvp_mode_that_passes():
    assert headline({"introduced": True, "native_only": True}) == "native_only"
    assert headline({"introduced": True, "native_only": False}) == "introduced"
    assert headline({"introduced": False, "native_only": False}) == "off"
    # a layer variant that passes is information, never the headline
    assert headline({"introduced": False, "native_only": False, "native_only-no-griis-no-griis_subnational": True}) == "off"
    with pytest.raises(ValueError, match="introduced"):
        headline({"native_only": False})
