import geopandas as gpd
import numpy as np
import polars as pl
import pytest
from shapely import MultiPolygon, box

from realms.grid import Grid, cell_of, locate
from realms.units import (
    MAJOR_GROUPS,
    Units,
    UnitSettings,
    coarse_id,
    impute,
    is_coarse,
    major_groups,
    network_units,
    part_regions,
    part_stats,
    removed_by_region,
)

# ------------------------------------------------------------------ a small world
# Fine cells (resolution 7) sit in a row on the equator at longitudes 0 to 5, each next to the one before. The coarse cells
# (resolution 5) are C1 over longitudes 0 to 2 and C2 over 3 to 5.
FINE = {"101": 0.0, "102": 1.0, "103": 2.0, "104": 3.0, "105": 4.0, "106": 5.0}
LAND, SEA = 100.0, 40.0
SMALL = UnitSettings(min_coverage=0.5, min_species=3, min_group_species=1, min_groups=1, impute_k=2)


def fine_grid(lons=None, resolution=7) -> Grid:
    lons = lons or FINE
    ids = list(lons)
    cells = pl.DataFrame(
        {
            "cell": ids,
            "lat": [0.0] * len(ids),
            "lon": [lons[c] for c in ids],
            "area_km2": [LAND + SEA] * len(ids),
            "land_km2": [LAND] * len(ids),
            "sea_km2": [SEA] * len(ids),
        }
    )
    adjacency = {c: set() for c in ids}
    for a, b in zip(ids, ids[1:], strict=False):
        adjacency[a].add(b)
        adjacency[b].add(a)
    return Grid(resolution, cells, adjacency)


def coarse_grid(tmp_path, squares, resolution=5):
    """The coarse grid and its file, from (id, centre lon, centre lat, geometry) rows."""
    path = tmp_path / f"ISEA3H-{resolution}.gpkg"
    gpd.GeoDataFrame(
        {
            "lat_c": [lat for _, _, lat, _ in squares],
            "lon_c": [lon for _, lon, _, _ in squares],
            "res": [float(resolution)] * len(squares),
            "geogrid_id": [i for i, *_ in squares],
        },
        geometry=[g for *_, g in squares],
        crs="EPSG:4326",
    ).to_file(path, driver="GPKG")
    cells = pl.DataFrame(
        {
            "cell": [i for i, *_ in squares],
            "lat": [lat for _, _, lat, _ in squares],
            "lon": [lon for _, lon, _, _ in squares],
            "area_km2": [1000.0] * len(squares),
            "land_km2": [500.0] * len(squares),
            "sea_km2": [500.0] * len(squares),
        }
    )
    return Grid(resolution, cells, {i: set() for i, *_ in squares}), path


def two_coarse(tmp_path):
    return coarse_grid(tmp_path, [("C1", 1.0, 0.0, box(-0.5, -1, 2.5, 1)), ("C2", 4.0, 0.0, box(2.5, -1, 5.5, 1))])


def presence(*parts):
    """`parts`: (part, {species: records}) as the cleaned presence has them."""
    rows = [(species, part, n) for part, counts in parts for species, n in counts.items()]
    return pl.DataFrame(rows, schema={"species_id": pl.String, "part": pl.String, "records": pl.Int64}, orient="row")


GROUPS = {"p": "plants", "i": "insects", "v": "land vertebrates", "f": "fish"}


def major(*species):
    """The major group of each species, from its first letter."""
    return pl.DataFrame({"species_id": list(species), "major_group": [GROUPS[s[0]] for s in species]}, schema={"species_id": pl.String, "major_group": pl.String})


ALL_SPECIES = [f"{g}{n}" for g in "pivf" for n in range(1, 10)]
MAJOR = major(*ALL_SPECIES)


def good(*species):
    """Five records of each species, which is enough for any part to be well sampled."""
    return {s: 5 for s in species}


def run(tmp_path, *parts, settings=SMALL, grid=None, coarse=None, gpkg=None):
    if coarse is None:
        coarse, gpkg = two_coarse(tmp_path)
    return network_units(presence(*parts), grid or fine_grid(), coarse, MAJOR, settings, coarse_gpkg=gpkg)


# ------------------------------------------------------------------ part_stats


def test_good_turing_coverage():
    rows = presence(("a:land", {"p1": 1, "p2": 1, "p3": 2, "p4": 6}), ("b:land", {"p1": 1}), ("c:sea", {"p1": 0}))
    stats = part_stats(rows, MAJOR)
    assert stats.columns == ["part", "kind", "records", "species", "singletons", "coverage", "groups_sampled"]
    a, b, c = (stats.filter(pl.col("part") == p).row(0, named=True) for p in ("a:land", "b:land", "c:sea"))
    assert (a["records"], a["species"], a["singletons"], a["kind"]) == (10, 4, 2, "land")
    assert a["coverage"] == pytest.approx(0.8)  # 1 - 2 / 10: singletons are species with exactly one record
    assert (b["records"], b["species"], b["singletons"], b["coverage"]) == (1, 1, 1, 0.0)
    # a part whose only row has no clean records is still a part, with nothing in it
    assert (c["records"], c["species"], c["singletons"], c["coverage"], c["kind"]) == (0, 0, 0, 0.0, "sea")


def test_clean_records_are_n_less_human_text_when_there_is_no_records_column():
    rows = pl.DataFrame(
        {"species_id": ["p1", "p2", "p3"], "part": ["a:land"] * 3, "n": [3, 2, 4], "n_human_text": [2, 0, 4]}
    )
    stats = part_stats(rows, MAJOR).row(0, named=True)
    # p1 has 1 clean record, p2 has 2 and p3 has none
    assert (stats["records"], stats["species"], stats["singletons"]) == (3, 2, 1)


def test_groups_sampled_counts_the_major_groups_with_enough_species():
    rows = presence(("a:land", {**good("p1", "p2", "p3"), **good("i1", "i2"), **good("v1"), "x1": 5}))  # x1 has no group
    settings = UnitSettings(min_group_species=2)
    assert part_stats(rows, MAJOR, settings).row(0, named=True)["groups_sampled"] == 2  # plants and insects, not the one vertebrate
    assert part_stats(rows, MAJOR, UnitSettings(min_group_species=1)).row(0, named=True)["groups_sampled"] == 3


# ------------------------------------------------------------------ coarsening


def test_a_sparse_part_coarsens_and_its_species_count_once(tmp_path):
    units = run(
        tmp_path,
        ("101:land", {"p1": 1, "p2": 5}),  # two species of the three that are needed, one of them a singleton
        ("102:land", {"p1": 1, "p3": 5}),  # p1 again: the same species, in a second part of C1
        ("103:land", good("p1", "p2", "p3")),  # not sparse
    )
    assert units.members["101:land"] == units.members["102:land"] == coarse_id(5, "C1", "land") == "r5:C1:land"
    assert units.members["103:land"] == "103:land"  # a part that is not sparse stays as it is
    merged = units.presence.filter(pl.col("unit") == "r5:C1:land").sort("species_id")
    assert merged.rows() == [("p1", "r5:C1:land", 2), ("p2", "r5:C1:land", 5), ("p3", "r5:C1:land", 5)]  # one row for p1, records summed
    stats = units.stats.filter(pl.col("unit") == "r5:C1:land").row(0, named=True)
    # p1 was a singleton in each part, and with two records is not one in the unit: this is what merging is for
    assert (stats["species"], stats["singletons"], stats["records"], stats["status"]) == (3, 0, 12, "network")
    assert units.excluded.is_empty()


def test_a_still_sparse_coarse_unit_is_excluded_and_imputed_by_its_kind(tmp_path):
    units = run(
        tmp_path,
        ("101:land", good("p1", "p2", "p3")),
        ("102:land", good("p1", "p2", "p3")),
        ("105:land", good("p1", "p2", "p3")),
        ("103:sea", good("p1", "p2", "p3")),
        ("104:sea", good("p1", "p2", "p3")),
        ("104:land", {"p1": 1}),  # the only sparse part of C2 and still sparse alone
    )
    assert units.excluded.select("part", "reason", "tried_unit").rows() == [("104:land", "sparse", "r5:C2:land")]
    assert "104:land" not in units.members
    # the realms of the units: the land parts either side of 104 disagree and the sea parts next to it would outvote both
    labels = {"101:land": 1, "102:land": 1, "105:land": 2, "103:sea": 7, "104:sea": 7}
    realms, imputed = impute(labels, units, fine_grid())
    # 104:land has 105 (one cell away) and 102 (two) as its nearest land parts, realms 2 and 1: a tie goes to the lowest id, and the sea is ignored
    assert realms["104:land"] == 1
    assert imputed == {"104:land"}
    assert realms["101:land"] == 1 and realms["105:land"] == 2 and realms["103:sea"] == 7  # network parts keep their unit's realm


def test_the_composition_gate_needs_two_major_groups(tmp_path):
    settings = UnitSettings(min_coverage=0.5, min_species=3, min_group_species=2, min_groups=2, impute_k=2)
    units = run(
        tmp_path,
        ("101:land", good("p1", "p2", "p3")),  # plants only
        ("102:land", {**good("p1", "p2"), **good("i1", "i2")}),  # plants and insects, two species each
        ("103:land", {**good("p1", "p2"), **good("i1")}),  # insects are one species short
        ("104:land", good("p1", "p2")),  # sparse, and with 105 makes three plants and nothing else
        ("105:land", good("p2", "p3")),
        settings=settings,
    )
    assert units.members == {"102:land": "102:land"}
    assert units.excluded.select("part", "reason", "tried_unit").sort("part").rows() == [
        ("101:land", "composition", "101:land"),
        ("103:land", "composition", "103:land"),
        ("104:land", "composition", "r5:C2:land"),  # the coarse unit passed the sparse test and failed this one
        ("105:land", "composition", "r5:C2:land"),
    ]
    assert units.stats.filter(pl.col("unit") == "102:land")["groups_sampled"].to_list() == [2]
    assert units.stats.filter(pl.col("unit") == "r5:C2:land").row(0, named=True)["status"] == "composition"


# ------------------------------------------------------------------ ids


def test_fine_and_coarse_ids_never_collide(tmp_path):
    coarse = coarse_id(5, "C1", "land")
    assert is_coarse(coarse) and not is_coarse("101:land") and not is_coarse("-101:sea")
    # read as a part id, a coarse unit has a "cell" that is not a cell of the fine grid
    parsed = pl.select(cell_of(pl.lit(coarse))).item()
    assert parsed == "r5:C1" and parsed not in fine_grid().adjacency
    # the resolution is the setting's, never a literal
    assert coarse_id(4, "C1", "sea") == "r4:C1:sea"

    units = run(tmp_path, ("101:land", {"p1": 1}), ("102:land", good("p1", "p2", "p3")), ("101:sea", {"p1": 1}))
    assert not any(is_coarse(p) for p in units.members) and not any(is_coarse(p) for p in units.excluded["part"])
    assert all(is_coarse(t) for t in units.excluded.filter(pl.col("part").is_in(["101:land", "101:sea"]))["tried_unit"])
    assert set(units.excluded["tried_unit"]) == {"r5:C1:land", "r5:C1:sea"}  # one coarse unit per kind

    # a fine cell whose name could be read as a coarse unit is refused, so the two can never meet
    clash = fine_grid({"r5:C1": 0.0, "102": 1.0})
    coarse_grid_, gpkg = two_coarse(tmp_path)
    with pytest.raises(ValueError, match="r5:C1"):
        network_units(presence(("r5:C1:land", {"p1": 1})), clash, coarse_grid_, MAJOR, SMALL, coarse_gpkg=gpkg)


def test_the_coarse_grid_must_be_the_one_the_settings_name(tmp_path):
    coarse, gpkg = coarse_grid(tmp_path, [("C1", 1.0, 0.0, box(-0.5, -1, 5.5, 1))], resolution=4)
    with pytest.raises(ValueError, match="resolution 4"):
        network_units(presence(("101:land", {"p1": 1})), fine_grid(), coarse, MAJOR, SMALL, coarse_gpkg=gpkg)
    # and with the setting changed it is accepted, and the id carries that resolution
    units = network_units(presence(("101:land", {"p1": 1})), fine_grid(), coarse, MAJOR, UnitSettings(coarse_resolution=4, min_species=3), coarse_gpkg=gpkg)
    assert units.excluded["tried_unit"].to_list() == ["r4:C1:land"]


# ------------------------------------------------------------------ adjacency and areas


def test_a_coarse_unit_is_adjacent_to_its_fine_neighbours(tmp_path):
    units = run(
        tmp_path,
        ("101:land", {"p1": 1, "p2": 5}),
        ("102:land", {"p1": 1, "p3": 5}),  # C1 merges these two
        ("102:sea", good("p1", "p2", "p3")),
        ("103:land", good("p1", "p2", "p3")),
        ("104:land", {"p1": 1, "p2": 5}),
        ("105:land", {"p1": 1, "p3": 5}),  # C2 merges these two
        ("106:land", good("p1", "p2", "p3")),
    )
    c1, c2 = "r5:C1:land", "r5:C2:land"
    assert units.adjacency[c1] == {"103:land", "102:sea"}  # the land part next to 102, and the sea part in 102's own cell
    assert units.adjacency["103:land"] == {c1, c2}
    assert units.adjacency[c2] == {"103:land", "106:land"}
    assert units.adjacency["102:sea"] == {c1}
    assert units.adjacency["106:land"] == {c2}
    for unit, near in units.adjacency.items():
        assert unit not in near  # no unit is its own neighbour
        assert all(unit in units.adjacency[other] for other in near)  # and the relation is symmetric
    assert set(units.adjacency) == set(units.areas) == set(units.members.values())


def test_unit_areas_sum_their_member_parts_by_kind(tmp_path):
    units = run(
        tmp_path,
        ("101:land", {"p1": 1, "p2": 5}),
        ("102:land", {"p1": 1, "p3": 5}),
        ("101:sea", {"p1": 1, "p2": 5}),
        ("102:sea", {"p1": 1, "p3": 5}),
        ("103:land", good("p1", "p2", "p3")),
    )
    assert units.areas == {"r5:C1:land": 2 * LAND, "r5:C1:sea": 2 * SEA, "103:land": LAND}


# ------------------------------------------------------------------ the antimeridian


def test_a_cell_across_the_antimeridian_finds_its_resolution_5_cell(tmp_path):
    # the file holds A in two pieces at 180 and -180, and B lies east of it
    a = MultiPolygon([box(170, -5, 180, 5), box(-180, -5, -170, 5)])
    coarse, gpkg = coarse_grid(tmp_path, [("A", 180.0, 0.0, a), ("B", -160.0, 0.0, box(-170, -5, -150, 5))])
    grid = fine_grid({"201": 175.0, "202": -175.0, "203": -165.0, "204": 179.9})
    units = network_units(presence(*((f"{c}:land", {"p1": 1}) for c in grid.adjacency)), grid, coarse, MAJOR, SMALL, coarse_gpkg=gpkg)
    assert dict(units.excluded.select("part", "tried_unit").iter_rows()) == {
        "201:land": "r5:A:land",
        "202:land": "r5:A:land",  # at -175, on the other side of the line from where A's polygon is drawn
        "203:land": "r5:B:land",
        "204:land": "r5:A:land",
    }


def test_locate_takes_the_lowest_index_on_a_shared_edge_and_the_nearest_outside():
    shapes = np.array([box(0, 0, 2, 2), box(2, 0, 4, 2), box(10, 0, 12, 2)], dtype=object)
    lon = np.array([1.0, 2.0, 3.0, 8.0, 12.5, -1.0])
    lat = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    # 2.0 lies on the edge of the first two. 8 is outside all and nearer to the third than the second (2 against 4)
    assert locate(lon, lat, shapes).tolist() == [0, 0, 1, 2, 2, 0]


# ------------------------------------------------------------------ every part is somewhere


def test_every_part_is_a_member_or_excluded_never_both(tmp_path):
    parts = [
        ("101:land", good("p1", "p2", "p3")),  # kept as it is
        ("102:land", {"p1": 1, "p2": 5}),  # sparse, joins C1
        ("103:land", {"p1": 1, "p3": 5}),  # sparse, joins C1, which passes
        ("104:land", {"p1": 1}),  # sparse and alone in C2
        ("101:sea", good("p1", "p2")),  # sparse sea, alone in C1's sea
        ("102:sea", {"p1": 0}),  # nothing clean at all
        ("106:sea", good("p1", "p2", "p3")),
    ]
    units = run(tmp_path, *parts)
    present = {p for p, _ in parts}
    members, excluded = set(units.members), set(units.excluded["part"])
    assert members | excluded == present
    assert not members & excluded
    assert units.excluded["part"].is_unique().all()
    # an excluded part is only in `excluded`
    in_units = set(units.members.values())
    assert not excluded & (set(units.presence["unit"]) | set(units.areas) | set(units.adjacency) | in_units)
    assert set(units.presence["unit"]) <= in_units
    assert set(units.excluded["reason"]) <= {"sparse", "composition"}
    assert units.excluded.columns[:3] == ["part", "reason", "tried_unit"]
    # the units' records are the records of the member parts, none lost and none counted twice
    kept = sum(n for p, counts in parts for n in counts.values() if p in members)
    assert units.presence["records"].sum() == kept


def test_a_part_without_a_cell_in_the_grid_is_refused(tmp_path):
    with pytest.raises(ValueError, match="999"):
        run(tmp_path, ("999:land", {"p1": 1}))


# ------------------------------------------------------------------ imputation


def made(members, excluded, grid_lons, k=2):
    """Units with the given member parts (all their own units) and excluded parts, for `impute` alone."""
    ex = pl.DataFrame(
        {"part": excluded, "reason": ["sparse"] * len(excluded), "tried_unit": excluded, "records": [1] * len(excluded), "species": [1] * len(excluded)}
    )
    units = Units(
        members={p: p for p in members},
        presence=pl.DataFrame(schema={"species_id": pl.String, "unit": pl.String, "records": pl.Int64}),
        areas={p: 1.0 for p in members},
        adjacency={p: set() for p in members},
        excluded=ex,
        stats=pl.DataFrame(),
        settings=UnitSettings(impute_k=k),
    )
    return units, fine_grid(grid_lons)


def test_imputation_takes_the_majority_of_the_nearest_k_and_ties_go_to_the_lowest_realm():
    lons = {"a": 0.0, "b": 1.0, "c": 2.0, "d": 3.0, "e": 20.0, "f": 21.0}
    units, grid = made(["b:land", "c:land", "d:land", "e:land", "f:land"], ["a:land"], lons, k=3)
    # a's nearest three are b, c, d (realms 5, 5, 3): the majority is 5, though 3 is the lowest id
    realms, imputed = impute({"b:land": 5, "c:land": 5, "d:land": 3, "e:land": 3, "f:land": 3}, units, grid)
    assert realms["a:land"] == 5 and imputed == {"a:land"}
    units, grid = made(["b:land", "c:land", "d:land", "e:land", "f:land"], ["a:land"], lons, k=2)
    realms, _ = impute({"b:land": 5, "c:land": 3, "d:land": 3, "e:land": 3, "f:land": 3}, units, grid)
    assert realms["a:land"] == 3  # b and c are the nearest two: one vote each, and 3 is lower
    # k reaching past the network takes everything there is
    units, grid = made(["b:land", "c:land"], ["a:land"], lons, k=6)
    assert impute({"b:land": 8, "c:land": 8}, units, grid)[0]["a:land"] == 8


def test_imputation_breaks_a_tie_in_distance_by_part_id():
    lons = {"x": 0.0, "m": -1.0, "n": 1.0, "far": 10.0}
    units, grid = made(["m:land", "n:land", "far:land"], ["x:land"], lons, k=1)
    # m and n are the same distance from x: the lower id, m, is the one that counts, and that is not the lower realm
    assert impute({"m:land": 9, "n:land": 2, "far:land": 2}, units, grid)[0]["x:land"] == 9


def test_imputation_uses_the_great_circle_not_the_degrees():
    # 80 degrees north, where one degree of longitude is about a sixth of one of latitude: p is 20 degrees of longitude from x
    # (3.5 degrees of arc) and q is 4 degrees of latitude from it (4 degrees of arc)
    grid = fine_grid({"x": 0.0, "p": 20.0, "q": 0.0})
    grid = Grid(7, grid.cells.with_columns(pl.Series("lat", [80.0, 80.0, 84.0])), grid.adjacency)
    units, _ = made(["p:land", "q:land"], ["x:land"], {"x": 0.0, "p": 20.0, "q": 0.0}, k=1)
    assert impute({"p:land": 1, "q:land": 2}, units, grid)[0]["x:land"] == 1


def test_imputation_by_kind_needs_a_network_part_of_that_kind():
    lons = {"a": 0.0, "b": 1.0}
    units, grid = made(["b:land"], ["a:sea"], lons, k=1)
    with pytest.raises(ValueError, match="sea"):
        impute({"b:land": 1}, units, grid)


def test_imputation_wants_a_label_for_every_unit_and_only_those():
    lons = {"a": 0.0, "b": 1.0, "c": 2.0}
    units, grid = made(["b:land", "c:land"], ["a:land"], lons, k=1)
    with pytest.raises(ValueError, match="c:land"):
        impute({"b:land": 1}, units, grid)
    with pytest.raises(ValueError, match="z:land"):
        impute({"b:land": 1, "c:land": 1, "z:land": 4}, units, grid)


def test_nothing_to_impute_still_labels_the_network_parts(tmp_path):
    units = run(tmp_path, ("101:land", good("p1", "p2", "p3")))
    realms, imputed = impute({"101:land": 4}, units, fine_grid())
    assert realms == {"101:land": 4} and imputed == set()


# ------------------------------------------------------------------ major groups


def taxa(rows):
    return pl.DataFrame(rows, schema={"species_id": pl.String, "phylum": pl.String, "class": pl.String}, orient="row")


def test_major_groups_follow_the_tax_group_and_the_class_of_chordates():
    t = taxa(
        [
            ("rose", "Tracheophyta", "Magnoliopsida"),
            ("beetle", "Arthropoda", "Insecta"),
            ("spider", "Arthropoda", "Arachnida"),
            ("snail", "Mollusca", "Gastropoda"),
            ("mushroom", "Basidiomycota", "Agaricomycetes"),
            ("worm", "Annelida", "Clitellata"),
            ("cod", "Chordata", "Teleostei"),
            ("owl", "Chordata", "Aves"),
            ("frog", "Chordata", "Amphibia"),
            ("seasquirt", "Chordata", "Ascidiacea"),
            ("lancelet", "Chordata", "Leptocardii"),
            ("lamprey", "Chordata", "Petromyzonti"),
            ("odd fish", "Chordata", "Teleostei"),  # grouped by COL as plain animals, which it is not
            ("classless", "Chordata", None),
            ("kelp", "Ochrophyta", "Phaeophyceae"),
            ("unplaced", None, None),
        ]
    )
    g = pl.DataFrame(
        {
            "species_id": ["rose", "beetle", "spider", "snail", "mushroom", "worm", "cod", "owl", "frog", "seasquirt", "lancelet", "lamprey", "odd fish", "classless", "kelp"],
            "tax_group": ["angiosperms", "coleoptera", "arachnids", "gastropods", "basidiomycetes", "otheranimals"] + ["chordates"] * 6 + ["animals", "chordates", "algae"],
        }
    )
    found = dict(major_groups(t, g).iter_rows())
    assert found == {
        "rose": "plants",
        "beetle": "insects",
        "spider": "other arthropods",
        "snail": "molluscs",
        "mushroom": "fungi",
        "worm": "other invertebrates",
        "cod": "fish",
        "owl": "land vertebrates",
        "frog": "land vertebrates",
        "seasquirt": "other invertebrates",
        "lancelet": "other invertebrates",
        "lamprey": "fish",
        "odd fish": "fish",
    }  # a chordate of no class, algae and a species with no group are not counted


def test_every_tax_group_in_the_brief_maps_and_an_unknown_one_is_refused():
    assert MAJOR_GROUPS["angiosperms"] == "plants" and MAJOR_GROUPS["otherinsects"] == "insects" and MAJOR_GROUPS["crustacean"] == "other arthropods"
    assert MAJOR_GROUPS["othermolluscs"] == "molluscs" and MAJOR_GROUPS["otherfungi"] == "fungi" and MAJOR_GROUPS["animals"] == "other invertebrates"
    assert not {"algae", "protists", "bacteria", "archaea", "viruses", "eukaryotes", "pseudofungi", "prokaryotes"} & set(MAJOR_GROUPS)
    with pytest.raises(ValueError, match="newgroup"):
        major_groups(taxa([("x", "P", "C")]), pl.DataFrame({"species_id": ["x"], "tax_group": ["newgroup"]}))


# ------------------------------------------------------------------ removed by region


def write_level3(path, regions):
    gpd.GeoDataFrame(
        {"LEVEL3_NAM": [c for c, *_ in regions], "LEVEL3_COD": [c for c, *_ in regions], "LEVEL2_COD": [10] * len(regions), "LEVEL1_COD": [l1 for _, l1, _ in regions]},
        geometry=[g for *_, g in regions],
        crs="EPSG:4326",
    ).to_file(path, driver="GeoJSON")
    return path


def test_part_regions_are_the_level_1_region_of_the_cell_centre(tmp_path):
    level3 = write_level3(
        tmp_path / "level3.geojson",
        [("AAA", 1, box(-0.5, -1, 1.5, 1)), ("BBB", 2, box(1.5, -1, 3.5, 1)), ("FJI", 6, MultiPolygon([box(176, -5, 180, 5), box(-180, -5, -176, 5)]))],
    )
    grid = fine_grid({"101": 0.0, "102": 1.0, "103": 2.0, "104": 10.0, "105": -178.0, "106": 178.0})
    regions = part_regions(grid, level3)
    assert regions.columns == ["part", "region"]
    found = dict(regions.iter_rows())
    assert found["101:land"] == found["102:land"] == "Europe"
    assert found["103:land"] == "Africa"
    assert found["104:land"] == "Africa"  # 10 degrees east of everything: a centre in no region goes to the nearest
    assert found["105:land"] == found["106:land"] == "Pacific"  # both sides of the antimeridian
    assert not any(p.endswith(":sea") for p in found)


def test_removed_by_region_counts_excluded_parts_and_records_by_reason(tmp_path):
    excluded = pl.DataFrame(
        {
            "part": ["a:land", "b:land", "c:land", "d:land", "e:sea", "f:sea", "g:land"],
            "reason": ["sparse", "sparse", "composition", "sparse", "sparse", "composition", "sparse"],
            "tried_unit": ["r5:X:land", "r5:X:land", "c:land", "r5:Y:land", "r5:Z:sea", "f:sea", "r5:W:land"],
            "records": [10, 5, 7, 3, 100, 40, 1],
            "species": [1, 1, 1, 1, 1, 1, 1],
        }
    )
    units = Units({}, pl.DataFrame(), {}, {}, excluded, pl.DataFrame())
    regions = pl.DataFrame({"part": ["a:land", "b:land", "c:land", "d:land"], "region": ["Europe", "Europe", "Europe", "Africa"]})
    out = removed_by_region(units, regions)
    assert out.columns == ["reason", "region", "parts", "records"]
    assert out.rows() == [
        ("composition", "Europe", 1, 7),
        ("composition", "sea", 1, 40),
        ("sparse", "Africa", 1, 3),
        ("sparse", "Europe", 2, 15),
        ("sparse", "sea", 1, 100),
        ("sparse", "unplaced", 1, 1),  # a land part no region was found for is said so
    ]
