"""The cleaning of occurrence counts, layers 1 to 4 and the way layer 5 is called, on hand-made tables. Layer 5's rules are in `test_outliers`.

The grid has cells named by numbers, so that "lowest cell id" can be told from the lowest id as text ("9" before "100"). A cell is
`(land_km2, sea_km2)`, and a part with less than `ZERO_KM2` of area is none.
"""

import numpy as np
import polars as pl
import pytest

from clean_helpers import GRID, counts, gone, kept, make_grid, refs, row, taxon
from realms.clean import LAYERS, RULES, SWITCHABLE, Settings, _Run, clean
from realms.clean_report import removal_summary
from realms.grid import ZERO_KM2, part_adjacency, parts
from realms.matrix import cell_taxa
from realms.nullmodel import excess_table
from realms.occurrences import COUNTS

ALL = frozenset(LAYERS)


# ---- settings

def test_settings_hold_the_design_values():
    s = Settings()
    assert (s.land_share, s.satellite_gap_km, s.satellite_record_share, s.satellite_cell_share) == (0.95, 500, 0.01, 0.05)
    assert (s.satellite_min_datasets, s.satellite_min_years, s.group_min_records) == (2, 3, 1000)
    assert (s.under_max_records, s.under_rate_factor, s.thin_max_records, s.thin_km) == (5, 0.01, 5, 1000)
    assert (s.core_share, s.heuristic_min_records, s.wcvp_mode) == (0.9, 20, "native_only")
    assert (s.side_conflict, s.rules_off) == ("environment", frozenset())
    with pytest.raises(AttributeError):
        s.land_share = 0.5


def test_settings_refuse_what_they_do_not_know():
    with pytest.raises(ValueError, match="wcvp_mode"):
        Settings(wcvp_mode="strict")
    with pytest.raises(ValueError, match="side_conflict"):
        Settings(side_conflict="records")
    with pytest.raises(ValueError, match="land_share"):
        Settings(land_share=0.5)  # a species could then be on both sides
    with pytest.raises(ValueError, match="zero_area_part"):
        Settings(rules_off=frozenset({"zero_area_part"}))
    assert Settings(rules_off=frozenset(SWITCHABLE)).rules_off == frozenset(SWITCHABLE)


def test_clean_refuses_what_it_does_not_know():
    with pytest.raises(ValueError, match="layers"):
        clean(counts(), refs(), GRID, layers=frozenset({6}))
    with pytest.raises(ValueError, match="columns"):
        clean(counts().drop("n_raw"), refs(), GRID)
    with pytest.raises(ValueError, match="nowhere"):
        clean(counts(row("s", "nowhere:land")), refs(), GRID)


# ---- zero-area parts

ZERO = make_grid({"10": (1000.0, 0.0), "11": (0.0, 1000.0), "20": (0.0, 0.0), "30": (ZERO_KM2 / 2, 1000.0)}, [("10", "11")])


@pytest.mark.parametrize("layers", [ALL, frozenset()])
def test_a_zero_area_part_folds_into_the_other_part_of_its_cell(layers):
    c = counts(
        row("s1", "10:sea", n=4, n_datasets=2, n_years=5),
        row("s1", "10:land", n=6, n_datasets=3, n_years=1),
        row("s2", "10:sea", n=3, country="BB"),
        row("s3", "30:land", n=7),
    )
    out = clean(c, refs(), ZERO, layers)
    assert kept(out, "n") == {("s1", "10:land"): 10, ("s2", "10:land"): 3, ("s3", "30:sea"): 7}
    merged = out.presence.filter(pl.col("species_id") == "s1").row(0, named=True)
    assert (merged["n_datasets"], merged["n_years"]) == (3, 5)  # distinct counts take the larger, as in load_presence
    zero = out.removed.filter(pl.col("rule") == "zero_area_part").sort("species_id")
    assert zero["layer"].to_list() == [0, 0, 0]
    assert zero["part"].to_list() == ["10:sea", "10:sea", "30:land"]
    assert zero["to_part"].to_list() == ["10:land", "10:land", "30:sea"]
    assert zero["records"].to_list() == [4, 3, 7]
    assert zero["countrycode"].to_list() == ["AA", "BB", "AA"]
    assert out.removed.height == 3  # nothing else is removed


def test_a_cell_with_only_zero_area_parts_keeps_its_records():
    out = clean(counts(row("s1", "20:land", n=5), row("s1", "20:sea", n=2)), refs(), ZERO)
    assert kept(out, "n") == {("s1", "20:land"): 5, ("s1", "20:sea"): 2}
    assert out.removed.height == 0


def test_a_zero_area_part_folds_into_a_part_that_holds_no_records_yet():
    out = clean(counts(row("s1", "11:land", n=5)), refs(), ZERO)  # 11 has no land
    assert kept(out, "n") == {("s1", "11:sea"): 5}


def test_the_zero_area_fold_is_not_a_layer_and_counts_only_records_that_count():
    c = counts(row("s1", "10:sea", n=0, n_raw=4), row("s2", "10:sea", n=4, n_raw=4))
    with_layer = clean(c, refs(), ZERO, ALL)
    assert gone(with_layer, rule="zero_area_part") == {("s2", "10:sea", "AA")}  # s1's records are no records at layer 1
    off = clean(c, refs(), ZERO, frozenset())  # layer 1 off: n_raw counts
    assert gone(off, rule="zero_area_part") == {("s1", "10:sea", "AA"), ("s2", "10:sea", "AA")}
    assert kept(off) == {("s1", "10:land"): 4, ("s2", "10:land"): 4}


def test_excess_table_runs_on_maps_built_from_cleaned_presence():
    grid = make_grid(
        {"1": (1000.0, 0.0), "2": (1000.0, 0.0), "3": (500.0, 500.0), "4": (1000.0, 0.0), "5": (1000.0, 0.0), "20": (0.0, 1000.0)},
        [("1", "2"), ("2", "3"), ("3", "4"), ("4", "5"), ("5", "20")],
    )
    rows = [row(s, p, n=3) for s in ("a", "b", "c") for p in ("1:land", "2:land", "3:land", "3:sea", "4:land", "5:land")]
    rows += [row(s, "20:land", n=3) for s in ("a", "b")]  # a sliver of land in a cell of sea
    c = counts(*rows)

    # after the fold no part of zero area is left, and the map's realms all have area
    cleaned = clean(c, refs(), grid)
    assert "20:land" not in set(cleaned.presence["part"]) and "20:sea" in set(cleaned.presence["part"])
    names = sorted(cleaned.presence["part"].unique())
    assert all(parts(grid, names)["area_km2"] >= ZERO_KM2)
    table = parts(grid, names)
    cell = dict(zip(table["part"], table["cell"], strict=True))
    labels = {p: 1 if cell[p] in {"4", "5", "20"} else 0 for p in names}
    ct = cell_taxa(cleaned.presence.select("species_id", pl.col("part").alias("cell")), None)
    result = excess_table({2: labels}, ct, part_adjacency(grid, names), dict(zip(table["part"], table["area_km2"], strict=True)), n_null=3)
    assert set(result) == {2} and isinstance(result[2], float)


# ---- layer 1

def test_layer_1_counts_the_clean_records_less_those_matching_human_text():
    c = counts(
        row("text", "10:land", n=5, n_raw=6, n_human_text=5),
        row("flags", "10:land", n=0, n_raw=3),
        row("some", "10:land", n=4, n_raw=7, n_human_text=1),
        row("none", "10:land", n=0, n_raw=0),
    )
    out = clean(c, refs(), GRID, frozenset({1}))
    assert kept(out) == {("some", "10:land"): 3}
    assert out.presence.row(0, named=True)["n_raw"] == 7  # the other counts are kept as they came
    assert gone(out, 1, "human_text") == {("text", "10:land", "AA")}
    assert gone(out, 1, "failed_flags") == {("flags", "10:land", "AA")}
    assert set(out.removed["layer"]) == {1}
    assert dict(zip(out.removed["species_id"], out.removed["records"], strict=True)) == {"text": 6, "flags": 3}  # all that the row held


def test_layer_1_off_counts_every_record():
    c = counts(row("text", "10:land", n=5, n_raw=6, n_human_text=5), row("flags", "10:land", n=0, n_raw=3), row("none", "10:land", n=0, n_raw=0))
    out = clean(c, refs(), GRID, frozenset({2, 3, 4}))
    assert kept(out) == {("text", "10:land"): 6, ("flags", "10:land"): 3}
    assert out.removed.height == 0


# ---- layer 2

def test_layer_2_removes_extinct_and_no_wild_range_species_each_by_its_own_rule():
    r = refs([taxon("ext", extinct=True), taxon("crop"), taxon("wild"), taxon("both", extinct=True)], no_wild_range={"crop", "both"})
    c = counts(*(row(s, "10:land", country=cc) for s in ("ext", "crop", "wild", "both") for cc in ("AA", "BB")))
    out = clean(c, r, GRID)
    assert set(kept(out)) == {("wild", "10:land")}
    assert gone(out, 2, "extinct") == {("ext", "10:land", "AA"), ("ext", "10:land", "BB"), ("both", "10:land", "AA"), ("both", "10:land", "BB")}
    assert gone(out, 2, "no_wild_range") == {("crop", "10:land", "AA"), ("crop", "10:land", "BB")}


def test_layer_2_rules_can_be_switched_off_one_at_a_time():
    r = refs([taxon("ext", extinct=True), taxon("crop")], no_wild_range={"crop"})
    c = counts(row("ext", "10:land"), row("crop", "10:land"))
    assert set(kept(clean(c, r, GRID, settings=Settings(rules_off=frozenset({"extinct"}))))) == {("ext", "10:land")}
    assert set(kept(clean(c, r, GRID, settings=Settings(rules_off=frozenset({"no_wild_range"}))))) == {("crop", "10:land")}
    assert set(kept(clean(c, r, GRID, frozenset({1, 3, 4})))) == {("ext", "10:land"), ("crop", "10:land")}


# ---- layer 3: WCVP

WCVP_GRID = make_grid(
    {"10": (1000.0, 0.0), "11": (500.0, 500.0), "13": (1000.0, 0.0), "14": (1000.0, 0.0), "15": (1000.0, 0.0)},
    [("10", "11"), ("11", "13"), ("13", "14")],
)
# 10 is in the native region N, 11 in the introduced region I, 13 in neither mentioned nor unmentioned alone: U, which WCVP says nothing
# about. 14 is in N and I at once (a border cell), and 15 is in no region at all, as the 80 specks of coast are.
REGIONS = [("10:land", "N"), ("11:land", "I"), ("13:land", "U"), ("14:land", "N"), ("14:land", "I"), ("15:land", None)]
PLANT_ROWS = [row("p", part) for part in ("10:land", "11:land", "11:sea", "13:land", "14:land", "15:land")]


def plant_refs(wcvp=(("p", "N", True), ("p", "I", False)), kingdom="Plantae"):
    return refs([taxon("p", kingdom=kingdom, cls="Magnoliopsida")], wcvp=wcvp, tdwg3=REGIONS)


def test_wcvp_modes():
    c = counts(*PLANT_ROWS)

    native_only = clean(c, plant_refs(), WCVP_GRID, settings=Settings(wcvp_mode="native_only"))
    # the introduced part, and the part WCVP does not mention. A border cell touching a native region counts as native.
    assert gone(native_only, 3, "wcvp_outside_native") == {("p", "11:land", "AA"), ("p", "13:land", "AA")}
    assert set(kept(native_only)) == {("p", "10:land"), ("p", "11:sea"), ("p", "14:land"), ("p", "15:land")}

    introduced = clean(c, plant_refs(), WCVP_GRID, settings=Settings(wcvp_mode="introduced"))
    assert gone(introduced, 3, "wcvp_introduced") == {("p", "11:land", "AA")}
    assert set(kept(introduced)) == {("p", "10:land"), ("p", "11:sea"), ("p", "13:land"), ("p", "14:land"), ("p", "15:land")}

    off = clean(c, plant_refs(), WCVP_GRID, settings=Settings(wcvp_mode="off"))
    assert off.removed.height == 0
    assert len(kept(off)) == len(PLANT_ROWS)


def test_the_default_wcvp_mode_is_native_only():
    assert set(gone(clean(counts(*PLANT_ROWS), plant_refs(), WCVP_GRID), 3)) == {("p", "11:land", "AA"), ("p", "13:land", "AA")}


def test_wcvp_never_touches_a_sea_part_or_a_land_part_without_a_region():
    out = clean(counts(*PLANT_ROWS), plant_refs(), WCVP_GRID)
    assert ("p", "11:sea") in kept(out) and ("p", "15:land") in kept(out)


def test_wcvp_applies_to_plants_only():
    c = counts(*PLANT_ROWS)
    for kingdom in ("Animalia", "Fungi", None):
        out = clean(c, plant_refs(kingdom=kingdom), WCVP_GRID)
        assert out.removed.height == 0, kingdom


def test_wcvp_leaves_a_plant_with_no_native_region_alone_unless_the_mode_is_introduced():
    only_introduced = (("p", "I", False),)
    c = counts(*PLANT_ROWS)
    assert clean(c, plant_refs(only_introduced), WCVP_GRID, settings=Settings(wcvp_mode="native_only")).removed.height == 0
    out = clean(c, plant_refs(only_introduced), WCVP_GRID, settings=Settings(wcvp_mode="introduced"))
    assert gone(out, 3, "wcvp_introduced") == {("p", "11:land", "AA"), ("p", "14:land", "AA")}  # 14 touches I and no native region now


def test_wcvp_leaves_a_plant_it_does_not_list_alone():
    for mode in ("native_only", "introduced"):
        assert clean(counts(*PLANT_ROWS), plant_refs(wcvp=()), WCVP_GRID, settings=Settings(wcvp_mode=mode)).removed.height == 0


def test_wcvp_removes_every_country_of_a_part():
    c = counts(row("p", "13:land", country="AA"), row("p", "13:land", country="BB"), row("p", "13:land", country=None), row("p", "10:land"))
    out = clean(c, plant_refs(), WCVP_GRID)
    assert gone(out, 3, "wcvp_outside_native") == {("p", "13:land", "AA"), ("p", "13:land", "BB"), ("p", "13:land", None)}


# ---- layer 3: GRIIS and WRiMS

def test_griis_removes_a_species_country_pair_from_every_part_of_that_country():
    r = refs([taxon("s1"), taxon("s2")], griis=[("s1", "AA", None)])
    c = counts(
        row("s1", "10:land"), row("s1", "11:sea"), row("s1", "10:land", country="BB"), row("s1", "10:land", country=None), row("s2", "10:land")
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis") == {("s1", "10:land", "AA"), ("s1", "11:sea", "AA")}
    assert set(kept(out)) == {("s1", "10:land"), ("s2", "10:land")}
    assert kept(out)[("s1", "10:land")] == 20  # BB and no country, which the list does not name


FR_CELLS = [("11", "FR", "Reunion", "centre")]


def test_a_sub_national_list_matches_on_cells_not_on_the_country_code():
    # an overseas place carries the sovereign's code in the list (FR), while GBIF's records say RE
    r = refs([taxon("s1")], griis=[("s1", "FR", "Reunion")], admin1=FR_CELLS)
    c = counts(row("s1", "11:land", country="RE"), row("s1", "11:sea", country="FR"), row("s1", "10:land", country="FR"), row("s1", "10:land", country="RE"))
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis_subnational") == {("s1", "11:land", "RE"), ("s1", "11:sea", "FR")}
    assert gone(out, 3, "griis") == set()
    assert set(kept(out)) == {("s1", "10:land")}  # outside the cells the list does not apply, whatever the country code says


def test_the_country_list_applies_to_the_rest_of_the_country_and_the_sub_national_list_to_its_cells():
    r = refs(
        [taxon("s1"), taxon("s2")],
        griis=[("s1", "FR", None), ("s2", "FR", "Reunion")],
        admin1=FR_CELLS,
    )
    c = counts(
        row("s1", "10:land", country="FR"),  # the rest of France: the country's list
        row("s1", "11:land", country="FR"),  # in Reunion's cells: Reunion's list rules, which does not name s1
        row("s2", "11:land", country="FR"),
        row("s2", "10:land", country="FR"),  # the country's list does not name s2
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis") == {("s1", "10:land", "FR")}
    assert gone(out, 3, "griis_subnational") == {("s2", "11:land", "FR")}
    assert set(kept(out)) == {("s1", "11:land"), ("s2", "10:land")}


def test_what_task_10_could_not_place_is_not_applied_and_is_counted():
    r = refs(
        [taxon("s1"), taxon("s2"), taxon("absent")],
        griis=[("s1", "FR", "Clipperton"), ("s2", "FR", "Clipperton")],
        wrims=[("s1", "a", "Open sea", None), ("absent", "a", "Open sea", None)],
        griis_no_cells=[("s1", "FR", "Clipperton"), ("s2", "FR", "Clipperton"), ("absent", "FR", "Clipperton")],
        wrims_no_country=[("s1", "a", "Open sea", None), ("absent", "a", "Open sea", None)],
    )
    c = counts(row("s1", "10:land", country="FR", n=5), row("s1", "11:sea", country="FR", n=2), row("s2", "10:land", country="FR", n=4))
    out = clean(c, r, GRID)
    assert out.removed.height == 0 and kept(out) == {("s1", "10:land"): 5, ("s1", "11:sea"): 2, ("s2", "10:land"): 4}
    # only pairs of species that have records are counted, with the records they have
    unplaced = out.unplaced.sort("rule", "species_id")
    assert unplaced.columns == ["rule", "species_id", "country", "place", "records"]
    assert unplaced.select("rule", "species_id", "country", "place", "records").rows() == [
        ("griis_subnational", "s1", "FR", "Clipperton", 7),
        ("griis_subnational", "s2", "FR", "Clipperton", 4),
        ("wrims", "s1", None, "Open sea", 7),
    ]
    summary = removal_summary(out, r).filter(pl.col("unplaced"))
    assert set(summary["rule"]) == {"unplaced_griis", "unplaced_wrims"}
    assert summary.filter(pl.col("rule") == "unplaced_griis")["rows"].sum() == 2  # pairs
    assert summary.filter(pl.col("rule") == "unplaced_wrims")["rows"].sum() == 1
    assert set(summary["layer"]) == {3}


def test_nothing_is_unplaced_for_a_rule_that_is_off():
    r = refs([taxon("s1")], griis_no_cells=[("s1", "FR", "Clipperton")], wrims_no_country=[("s1", "a", "Open sea", None)])
    c = counts(row("s1", "10:land"))
    assert clean(c, r, GRID).unplaced.height == 2
    assert clean(c, r, GRID, frozenset({1, 2, 4})).unplaced.height == 0
    assert clean(c, r, GRID, settings=Settings(rules_off=frozenset({"wrims"}))).unplaced["rule"].to_list() == ["griis_subnational"]


def test_wrims_removes_a_species_from_the_sea_parts_of_the_countries_it_is_introduced_in():
    r = refs([taxon("s1"), taxon("s2")], wrims=[("s1", "a", "Italian part of the Tyrrhenian Sea", "IT")])
    c = counts(
        row("s1", "11:sea", country="IT"),
        row("s1", "11:land", country="IT"),  # on land: WRiMS lists the sea
        row("s1", "11:sea", country="FR"),
        row("s1", "12:sea", country=None),
        row("s2", "11:sea", country="IT"),
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "wrims") == {("s1", "11:sea", "IT")}
    assert set(kept(out)) == {("s1", "11:land"), ("s1", "11:sea"), ("s1", "12:sea"), ("s2", "11:sea")}


def test_layer_3_rules_can_be_switched_off_one_at_a_time():
    r = refs(
        [taxon("g"), taxon("n"), taxon("w")],
        griis=[("g", "AA", None), ("n", "FR", "Reunion")],
        admin1=FR_CELLS,
        wrims=[("w", "a", "x", "AA")],
    )
    c = counts(row("g", "10:land"), row("n", "11:land", country="RE"), row("w", "11:sea"))
    assert len(clean(c, r, GRID).presence) == 0
    for rule, species in (("griis", "g"), ("griis_subnational", "n"), ("wrims", "w")):
        out = clean(c, r, GRID, settings=Settings(rules_off=frozenset({rule})))
        assert set(out.presence["species_id"]) == {species}, rule
    assert len(clean(c, r, GRID, frozenset({1, 2, 4})).presence) == 3  # layer 3 off


def test_a_row_the_alien_lists_and_wcvp_both_catch_goes_to_the_lists():
    # WCVP comes last in layer 3, so what it removes is what no alien list had already removed
    r = refs([taxon("p", kingdom="Plantae")], griis=[("p", "AA", None)], wcvp=[("p", "N", True)], tdwg3=REGIONS)
    c = counts(row("p", "13:land", country="AA"), row("p", "13:land", country="BB"))
    out = clean(c, r, WCVP_GRID)
    assert gone(out, 3, "griis") == {("p", "13:land", "AA")}
    assert gone(out, 3, "wcvp_outside_native") == {("p", "13:land", "BB")}


# ---- layer 4

L4 = make_grid(
    {"9": (1000.0, 0.0), "100": (1000.0, 0.0), "12": (0.0, 1000.0), "11": (500.0, 500.0), "500": (1000.0, 0.0), "13": (1000.0, 0.0)},
    [("9", "12"), ("12", "100"), ("100", "13"), ("11", "13")],
)
LAND = taxon("land", environment="terrestrial")
SEA = taxon("sea", environment="marine")


def test_a_land_species_at_sea_folds_into_the_land_part_of_its_own_cell():
    r = refs([LAND])
    c = counts(row("land", "11:land", n=6), row("land", "11:sea", n=2, country="BB"), row("land", "9:land", n=5))
    out = clean(c, r, L4)
    assert kept(out) == {("land", "11:land"): 8, ("land", "9:land"): 5}
    assert gone(out, 4, "land_sea_fold") == {("land", "11:sea", "BB")}
    folded = out.removed.filter(pl.col("rule") == "land_sea_fold").row(0, named=True)
    assert (folded["to_part"], folded["records"]) == ("11:land", 2)
    assert out.presence.filter(pl.col("part") == "11:land").row(0, named=True)["n"] == 8


def test_a_sea_species_on_land_folds_into_the_sea_part_of_its_own_cell():
    out = clean(counts(row("sea", "11:land", n=3), row("sea", "11:sea", n=4)), refs([SEA]), L4)
    assert kept(out) == {("sea", "11:sea"): 7}
    assert gone(out, 4, "land_sea_fold") == {("sea", "11:land", "AA")}


def test_with_no_such_part_in_its_cell_a_row_folds_into_the_neighbour_with_most_of_the_records():
    # 12 is open sea, next to 9 and 100. The land species holds 4 records at 9 and 6 at 100
    c = counts(row("land", "12:sea", n=2), row("land", "9:land", n=4), row("land", "100:land", n=6))
    out = clean(c, refs([LAND]), L4)
    assert kept(out) == {("land", "9:land"): 4, ("land", "100:land"): 8}
    assert gone(out, 4, "land_sea_neighbour") == {("land", "12:sea", "AA")}
    assert out.removed.row(0, named=True)["to_part"] == "100:land"


def test_a_tie_between_neighbours_goes_to_the_lowest_cell_id_as_a_number():
    c = counts(row("land", "12:sea", n=2), row("land", "9:land", n=4), row("land", "100:land", n=4))
    out = clean(c, refs([LAND]), L4)
    assert kept(out) == {("land", "9:land"): 6, ("land", "100:land"): 4}  # 9 is below 100, though "100" sorts before "9" as text


def test_the_neighbour_a_row_folds_into_is_a_part_of_the_right_side_that_holds_records():
    # 100 is next to 12 and has land, but the species has no record there. 11's sea part is the other kind than the species'.
    c = counts(row("land", "12:sea", n=2), row("land", "11:sea", n=1), row("land", "13:land", n=50))
    out = clean(c, refs([LAND]), L4)
    assert gone(out, 4, "land_sea_dropped") == {("land", "12:sea", "AA")}  # 9 and 100 have land and no records, 13 is no neighbour of 12
    assert kept(out) == {("land", "11:land"): 1, ("land", "13:land"): 50}  # the 11:sea record has land in its own cell
    assert out.removed.filter(pl.col("rule") == "land_sea_dropped")["to_part"].null_count() == 1
    assert out.removed.filter(pl.col("rule") == "land_sea_dropped")["records"].to_list() == [2]


def test_a_row_with_no_neighbour_to_fold_into_is_dropped():
    # 12 is open sea, and the species has records far away at 500 but none at 9 or 100, which are next to it
    out = clean(counts(row("land", "12:sea", n=2), row("land", "500:land", n=10)), refs([LAND]), L4)
    assert kept(out) == {("land", "500:land"): 10}
    assert gone(out, 4, "land_sea_dropped") == {("land", "12:sea", "AA")}
    assert out.removed["presence_lost"].to_list() == [True]


def test_a_sea_species_folds_into_the_neighbour_sea_part_that_holds_records():
    c = counts(row("sea", "13:land", n=2), row("sea", "11:sea", n=5))  # 13 is land only and next to 11, which is coastal
    out = clean(c, refs([SEA]), L4)
    assert kept(out) == {("sea", "11:sea"): 7}
    assert gone(out, 4, "land_sea_neighbour") == {("sea", "13:land", "AA")}


@pytest.mark.parametrize("environment, cls", [("freshwater,terrestrial", "Insecta"), ("brackish,marine", "Insecta"), ("terrestrial", "Aves"), (None, "Aves")])
def test_freshwater_brackish_and_all_birds_are_exempt(environment, cls):
    r = refs([taxon("x", environment=environment, cls=cls)])
    c = counts(row("x", "11:land", n=98), row("x", "11:sea", n=1), row("x", "12:sea", n=1))
    out = clean(c, r, L4)
    assert out.removed.height == 0 and len(kept(out)) == 3


def test_exactly_terrestrial_is_enough_whatever_the_records_say():
    c = counts(row("land", "11:land", n=5), row("land", "11:sea", n=5))
    out = clean(c, refs([LAND]), L4)
    assert kept(out) == {("land", "11:land"): 10}
    # "terrestrial" plus another environment is not "exactly terrestrial": the records decide
    r = refs([taxon("land", environment="marine,terrestrial")])
    assert len(kept(clean(c, r, L4))) == 2


def test_at_least_95_percent_of_records_on_one_side_makes_a_species_of_that_side():
    r = refs([taxon("land"), taxon("near"), taxon("sea")])
    c = counts(
        row("land", "11:land", n=19), row("land", "11:sea", n=1),  # exactly 95%
        row("near", "11:land", n=18), row("near", "11:sea", n=2),  # 90%
        row("sea", "11:sea", n=95), row("sea", "11:land", n=5),
    )
    out = clean(c, r, L4)
    assert kept(out) == {("land", "11:land"): 20, ("near", "11:land"): 18, ("near", "11:sea"): 2, ("sea", "11:sea"): 100}


def test_the_share_is_a_setting():
    c = counts(row("near", "11:land", n=18), row("near", "11:sea", n=2))
    out = clean(c, refs([taxon("near")]), L4, settings=Settings(land_share=0.9))
    assert kept(out) == {("near", "11:land"): 20}


def test_where_the_environment_and_the_records_disagree_the_environment_decides():
    # exactly terrestrial, but 99% of its records are in the sea. COL says land, and the part label is only GADM's
    r = refs([taxon("odd", environment="terrestrial"), taxon("shore", environment="marine")])
    c = counts(row("odd", "11:land", n=1), row("odd", "11:sea", n=99), row("shore", "11:land", n=99), row("shore", "11:sea", n=1))
    out = clean(c, r, L4)
    assert kept(out) == {("odd", "11:land"): 100, ("shore", "11:sea"): 100}
    assert gone(out, 4, "land_sea_fold") == {("odd", "11:sea", "AA"), ("shore", "11:land", "AA")}


def test_a_species_where_the_two_disagree_can_be_left_alone_instead():
    r = refs([taxon("odd", environment="terrestrial"), taxon("shore", environment="marine")])
    c = counts(row("odd", "11:land", n=1), row("odd", "11:sea", n=99), row("shore", "11:land", n=99), row("shore", "11:sea", n=1))
    out = clean(c, r, L4, settings=Settings(side_conflict="leave"))
    assert out.removed.height == 0 and len(kept(out)) == 4


def test_the_environment_decides_a_conflict_but_a_freshwater_species_is_still_exempt():
    r = refs([taxon("fw", environment="freshwater,terrestrial")])
    out = clean(counts(row("fw", "11:land", n=1), row("fw", "11:sea", n=99)), r, L4)
    assert out.removed.height == 0


def test_a_species_not_in_the_taxa_is_decided_by_its_records():
    c = counts(row("unknown", "11:land", n=99), row("unknown", "11:sea", n=1))
    assert kept(clean(c, refs(), L4)) == {("unknown", "11:land"): 100}


def test_layer_4_shares_are_taken_after_the_layers_before_it():
    # 90 of the 100 records are removed by GRIIS in AA, leaving 5 on land and 5 at sea... which is no side
    r = refs([taxon("x")], griis=[("x", "AA", None)])
    c = counts(row("x", "11:land", n=90, country="AA"), row("x", "11:land", n=5, country="BB"), row("x", "11:sea", n=5, country="BB"))
    assert len(kept(clean(c, r, L4))) == 2


# ---- the layers together

def test_a_narrow_endemic_survives_every_layer():
    endemic = taxon("endemic", kingdom="Plantae", cls="Magnoliopsida", environment="terrestrial")
    r = refs([endemic, taxon("other")], wcvp=[("endemic", "N", True)], tdwg3=[("10:land", "N")], griis=[("other", "AA", None)], no_wild_range={"other"})
    c = counts(row("endemic", "10:land", n=3), row("other", "10:land", n=50))
    for mode in ("native_only", "introduced", "off"):
        out = clean(c, r, GRID, ALL, Settings(wcvp_mode=mode))
        assert kept(out) == {("endemic", "10:land"): 3}, mode
        assert not (set(out.removed["species_id"]) & {"endemic"})


def test_switching_a_layer_off_restores_its_rows_and_only_its_rows():
    r = refs(
        [taxon("flag"), taxon("ext", extinct=True), taxon("alien"), taxon("land", environment="terrestrial"), taxon("ok")],
        griis=[("alien", "AA", None)],
    )
    c = counts(
        row("flag", "10:land", n=0, n_raw=4),
        row("ext", "10:land"),
        row("alien", "10:land"),
        row("land", "11:land", n=5), row("land", "11:sea", n=5),
        row("ok", "10:land"),
    )
    everything = {("ok", "10:land")}
    assert set(kept(clean(c, r, GRID, ALL))) == everything | {("land", "11:land")}
    assert set(clean(c, r, GRID, ALL).removed["layer"]) == {1, 2, 3, 4}
    # layer k off: what layer k removed returns, and the rest stays removed
    back = {1: ("flag", "10:land"), 2: ("ext", "10:land"), 3: ("alien", "10:land")}
    for layer, comeback in back.items():
        out = clean(c, r, GRID, ALL - {layer})
        assert set(kept(out)) == everything | {("land", "11:land"), comeback}, layer
        assert layer not in set(out.removed["layer"])
    out = clean(c, r, GRID, ALL - {4})
    assert set(kept(out)) == everything | {("land", "11:land"), ("land", "11:sea")}


def test_presence_sums_the_countries_and_takes_the_larger_distinct_counts():
    c = counts(row("s", "10:land", country="AA", n=3, n_datasets=2, n_years=1), row("s", "10:land", country="BB", n=4, n_datasets=1, n_years=7), row("s", "13:land", country=None, n=1))
    out = clean(c, refs([taxon("s")]), GRID)
    assert out.presence.columns == ["species_id", "part", *COUNTS, "records"]
    both = out.presence.filter(pl.col("part") == "10:land").row(0, named=True)
    assert (both["n"], both["n_raw"], both["records"], both["n_datasets"], both["n_years"]) == (7, 7, 7, 2, 7)
    assert out.presence.height == 2


def test_clean_does_not_change_its_input():
    c = counts(row("s", "10:sea", n=3), row("t", "10:land", n=2, n_human_text=2))
    before = c.clone()
    clean(c, refs([taxon("s")]), ZERO)
    assert c.equals(before)


def test_a_part_is_removed_once_by_one_rule():
    r = refs([taxon("both", extinct=True)], no_wild_range={"both"}, griis=[("both", "AA", None)])
    out = clean(counts(row("both", "10:land")), r, GRID)
    assert out.removed.height == 1 and out.removed["rule"].to_list() == ["extinct"]




def test_a_rule_that_is_not_in_rules_raises_before_anything_changes():
    rows = counts(row("s", "10:land")).with_columns(pl.col("n").alias("records"))
    run = _Run(rows)
    with pytest.raises(ValueError, match="no_such_rule"):
        run.drop(rows, "no_such_rule")
    assert run.rows.equals(rows) and not run.removed
    moves = pl.DataFrame({"part": ["10:land"], "to_part": ["11:land"], "rule": ["no_such_rule"]})
    with pytest.raises(ValueError, match="no_such_rule"):
        run.move(moves, ["part"])
    assert run.rows.equals(rows)


# ---- conservation: records are never lost without a rule

# 1 and 7 are land only, 4 has a sea of 0.4 m2 and 5 a land of 0.5 m2 (zero-area parts), 2 and 6 are coastal, 3 and 8 are open sea. A ring, with a chord.
RANDOM_CELLS = {
    "1": (1000.0, 0.0), "2": (500.0, 500.0), "3": (0.0, 1000.0), "4": (1000.0, 0.0004),
    "5": (0.0005, 1000.0), "6": (300.0, 700.0), "7": (1000.0, 0.0), "8": (0.0, 1000.0),
}
RANDOM_EDGES = [("1", "2"), ("2", "3"), ("3", "4"), ("4", "5"), ("5", "6"), ("6", "7"), ("7", "8"), ("8", "1"), ("2", "6")]
# 15 degrees of longitude apart, so that cells are over 1,000 km apart and the neighbour spacing, the median of the edges, is 15 degrees
RANDOM_CENTRES = {cell: (0.0, 15.0 * i) for i, cell in enumerate(RANDOM_CELLS)}
RANDOM_GRID = make_grid(RANDOM_CELLS, RANDOM_EDGES, RANDOM_CENTRES)
# thresholds of layer 5 that the random tables of 8 cells and a few hundred records left by layers 1 to 4 can meet: the design's would need 20 cells and 1,000 records of a group
LOOSE_5 = {"satellite_record_share": 0.3, "satellite_cell_share": 0.6, "group_min_records": 60, "under_max_records": 20, "under_rate_factor": 0.5}
RANDOM_PARTS = [f"{cell}:{kind}" for cell in RANDOM_CELLS for kind in ("land", "sea")]
# the land parts with area, and the TDWG regions they touch: 7 is in none, and 5's land has no area
RANDOM_TDWG3 = [("1:land", "R1"), ("2:land", "R1"), ("2:land", "R2"), ("4:land", "R2"), ("6:land", "R3"), ("7:land", None)]
RANDOM_ADMIN1 = [("2", "AA", "Sub", "centre"), ("6", "AA", "Sub", "overlap")]
COUNTRIES = ["AA", "BB", "CC", None]


def random_case(seed: int):
    """A small random table of country counts, and references that have a hit for every rule of layers 2 to 5 (with `LOOSE_5` for layer 5)."""
    rng = np.random.default_rng(seed)
    ids = [f"s{i}" for i in range(40)]
    kingdoms = ["Plantae", "Animalia", "Fungi"]
    classes = ["Aves", "Insecta", "Magnoliopsida"]
    environments = [None, "terrestrial", "marine", "freshwater", "brackish,marine", "marine,terrestrial"]
    taxa, wcvp, griis, wrims, no_wild_range, griis_no_cells, wrims_no_country, rows = [], [], [], [], [], [], [], []
    for species in ids:
        kingdom = kingdoms[rng.choice(3, p=[0.5, 0.4, 0.1])]
        taxa.append(taxon(species, kingdom, classes[rng.choice(3, p=[0.1, 0.4, 0.5])], environments[rng.integers(len(environments))], bool(rng.random() < 0.1)))
        if kingdom == "Plantae" and rng.random() < 0.8:
            for region in ("R1", "R2", "R3", "R4"):
                draw = rng.random()
                if draw < 0.3:
                    wcvp.append((species, region, True))
                elif draw < 0.5:
                    wcvp.append((species, region, False))
        for country in ("AA", "BB"):
            if rng.random() < 0.15:
                griis.append((species, country, None))
        if rng.random() < 0.1:
            griis.append((species, "AA", "Sub"))
        if rng.random() < 0.05:
            griis_no_cells.append((species, "AA", "Nowhere"))
        if rng.random() < 0.15:
            wrims.append((species, "a", "an area", COUNTRIES[rng.integers(3)]))
        if rng.random() < 0.05:
            wrims_no_country.append((species, "b", "open sea", None))
        for _ in range(rng.integers(1, 7)):
            n = 0 if rng.random() < 0.1 else int(rng.integers(1, 60))
            human = n if rng.random() < 0.15 else int(rng.integers(0, n + 1))
            rows.append(row(species, RANDOM_PARTS[rng.integers(len(RANDOM_PARTS))], COUNTRIES[rng.integers(len(COUNTRIES))], n=n, n_raw=n + int(rng.integers(0, 6)), n_human_text=human))
    no_wild_range = set(rng.choice(ids, 3, replace=False))
    for i in range(5):  # thin species: 2 to 4 records in as many parts, which are 15 degrees or more apart if they are in different cells
        species = f"t{i}"
        taxa.append(taxon(species, "Animalia", "Insecta", None, False))
        for _ in range(rng.integers(2, 4)):
            rows.append(row(species, RANDOM_PARTS[rng.integers(len(RANDOM_PARTS))], COUNTRIES[rng.integers(len(COUNTRIES))], n=1))
    for i in range(3):  # a core of 80 records in one cell and a stray record in another, where the group of every species is `Insects`
        species = f"z{i}"
        taxa.append(taxon(species, "Animalia", "Insecta", None, False))
        at = int(rng.integers(len(RANDOM_PARTS)))
        rows.append(row(species, RANDOM_PARTS[at], "AA", n=80))
        rows.append(row(species, RANDOM_PARTS[(at + 2 * int(rng.integers(1, 4))) % len(RANDOM_PARTS)], "AA", n=1))
    return counts(*rows), refs(taxa, wcvp=wcvp, griis=griis, wrims=wrims, tdwg3=RANDOM_TDWG3, admin1=RANDOM_ADMIN1, griis_no_cells=griis_no_cells, wrims_no_country=wrims_no_country, no_wild_range=no_wild_range)


@pytest.mark.parametrize("heuristics", [{}, LOOSE_5], ids=["design-layer5", "loose-layer5"])
@pytest.mark.parametrize("layer_1", [True, False], ids=["layer1", "no-layer1"])
@pytest.mark.parametrize("side_conflict", ["environment", "leave"])
@pytest.mark.parametrize("wcvp_mode", ["native_only", "introduced", "off"])
@pytest.mark.parametrize("seed", range(4))
def test_records_are_conserved(seed, wcvp_mode, side_conflict, layer_1, heuristics):
    """What is left plus what layers 2 to 5 removed is what came in. Folds move records and the zero-area fold loses none."""
    c, r = random_case(seed)
    layers = ALL if layer_1 else ALL - {1}
    out = clean(c, r, RANDOM_GRID, layers, Settings(wcvp_mode=wcvp_mode, side_conflict=side_conflict, **heuristics))
    came_in = int((c["n"] - c["n_human_text"]).sum() if layer_1 else c["n_raw"].sum())
    lost = out.removed.filter(pl.col("layer") >= 2, pl.col("to_part").is_null())["records"].sum()
    assert int(out.presence["records"].sum()) + int(lost) == came_in
    # a removal has no `to_part` and a move has one, and the layer of every row is that of its rule
    assert out.removed.filter(pl.col("layer") == 0)["to_part"].null_count() == 0
    assert all(RULES[rule] == layer for rule, layer in out.removed.select("rule", "layer").unique().iter_rows())
    assert out.presence["records"].min() > 0
    assert not out.presence.select("species_id", "part").is_duplicated().any()


def test_the_random_tables_reach_every_rule():
    reached = set()
    for seed in range(4):
        for mode in ("native_only", "introduced"):
            c, r = random_case(seed)
            reached |= set(clean(c, r, RANDOM_GRID, ALL, Settings(wcvp_mode=mode, **LOOSE_5)).removed["rule"])
    assert reached == set(RULES)
