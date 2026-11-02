"""The cleaning of occurrence counts, layers 1 to 4 and the way layer 5 is called, on hand-made tables. Layer 5's rules are in `test_outliers`.

The grid has cells named by numbers, so that "lowest cell id" can be told from the lowest id as text ("9" before "100"). A cell is
`(land_km2, sea_km2)`, and a part with less than `ZERO_KM2` of area is none.
"""

from dataclasses import replace
from pathlib import Path

import numpy as np
import polars as pl
import pytest

import realms.clean as clean_module
from clean_helpers import DESIGN_5, GRID, counts, gone, kept, make_grid, refs, row, taxon
from realms.clean import (
    KEY,
    LAYERS,
    RULES,
    SWITCHABLE,
    LIST_CODES_SCHEMA,
    Settings,
    _layer_3,
    _part_table,
    _Run,
    clean,
    code_rows,
    derive_list_codes,
    frozen_list_codes,
    list_code_shares,
    read_list_codes,
)
from realms.clean_report import removal_summary
from realms.grid import ZERO_KM2, part_adjacency, parts
from realms.matrix import cell_taxa
from realms.nullmodel import excess_table
from realms.occurrences import COUNTS

ALL = frozenset(LAYERS)


# ---- settings

def test_default_settings_are_setting_132():
    # Markus chose setting 132 of the calibration's 144 on 2026-11-01: record share 2%, cell share 20%, rate factor 0.5%, group minimum 300
    s = Settings()
    assert (s.satellite_record_share, s.satellite_cell_share, s.under_rate_factor, s.group_min_records) == (0.02, 0.20, 0.005, 300)
    assert s != DESIGN_5
    assert (DESIGN_5.satellite_record_share, DESIGN_5.satellite_cell_share, DESIGN_5.under_rate_factor, DESIGN_5.group_min_records) == (0.01, 0.05, 0.01, 1000)
    assert replace(s, **{k: getattr(DESIGN_5, k) for k in ("satellite_record_share", "satellite_cell_share", "under_rate_factor", "group_min_records")}) == DESIGN_5


def test_settings_hold_the_design_values_but_for_the_layer_5_thresholds():
    s = Settings()
    assert (s.land_share, s.satellite_gap_km) == (0.95, 500)
    assert (s.satellite_min_datasets, s.satellite_min_years) == (2, 3)
    assert (s.under_max_records, s.thin_max_records, s.thin_km) == (5, 5, 1000)
    assert (s.core_share, s.heuristic_min_records, s.wcvp_mode) == (0.9, 20, "native_only")
    assert (s.side_conflict, s.rules_off) == ("environment", frozenset())
    assert s.griis_without_status == "skip"  # a GRIIS list with no distributions is not applied
    with pytest.raises(AttributeError):
        s.land_share = 0.5


def test_settings_refuse_what_they_do_not_know():
    with pytest.raises(ValueError, match="wcvp_mode"):
        Settings(wcvp_mode="strict")
    with pytest.raises(ValueError, match="side_conflict"):
        Settings(side_conflict="records")
    with pytest.raises(ValueError, match="griis_without_status"):
        Settings(griis_without_status="guess")
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


def test_a_sub_national_list_applies_to_the_code_gbif_uses_in_its_cells_and_not_to_its_own():
    # an overseas place carries the sovereign's code in the list (FR), while GBIF's records say RE: RE has most rows in Reunion's cell (2 against MU's 1, and no row
    # says FR there), so the list applies to RE rows, and the row of a neighbour country in the same cell is spared
    r = refs([taxon("s1"), taxon("s2")], griis=[("s1", "FR", "Reunion")], admin1=FR_CELLS)
    c = counts(
        row("s1", "11:land", country="RE"),
        row("s2", "11:land", country="RE"),
        row("s1", "11:sea", country="MU"),
        row("s1", "10:land", country="FR"),
        row("s1", "10:land", country="RE"),
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis_subnational") == {("s1", "11:land", "RE")}
    assert gone(out, 3, "griis") == set()
    assert set(kept(out)) == {("s1", "11:sea"), ("s1", "10:land"), ("s2", "11:land")}  # outside the cells, or another code, the list does not apply


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


def test_a_country_list_applies_to_the_rest_of_the_country_and_a_sub_national_list_to_every_cell_it_overlaps():
    # the cells of Reunion's list are 11 (by centre) and 12 (a cell that overlaps it with its centre at sea): the country list leaves both to Reunion's list
    # and keeps applying to the rest of France, which is 10 and 13
    r = refs(
        [taxon("s1"), taxon("s2")],
        griis=[("s1", "FR", None), ("s2", "FR", "Reunion")],
        admin1=[("11", "FR", "Reunion", "centre"), ("12", "FR", "Reunion", "overlap")],
    )
    c = counts(
        row("s1", "10:land", country="FR"),
        row("s1", "13:land", country="FR"),
        row("s1", "11:land", country="FR"),  # in Reunion's cells, where its list rules, and it does not name s1
        row("s1", "12:sea", country="FR"),
        row("s2", "11:land", country="FR"),
        row("s2", "12:sea", country="FR"),
        row("s2", "10:land", country="FR"),  # the country's list does not name s2
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis") == {("s1", "10:land", "FR"), ("s1", "13:land", "FR")}
    assert gone(out, 3, "griis_subnational") == {("s2", "11:land", "FR"), ("s2", "12:sea", "FR")}
    assert set(kept(out)) == {("s1", "11:land"), ("s1", "12:sea"), ("s2", "10:land")}


US_CELLS = [("11", "US", "United", "overlap"), ("12", "US", "United", "overlap")]


def test_subnational_list_spares_the_neighbour_country_in_a_border_cell():
    # `US` has most rows in the cells (3 against CA's 2), so only US rows count. A CA row in a cell that overlaps a US state is Canada's, and Canada's lists are not this one
    r = refs([taxon("s1"), taxon("s2")], griis=[("s1", "US", "United")], admin1=US_CELLS)
    c = counts(
        row("s1", "11:land", country="US"),
        row("s1", "11:sea", country="US"),
        row("s2", "11:land", country="US"),
        row("s1", "11:land", country="CA"),
        row("s1", "12:sea", country="CA"),
        row("s1", "10:land", country="US"),  # a US row outside the list's cells
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis_subnational") == {("s1", "11:land", "US"), ("s1", "11:sea", "US")}
    assert kept(out) == {("s1", "10:land"): 10, ("s1", "11:land"): 10, ("s1", "12:sea"): 10, ("s2", "11:land"): 10}


def test_a_list_applies_to_the_code_gbif_uses_where_its_own_never_occurs():
    # the Canary Islands are `IC` in GRIIS and `ES` in GBIF, whose records never carry IC: ES has the rows there, so the list applies to ES rows. A row with no
    # country code is no candidate and is judged by no sub-national list
    r = refs([taxon("s1")], griis=[("s1", "IC", "Tenerife")], admin1=[("11", "IC", "Tenerife", "overlap")])
    c = counts(row("s1", "11:land", country="ES"), row("s1", "11:sea", country=None), row("s1", "10:land", country="ES"))
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis_subnational") == {("s1", "11:land", "ES")}
    assert set(kept(out)) == {("s1", "11:sea"), ("s1", "10:land")}


def test_a_list_whose_cells_are_dominated_by_another_code_applies_to_that_codes_rows_only():
    # Svalbard is a list of `NO`, and GBIF codes its records `SJ`: 12 rows of SJ, 1 of NO, 1 with no code and 1 of RU in the list's cell, so NO has 6.7% of the rows, under the
    # 10% that a list keeps its own code at. The list applies to the SJ rows. The NO row of the cell is Norway's, so Norway's own list takes it (it is not left to Svalbard's
    # list, which is not for NO rows), and the others are spared
    r = refs([taxon("s1"), taxon("s2"), taxon("s3")], griis=[("s1", "NO", "Svalbard"), ("s1", "NO", None)], admin1=[("11", "NO", "Svalbard", "overlap")])
    c = counts(
        row("s1", "11:land", country="SJ"),
        row("s2", "11:land", country="SJ"),
        row("s3", "11:sea", country="SJ"),
        *crowd(9, "SJ"),
        row("s1", "11:sea", country="NO"),
        row("s1", "11:land", country=None),
        row("s1", "11:land", country="RU"),
        row("s1", "10:land", country="NO"),
    )
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis_subnational") == {("s1", "11:land", "SJ")}
    assert gone(out, 3, "griis") == {("s1", "11:sea", "NO"), ("s1", "10:land", "NO")}
    assert {k for k in kept(out) if k[0].startswith("s")} == {("s1", "11:land"), ("s2", "11:land"), ("s3", "11:sea")}
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("NO", "Svalbard", "SJ")]


def test_a_list_keeps_its_own_code_where_gbif_uses_it_materially_and_not_below_a_tenth_of_the_rows():
    # Saba is a list of `BQ`, and its cell is mostly the rows of a bigger neighbour, `MF`. The list's own code needs 10% of the rows in the cell to be the one it applies to
    r = refs([taxon("s1")], griis=[("s1", "BQ", "Saba")], admin1=[("11", "BQ", "Saba", "overlap")])
    mine, theirs = row("s1", "11:land", country="BQ"), row("s1", "11:sea", country="MF")
    # 1 BQ row of 9 (11.1%): the list is for BQ rows, and the MF rows around them are spared
    c = counts(mine, theirs, *crowd(7, "MF"))
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("BQ", "Saba", "BQ")]
    assert gone(clean(c, r, GRID), 3, "griis_subnational") == {("s1", "11:land", "BQ")}
    # 1 of 10 is exactly 10%: it still is
    c = counts(mine, theirs, *crowd(8, "MF"))
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("BQ", "Saba", "BQ")]
    assert gone(clean(c, r, GRID), 3, "griis_subnational") == {("s1", "11:land", "BQ")}
    # 1 of 11 is 9.1%: the code with most rows, MF, takes over, and the BQ row is spared
    c = counts(mine, theirs, *crowd(9, "MF"))
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("BQ", "Saba", "MF")]
    assert gone(clean(c, r, GRID), 3, "griis_subnational") == {("s1", "11:sea", "MF")}
    # the rows with no code count in the denominator: 1 of 10 rows, 8 of them with no code, is 10%
    c = counts(mine, theirs, *crowd(8, None))
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("BQ", "Saba", "BQ")]
    # exactly 10% is not lost to float arithmetic (0.1 * 30 is 3.0000000000000004), and it is a share of the rows and not of those with a code
    seen = pl.DataFrame({"cell": ["11", "11"], "countrycode": ["BQ", "MF"], "rows": [3, 27]})
    assert derive_list_codes(r.admin1_cells, seen).rows() == [("BQ", "Saba", "BQ")]
    assert derive_list_codes(r.admin1_cells, seen.with_columns(pl.when(pl.col("countrycode") == "BQ").then(2).otherwise(28).alias("rows"))).rows() == [("BQ", "Saba", "MF")]


def test_a_tie_between_codes_goes_to_the_lists_own_code_and_then_to_the_lower_code():
    # one IC row and one ES row in the cell, and 20 rows with no code that make IC 4.5% of the rows, under a tenth: the codes with most rows are IC and ES, and the list's own
    # code wins the tie
    r = refs([taxon("s1")], griis=[("s1", "IC", "Tenerife")], admin1=[("11", "IC", "Tenerife", "overlap")])
    c = counts(row("s1", "11:land", country="IC"), row("s1", "11:sea", country="ES"), *crowd(20, None))
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("IC", "Tenerife", "IC")]
    assert gone(clean(c, r, GRID), 3, "griis_subnational") == {("s1", "11:land", "IC")}
    # with no row of the list's own code the lower of the tied codes wins: MU before RE
    r = refs([taxon("s1")], griis=[("s1", "FR", "Reunion")], admin1=FR_CELLS)
    c = counts(row("s1", "11:land", country="RE"), row("s1", "11:sea", country="MU"), row("s1", "10:land", country="FR"), row("s1", "13:land", country="FR"))
    assert gone(clean(c, r, GRID), 3, "griis_subnational") == {("s1", "11:sea", "MU")}
    # and FR rows elsewhere do not count: only the rows in the list's cells do
    assert derive_list_codes(r.admin1_cells, code_rows(c)).rows() == [("FR", "Reunion", "MU")]
    # a code with more rows wins over the list's own where the own code has under a tenth: 12 rows of RE against 1 of FR
    c = counts(row("s1", "11:land", country="RE"), row("s1", "12:sea", country="RE"), row("s1", "11:sea", country="FR"), *crowd(10, "RE"))
    r2 = refs([taxon("s1")], griis=[("s1", "FR", "Reunion")], admin1=[("11", "FR", "Reunion", "centre"), ("12", "FR", "Reunion", "overlap")])
    assert gone(clean(c, r2, GRID), 3, "griis_subnational") == {("s1", "11:land", "RE"), ("s1", "12:sea", "RE")}


# ---- the frozen country codes of the sub-national lists

IC = refs([taxon("s1")], griis=[("s1", "IC", "Tenerife")], admin1=[("11", "IC", "Tenerife", "overlap")])


def crowd(n: int, code, part: str = "11:land", prefix: str = "x") -> list[dict]:
    """`n` rows of `code` at `part`, each of a species of its own that no list names: they only count in the cell."""
    return [row(f"{prefix}{i}", part, country=code) for i in range(n)]


def frozen(*rows_) -> pl.DataFrame:
    """A frozen table: (country, subdivision, code) rows."""
    return pl.DataFrame(list(rows_), schema=LIST_CODES_SCHEMA, orient="row")


def test_a_frozen_list_codes_table_overrides_what_the_counts_say():
    # no row of the counts carries `IC`, so the counts say ES and the list applies to the ES row
    c = counts(row("s1", "11:land", country="ES"), row("s1", "10:land", country="ES"))
    assert gone(clean(c, IC, GRID), 3, "griis_subnational") == {("s1", "11:land", "ES")}
    # frozen as IC, the list applies to IC rows: there is none, so nothing is taken
    out = clean(c, IC, GRID, list_codes=frozen(("IC", "Tenerife", "IC")))
    assert gone(out, 3, "griis_subnational") == set()
    assert set(kept(out)) == {("s1", "11:land"), ("s1", "10:land")}
    # and the other way round: IC has 1 of 13 rows in the cell, under a tenth, so the counts say ES and the IC row stays, but frozen as IC it is the IC row that goes
    c = counts(row("s1", "11:land", country="IC"), row("s1", "11:sea", country="ES"), *crowd(11, "ES"))
    assert gone(clean(c, IC, GRID), 3, "griis_subnational") == {("s1", "11:sea", "ES")}
    assert gone(clean(c, IC, GRID, list_codes=frozen(("IC", "Tenerife", "IC"))), 3, "griis_subnational") == {("s1", "11:land", "IC")}


def test_a_frozen_table_that_agrees_with_the_counts_changes_nothing():
    c = counts(row("s1", "11:land", country="IC"), row("s1", "11:land", country="ES"), row("s1", "10:land", country="FR"))
    derived = derive_list_codes(IC.admin1_cells, code_rows(c))
    assert derived.rows() == [("IC", "Tenerife", "IC")]  # the list's own code, with half the rows
    same = clean(c, IC, GRID)
    assert clean(c, IC, GRID, list_codes=derived).removed.equals(same.removed)


def test_the_frozen_table_decides_every_pass_of_layer_3_not_only_the_first(monkeypatch):
    # a land species has a row in the open sea of cell 12 and one on the land of cell 13: layer 4 moves the first to 13:land, and layer 3 judges it there again
    r = refs([taxon("s1", environment="terrestrial")], griis=[("s1", "IC", "Tenerife")], admin1=[("13", "IC", "Tenerife", "overlap")])
    c = counts(row("s1", "12:sea", country="ES"), row("s1", "13:land", country="ES"))
    table = frozen(("IC", "Tenerife", "IC"))
    calls = []
    real = clean_module._layer_3

    def spy(*args, **kwargs):
        calls.append((kwargs.get("only") is not None, kwargs.get("list_codes")))
        return real(*args, **kwargs)

    monkeypatch.setattr(clean_module, "_layer_3", spy)
    out = clean(c, r, GRID, frozenset({1, 2, 3, 4}), list_codes=table)
    assert [first for first, _ in calls] == [False, True]  # the first pass, and the pass over the rows layer 4 moved
    assert all(t is not None and t.equals(table) for _, t in calls)  # the frozen table, both times, and not the counts' own decision
    assert gone(out, 4, "land_sea_neighbour") == {("s1", "12:sea", "ES")} and gone(out, 3) == set()  # applied to IC rows: no ES row is taken
    calls.clear()
    assert gone(clean(c, r, GRID, frozenset({1, 2, 3, 4})), 3, "griis_subnational") != set()  # derived, the code is ES and the rows are taken


def test_a_frozen_table_is_not_read_when_layer_3_is_off():
    c = counts(row("s1", "11:land", country="ES"))
    broken = frozen(("IC", "Tenerife", "IC"), ("IC", "Tenerife", "ES"))  # twice, which layer 3 would refuse
    assert clean(c, IC, GRID, frozenset({1, 2, 4, 5}), list_codes=broken).presence.height == 1


def test_a_frozen_table_that_misses_an_applied_list_is_refused():
    c = counts(row("s1", "11:land", country="ES"))
    with pytest.raises(ValueError, match=r"IC.*Tenerife"):
        clean(c, IC, GRID, list_codes=frozen(("IC", "Gomera", "ES")))
    with pytest.raises(ValueError, match=r"IC.*Tenerife"):
        clean(c, IC, GRID, list_codes=frozen())
    # a place of the table that no list names does no harm, and a list that is not applied needs no decision (under `skip` a name-only list is not applied)
    name_only = refs([taxon("s1")], griis_name_only=[("s1", "IC", "Tenerife")], admin1=[("11", "IC", "Tenerife", "overlap")])
    assert clean(c, name_only, GRID, list_codes=frozen(("FR", "Reunion", "RE"))).presence.height == 1
    with pytest.raises(ValueError, match=r"IC.*Tenerife"):
        clean(c, name_only, GRID, settings=Settings(griis_without_status="apply"), list_codes=frozen(("FR", "Reunion", "RE")))


def test_a_frozen_table_of_the_wrong_shape_is_refused():
    c = counts(row("s1", "11:land", country="ES"))
    with pytest.raises(ValueError, match="columns"):
        clean(c, IC, GRID, list_codes=frozen(("IC", "Tenerife", "ES")).drop("code"))
    with pytest.raises(ValueError, match="code is"):
        clean(c, IC, GRID, list_codes=pl.DataFrame({"country": ["IC"], "subdivision": ["Tenerife"], "code": [1]}))
    with pytest.raises(ValueError, match="twice"):
        clean(c, IC, GRID, list_codes=frozen(("IC", "Tenerife", "IC"), ("IC", "Tenerife", "ES")))
    with pytest.raises(ValueError, match="no value"):
        clean(c, IC, GRID, list_codes=frozen(("IC", "Tenerife", None)))
    with pytest.raises(ValueError, match="no value"):
        clean(c, IC, GRID, list_codes=frozen(("IC", "Tenerife", "")))


def test_the_code_rows_count_the_rows_of_each_code_in_each_cell_with_the_rows_that_have_none():
    c = counts(row("s1", "11:land", country="RE"), row("s2", "11:land", country="RE"), row("s1", "11:sea", country="US"), row("s1", "13:land", country="ES"), row("s1", "13:sea", country=None))
    seen = code_rows(c)
    assert seen.columns == ["cell", "countrycode", "rows"]
    assert seen.sort("cell", "countrycode", nulls_last=True).rows() == [("11", "RE", 2), ("11", "US", 1), ("13", "ES", 1), ("13", None, 1)]  # rows, over parts and species; a row with no code is counted


def test_the_derivation_gives_every_place_of_the_lists_its_own_code_where_it_has_a_tenth_of_the_rows_and_else_the_code_with_most():
    cells = pl.DataFrame(
        [
            ("11", "FR", "Reunion", "centre"),
            ("12", "FR", "Reunion", "overlap"),
            ("15", "US", "Texas", "centre"),
            ("13", "IC", "Tenerife", "overlap"),
            ("14", "XX", "Empty", "overlap"),
            ("16", "NO", "Svalbard", "overlap"),
            ("17", "BQ", "Saba", "overlap"),
        ],
        schema={"cell": pl.String, "country": pl.String, "subdivision": pl.String, "via": pl.String},
        orient="row",
    )
    c = counts(
        row("s1", "11:land", country="RE"),
        row("s2", "11:land", country="RE"),
        row("s1", "11:sea", country="US"),
        row("s3", "12:land", country="MU"),
        row("s1", "10:land", country="FR"),  # not in a cell of a list
        row("s1", "13:land", country="ES"),
        row("s1", "13:sea", country=None),
        row("s1", "15:land", country="US"),
        row("s2", "15:land", country="US"),
        *crowd(12, "SJ", "16:land"),
        row("s1", "16:sea", country="NO"),
        *crowd(8, "MF", "17:land"),
        row("s1", "17:sea", country="BQ"),
    )
    seen = code_rows(c)
    table = derive_list_codes(cells, seen)
    assert table.schema == pl.Schema(LIST_CODES_SCHEMA)
    # Reunion: no row of FR, and RE has 2 of the 4; Tenerife: ES has 1 of 2 (the row with no code is no candidate and no row says IC); Texas is all US; the place with no row has
    # its list's own code; Svalbard: NO has 1 of 13 rows, under a tenth, so SJ; Saba: BQ has 1 of 9, over a tenth, so it keeps BQ though MF has 8
    assert table.rows() == [
        ("BQ", "Saba", "BQ"),
        ("FR", "Reunion", "RE"),
        ("IC", "Tenerife", "ES"),
        ("NO", "Svalbard", "SJ"),
        ("US", "Texas", "US"),
        ("XX", "Empty", "XX"),
    ]
    assert derive_list_codes(cells.head(0), seen).rows() == []
    shares = list_code_shares(cells, seen, table)
    assert shares.columns == ["country", "subdivision", "code", "rows", "with_code", "share", "own_rows", "own_share", "next_code", "next_rows", "low"]
    got = {(r["country"], r["subdivision"]): r for r in shares.iter_rows(named=True)}
    expected = {  # code, rows, with_code, own_rows, next_code, next_rows, low
        ("BQ", "Saba"): ("BQ", 9, 1, 1, "MF", 8, True),
        ("FR", "Reunion"): ("RE", 4, 2, 0, "MU", 1, True),  # the next code is the lower of MU and US, which have one row each
        ("IC", "Tenerife"): ("ES", 2, 1, 0, None, None, True),
        ("NO", "Svalbard"): ("SJ", 13, 12, 1, "NO", 1, False),
        ("US", "Texas"): ("US", 2, 2, 2, None, None, False),
        ("XX", "Empty"): ("XX", 0, 0, 0, None, None, False),
    }
    assert set(got) == set(expected)
    for key, (code, rows_, with_code, own_rows, next_code, next_rows, low) in expected.items():
        r = got[key]
        assert (r["code"], r["rows"], r["with_code"], r["own_rows"], r["next_code"], r["next_rows"], r["low"]) == (code, rows_, with_code, own_rows, next_code, next_rows, low), key
        assert r["share"] == (with_code / rows_ if rows_ else None) and r["own_share"] == (own_rows / rows_ if rows_ else None), key


def test_a_frozen_table_is_read_from_its_file_and_written_back_the_same(tmp_path):
    table = frozen(("US", "Texas", "US"), ("FR", "Reunion", "RE"))
    path = tmp_path / "griis_uses_code-res7.csv"
    table.sort("country", "subdivision").write_csv(path)
    assert path.read_text().splitlines() == ["country,subdivision,code", "FR,Reunion,RE", "US,Texas,US"]
    assert read_list_codes(path).equals(table.sort("country", "subdivision"))
    assert frozen_list_codes(7, tmp_path).equals(read_list_codes(path))
    with pytest.raises(FileNotFoundError, match=r"griis_codes\.py"):
        frozen_list_codes(6, tmp_path)  # the message says what to run
    path.write_text("country,subdivision,uses_code\nFR,Reunion,true\n")  # the first form of the table, which said whether GBIF used the list's code
    with pytest.raises(ValueError, match="columns"):
        read_list_codes(path)
    path.write_text("country,subdivision,code\nFR,Reunion,RE\nFR,Reunion,GP\n")
    with pytest.raises(ValueError, match="twice"):
        read_list_codes(path)


ROOT = Path(__file__).parent.parent
FROZEN_RESOLUTIONS = {7: "0020744-261023074918325", 6: "0019721-261023074918325"}


def derived_from_half_a(resolution: int) -> pl.DataFrame:
    """The table as `scripts/griis_codes.py` derives it, from half A's counts, or a skip where the data is not on this machine."""
    import importlib.util

    key = FROZEN_RESOLUTIONS[resolution]
    needed = [
        ROOT / "data" / "gbif" / key / "_COMPLETE",
        ROOT / "data" / "split" / "halves-315557.parquet",
        ROOT / "data" / "reference" / "manifest.json",
        ROOT / "data" / "reference" / "grids" / f"ISEA3H-{resolution}.gpkg",
    ]
    if missing := [p for p in needed if not p.exists()]:
        pytest.skip(f"{missing[0]} is not there")
    spec = importlib.util.spec_from_file_location("griis_codes_script", ROOT / "scripts" / "griis_codes.py")
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    return script.derive(ROOT / "data" / "gbif" / key, resolution, ROOT / "data" / "reference", ROOT / "data" / "split" / "halves-315557.parquet")


def test_the_frozen_res7_table_equals_the_derivation_from_half_a():
    derived = derived_from_half_a(7)
    frozen_file = read_list_codes(ROOT / "reference" / "griis_uses_code-res7.csv")
    assert frozen_file.equals(derived) and frozen_file.height == 24  # every place of the sub-national lists has its code
    other = {(c, d): k for c, d, k in frozen_file.filter(pl.col("code") != pl.col("country")).iter_rows()}
    # the lists that apply to another code than their own country's, which is where the own code has under a tenth of the rows in the list's cells: Svalbard's records are `SJ`
    # and the French overseas places have codes of their own. The Dutch Caribbean lists keep BQ, which has 17.6% of Saba's and St Eustatius's cells and 45.7% of Bonaire's
    assert other == {("FR", "Guadeloupe"): "GP", ("FR", "Martinique"): "MQ", ("FR", "Reunion"): "RE", ("FR", "Saint-Pierre"): "PM", ("NO", "Svalbard"): "SJ"}
    assert set(frozen_file.filter(pl.col("code") == pl.col("country"))["country"]) == {"BQ", "FM", "KM", "PT", "SH", "TO", "US"}


def test_the_frozen_res6_table_equals_the_derivation_from_half_a():
    derived = derived_from_half_a(6)
    assert read_list_codes(ROOT / "reference" / "griis_uses_code-res6.csv").equals(derived)
    # the cells at resolution 6 are bigger, so more of the codes in them are a neighbour's: the res 6 table is not the res 7 one, and res 7 is the analysis resolution. Saba has
    # 9.9% of its cells' rows in `BQ` there and St Eustatius 11.0%, either side of the tenth that a list keeps its own code at, so that Saba goes to the dominant `VI`
    at_7 = {(c, d): k for c, d, k in read_list_codes(ROOT / "reference" / "griis_uses_code-res7.csv").iter_rows()}
    differs = {(c, d): k for c, d, k in derived.iter_rows() if at_7[(c, d)] != k}
    assert differs == {("BQ", "Saba"): "VI", ("FR", "Martinique"): "GP"}


def test_the_code_gbif_uses_is_read_from_the_counts_and_not_from_what_is_left_of_them():
    # the only IC row has no clean record: layer 1 takes it out before layer 3, but GBIF did use the code
    r = refs([taxon("s1")], griis=[("s1", "IC", "Tenerife")], admin1=[("11", "IC", "Tenerife", "overlap")])
    c = counts(row("s1", "11:land", country="IC", n=3, n_human_text=3), row("s1", "11:land", country="ES"))
    assert gone(clean(c, r, GRID), 3, "griis_subnational") == set()


def test_a_sub_national_list_that_is_not_applied_exempts_nothing_from_the_country_list():
    # France's list has a table and Reunion's has none: under `skip` Reunion's list is not applied, so its cell is not left to it and the country list
    # applies there, as everywhere in France. Under `apply` Reunion's list is applied and the cell is its own
    r = refs([taxon("s1"), taxon("s2")], griis=[("s1", "FR", None)], griis_name_only=[("s2", "FR", "Reunion")], admin1=FR_CELLS)
    c = counts(row("s1", "10:land", country="FR"), row("s1", "11:land", country="FR"), row("s2", "11:land", country="FR"))
    skipped = clean(c, r, GRID)
    assert gone(skipped, 3, "griis") == {("s1", "10:land", "FR"), ("s1", "11:land", "FR")}
    assert gone(skipped, 3, "griis_subnational") == set()
    assert set(kept(skipped)) == {("s2", "11:land")}
    applied = clean(c, r, GRID, settings=Settings(griis_without_status="apply"))
    assert gone(applied, 3, "griis") == {("s1", "10:land", "FR")}
    assert gone(applied, 3, "griis_subnational") == {("s2", "11:land", "FR")}
    assert set(kept(applied)) == {("s1", "11:land")}


def test_a_cell_shared_by_two_subdivisions_of_a_country_takes_both_lists():
    # cell 11 overlaps both of the United States' subdivisions here (a border cell of two states), cell 12 only the second
    cells = [("11", "US", "A", "overlap"), ("11", "US", "B", "overlap"), ("12", "US", "B", "overlap")]
    r = refs([taxon("s1"), taxon("s2")], griis=[("s1", "US", "A"), ("s2", "US", "B")], admin1=cells)
    c = counts(row("s1", "11:land", country="US"), row("s2", "11:land", country="US"), row("s1", "12:sea", country="US"), row("s2", "12:sea", country="US"))
    out = clean(c, r, GRID)
    assert gone(out, 3, "griis_subnational") == {("s1", "11:land", "US"), ("s2", "11:land", "US"), ("s2", "12:sea", "US")}
    assert set(kept(out)) == {("s1", "12:sea")}  # cell 12 is not A's, and s1 is on A's list only


def test_a_list_with_no_distributions_is_applied_only_where_asked():
    # the name-only pairs of a list whose export has no distribution table: not applied by default, applied name-only under `apply`
    r = refs(
        [taxon("s1"), taxon("s2")],
        griis=[("s1", "AA", None)],
        griis_name_only=[("s2", "AA", None), ("s2", "FR", "Reunion")],
        admin1=FR_CELLS,
    )
    c = counts(row("s1", "10:land"), row("s2", "10:land"), row("s2", "11:land", country="FR"))
    skipped = clean(c, r, GRID)
    assert gone(skipped, 3) == {("s1", "10:land", "AA")}
    assert set(kept(skipped)) == {("s2", "10:land"), ("s2", "11:land")}
    applied = clean(c, r, GRID, settings=Settings(griis_without_status="apply"))
    assert gone(applied, 3, "griis") == {("s1", "10:land", "AA"), ("s2", "10:land", "AA")}
    assert gone(applied, 3, "griis_subnational") == {("s2", "11:land", "FR")}


def test_what_a_list_with_no_distributions_could_not_place_is_unplaced_only_where_it_is_applied():
    r = refs(
        [taxon("s1"), taxon("s2")],
        griis=[("s1", "FR", "Clipperton")],
        griis_no_cells=[("s1", "FR", "Clipperton")],
        griis_name_only=[("s2", "FR", "Clipperton")],
        griis_name_only_no_cells=[("s2", "FR", "Clipperton")],
    )
    c = counts(row("s1", "10:land", country="FR", n=5), row("s2", "10:land", country="FR", n=4))
    assert clean(c, r, GRID).unplaced.select("species_id", "records").rows() == [("s1", 5)]
    both = clean(c, r, GRID, settings=Settings(griis_without_status="apply")).unplaced
    assert both.select("species_id", "records").rows() == [("s1", 5), ("s2", 4)]


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
    # RE has 2 rows in Reunion's cell against AA's 1, so Reunion's list applies to RE
    c = counts(row("g", "10:land"), row("n", "11:land", country="RE"), row("n", "11:sea", country="RE"), row("w", "11:sea"))
    assert len(clean(c, r, GRID).presence) == 0
    for rule, species in (("griis", "g"), ("griis_subnational", "n"), ("wrims", "w")):
        out = clean(c, r, GRID, settings=Settings(rules_off=frozenset({rule})))
        assert set(out.presence["species_id"]) == {species}, rule
    assert len(clean(c, r, GRID, frozenset({1, 2, 4})).presence) == 4  # layer 3 off


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


# ---- layer 3 judges again the rows that layer 4 moved, where they went

CRAB = refs([taxon("crab", environment="marine")], wrims=[("crab", "1", "Aa part of the Sea", "AA")])
HERB = refs([taxon("herb", kingdom="Plantae", cls="Magnoliopsida", environment="terrestrial")], wcvp=[("herb", "N", True), ("herb", "X", False)], tdwg3=[("10:land", "N"), ("11:land", "X")])
HERB_ROWS = counts(row("herb", "10:land", n=100), row("herb", "11:land", n=20), row("herb", "11:sea", n=3))


def test_a_fold_does_not_put_back_what_wrims_took_from_the_part_it_folds_into():
    # a marine crab introduced in AA: GBIF puts 50 AA records in the sea part of coastal cell 11 and 5 in its land part. WRiMS takes the sea part,
    # and layer 4 folds the land part's 5 records into it, where WRiMS takes them too
    out = clean(counts(row("crab", "11:sea", n=50), row("crab", "11:land", n=5)), CRAB, GRID, frozenset({1, 2, 3, 4}), Settings(wcvp_mode="off"))
    assert kept(out) == {}
    log = out.removed.select("layer", "rule", "part", "records", "to_part", "presence_lost").sort("layer", "records", descending=[False, True])
    assert log.rows() == [(3, "wrims", "11:sea", 50, None, True), (3, "wrims", "11:sea", 5, None, True), (4, "land_sea_fold", "11:land", 5, "11:sea", True)]


def test_a_fold_does_not_put_back_what_wcvp_took_from_the_part_it_folds_into():
    # a terrestrial plant native in N (cell 10) and introduced in X (cell 11), recorded in coastal cell 11 on land and, by GADM's coast, in its sea part
    for mode, rule in (("native_only", "wcvp_outside_native"), ("introduced", "wcvp_introduced")):
        out = clean(HERB_ROWS, HERB, GRID, frozenset({1, 2, 3, 4}), Settings(wcvp_mode=mode))
        assert kept(out) == {("herb", "10:land"): 100}, mode
        assert sorted(out.removed.filter(pl.col("layer") == 3).select("rule", "part", "records").rows()) == [(rule, "11:land", 3), (rule, "11:land", 20)], mode
        assert gone(out, 4, "land_sea_fold") == {("herb", "11:sea", "AA")}, mode
    # with WCVP off the fold stands
    assert kept(clean(HERB_ROWS, HERB, GRID, frozenset({1, 2, 3, 4}), Settings(wcvp_mode="off"))) == {("herb", "10:land"): 100, ("herb", "11:land"): 23}


def test_a_row_folded_into_the_neighbouring_cell_is_judged_where_it_went():
    # 13 is land only and next to 11. The crab's 2 AA records on 13's land fold into 11's sea part, where its BB records are, and WRiMS takes them
    # there: it lists the crab in AA
    c = counts(row("crab", "13:land", n=2), row("crab", "11:sea", n=10, country="BB"))
    out = clean(c, CRAB, L4, frozenset({1, 2, 3, 4}))
    assert kept(out) == {("crab", "11:sea"): 10}
    assert gone(out, 4, "land_sea_neighbour") == {("crab", "13:land", "AA")}
    assert out.removed.filter(pl.col("rule") == "land_sea_neighbour")["to_part"].to_list() == ["11:sea"]
    assert out.removed.filter(pl.col("layer") == 3).select("rule", "part", "countrycode", "records").rows() == [("wrims", "11:sea", "AA", 2)]
    assert kept(clean(c, CRAB, L4, frozenset({1, 2, 4}))) == {("crab", "11:sea"): 12}  # layer 3 off: the fold stands


def test_a_fold_brings_back_what_layer_3_does_not_judge_out_where_it_went():
    # WRiMS lists the crab in AA and not in BB: the BB records folded into the sea part are the sea part's, though WRiMS emptied it of its AA records.
    # `presence_lost` says what the species had left right after the rule, and so it stays true for WRiMS's row
    c = counts(row("crab", "11:sea", n=50), row("crab", "11:land", n=5, country="BB"))
    out = clean(c, CRAB, GRID, frozenset({1, 2, 3, 4}), Settings(wcvp_mode="off"))
    assert kept(out) == {("crab", "11:sea"): 5}
    assert out.removed.filter(pl.col("layer") == 3).select("rule", "countrycode", "records", "presence_lost").rows() == [("wrims", "AA", 50, True)]


def test_layer_3_judges_the_moved_rows_with_the_same_settings():
    c = counts(row("crab", "11:sea", n=50), row("crab", "11:land", n=5))
    off = clean(c, CRAB, GRID, frozenset({1, 2, 3, 4}), Settings(rules_off=frozenset({"wrims"})))
    assert kept(off) == {("crab", "11:sea"): 55} and not gone(off, 3)
    # a sub-national list with no distribution table, which only `apply` applies. The land species' AA row at open sea (12) folds into the neighbour's
    # land part (100), in the list's cell, where it already has BB records. AA and BB have a row each in cell 100, a tie that goes to the list's own code, so the list applies to AA rows.
    r = refs(
        [taxon("land", environment="terrestrial"), taxon("other")],
        griis_name_only=[("land", "AA", "Sub")],
        admin1=[("100", "AA", "Sub", "overlap")],
    )
    c = counts(row("land", "12:sea", n=2), row("land", "100:land", n=8, country="BB"), row("other", "100:land", n=1))
    skip = clean(c, r, L4, frozenset({1, 2, 3, 4}), Settings(griis_without_status="skip"))
    assert kept(skip)[("land", "100:land")] == 10 and not gone(skip, 3)
    applied = clean(c, r, L4, frozenset({1, 2, 3, 4}), Settings(griis_without_status="apply"))
    assert kept(applied)[("land", "100:land")] == 8
    assert applied.removed.filter(pl.col("layer") == 3).select("rule", "part", "countrycode", "records").rows() == [("griis_subnational", "100:land", "AA", 2)]


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
# thresholds of layer 5 that the random tables of 8 cells and a few hundred records left by layers 1 to 4 can meet, which the defaults (setting 132) and the design's need not
LOOSE_5 = {"satellite_record_share": 0.3, "satellite_cell_share": 0.6, "group_min_records": 60, "under_max_records": 20, "under_rate_factor": 0.5}
RANDOM_PARTS = [f"{cell}:{kind}" for cell in RANDOM_CELLS for kind in ("land", "sea")]
# the land parts with area, and the TDWG regions they touch: 7 is in none, and 5's land has no area
RANDOM_TDWG3 = [("1:land", "R1"), ("2:land", "R1"), ("2:land", "R2"), ("4:land", "R2"), ("6:land", "R3"), ("7:land", None)]
# `Isle` is a list of a country, XX, that no row carries (GBIF does not use its code), so it applies to the code with most rows in its cells; it shares cell 6 with `Sub`,
# whose code AA is used, so that cell takes both lists. `Sub` also reaches cell 5, which has no land, so that a land species' AA row there, which the country's
# list leaves to `Sub`, folds into a neighbour that is the country list's
RANDOM_ADMIN1 = [("2", "AA", "Sub", "centre"), ("6", "AA", "Sub", "overlap"), ("5", "AA", "Sub", "overlap"), ("4", "XX", "Isle", "overlap"), ("6", "XX", "Isle", "overlap")]
# the parts of `Sub`'s cells that have area: both parts of the coastal cells 2 and 6, and cell 5's sea
SUB_PARTS = ["2:land", "2:sea", "6:land", "6:sea", "5:sea"]
# the parts of the shore species: both parts of the coastal cells, the open sea next to them, cell 5's sea and the land of cell 4 next to it
SHORE_PARTS = ["2:land", "2:sea", "6:land", "6:sea", "3:sea", "8:sea", "7:land", "5:sea", "4:land"]
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
    # the lists of layer 3's GRIIS rules that the draws above leave out, from a stream of their own so that those draws stay what they were: the cells-alone
    # list `Isle`, and the lists with no distribution table, which only `griis_without_status="apply"` applies
    extra = np.random.default_rng(seed + 1000)
    name_only, name_only_no_cells = [], []
    for species in ids:
        if extra.random() < 0.1:
            griis.append((species, "XX", "Isle"))
        for country, subdivision in (("AA", None), ("BB", None), ("AA", "Sub"), ("XX", "Isle")):
            if extra.random() < 0.15:
                name_only.append((species, country, subdivision))
        if extra.random() < 0.05:
            name_only_no_cells.append((species, "BB", "Nowhere"))
    # shore species, which layer 4 folds across the coast and into a neighbour, into parts and cells that layer 3 judges: from a stream of their own again
    shore = np.random.default_rng(seed + 2000)
    for i in range(12):
        species = f"c{i}"
        plant = bool(shore.random() < 0.5)
        environment = "terrestrial" if plant or shore.random() < 0.5 else "marine"
        taxa.append(taxon(species, "Plantae" if plant else "Animalia", "Magnoliopsida" if plant else "Insecta", environment, False))
        if plant:
            wcvp.append((species, "R3", True))  # native only in cell 6's land, so 2:land and 4:land are outside
            if shore.random() < 0.5:
                wcvp.append((species, "R1", False))
        if shore.random() < 0.4:
            wrims.append((species, "a", "an area", COUNTRIES[shore.integers(2)]))
        for country, subdivision in (("AA", "Sub"), ("AA", None), ("BB", None)):
            if shore.random() < 0.5:
                griis.append((species, country, subdivision))
        home = COUNTRIES[shore.integers(2)]  # most of a species' rows are of one country, so that a fold meets a row of its own country
        for part in SHORE_PARTS:
            if shore.random() < 0.6:
                country = home if shore.random() < 0.7 else COUNTRIES[shore.integers(3)]
                rows.append(row(species, part, country, n=int(shore.integers(1, 30))))
    # rows of AA in the cells of `Sub`, from a stream of their own, so that AA has most rows there in every seed and is the code that `Sub`'s list applies to: the
    # country list leaves those rows to `Sub`, and layer 4 folds some of them into places where the country list takes them
    sub = np.random.default_rng(seed + 3000)
    have = {(r["species_id"], r["countrycode"], r["part"]) for r in rows}
    for species in ids:
        for _ in range(2):
            part = SUB_PARTS[sub.integers(len(SUB_PARTS))]
            if sub.random() < 0.8 and (species, "AA", part) not in have:
                have.add((species, "AA", part))
                rows.append(row(species, part, "AA", n=int(sub.integers(1, 30))))
    # and rows of BB in cell 4, which only `Isle` has, so that `Isle`'s effective code is BB and not the AA of cell 6: the rows of AA in cell 4 are then the country list's
    isle = np.random.default_rng(seed + 4000)
    for species in ids:
        for part in ("4:land", "4:sea"):
            if isle.random() < 0.9 and (species, "BB", part) not in have:
                have.add((species, "BB", part))
                rows.append(row(species, part, "BB", n=int(isle.integers(1, 30))))
    return counts(*rows), refs(
        taxa, wcvp=wcvp, griis=griis, wrims=wrims, tdwg3=RANDOM_TDWG3, admin1=RANDOM_ADMIN1, griis_no_cells=griis_no_cells, wrims_no_country=wrims_no_country,
        no_wild_range=no_wild_range, griis_name_only=name_only, griis_name_only_no_cells=name_only_no_cells,
    )


@pytest.mark.parametrize("heuristics", [{}, LOOSE_5], ids=["default-layer5", "loose-layer5"])
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


@pytest.mark.parametrize("griis_without_status", ["skip", "apply"])
@pytest.mark.parametrize("wcvp_mode", ["native_only", "off"])
@pytest.mark.parametrize("seed", range(6))
def test_records_are_conserved_with_every_branch_of_the_griis_rules(seed, wcvp_mode, griis_without_status):
    """The same sum with the lists that have no table applied or not, and with the sub-national list whose own code no row carries (`XX`), which shares a cell
    with one that applies to its own (`AA`). The random world has both in every seed."""
    c, r = random_case(seed)
    out = clean(c, r, RANDOM_GRID, ALL, Settings(wcvp_mode=wcvp_mode, griis_without_status=griis_without_status, **LOOSE_5))
    lost = out.removed.filter(pl.col("layer") >= 2, pl.col("to_part").is_null())["records"].sum()
    assert int(out.presence["records"].sum()) + int(lost) == int((c["n"] - c["n_human_text"]).sum())
    assert not out.presence.select("species_id", "part").is_duplicated().any()
    # layer 3 takes a row twice only where a fold of layer 4 brought records back to a row it had emptied, and then by the same rule: once in each pass
    taken = out.removed.filter(pl.col("layer") == 3).group_by("species_id", "part", "countrycode").agg(pl.len().alias("times"), pl.col("rule").n_unique().alias("rules"))
    twice = taken.filter(pl.col("times") > 1)
    assert taken["times"].max() is None or taken["times"].max() <= 2
    assert (twice["rules"] == 1).all()
    folded_to = out.removed.filter(pl.col("layer") == 4, pl.col("to_part").is_not_null()).select("species_id", pl.col("to_part").alias("part"), "countrycode")
    assert twice.join(folded_to, on=["species_id", "part", "countrycode"], how="anti", nulls_equal=True).is_empty()


def end_rows(c: pl.DataFrame, out, layer_1: bool = True) -> pl.DataFrame:
    """The country rows that `clean` leaves (`KEY` and `records`), made from its input and its log alone: what each (species, country, part) came in with,
    less the records the log took or moved from it, plus those it moved there. Layer 1's rows are left out, as they log `n_raw` and take no records."""
    records = (pl.col("n") - pl.col("n_human_text")) if layer_1 else pl.col("n_raw")
    came = c.group_by(KEY).agg(records.sum().alias("records"))
    log = out.removed.filter(pl.col("layer") != 1)
    took = log.select(*KEY, (-pl.col("records")).alias("records"))
    brought = log.drop_nulls("to_part").select("species_id", "countrycode", pl.col("to_part").alias("part"), "records").select(*KEY, "records")
    every = pl.concat([came, took, brought]).group_by(KEY).agg(pl.col("records").sum())
    assert every["records"].min() >= 0, "the log takes more from a row than it holds"
    return every.filter(pl.col("records") > 0)


@pytest.mark.parametrize("layers", [ALL, frozenset({1, 2, 3, 4})], ids=["layers1-5", "layers1-4"])
@pytest.mark.parametrize("griis_without_status", ["skip", "apply"])
@pytest.mark.parametrize("wcvp_mode", ["native_only", "introduced", "off"])
@pytest.mark.parametrize("seed", range(6))
def test_layer_3_would_take_nothing_of_what_clean_leaves(seed, wcvp_mode, griis_without_status, layers):
    """Layer 3 judges a row by where it is, so no fold of layer 4 may leave a row there that its rules take. The rows `clean` leaves are made from the input
    and the log (`end_rows`), which also checks that the log accounts for every (species, country, part), and layer 3 run on them again takes nothing."""
    c, r = random_case(seed)
    settings = Settings(wcvp_mode=wcvp_mode, griis_without_status=griis_without_status, **LOOSE_5)
    out = clean(c, r, RANDOM_GRID, layers, settings)
    end = end_rows(c, out)
    by_part = end.group_by("species_id", "part").agg(pl.col("records").sum()).sort("species_id", "part")
    assert by_part.equals(out.presence.select("species_id", "part", "records"))
    run = _Run(end)
    _layer_3(run, r, _part_table(RANDOM_GRID), settings, code_rows(c))
    assert not run.removed and run.rows.equals(end)


def test_the_random_world_reaches_layer_3_after_the_folds():
    """The rows that layer 4 moved and layer 3 took where they went, over the seeds: for each rule that the random world reaches that way, and the
    rows it took twice, once before the fold and once after. `griis_subnational` is reached by the hand-made test of the neighbour above."""
    after, twice = set(), set()
    for seed in range(6):
        for mode in ("native_only", "introduced", "off"):
            for status in ("skip", "apply"):
                c, r = random_case(seed)
                removed = clean(c, r, RANDOM_GRID, ALL, Settings(wcvp_mode=mode, griis_without_status=status, **LOOSE_5)).removed
                folded_to = removed.filter(pl.col("layer") == 4, pl.col("to_part").is_not_null()).select("species_id", pl.col("to_part").alias("part"), "countrycode")
                taken = removed.filter(pl.col("layer") == 3)
                after |= set(taken.join(folded_to, on=["species_id", "part", "countrycode"], how="semi", nulls_equal=True)["rule"])
                twice |= set(taken.filter(pl.struct("species_id", "part", "countrycode").is_duplicated())["rule"])
    assert after >= {"griis", "wrims", "wcvp_outside_native", "wcvp_introduced"}
    assert twice >= {"wrims", "wcvp_outside_native", "wcvp_introduced"}


def test_the_random_world_reaches_a_list_that_applies_to_a_code_not_its_own_and_the_switch():
    by_alone = by_skip = by_apply = 0
    for seed in range(6):
        c, r = random_case(seed)
        skipped = clean(c, r, RANDOM_GRID, ALL, Settings(wcvp_mode="off", **LOOSE_5))
        applied = clean(c, r, RANDOM_GRID, ALL, Settings(wcvp_mode="off", griis_without_status="apply", **LOOSE_5))
        # a row in cell 4, which only `Isle` has, is taken where it carries `Isle`'s effective code (no row says XX)
        alone = skipped.removed.filter(pl.col("rule") == "griis_subnational", pl.col("part").str.starts_with("4:"))
        by_alone += alone.height
        by_skip += skipped.removed.filter(pl.col("layer") == 3, pl.col("rule").str.starts_with("griis")).height
        by_apply += applied.removed.filter(pl.col("layer") == 3, pl.col("rule").str.starts_with("griis")).height
        assert not c.filter(pl.col("countrycode") == "XX").height
    assert by_alone > 0 and by_apply > by_skip > 0


def test_the_random_tables_reach_every_rule():
    reached = set()
    for seed in range(4):
        for mode in ("native_only", "introduced"):
            c, r = random_case(seed)
            reached |= set(clean(c, r, RANDOM_GRID, ALL, Settings(wcvp_mode=mode, **LOOSE_5)).removed["rule"])
    assert reached == set(RULES)
