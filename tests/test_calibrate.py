"""The calibration of the layer 5 thresholds on WCVP plants, and the known cases, on hand-made tables.

The world is a chain of 60 cells along the equator, 3 degrees (334 km) apart, so that cells are of one patch when their centres are at most
834 km apart (as in `test_outliers`), and four cells off the chain, where a rare record can sit. Every plant is native in `NAT`, the region of
the whole chain, and `p2` is native in `ISL` too. The cells off the chain are `FAR` (100, 102, 103) and `ISL` (101).

- `p1` has a satellite of 1 record in cell 100, outside its native regions.
- `p2` has one in cell 101, where it is native: a native presence that the heuristic would lose.
- `p3` has one of 8 records in cell 102, outside.
- `p4` has 3 records in cell 103, outside, 400 km from the chain and so joined to its range, in a cell where `heavy` has 2,000 records of the same group.

So, at the design's thresholds, layer 5 takes all four, and at other settings it takes some of them (see `test_the_settings_decide_what_is_taken`).
"""

import importlib.util
import itertools
import json
import re
import zipfile
from dataclasses import replace
from pathlib import Path

import polars as pl
import pytest

from clean_helpers import GRID, counts, make_grid, refs, row, taxon
from realms.calibrate import (
    CASES,
    GRID as THRESHOLDS,
    NATIVE_LOSS_LIMIT,
    Baseline,
    Case,
    area_cells,
    baseline,
    calibrate,
    calibration_counts,
    calibration_species,
    choose,
    evaluate_case,
    known_cases,
    native_status,
    pipeline_tally,
    rates,
    region_loss,
    resolve_species,
    score,
    settings_grid,
    tally,
    tdwg_regions,
    top_regions,
)
from realms.clean import Settings, clean, outlier_rows
from realms.occurrences import COLUMNS, COUNTS as COUNT_COLUMNS, SealedHalf, ingest, load_country_counts

ROOT = Path(__file__).parent.parent
FAR, ISLAND = ("100", "102", "103"), "101"
OFF_CHAIN = {"100": (50.0, 0.0), "101": (50.0, 40.0), "102": (50.0, 80.0), "103": (3.0, 10.0)}
CHAIN = [str(i) for i in range(1, 61)]
PLANTS = ("p1", "p2", "p3", "p4")


def world():
    cells = {c: (500.0, 500.0) for c in CHAIN + list(OFF_CHAIN)}
    centres = {str(i): (0.0, (i - 1) * 3.0) for i in range(1, 61)} | OFF_CHAIN
    return make_grid(cells, [(str(i), str(i + 1)) for i in range(1, 60)], centres)


def land(species, cells: dict, **given):
    return [row(species, f"{cell}:land", n=n, **given) for cell, n in cells.items()]


def range_of(species, extra: dict):
    """1,000 records over the first 45 cells of the chain (22 each, 32 in the first), and `extra`."""
    return land(species, {str(i): 22 for i in range(1, 46)} | {"1": 32} | extra)


COUNTS = counts(
    *range_of("p1", {"100": 1}),
    *range_of("p2", {"101": 1}),
    *range_of("p3", {"102": 8}),
    *range_of("p4", {"103": 3}),
    *land("heavy", {"103": 2000}),
    *land("bird", {"1": 500, "2": 500}),
)
TAXA = [taxon(s, "Plantae", "Magnoliopsida") for s in (*PLANTS, "heavy")] + [taxon("bird", "Animalia", "Aves")]
GROUPS = [(s, "angiosperms") for s in (*PLANTS, "heavy")] + [("bird", "birds")]
WCVP = [(s, "NAT", True) for s in PLANTS] + [("p2", "ISL", True), ("p1", "FAR", False)]
TDWG = [(f"{c}:land", "NAT") for c in CHAIN] + [(f"{c}:land", "FAR") for c in FAR] + [(f"{ISLAND}:land", "ISL")]
REFS = refs(TAXA, GROUPS, WCVP, tdwg3=TDWG)
WORLD = world()
LAYERS_WITHOUT_3 = frozenset({1, 2, 4, 5})
BEFORE_5 = frozenset({1, 2, 4})


def status_table(*rows_) -> pl.DataFrame:
    return pl.DataFrame(list(rows_), schema={"species_id": pl.String, "part": pl.String, "status": pl.String}, orient="row")


def pairs(*rows_) -> pl.DataFrame:
    return pl.DataFrame(list(rows_), schema={"species_id": pl.String, "part": pl.String}, orient="row")


# ---- native_status

def test_native_status_calls_a_part_native_when_it_touches_a_native_region_and_outside_when_it_touches_none():
    wcvp = pl.DataFrame(
        [("s1", "R1", True), ("s1", "R2", False), ("s2", "R2", True), ("s4", "R1", False)],
        schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean},
        orient="row",
    )
    parts_ = pl.DataFrame(
        [("1:land", "R1"), ("2:land", "R2"), ("3:land", "R1"), ("3:land", "R2"), ("4:land", None), ("6:land", "R3")],
        schema={"part": pl.String, "tdwg3": pl.String},
        orient="row",
    )
    present = pairs(*[(s, p) for s in ("s1", "s2", "s3", "s4") for p in ("1:land", "2:land", "3:land", "4:land", "5:sea", "6:land")])
    out = native_status(present, wcvp, parts_)
    assert out.columns == ["species_id", "part", "status"]
    got = {(s, p): v for s, p, v in out.iter_rows()}
    assert got == {
        ("s1", "1:land"): "native",
        ("s1", "2:land"): "outside",  # introduced there
        ("s1", "3:land"): "native",  # touches a native and an introduced region
        ("s1", "6:land"): "outside",  # a region WCVP does not list the species in
        ("s2", "1:land"): "outside",
        ("s2", "2:land"): "native",
        ("s2", "3:land"): "native",
        ("s2", "6:land"): "outside",
    }
    # s3 is not in WCVP and s4 has no native region, so neither has a status; a part with no region or on the sea has none either
    assert not {"s3", "s4"} & set(out["species_id"])
    assert not {"4:land", "5:sea"} & set(out["part"])


def test_native_status_counts_a_presence_once_and_takes_a_table_with_counts():
    wcvp = pl.DataFrame([("s1", "R1", True)], schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean}, orient="row")
    parts_ = pl.DataFrame([("1:land", "R1")], schema={"part": pl.String, "tdwg3": pl.String}, orient="row")
    present = counts(row("s1", "1:land", "AA"), row("s1", "1:land", "BB"))  # a row per country
    assert native_status(present, wcvp, parts_).height == 1
    assert native_status(present.head(0), wcvp, parts_).columns == ["species_id", "part", "status"]


def test_native_status_agrees_with_layer_3_in_native_only_mode():
    out = clean(COUNTS, REFS, WORLD, frozenset({3}), Settings(wcvp_mode="native_only"))
    taken = {(s, p) for s, p in out.removed.filter(pl.col("rule") == "wcvp_outside_native").select("species_id", "part").iter_rows()}
    status = native_status(COUNTS, REFS.wcvp, REFS.tdwg3_parts)
    assert taken == {(s, p) for s, p, v in status.iter_rows() if v == "outside"}
    assert taken == {("p1", "100:land"), ("p3", "102:land"), ("p4", "103:land")}


# ---- rates

def test_rates_are_the_share_of_outside_presences_removed_and_of_native_presences_lost():
    status = status_table(*[("a", f"{i}:land", "outside") for i in range(4)], *[("b", f"{i}:land", "native") for i in range(6)])
    removed = pairs(("a", "0:land"), ("a", "1:land"), ("b", "0:land"), ("c", "9:land"), ("a", "0:land"))  # c has no status, and a repeats
    assert rates(status, removed) == (0.5, 1 / 6)
    assert tally(status, removed) == {"outside": 4, "outside_removed": 2, "native": 6, "native_lost": 1}
    assert rates(status, removed.head(0)) == (0.0, 0.0)


def test_rates_count_presences_and_not_records():
    status = status_table(("a", "0:land", "native"), ("a", "1:land", "native"))
    # the same (species, part) twice, as two rows of countries, is one presence
    assert tally(status, pairs(("a", "0:land"), ("a", "0:land")))["native_lost"] == 1


def test_rates_are_not_a_number_where_a_class_is_empty():
    only_native = status_table(("a", "0:land", "native"))
    outside, native = rates(only_native, pairs(("a", "0:land")))
    assert outside != outside and native == 1.0
    assert tally(status_table(), pairs()) == {"outside": 0, "outside_removed": 0, "native": 0, "native_lost": 0}


# ---- the settings

def test_the_grid_has_81_settings_with_the_designs_values():
    grid = settings_grid()
    assert len(grid) == 81 == len(set(grid))
    assert THRESHOLDS == {
        "satellite_record_share": (0.005, 0.01, 0.02),
        "satellite_cell_share": (0.02, 0.05, 0.10),
        "under_rate_factor": (0.005, 0.01, 0.02),
        "group_min_records": (300, 1000, 3000),
    }
    assert {(s.satellite_record_share, s.satellite_cell_share, s.under_rate_factor, s.group_min_records) for s in grid} == set(itertools.product(*THRESHOLDS.values()))
    assert grid[0] == Settings(satellite_record_share=0.005, satellite_cell_share=0.02, under_rate_factor=0.005, group_min_records=300)
    assert grid[1].group_min_records == 1000 and grid[3].under_rate_factor == 0.01  # the last varies fastest
    assert grid[40] == Settings()  # the design's own setting is the middle of the grid


def test_the_grid_keeps_the_other_settings_of_its_base():
    base = Settings(heuristic_min_records=30, wcvp_mode="introduced", land_share=0.9)
    assert all((s.heuristic_min_records, s.wcvp_mode, s.land_share) == (30, "introduced", 0.9) for s in settings_grid(base))


# ---- choosing a setting

def table(*rows_, outside=100, native=1000) -> pl.DataFrame:
    """Settings with (outside_removed, native_lost), in grid order."""
    return pl.DataFrame(
        [{"setting": i, "outside": outside, "outside_removed": o, "native": native, "native_lost": n} for i, (o, n) in enumerate(rows_)]
    )


def test_choose_takes_the_most_outside_removed_among_the_settings_within_the_limit():
    assert NATIVE_LOSS_LIMIT == 0.05
    assert choose(table((10, 5), (60, 49), (90, 51), (50, 0))) == (1, True)  # 90 loses 5.1% and is out


def test_choose_breaks_a_tie_by_lower_native_loss_and_then_by_the_order_of_the_grid():
    assert choose(table((60, 40), (60, 20), (60, 20), (50, 0))) == (1, True)
    assert choose(table((60, 20), (60, 20))) == (0, True)


def test_a_native_loss_of_exactly_the_limit_is_within_it_whatever_the_float():
    # 3 of 60 is 5% and 0.05 * 60 is 3.0000000000000004, but 0.1 * 3 style noise must not push it out either
    assert choose(table((10, 3), (5, 0), outside=10, native=60)) == (0, True)
    assert choose(table((10, 150), (5, 0), outside=10, native=3000)) == (0, True)
    assert choose(table((10, 151), (5, 0), outside=10, native=3000)) == (1, True)


def test_choose_reports_when_no_setting_is_within_the_limit_and_takes_the_lowest_native_loss():
    assert choose(table((90, 70), (30, 60), (95, 60), (10, 80))) == (1, False)  # 60 twice: the lower index
    assert choose(table((90, 70), (30, 60), (95, 60), (10, 80)), limit=0.07) == (2, True)  # 70 is within 7%, and 95 is the most removed


def test_choose_refuses_an_empty_table():
    with pytest.raises(ValueError, match="no settings"):
        choose(table())


# ---- the plants that are calibrated

def test_the_calibration_runs_on_the_plants_wcvp_gives_a_native_region_and_on_the_groups_of_those():
    taxa = [
        taxon("plant", "Plantae"),
        taxon("mate", "Plantae"),  # no WCVP, but of the plant's group: its records count in the group's total
        taxon("alga", "Chromista"),  # likewise, of a group that spans kingdoms
        taxon("bird", "Animalia", "Aves"),
        taxon("introduced_only", "Plantae"),  # WCVP lists it but with no native region: no status, but not a calibration plant
    ]
    groups = [("plant", "angiosperms"), ("mate", "angiosperms"), ("alga", "angiosperms"), ("bird", "birds"), ("introduced_only", "weeds")]
    wcvp = [("plant", "R1", True), ("introduced_only", "R1", False)]
    r = refs(taxa, groups, wcvp)
    assert set(calibration_species(r)) == {"plant", "mate", "alga"}
    c = counts(*[row(s, "1:land") for s in ("plant", "mate", "alga", "bird", "introduced_only", "unknown")])
    assert set(calibration_counts(c, r)["species_id"]) == {"plant", "mate", "alga"}


def test_a_plant_with_no_group_is_calibrated_on_its_own():
    r = refs([taxon("plant", "Plantae")], groups=[], wcvp=[("plant", "R1", True)])
    assert set(calibration_species(r)) == {"plant"}


# ---- calibrate

def reference(setting: Settings) -> tuple[dict, dict]:
    """What `clean` itself says: the rates and counts from layers 1, 2 and 4, and layers 1, 2, 4 and 5 at `setting`, on the plants."""
    plants = calibration_counts(COUNTS, REFS)
    before = clean(plants, REFS, WORLD, BEFORE_5, setting)
    after = clean(plants, REFS, WORLD, LAYERS_WITHOUT_3, setting)
    removed = before.presence.join(after.presence, on=["species_id", "part"], how="anti")
    status = native_status(before.presence, REFS.wcvp, REFS.tdwg3_parts)
    return tally(status, removed), dict(zip(("outside_rate", "native_loss"), rates(status, removed), strict=True))


def test_calibrate_gives_one_row_per_setting_with_the_four_thresholds_and_the_counts():
    grid = settings_grid()[:5]
    out = calibrate(COUNTS, grid, REFS, WORLD)
    assert out.height == 5
    assert out.columns == [
        "setting", "satellite_record_share", "satellite_cell_share", "under_rate_factor", "group_min_records",
        "outside", "outside_removed", "native", "native_lost", "outside_rate", "native_loss",
    ]  # fmt: skip
    assert out["setting"].to_list() == [0, 1, 2, 3, 4]
    assert out["group_min_records"].to_list() == [s.group_min_records for s in grid]


def test_calibrate_equals_cleaning_with_every_setting_of_the_grid():
    grid = settings_grid()
    out = calibrate(COUNTS, grid, REFS, WORLD)
    assert out.height == 81
    for setting, got in zip(grid, out.iter_rows(named=True), strict=True):
        counted, shares = reference(setting)
        assert {k: got[k] for k in counted} == counted
        assert got["outside_rate"] == shares["outside_rate"] and got["native_loss"] == shares["native_loss"]


def test_the_settings_decide_what_is_taken():
    grid = settings_grid()
    out = calibrate(COUNTS, grid, REFS, WORLD)
    # 3 outside presences (p1, p3, p4 off the chain) and 181 native ones (45 for each of 4 plants, and p2 in the island)
    by = {(s.satellite_record_share, s.satellite_cell_share, s.under_rate_factor, s.group_min_records): r for s, r in zip(grid, out.iter_rows(named=True), strict=True)}
    design = by[(0.01, 0.05, 0.01, 1000)]
    assert (design["outside"], design["outside_removed"], design["native"], design["native_lost"]) == (3, 3, 181, 1)
    # p3's 8 records are 0.79% of its 1,008: kept by a record share of 0.5%
    assert by[(0.005, 0.05, 0.01, 1000)]["outside_removed"] == 2
    # a patch of 1 of 46 cells is 2.2%, not under 2%: no satellite goes, and p4's under-reported cell goes alone
    assert (by[(0.01, 0.02, 0.01, 1000)]["outside_removed"], by[(0.01, 0.02, 0.01, 1000)]["native_lost"]) == (1, 0)
    # the group has 2,003 records in cell 103: not 3,000
    assert by[(0.01, 0.05, 0.01, 3000)]["outside_removed"] == 2
    # a rate of 0.15% is not under 0.5% of p4's core rate of 25%
    assert by[(0.01, 0.02, 0.005, 1000)]["outside_removed"] == 0
    assert design["outside_rate"] == 1.0 and design["native_loss"] == pytest.approx(1 / 181)


def test_calibrate_can_take_the_baseline_a_caller_made():
    start = baseline(COUNTS, REFS, WORLD)
    assert isinstance(start, Baseline) and start.status.columns == ["species_id", "part", "status"]
    grid = settings_grid()[::20]
    assert calibrate(COUNTS, grid, REFS, WORLD, start).equals(calibrate(COUNTS, grid, REFS, WORLD))


def test_the_pipeline_tally_counts_what_every_layer_took():
    raw = pairs(("p1", "1:land"), ("p1", "100:land"), ("p1", "2:land"), ("p1", "2:sea"), ("zoo", "1:land"))
    final = pairs(("p1", "1:land"))
    # p1 was native at 1 and 2 and outside at 100, and had a sea part with no status: the whole run took one native and one outside presence
    wcvp = pl.DataFrame([("p1", "NAT", True)], schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean}, orient="row")
    parts_ = pl.DataFrame([("1:land", "NAT"), ("2:land", "NAT"), ("100:land", "FAR")], schema={"part": pl.String, "tdwg3": pl.String}, orient="row")
    assert pipeline_tally(raw, final, wcvp, parts_) == {"outside": 1, "outside_removed": 1, "native": 2, "native_lost": 1}


def test_calibrate_leaves_layer_3_off():
    # with WCVP's layer on, the outside presences would be gone before layer 5 ever saw them
    out = calibrate(COUNTS, [Settings()], REFS, WORLD)
    assert out["outside"].to_list() == [3]


def test_calibrate_refuses_a_grid_that_differs_in_what_layers_1_2_and_4_use():
    with pytest.raises(ValueError, match="land_share"):
        calibrate(COUNTS, [Settings(), Settings(land_share=0.9)], REFS, WORLD)
    with pytest.raises(ValueError, match="no settings"):
        calibrate(COUNTS, [], REFS, WORLD)


def test_calibrate_refuses_a_table_with_no_outside_or_no_native_presences():
    only_native = refs(TAXA, GROUPS, [(s, "NAT", True) for s in PLANTS] + [("p2", "ISL", True), ("p1", "FAR", True), ("p3", "FAR", True), ("p4", "FAR", True)], tdwg3=TDWG)
    with pytest.raises(ValueError, match="outside"):
        calibrate(COUNTS, [Settings()], only_native, WORLD)


def test_outlier_rows_are_what_clean_logs_for_layer_5():
    plants = calibration_counts(COUNTS, REFS)
    before = clean(plants, REFS, WORLD, BEFORE_5)
    after = clean(plants, REFS, WORLD, LAYERS_WITHOUT_3)
    rows_ = outlier_rows(before.presence, REFS, WORLD, Settings())
    assert {"species_id", "part", "rule", "records"} <= set(rows_.columns)
    logged = after.removed.filter(pl.col("layer") == 5)
    assert set(rows_.select("species_id", "part", "rule").iter_rows()) == set(logged.select("species_id", "part", "rule").iter_rows())
    assert rows_.height == 4


# ---- scoring one setting (the hold-out)

def test_score_is_the_calibration_of_one_setting_as_a_dict():
    setting = Settings()
    got = score(COUNTS, setting, REFS, WORLD)
    row_ = calibrate(COUNTS, [setting], REFS, WORLD).row(0, named=True)
    assert got == {k: row_[k] for k in ("outside", "outside_removed", "native", "native_lost", "outside_rate", "native_loss")}
    assert got["outside_removed"] == 3 and got["native_lost"] == 1


def test_scoring_names_no_species():
    got = score(COUNTS, Settings(), REFS, WORLD)
    assert all(isinstance(v, int | float) for v in got.values())


def test_the_one_unsealed_call_is_in_the_calibration_script_and_nowhere_else():
    """`DECISIONS.md` lists the unsealed call sites. A new one has to be listed, and this fails until it is."""
    found = {}
    for path in [*(ROOT / "src").rglob("*.py"), *(ROOT / "scripts").glob("*.py")]:
        if path.name == "occurrences.py":
            continue  # where the seal is
        n = len(re.findall(r"unsealed\s*=\s*True", path.read_text()))
        if n:
            found[str(path.relative_to(ROOT))] = n
    assert found == {"scripts/calibrate.py": 1}


# ---- regions

def test_regions_are_read_from_the_level_3_file_without_its_polygons(tmp_path):
    features = [
        {"type": "Feature", "properties": {"LEVEL3_NAM": "Iceland", "LEVEL3_COD": "ICE", "LEVEL2_COD": 10, "LEVEL1_COD": 1}, "geometry": None},
        {"type": "Feature", "properties": {"LEVEL3_NAM": "Madagascar", "LEVEL3_COD": "MDG", "LEVEL2_COD": 29, "LEVEL1_COD": 2}, "geometry": None},
        {"type": "Feature", "properties": {"LEVEL3_NAM": "Alaska", "LEVEL3_COD": "ASK", "LEVEL2_COD": 71, "LEVEL1_COD": 7}, "geometry": None},
    ]
    path = tmp_path / "level3.geojson"
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    regions = tdwg_regions(path)
    assert regions.columns == ["tdwg3", "name", "level2", "level1"]
    assert regions.filter(pl.col("tdwg3") == "MDG").row(0) == ("MDG", "Madagascar", 29, 2)
    assert regions["level1"].dtype == pl.Int64


REGIONS = pl.DataFrame(
    [("ICE", "Iceland", 10, 1), ("NOR", "Norway", 10, 1), ("MDG", "Madagascar", 29, 2), ("ASK", "Alaska", 71, 7), ("HAW", "Hawaii", 61, 6)],
    schema={"tdwg3": pl.String, "name": pl.String, "level2": pl.Int64, "level1": pl.Int64},
    orient="row",
)
PARTS = pl.DataFrame(
    [("1:land", "ICE"), ("2:land", "NOR"), ("2:land", "ICE"), ("3:land", "MDG"), ("4:land", "ASK"), ("5:land", "HAW"), ("6:land", None)],
    schema={"part": pl.String, "tdwg3": pl.String},
    orient="row",
)


def test_an_area_is_the_cells_whose_land_touches_a_region_of_a_level_or_a_code():
    assert area_cells(PARTS, REGIONS, "level1:1") == {"1", "2"}
    assert area_cells(PARTS, REGIONS, "tdwg3:MDG") == {"3"}
    assert area_cells(PARTS, REGIONS, "level1:7") == {"4"}
    with pytest.raises(ValueError, match="area"):
        area_cells(PARTS, REGIONS, "continent:1")
    with pytest.raises(ValueError, match="no region"):
        area_cells(PARTS, REGIONS, "tdwg3:XXX")


# ---- native loss per region

def test_region_loss_counts_native_presences_per_region_the_species_is_native_in():
    wcvp = pl.DataFrame(
        [("a", "ICE", True), ("a", "NOR", True), ("b", "NOR", True), ("b", "ICE", False)],
        schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean},
        orient="row",
    )
    parts_ = pl.DataFrame([("1:land", "ICE"), ("2:land", "NOR"), ("3:land", "NOR"), ("3:land", "ICE")], schema={"part": pl.String, "tdwg3": pl.String}, orient="row")
    status = status_table(
        ("a", "1:land", "native"), ("a", "2:land", "native"), ("a", "3:land", "native"), ("b", "2:land", "native"), ("b", "1:land", "outside"), ("b", "3:land", "native")
    )
    removed = pairs(("a", "1:land"), ("a", "3:land"), ("b", "1:land"))
    out = region_loss(status, removed, wcvp, parts_, REGIONS).sort("tdwg3")
    assert out.columns == ["tdwg3", "name", "level1", "parts", "native", "lost", "loss"]
    got = {r["tdwg3"]: r for r in out.iter_rows(named=True)}
    # ICE: a is native there at parts 1 and 3 (both gone); b is introduced there, and its outside part does not count
    assert (got["ICE"]["native"], got["ICE"]["lost"], got["ICE"]["loss"]) == (2, 2, 1.0)
    # NOR: a at 2 and 3 (3 gone), b at 2 and 3 (none gone)
    assert (got["NOR"]["native"], got["NOR"]["lost"], got["NOR"]["loss"]) == (4, 1, 0.25)
    assert (got["ICE"]["name"], got["ICE"]["level1"], got["ICE"]["parts"]) == ("Iceland", 1, 2)
    assert got["NOR"]["parts"] == 2
    assert out.sort("loss", descending=True)["tdwg3"].to_list() == ["ICE", "NOR"]  # the order of the result is by loss
    assert region_loss(status, removed, wcvp, parts_, REGIONS)["tdwg3"].to_list() == ["ICE", "NOR"]
    assert top_regions(region_loss(status, removed, wcvp, parts_, REGIONS), min_native=3, n=5)["tdwg3"].to_list() == ["NOR"]  # ICE has 2 native presences
    assert top_regions(region_loss(status, removed, wcvp, parts_, REGIONS), min_native=1, n=1)["tdwg3"].to_list() == ["ICE"]


# ---- the known cases

def pres(*cells, species="s") -> pl.DataFrame:
    return pairs(*[(species, f"{c}:land") for c in cells])


def case(kind, **kw) -> Case:
    return Case("c", "a claim", kind, name="Sp", kingdom="Animalia", **kw)


def lost_rows(*rows_) -> pl.DataFrame:
    return pl.DataFrame(list(rows_), schema={"rule": pl.String, "species_id": pl.String, "part": pl.String, "presence_lost": pl.Boolean}, orient="row")


NONE = lost_rows()


def test_an_absent_species_has_no_presence_in_the_area_and_keeps_some_elsewhere():
    c = case("absent", area="level1:1")
    ok = evaluate_case(c, ["s"], pres("1", "2", "3", "4"), pres("3", "4"), NONE, {"1", "2"})
    assert (ok["area_before"], ok["area_after"], ok["presences_before"], ok["presences_after"], ok["passed"]) == (2, 0, 4, 2, True)
    bad = evaluate_case(c, ["s"], pres("1", "2", "3", "4"), pres("2", "3", "4"), NONE, {"1", "2"})
    assert (bad["area_after"], bad["passed"]) == (1, False)
    gone_all = evaluate_case(c, ["s"], pres("1", "2"), pres(), NONE, {"1", "2"})
    assert gone_all["passed"] is False  # a species that vanished is not a species kept out of Europe


def test_an_absent_species_reports_what_removed_it_in_the_area():
    c = case("absent", area="level1:7")
    removed = lost_rows(("griis", "s", "4:land", True), ("griis", "s", "4:land", True), ("griis_subnational", "s", "5:land", True), ("griis", "s", "9:land", True), ("under_reported", "s", "4:land", False), ("griis", "t", "4:land", True))
    out = evaluate_case(c, ["s"], pres("4", "5", "9"), pres("9"), removed, {"4", "5"})
    assert out["removed_by_rule"] == {"griis": 1, "griis_subnational": 1}  # in the area, once per presence, for the case's species only
    assert out["passed"] is True


def test_a_case_says_where_the_presences_that_break_it_are():
    only = case("only_in", area="tdwg3:MDG")
    out = evaluate_case(only, ["s"], pres("3", "1", "2", "6", "9"), pres("3", "1", "2", "6", "9"), NONE, {"3"}, None, PARTS, REGIONS)
    # 1:land touches Iceland, 2:land Norway and Iceland, 6:land and 9:land no region at all
    assert out["leftover_regions"] == [{"region": "(no region)", "presences": 2}, {"region": "Iceland", "presences": 2}, {"region": "Norway", "presences": 1}]
    absent = case("absent", area="level1:7")
    out = evaluate_case(absent, ["s"], pres("4", "5", "9"), pres("4", "9"), NONE, {"4", "5"}, None, PARTS, REGIONS)
    assert out["leftover_regions"] == [{"region": "Alaska", "presences": 1}]  # only the presences in the area
    assert evaluate_case(absent, ["s"], pres("4", "9"), pres("9"), NONE, {"4"}, None, PARTS, REGIONS)["leftover_regions"] == []
    wcvp = pl.DataFrame([("s", "ICE", True)], schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean}, orient="row")
    out = evaluate_case(case("native_only"), ["s"], pres("1", "3", "4"), pres("1", "3", "4"), NONE, set(), wcvp, PARTS, REGIONS)
    assert out["leftover_regions"] == [{"region": "Alaska", "presences": 1}, {"region": "Madagascar", "presences": 1}]  # the outside ones
    assert "leftover_regions" not in evaluate_case(only, ["s"], pres("3"), pres("3"), NONE, {"3"}, None, PARTS)  # no regions given
    kept = evaluate_case(case("kept", area="tdwg3:HAW"), ["s"], pres("5"), pres("5"), NONE, {"5"}, None, PARTS, REGIONS)
    assert kept["leftover_regions"] == []


def test_a_species_kept_loses_at_most_the_native_loss_limit_of_its_presences_in_the_area():
    c = case("kept", area="tdwg3:HAW")
    cells = [str(i) for i in range(20)]
    area = set(cells)
    ok = evaluate_case(c, ["s"], pres(*cells, "x"), pres(*cells[1:], "x"), NONE, area)  # 1 of 20 is 5%
    assert (ok["area_before"], ok["area_after"], ok["passed"]) == (20, 19, True)
    bad = evaluate_case(c, ["s"], pres(*cells, "x"), pres(*cells[2:], "x"), NONE, area)
    assert (bad["area_after"], bad["passed"]) == (18, False)
    assert evaluate_case(c, ["s"], pres("x"), pres("x"), NONE, area)["passed"] is False  # nothing there to keep


def test_an_endemic_has_presences_only_in_its_area_and_some_of_them():
    c = case("only_in", area="tdwg3:MDG")
    ok = evaluate_case(c, ["s"], pres("3", "7", "8"), pres("3"), NONE, {"3"})
    assert (ok["area_after"], ok["outside_before"], ok["outside_after"], ok["passed"]) == (1, 2, 0, True)
    assert evaluate_case(c, ["s"], pres("3", "7"), pres("3", "7"), NONE, {"3"})["passed"] is False
    assert evaluate_case(c, ["s"], pres("3", "7"), pres("7"), NONE, {"3"})["passed"] is False  # nothing left on the island
    # what removed the presences outside is reported, and what removed those inside is not
    removed = lost_rows(("detached_satellite", "s", "7:land", True), ("under_reported", "s", "3:land", True))
    assert evaluate_case(c, ["s"], pres("3", "7"), pres("3"), removed, {"3"})["removed_by_rule"] == {"detached_satellite": 1}


def test_a_species_native_only_has_no_presence_outside_its_native_regions_and_keeps_some_inside():
    wcvp = pl.DataFrame([("s", "ICE", True)], schema={"species_id": pl.String, "tdwg3": pl.String, "native": pl.Boolean}, orient="row")
    c = case("native_only")
    ok = evaluate_case(c, ["s"], pres("1", "3", "4"), pres("1"), NONE, set(), wcvp, PARTS)
    assert (ok["native_before"], ok["outside_before"], ok["native_after"], ok["outside_after"], ok["passed"]) == (1, 2, 1, 0, True)
    assert evaluate_case(c, ["s"], pres("1", "3"), pres("1", "3"), NONE, set(), wcvp, PARTS)["passed"] is False
    assert evaluate_case(c, ["s"], pres("1", "3"), pres("3"), NONE, set(), wcvp, PARTS)["passed"] is False
    assert evaluate_case(c, ["s"], pres("1", "3"), pres(), NONE, set(), wcvp, PARTS)["passed"] is False


def test_a_case_with_several_species_counts_each_species_presence():
    c = Case("lemurs", "stays on the island", "only_in", family="F", kingdom="Animalia", area="tdwg3:MDG")
    before = pairs(("a", "3:land"), ("b", "3:land"), ("b", "7:land"))
    out = evaluate_case(c, ["a", "b"], before, before.filter(pl.col("part") == "3:land"), NONE, {"3"})
    assert (out["species"], out["presences_before"], out["presences_after"], out["outside_before"], out["passed"]) == (2, 3, 2, 1, True)


def test_the_designs_five_cases():
    by = {c.key: c for c in CASES}
    assert list(by) == ["tiger", "potato", "starling", "hawaiian_bird", "lemurs"]
    assert (by["tiger"].name, by["tiger"].kind, by["tiger"].area) == ("Panthera tigris", "absent", "level1:1")
    assert (by["potato"].name, by["potato"].kind) == ("Solanum tuberosum", "native_only")
    assert (by["starling"].name, by["starling"].kind, by["starling"].area) == ("Sturnus vulgaris", "absent", "level1:7")
    assert (by["hawaiian_bird"].name, by["hawaiian_bird"].kind, by["hawaiian_bird"].area) == ("Drepanis coccinea", "kept", "tdwg3:HAW")
    assert (by["lemurs"].family, by["lemurs"].kind, by["lemurs"].area) == ("Lemuridae", "only_in", "tdwg3:MDG")


INDEX = pl.DataFrame(
    [("Panthera tigris", "Animalia", "T1", "accepted"), ("Drepanis coccinea", "Animalia", "H1", "synonym"), ("Vestiaria coccinea", "Animalia", "H1", "accepted"), ("Sturnus vulgaris", "Animalia", "S1", "accepted"), ("Solanum tuberosum", "Plantae", "P1", "accepted")],
    schema={"name": pl.String, "kingdom": pl.String, "species_id": pl.String, "status": pl.String},
    orient="row",
)
CASE_TAXA = pl.DataFrame(
    [("T1", "Felidae", "Animalia"), ("H1", "Fringillidae", "Animalia"), ("S1", "Sturnidae", "Animalia"), ("P1", "Solanaceae", "Plantae"), ("L1", "Lemuridae", "Animalia"), ("L2", "Lemuridae", "Animalia"), ("L3", "Indriidae", "Animalia")],
    schema={"species_id": pl.String, "family": pl.String, "kingdom": pl.String},
    orient="row",
)


def test_a_case_resolves_to_its_species_by_name_through_a_synonym_or_to_a_family():
    by = {c.key: c for c in CASES}
    assert resolve_species(by["hawaiian_bird"], INDEX, CASE_TAXA) == ["H1"]  # COL's accepted name is Vestiaria coccinea
    assert resolve_species(by["lemurs"], INDEX, CASE_TAXA) == ["L1", "L2"]
    assert resolve_species(by["potato"], INDEX, CASE_TAXA) == ["P1"]
    with pytest.raises(ValueError, match="not found"):
        resolve_species(Case("x", "y", "absent", name="Nothing here", kingdom="Animalia", area="level1:1"), INDEX, CASE_TAXA)
    with pytest.raises(ValueError, match="no species"):
        resolve_species(Case("x", "y", "absent", family="Nothingidae", kingdom="Animalia", area="level1:1"), INDEX, CASE_TAXA)


def test_known_cases_run_those_of_half_a_and_name_only_the_family_and_half_of_the_others():
    ids = ["T1", "H1", "S1", "P1", "L1", "L2", "L3"]
    families = ["Felidae", "Fringillidae", "Sturnidae", "Solanaceae", "Lemuridae", "Lemuridae", "Indriidae"]
    kingdoms = ["Plantae" if s == "P1" else "Animalia" for s in ids]
    split = pl.DataFrame(
        {"species_id": ids, "lineage": [f"family:{k}:{f}" for k, f in zip(kingdoms, families, strict=True)], "half": ["B"] + ["A"] * 6}
    )
    taxa = pl.DataFrame({"species_id": ids, "family": families, "kingdom": kingdoms})
    parts_ = [("1:land", "NOR"), ("3:land", "MDG"), ("4:land", "ASK"), ("5:land", "HAW"), ("6:land", "ICE")]
    r = replace(refs(wcvp=[("P1", "ICE", True)], tdwg3=parts_), taxa=taxa)
    a = counts(row("H1", "5:land"), row("S1", "4:land"), row("S1", "1:land"), row("P1", "6:land"), row("L1", "3:land"), row("L2", "3:land"))
    cleaned = clean(a, r, make_grid({str(i): (500.0, 500.0) for i in range(1, 7)}), frozenset({1}))
    out = known_cases(CASES, INDEX, split, a, cleaned, r, REGIONS)
    by = {o["case"]: o for o in out}
    assert list(by) == [c.key for c in CASES]
    # the tiger's family is in half B: nothing is run, and nothing is said but the family and the half
    assert by["tiger"] == {"case": "tiger", "claim": "no tiger presence in Europe", "family": "Felidae", "half": "B", "status": "waiting for step 9"}
    assert all(by[k]["status"] == "run" for k in ("potato", "starling", "hawaiian_bird", "lemurs"))
    assert by["potato"]["passed"] is True
    assert by["lemurs"]["passed"] is True and by["lemurs"]["half"] == "A" and by["lemurs"]["species"] == 2
    assert by["hawaiian_bird"]["passed"] is True
    assert by["starling"]["passed"] is False  # nothing removed it from Alaska: layers 3 and 5 were not run here
    assert by["starling"]["area_before"] == 1 and by["starling"]["area_after"] == 1
    assert by["starling"]["leftover_regions"] == [{"region": "Alaska", "presences": 1}]


def test_a_family_that_straddles_the_halves_is_an_error():
    split = pl.DataFrame({"species_id": ["L1", "L2"], "lineage": ["family:Animalia:Lemuridae"] * 2, "half": ["A", "B"]})
    r = replace(refs(), taxa=CASE_TAXA.filter(pl.col("family") == "Lemuridae"))
    by = {c.key: c for c in CASES}
    with pytest.raises(ValueError, match="halves"):
        known_cases([by["lemurs"]], INDEX, split, counts(), clean(counts(), r, GRID, frozenset()), r, REGIONS)


# ---- the script

def load_script():
    spec = importlib.util.spec_from_file_location("calibrate_script", ROOT / "scripts" / "calibrate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_script_refuses_a_second_scoring_of_half_b(tmp_path, monkeypatch):
    script = load_script()
    done = tmp_path / "res7"
    done.mkdir()
    (done / "held_out.json").write_text("{}")
    monkeypatch.setattr("sys.argv", ["calibrate.py", "--counts", str(tmp_path / "none"), "--resolution", "7", "--results", str(tmp_path), "--score-half-b"])
    with pytest.raises(SystemExit, match="already"):
        script.main()


def test_score_held_out_reads_the_plants_of_half_b_and_no_other_species(tmp_path, monkeypatch):
    """On a made-up download, with every species in half B but one: only the calibration plants of B are asked for, and the score is `score`'s."""
    script = load_script()
    lines = ["\t".join(COLUMNS)]
    for r in COUNTS.iter_rows(named=True):
        cell, _, kind = r["part"].rpartition(":")
        lines.append("\t".join(str(v) for v in (r["species_id"], r["countrycode"], cell, kind, *(r[c] for c in COUNT_COLUMNS))))
    lines.append("\t".join(str(v) for v in ("elsewhere", "AA", "1", "land", *([5] * len(COUNT_COLUMNS)))))
    with zipfile.ZipFile(tmp_path / "d.zip", "w") as archive:
        archive.writestr("d.csv", "\n".join(lines) + "\n")
    out = tmp_path / "out"
    ingest(tmp_path / "d.zip", out)
    split = pl.DataFrame({"species_id": [*PLANTS, "heavy", "bird", "elsewhere"], "lineage": ["f"] * 7, "half": ["B"] * 6 + ["A"]})
    calls = []

    def spy(out_dir, split_, half, **kwargs):
        calls.append((sorted(split_["species_id"]), half, kwargs))
        return load_country_counts(out_dir, split_, half, **kwargs)

    monkeypatch.setattr(script, "load_country_counts", spy)
    got = script.score_held_out(out, split, REFS, WORLD, Settings())
    assert calls == [(["heavy", "p1", "p2", "p3", "p4"], "B", {"unsealed": True})]  # the bird of half B is never loaded
    assert got == score(COUNTS, Settings(), REFS, WORLD)
    with pytest.raises(SealedHalf):
        load_country_counts(out, split, "B")  # and the seal itself is as it was


def test_the_script_has_the_flag_and_the_resolution_and_the_counts():
    text = (ROOT / "scripts" / "calibrate.py").read_text()
    for flag in ("--counts", "--resolution", "--score-half-b"):
        assert flag in text
