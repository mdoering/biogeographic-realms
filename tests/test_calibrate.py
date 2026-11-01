"""The calibration of the layer 5 thresholds on WCVP plants, and the known cases, on hand-made tables.

The world is a chain of 60 cells along the equator, 3 degrees (334 km) apart, so that cells are of one patch when their centres are at most
834 km apart (as in `test_outliers`), and four cells off the chain, where a rare record can sit. Every plant is native in `NAT`, the region of
the whole chain, and `p2` is native in `ISL` too. The cells off the chain are `FAR` (100, 102, 103) and `ISL` (101).

- `p1` has a satellite of 1 record in cell 100, outside its native regions.
- `p2` has one in cell 101, where it is native: a native presence that the heuristic would lose.
- `p3` has one of 8 records in cell 102, outside, 5 of them of one country and 3 of another.
- `p4` has 3 records in cell 103, outside, 400 km from the chain and so joined to its range, in a cell where `heavy` has 2,000 records of the same group.

So, at the design's thresholds, layer 5 takes all four, and at other settings it takes some of them (see `test_the_settings_decide_what_is_taken`).
"""

import ast
import importlib.util
import itertools
import json
import re
import zipfile
from dataclasses import replace
from datetime import datetime
from pathlib import Path

import polars as pl
import pytest

from clean_helpers import GRID, counts, make_grid, refs, row, taxon
from realms.calibrate import (
    CASES,
    GRID as THRESHOLDS,
    EARLY_RULES,
    MAX_MARGINAL,
    NATIVE_LOSS_LIMIT,
    REGION_SCHEMA,
    SCORE,
    Baseline,
    Case,
    area_cells,
    baseline,
    calibrate,
    calibration_counts,
    calibration_species,
    choose,
    evaluate_case,
    frontier,
    griis_cross_check,
    kingdom_counts,
    knee,
    known_cases,
    native_status,
    rates,
    readings,
    region_loss,
    resolve_species,
    score,
    settings_grid,
    tally,
    tdwg_regions,
    top_regions,
    within_limit,
)
from realms.clean import RULES, SWITCHABLE, Settings, clean, outlier_rows
from realms.grid import Grid
from realms.occurrences import COLUMNS, COUNTS as COUNT_COLUMNS, SealedHalf, ingest, load_country_counts, mark_complete

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
    row("p2", "7:land", "BB", n=2),  # a part of the range that has the rows of a second country
    *range_of("p3", {"102": 5}),
    row("p3", "102:land", "BB", n=3),  # the satellite's 8 records are of two countries, so the presence has two rows to take
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


def plantae(*species) -> pl.DataFrame:
    """`species_id` and `kingdom` of `species`, all of them Plantae."""
    return pl.DataFrame({"species_id": list(species), "kingdom": ["Plantae"] * len(species)})


def pairs(*rows_) -> pl.DataFrame:
    return pl.DataFrame(list(rows_), schema={"species_id": pl.String, "part": pl.String}, orient="row")


# ---- the cross-check against GRIIS

def cross_check_world(**given):
    """The world's counts and references with three more animals, `a1` to `a3`: each has a range of 45 cells, `a1` and `a3` a satellite of 1 record
    (a1 two of them, in cells 100 and 102, a3 one in cell 100) and `a2` 3 records in 103 that layer 5 leaves. GRIIS lists `a1` as introduced in AA, the country of all rows."""
    animals = [*range_of("a1", {"100": 1, "102": 1}), *range_of("a2", {"103": 3}), *range_of("a3", {"100": 1})]
    taxa = TAXA + [taxon(s, "Animalia", "Insecta") for s in ("a1", "a2", "a3")]
    groups = GROUPS + [(s, "insects") for s in ("a1", "a2", "a3")]
    return pl.concat([COUNTS, counts(*animals)]), refs(taxa, groups, WCVP, tdwg3=TDWG, griis=[("a1", "AA", None)], **given)


def test_the_cross_check_sets_what_layer_5_takes_of_listed_animals_against_the_other_animals():
    c, r = cross_check_world()
    found = griis_cross_check(c, r, WORLD)
    listed, other = found["listed"], found["other"]
    # a1's 45 cells of range and its two satellites; the other animals are a2, a3 and the bird, and layer 5 takes both satellites of a1 and the one of a3
    assert (listed["pairs"], listed["pairs_removed"], listed["presences"], listed["presences_removed"], listed["species"]) == (47, 2, 47, 2, 1)
    assert (listed["records"], listed["records_removed"]) == (1002, 2)
    assert (other["pairs"], other["pairs_removed"], other["presences"], other["presences_removed"], other["species"]) == (94, 1, 94, 1, 3)
    assert (other["records"], other["records_removed"]) == (3004, 1)
    assert listed["pairs_share"] == pytest.approx(2 / 47) and other["presences_share"] == pytest.approx(1 / 94) and listed["records_share"] == pytest.approx(2 / 1002)
    assert found["animal_presences"] == 141  # the plants are not animals


def test_the_cross_check_is_what_clean_does_with_layer_3_off():
    c, r = cross_check_world()
    full = clean(c, r, WORLD, frozenset({1, 2, 4, 5}))
    taken = full.removed.filter(pl.col("layer") == 5, pl.col("species_id").is_in(["a1", "a2", "a3", "bird"]))
    found = griis_cross_check(c, r, WORLD)
    assert found["listed"]["presences_removed"] + found["other"]["presences_removed"] == taken.select("species_id", "part").unique().height
    # at a looser setting layer 5 takes fewer, and the cross-check follows the setting it is given
    mild = replace(Settings(), satellite_record_share=0.0005)
    assert griis_cross_check(c, r, WORLD, mild)["listed"]["presences_removed"] <= found["listed"]["presences_removed"]


def test_the_cross_check_lists_only_what_a_list_with_distributions_names_as_introduced():
    # a2 is on a list with no distribution table: it is not listed as introduced whatever `griis_without_status` says, and so stays among the others
    c, r = cross_check_world(griis_name_only=[("a2", "AA", None)])
    for setting in (Settings(), Settings(griis_without_status="apply")):
        found = griis_cross_check(c, r, WORLD, setting)
        assert found["listed"]["species"] == 1 and found["other"]["species"] == 3


def test_an_animal_taken_by_the_sub_national_list_is_listed_too():
    c, r = cross_check_world()
    r = replace(r, griis=pl.concat([r.griis, pl.DataFrame({"species_id": ["a3"], "country": ["AA"], "subdivision": ["Sub"]})]), admin1_cells=pl.DataFrame([("100", "AA", "Sub", "overlap")], schema=r.admin1_cells.schema, orient="row"))
    found = griis_cross_check(c, r, WORLD)
    # cell 100 is Sub's. There a3 is listed, and a1 is not, as AA's list leaves Sub's cells to Sub's: a1 keeps its other 46 pairs and a3 has 1
    assert found["listed"]["species"] == 2 and found["listed"]["pairs"] == 47 and found["listed"]["pairs_removed"] == 2


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
    out = native_status(present, wcvp, parts_, plantae("s1", "s2", "s3", "s4"))
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
    assert native_status(present, wcvp, parts_, plantae("s1")).height == 1
    assert native_status(present.head(0), wcvp, parts_, plantae("s1")).columns == ["species_id", "part", "status"]


def test_native_status_agrees_with_layer_3_in_native_only_mode():
    out = clean(COUNTS, REFS, WORLD, frozenset({3}), Settings(wcvp_mode="native_only"))
    taken = {(s, p) for s, p in out.removed.filter(pl.col("rule") == "wcvp_outside_native").select("species_id", "part").iter_rows()}
    status = native_status(COUNTS, REFS.wcvp, REFS.tdwg3_parts, REFS.taxa)
    assert taken == {(s, p) for s, p, v in status.iter_rows() if v == "outside"}
    assert taken == {("p1", "100:land"), ("p3", "102:land"), ("p4", "103:land")}


def test_native_status_counts_plants_only_as_layer_3_does():
    """WCVP lists species that COL puts in no kingdom or another one (5 of the WCVP natives of the release have none). Layer 3 leaves them alone, and so
    `native_status` gives them no status: what the calibration counts as outside is what layer 3 would take."""
    stray = [taxon(s, kingdom, "Magnoliopsida") for s, kingdom in (("p1", None), ("p3", "Fungi"))]
    r = replace(REFS, taxa=pl.DataFrame([t for t in TAXA if t["species_id"] not in ("p1", "p3")] + stray, schema=REFS.taxa.schema))
    out = clean(COUNTS, r, WORLD, frozenset({3}), Settings(wcvp_mode="native_only"))
    taken = {(s, p) for s, p in out.removed.filter(pl.col("rule") == "wcvp_outside_native").select("species_id", "part").iter_rows()}
    status = native_status(COUNTS, r.wcvp, r.tdwg3_parts, r.taxa)
    assert taken == {(s, p) for s, p, v in status.iter_rows() if v == "outside"} == {("p4", "103:land")}
    assert not {"p1", "p3"} & set(status["species_id"])


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

def test_the_grid_has_144_settings_in_a_stable_order():
    grid = settings_grid()
    assert len(grid) == 144 == len(set(grid))
    assert THRESHOLDS == {
        "satellite_record_share": (0.005, 0.01, 0.02),
        "satellite_cell_share": (0.02, 0.05, 0.10, 0.20),
        "under_rate_factor": (0.005, 0.01, 0.015, 0.02),
        "group_min_records": (300, 1000, 3000),
    }
    assert [(s.satellite_record_share, s.satellite_cell_share, s.under_rate_factor, s.group_min_records) for s in grid] == list(itertools.product(*THRESHOLDS.values()))
    assert grid == settings_grid()  # the same order every time: a setting's number means the same thing in every table, report and scoring
    assert grid[0] == Settings(satellite_record_share=0.005, satellite_cell_share=0.02, under_rate_factor=0.005, group_min_records=300)
    assert grid[143] == Settings(satellite_record_share=0.02, satellite_cell_share=0.20, under_rate_factor=0.02, group_min_records=3000)
    assert grid[1].group_min_records == 1000 and grid[3].under_rate_factor == 0.01  # the last varies fastest
    assert grid[12].satellite_cell_share == 0.05 and grid[48].satellite_record_share == 0.01  # and the first slowest: 48 settings to a record share, 12 to a cell share
    assert grid[64] == Settings() and grid.index(Settings()) == 64  # the design's setting is the second value of every threshold: 1 * 48 + 1 * 12 + 1 * 3 + 1


def test_the_grid_keeps_the_other_settings_of_its_base():
    base = Settings(heuristic_min_records=30, wcvp_mode="introduced", land_share=0.9)
    assert all((s.heuristic_min_records, s.wcvp_mode, s.land_share) == (30, "introduced", 0.9) for s in settings_grid(base))


# ---- choosing a setting

def table(*rows_, outside=100, native=1000, run=None) -> pl.DataFrame:
    """Settings with (outside_removed, native_lost) by layer 5's reading, and by the whole run's if `run` gives them, in grid order."""
    run = rows_ if run is None else run
    return pl.DataFrame(
        [
            {
                "setting": i,
                "outside": outside,
                "outside_removed": o,
                "native": native,
                "native_lost": n,
                "run_outside": outside,
                "run_outside_removed": run_o,
                "run_native": native,
                "run_native_lost": run_n,
            }
            for i, ((o, n), (run_o, run_n)) in enumerate(zip(rows_, run, strict=True))
        ]
    )


OLD = {"rule": "most_removed"}  # the rule of before the knee: the most outside presences removed among the settings within the limit


def test_the_old_rule_takes_the_most_outside_removed_among_the_settings_within_the_limit():
    assert NATIVE_LOSS_LIMIT == 0.05
    assert choose(table((10, 5), (60, 49), (90, 51), (50, 0)), **OLD) == (1, True)  # 90 loses 5.1% and is out


def test_the_old_rule_breaks_a_tie_by_lower_native_loss_and_then_by_the_order_of_the_grid():
    assert choose(table((60, 40), (60, 20), (60, 20), (50, 0)), **OLD) == (1, True)
    assert choose(table((60, 20), (60, 20)), **OLD) == (0, True)


def test_a_native_loss_of_exactly_the_limit_is_within_it_whatever_the_float():
    # 3 of 60 is 5% and 0.05 * 60 is 3.0000000000000004, but 0.1 * 3 style noise must not push it out either
    for rule in ("most_removed", "knee"):
        assert choose(table((10, 3), (5, 0), outside=10, native=60), rule=rule, max_marginal=1e6) == (0, True)
        assert choose(table((10, 150), (5, 0), outside=10, native=3000), rule=rule, max_marginal=1e6) == (0, True)
        assert choose(table((10, 151), (5, 0), outside=10, native=3000), rule=rule, max_marginal=1e6) == (1, True)


def test_the_old_rule_reports_when_no_setting_is_within_the_limit_and_takes_the_lowest_native_loss():
    assert choose(table((90, 70), (30, 60), (95, 60), (10, 80)), **OLD) == (1, False)  # 60 twice: the lower index
    assert choose(table((90, 70), (30, 60), (95, 60), (10, 80)), limit=0.07, **OLD) == (2, True)  # 70 is within 7%, and 95 is the most removed


def test_choose_reads_the_whole_run_by_default_and_layer_5_alone_on_request():
    # layer 5 alone would take setting 0 (49 of 1,000 lost), but over the whole run it loses 7% and is out, so setting 1 is chosen
    t = table((90, 49), (60, 10), run=[(50, 70), (40, 20)])
    for rule in ("most_removed", "knee"):
        assert choose(t, rule=rule) == (1, True)
        assert choose(t, reading="run", rule=rule) == (1, True)
        assert choose(t, reading="layer5", rule=rule) == (0, True)
    assert t.filter(within_limit(reading="run"))["setting"].to_list() == [1]
    assert t.filter(within_limit(reading="layer5"))["setting"].to_list() == [0, 1]
    # and where no setting is within the limit of the whole run, the lowest whole-run native loss is taken
    assert choose(table((90, 10), (60, 10), run=[(50, 90), (40, 70)])) == (1, False)
    assert choose(table((90, 10), (60, 10), run=[(50, 90), (40, 70)]), **OLD) == (1, False)
    with pytest.raises(ValueError, match="reading"):
        choose(t, reading="records")
    with pytest.raises(ValueError, match="reading"):
        within_limit(reading="records")
    with pytest.raises(ValueError, match="reading"):
        frontier(t, reading="records")


def test_choose_refuses_an_empty_table_and_an_unknown_rule():
    with pytest.raises(ValueError, match="no settings"):
        choose(table())
    with pytest.raises(ValueError, match="no settings"):
        frontier(table())
    with pytest.raises(ValueError, match="rule"):
        choose(table((10, 5)), rule="max")


# ---- the knee rule

# the settings of this table, as (outside removed, native lost): the frontier is 1, 2 and 5
FRONT = table((50, 10), (40, 5), (60, 10), (45, 20), (60, 30), (80, 40), (30, 5), (40, 5), native=10_000)


def test_frontier_keeps_only_settings_that_remove_more_for_their_loss():
    steps = frontier(FRONT)
    # sorted by native lost, then outside removed (highest first), then the setting: 1 and 7 are alike (the lower index stays), 6 loses as much for less, 0 loses as much as 2 for less,
    # 3 and 4 lose more than 2 and remove no more, and 5 removes more than any before it
    assert steps["setting"].to_list() == [1, 2, 5]
    assert steps["native_lost"].to_list() == [5, 10, 40] and steps["outside_removed"].to_list() == [40, 60, 80]
    assert steps["native"].to_list() == [10_000] * 3
    # a step's marginal cost is the native presences it loses for each outside presence it removes more, and the first row has none
    assert steps["extra_outside"].to_list() == [None, 20, 20] and steps["extra_native"].to_list() == [None, 5, 30]
    assert steps["marginal"].to_list() == [None, 0.25, 1.5]
    # no setting of the table is better than a frontier one on both counts: removes more and loses less
    for row in FRONT.iter_rows(named=True):
        assert any(f["native_lost"] <= row["native_lost"] and f["outside_removed"] >= row["outside_removed"] for f in steps.iter_rows(named=True))
    # the order of the table does not matter, and a reading is the columns it names
    assert frontier(FRONT.reverse())["setting"].to_list() == [1, 2, 5]
    other = table((1, 1), (2, 2), (3, 3), run=[(7, 1), (6, 2), (8, 3)])  # by layer 5 the frontier is all three, over the whole run 1 removes less than 0 for more
    assert frontier(other, reading="layer5")["setting"].to_list() == [0, 1, 2]
    assert frontier(other)["setting"].to_list() == [0, 2]


def test_a_frontier_of_one_setting_is_that_setting_and_has_no_step():
    steps = frontier(table((7, 3)))
    assert steps["setting"].to_list() == [0] and steps["marginal"].to_list() == [None]
    assert choose(table((7, 3))) == (0, True)


def test_knee_rule_stops_at_the_first_step_dearer_than_max_marginal():
    # four frontier settings, 100 more outside removed at each step but the last, which removes 700 more: the native lost by each step is 20, 80 and 1
    t = table((100, 10), (200, 30), (300, 110), (1000, 111), native=10_000)
    assert frontier(t)["marginal"].to_list() == [None, 0.2, 0.8, pytest.approx(1 / 700)]
    assert MAX_MARGINAL == 4.0
    assert choose(t) == (3, True)  # at the default every step is dear enough to pay
    assert choose(t, max_marginal=0.5) == (1, True)  # 0.8 is over: stop there, though the step after it is the cheapest of all
    assert choose(t, max_marginal=0.1) == (0, True)  # even the first step is over
    assert choose(t, max_marginal=0.2) == (1, True)  # a step that costs exactly max_marginal is accepted: at most
    found = knee(t, max_marginal=0.5)
    assert (found.setting, found.met) == (1, True)
    assert found.steps["accepted"].to_list() == [True, True, False, False]  # the walk stopped at the first step it refused
    assert found.steps["affordable"].to_list() == [True, True, False, True]
    # the same at the default: 40 natives for 10 outside is 4.0, which is accepted, and 41 is not
    assert choose(table((100, 10), (110, 50), native=10_000)) == (1, True)
    assert choose(table((100, 10), (110, 51), native=10_000)) == (0, True)
    assert choose(table((100, 10), (110, 50), (120, 91), (500, 92), native=10_000)) == (1, True)


def test_knee_rule_never_passes_the_native_loss_limit():
    # every step is cheap, but the third loses 51 of 1,000 natives, past the 5%, which ends the walk there
    t = table((100, 10), (200, 40), (300, 51), (400, 52), native=1000)
    assert max(frontier(t)["marginal"].drop_nulls()) < 1
    assert choose(t) == (1, True)
    assert choose(t, limit=0.051) == (2, True)  # at 5.1% the third is within it, and the fourth, at 5.2%, is not
    assert choose(t, limit=0.052) == (3, True)  # and at 5.2% all four are, and the last step is cheap too
    assert choose(table((100, 10), (200, 40), (300, 50), native=1000)) == (2, True)  # exactly 5% is within it
    assert knee(t).steps["within_limit"].to_list() == [True, True, False, False]
    assert knee(t).steps["accepted"].to_list() == [True, True, False, False]
    # a cheap step does not rescue a setting past the limit, whatever max_marginal is
    assert choose(t, max_marginal=1e9) == (1, True)
    # where even the setting that loses least is past the limit, it is the one taken, and the caller is told
    assert choose(table((100, 60), (200, 70), native=1000)) == (0, False)
    assert choose(table((100, 60), (200, 70), native=1000), limit=0.07) == (1, True)
    assert knee(table((100, 60), (200, 70), native=1000)).steps["accepted"].to_list() == [False, False]
    # and by the whole run's reading, not layer 5's
    t = table((100, 10), (200, 20), run=[(100, 10), (200, 60)], native=1000)
    assert choose(t, reading="layer5") == (1, True) and choose(t, reading="run") == (0, True)


def test_the_knee_rule_chooses_a_setting_the_old_rule_would_not():
    # the old rule takes the most removed within the limit, however much each extra costs: here 29 natives for 2 outside, a step that the knee rule refuses
    t = table((100, 10), (120, 20), (122, 49), native=1000)
    assert choose(t, **OLD) == (2, True)
    assert choose(t) == (1, True)


FIXTURE = ROOT / "tests" / "fixtures" / "res7_grid81_settings.csv"


def test_knee_rule_on_the_resolution_7_table_of_81_chooses_75():
    # the aggregate counts of half A's dry run at resolution 7 over the 81 settings of the first grid, with the thresholds, and no taxon or cell
    t = pl.read_csv(FIXTURE)
    assert t.height == 81 and t["setting"].to_list() == list(range(81))
    steps = knee(t)
    assert (steps.setting, steps.met) == (75, True) == choose(t)
    chosen = t.filter(pl.col("setting") == 75).row(0, named=True)
    assert (chosen["satellite_record_share"], chosen["satellite_cell_share"], chosen["under_rate_factor"], chosen["group_min_records"]) == (0.02, 0.10, 0.01, 300)
    walk = steps.steps
    ids = walk["setting"].to_list()
    at = ids.index(75)
    assert ids[at + 1] == 8  # the step after it is to setting 8, which the rule refuses
    marginals = walk["marginal"].to_list()
    assert all(m <= 3.4 for m in marginals[1 : at + 1])  # every step up to 75 costs at most 3.4 natives for each outside presence
    assert marginals[at + 1] == pytest.approx(9.1, abs=0.05)  # and the next one 9.1
    assert walk["accepted"].to_list() == [True] * (at + 1) + [False] * (len(ids) - at - 1)
    # the frontier is a few settings of the 81, and the rule chose in it
    assert 1 < walk.height < 81
    # the old rule, which took the most removed within 5%, took setting 7: it removes more outside presences for 5% of the native ones, at about 9 natives each
    assert choose(t, **OLD) == (7, True)
    assert t.filter(within_limit()).height == 58  # 58 of the 81 settings are within the limit of the whole run, as the dry run said


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

def reference(setting: Settings) -> dict:
    """What `clean` itself says, on the plants: the counts of layer 5 alone and of the whole run counting drops only, from its own removal log."""
    plants = calibration_counts(COUNTS, REFS)
    before = clean(plants, REFS, WORLD, BEFORE_5, setting)
    after = clean(plants, REFS, WORLD, LAYERS_WITHOUT_3, setting)
    removed = before.presence.join(after.presence, on=["species_id", "part"], how="anti")
    status = native_status(before.presence, REFS.wcvp, REFS.tdwg3_parts, REFS.taxa)
    raw = plants.filter(pl.col("n_raw") > 0).select("species_id", "part").unique()
    moved = after.removed.filter(pl.col("to_part").is_not_null()).select("species_id", "part").unique()
    lost = raw.join(after.presence.select("species_id", "part"), on=["species_id", "part"], how="anti").join(moved, on=["species_id", "part"], how="anti")
    return {
        **tally(status, removed),
        **{f"run_{k}": v for k, v in tally(native_status(raw, REFS.wcvp, REFS.tdwg3_parts, REFS.taxa), lost).items()},
    }


def test_calibrate_gives_one_row_per_setting_with_the_four_thresholds_and_the_counts():
    grid = settings_grid()[:5]
    out = calibrate(COUNTS, grid, REFS, WORLD)
    assert out.height == 5
    assert out.columns == [
        "setting", "satellite_record_share", "satellite_cell_share", "under_rate_factor", "group_min_records",
        "outside", "outside_removed", "native", "native_lost", "outside_rate", "native_loss",
        "run_outside", "run_outside_removed", "run_native", "run_native_lost", "run_outside_rate", "run_native_loss",
    ]  # fmt: skip
    assert out["setting"].to_list() == [0, 1, 2, 3, 4]
    assert out["group_min_records"].to_list() == [s.group_min_records for s in grid]


def test_calibrate_equals_cleaning_with_every_setting_of_the_grid():
    """Including the rows of two countries for one part of a plant, a satellite's and a range's: `clean` judges rows by country, `calibrate` presences."""
    grid = settings_grid()
    out = calibrate(COUNTS, grid, REFS, WORLD)
    assert out.height == len(grid) == 144
    for setting, got in zip(grid, out.iter_rows(named=True), strict=True):
        counted = reference(setting)
        assert {k: got[k] for k in counted} == counted
        assert got["outside_rate"] == counted["outside_removed"] / counted["outside"]
        assert got["native_loss"] == counted["native_lost"] / counted["native"]
        assert got["run_outside_rate"] == counted["run_outside_removed"] / counted["run_outside"]
        assert got["run_native_loss"] == counted["run_native_lost"] / counted["run_native"]


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
    # no layer takes or folds anything else here, so the whole run (over the raw presences) says the same
    assert {k: design[f"run_{k}"] for k in ("outside", "outside_removed", "native", "native_lost")} == {k: design[k] for k in ("outside", "outside_removed", "native", "native_lost")}


def test_calibrate_can_take_the_baseline_a_caller_made():
    start = baseline(COUNTS, REFS, WORLD)
    assert isinstance(start, Baseline) and start.status.columns == ["species_id", "part", "status"]
    grid = settings_grid()[::20]
    assert calibrate(COUNTS, grid, REFS, WORLD, start).equals(calibrate(COUNTS, grid, REFS, WORLD))


def test_a_baseline_records_its_settings_and_calibrate_refuses_one_that_does_not_match_the_grid():
    start = baseline(COUNTS, REFS, WORLD, Settings(land_share=0.9))
    assert start.settings == Settings(land_share=0.9)
    with pytest.raises(ValueError, match="baseline"):
        calibrate(COUNTS, [Settings()], REFS, WORLD, start)
    with pytest.raises(ValueError, match="baseline"):
        calibrate(COUNTS, [Settings()], REFS, WORLD, baseline(COUNTS, REFS, WORLD, Settings(rules_off={"extinct"})))
    # what layer 5 reads is not in the baseline, so a different threshold or a layer 5 rule off does not matter
    other = baseline(COUNTS, REFS, WORLD, Settings(satellite_cell_share=0.1, rules_off={"thin_long_range"}))
    assert calibrate(COUNTS, [Settings()], REFS, WORLD, other).equals(calibrate(COUNTS, [Settings()], REFS, WORLD))


def test_the_rules_layers_1_2_and_4_read_come_from_rules():
    assert EARLY_RULES == {"extinct", "no_wild_range"}
    assert EARLY_RULES == {rule for rule in SWITCHABLE if RULES[rule] in (1, 2, 4)}
    assert {RULES[rule] for rule in SWITCHABLE} <= {2, 3, 5}  # so that no switchable rule of another layer goes unchecked
    calibrate(COUNTS, [Settings(), Settings(rules_off={"thin_long_range", "under_reported"})], REFS, WORLD)  # layer 5 rules may differ
    for rule in EARLY_RULES:
        with pytest.raises(ValueError, match="rules_off"):
            calibrate(COUNTS, [Settings(), Settings(rules_off={rule})], REFS, WORLD)


KELP_COUNTS = counts(*COUNTS.iter_rows(named=True), *land("kelp", {"1": 30, "2": 30}), row("kelp", "3:land", n=0, n_raw=5))
KELP_REFS = refs(
    [*TAXA, taxon("kelp", "Plantae", "Phaeophyceae", "marine")], [*GROUPS, ("kelp", "angiosperms")], [*WCVP, ("kelp", "NAT", True)], tdwg3=TDWG
)


def test_in_the_drops_only_reading_a_presence_a_fold_moved_is_moved_and_not_lost():
    """`kelp` is a marine species native in NAT, with records on land parts. Layer 4 folds its land parts 1 and 2 into the sea parts of their cells, and layer 1
    drops its land part 3, whose records are all flagged. Both ways the land parts are gone, and only the dropped one is lost."""
    start = baseline(KELP_COUNTS, KELP_REFS, WORLD)
    assert set(start.moved.iter_rows()) == {("kelp", "1:land"), ("kelp", "2:land")}
    assert set(start.absent.iter_rows()) == {("kelp", "1:land"), ("kelp", "2:land"), ("kelp", "3:land")}
    got = readings(start, KELP_REFS, WORLD, Settings())
    assert set(got) == {"layer5", "run_all", "run_drops"}
    # layer 5 judges what layers 1, 2 and 4 leave, and the kelp has no land presence then: the 181 native presences of the plants, one of them lost
    assert (got["layer5"]["native"], got["layer5"]["native_lost"], got["layer5"]["outside"], got["layer5"]["outside_removed"]) == (181, 1, 3, 3)
    # the raw presences hold the kelp's 3 land parts as well: 184 native ones
    assert got["run_all"]["native"] == got["run_drops"]["native"] == 184
    assert got["run_all"]["native_lost"] == 4  # p2's satellite, the part layer 1 dropped, and the two the fold moved
    assert got["run_drops"]["native_lost"] == 2  # p2's satellite and the dropped part: the moved ones are moved, not lost
    assert got["run_drops"]["native_loss"] == pytest.approx(2 / 184) and got["run_all"]["native_loss"] == pytest.approx(4 / 184)
    assert got["run_all"]["outside_removed"] == got["run_drops"]["outside_removed"] == 3
    assert got["run_drops"]["outside_rate"] == 1.0
    out = calibrate(KELP_COUNTS, [Settings()], KELP_REFS, WORLD).row(0, named=True)
    assert (out["run_native"], out["run_native_lost"], out["run_outside"], out["run_outside_removed"]) == (184, 2, 3, 3)
    assert (out["native"], out["native_lost"]) == (181, 1)


def test_a_part_that_only_received_a_fold_is_no_raw_presence_and_is_not_counted():
    """The sea parts the kelp's rows moved into are in what layers 1, 2 and 4 leave, but they are not raw presences, so no reading counts them."""
    start = baseline(KELP_COUNTS, KELP_REFS, WORLD)
    assert ("kelp", "1:sea") not in set(start.raw_status.iter_rows()) and ("kelp", "1:sea") not in set(start.absent.iter_rows())
    assert ("kelp", "1:sea") in set(start.presence.select("species_id", "part").iter_rows())


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
    assert rows_.height == 4  # a presence each
    assert logged.height == 5  # and `clean` logs the rows: the satellite of two countries has two, which go with the one presence


# ---- scoring one setting (the hold-out)

def test_score_is_the_calibration_of_one_setting_as_a_dict():
    setting = Settings()
    got = score(COUNTS, setting, REFS, WORLD)
    row_ = calibrate(COUNTS, [setting], REFS, WORLD).row(0, named=True)
    assert got == {k: row_[k] for k in SCORE}
    assert got["outside_removed"] == 3 and got["native_lost"] == 1 and got["run_native_lost"] == 1


def test_scoring_names_no_species():
    got = score(COUNTS, Settings(), REFS, WORLD)
    assert set(got) == set(SCORE) and all(isinstance(v, int | float) for v in got.values())


def test_kingdom_counts_say_how_many_species_are_not_plants_from_the_taxonomy_alone():
    taxa = pl.DataFrame({"species_id": ["a", "b", "c", "e"], "kingdom": ["Plantae", "Fungi", None, "Plantae"]})
    assert kingdom_counts(["a", "b", "c", "d", "a"], taxa) == {"species": 4, "plantae": 1, "not_plantae": 3}  # no kingdom, and not in taxa, are not Plantae
    assert kingdom_counts([], taxa) == {"species": 0, "plantae": 0, "not_plantae": 0}


SEALED_LOADERS = {"load_presence", "load_country_counts", "_select"}


def _called(func: ast.expr) -> str | None:
    return func.id if isinstance(func, ast.Name) else func.attr if isinstance(func, ast.Attribute) else None


def unsealing_calls(source: str, where: str) -> list[tuple[str, int, str, str]]:
    """The places in `source` (a module at `where`, relative to the repository) that can open half B: a call passing `unsealed=` with anything but the
    constant `False`, `**kwargs` into a sealed loader, a call of `_select` outside `occurrences.py` (it takes `unsealed` by position), a call of `_select`
    there that does not hand on its caller's `unsealed`, and the string `"unsealed"` (a dict of keywords, a `getattr`). (where, line, function, how)."""
    found = []
    inside = where == "src/realms/occurrences.py"
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Constant) and node.value == "unsealed":
            found.append((where, node.lineno, "", "the string"))
        if not isinstance(node, ast.Call):
            continue
        name = _called(node.func)
        for keyword in node.keywords:
            if keyword.arg == "unsealed" and not (isinstance(keyword.value, ast.Constant) and keyword.value.value is False):
                found.append((where, node.lineno, name, "unsealed="))
            if keyword.arg is None and name in SEALED_LOADERS:
                found.append((where, node.lineno, name, "**"))
        if name == "_select":
            if not inside:
                found.append((where, node.lineno, name, "_select"))
            elif not (len(node.args) > 3 and isinstance(node.args[3], ast.Name) and node.args[3].id == "unsealed"):
                found.append((where, node.lineno, name, "_select"))
    return found


def test_the_scan_finds_every_way_to_unseal():
    """The scan itself, on code that unseals in each way it looks for, and on code that does not."""
    caught = [
        'load_presence(d, s, "B", unsealed=True)',
        'load_country_counts(d, s, "B", unsealed=flag)',
        "load_country_counts(**keywords)",
        'occurrences.load_presence(d, s, half, *rest, **{"unsealed": True})',
        'partial(load_presence, unsealed=1)(d, s, "B")',
        '_select(d, s, "B", True, ("species_id",))',
        'getattr(module, "load_presence")(d, s, "B", **dict(unsealed=True))',
    ]
    for code in caught:
        assert unsealing_calls(code, "scripts/x.py"), code
    for code in ('load_presence(d, s, "A")', 'load_country_counts(d, s, "A", unsealed=False)', "print(unsealed)"):
        assert not unsealing_calls(code, "scripts/x.py"), code
    assert unsealing_calls('_select(d, s, "B", True, ())', "src/realms/occurrences.py")
    assert not unsealing_calls('_select(d, s, half, unsealed, ())', "src/realms/occurrences.py")


def test_the_one_unsealed_call_is_in_the_calibration_script_and_nowhere_else():
    """`DECISIONS.md` lists the unsealed call sites. A new one has to be listed, and this fails until it is. Scanned by the syntax tree of every module
    under `src/` and `scripts/`, `occurrences.py` included, so neither a comment nor a line break gets past it."""
    found = []
    for path in [*(ROOT / "src").rglob("*.py"), *(ROOT / "scripts").rglob("*.py")]:
        found += unsealing_calls(path.read_text(), path.relative_to(ROOT).as_posix())
    assert [(where, name, how) for where, _, name, how in found] == [("scripts/calibrate.py", "load_country_counts", "unsealed=")]
    line = found[0][1]
    held_out = [n for n in ast.walk(ast.parse((ROOT / "scripts" / "calibrate.py").read_text())) if isinstance(n, ast.FunctionDef) and n.name == "score_held_out"][0]
    assert held_out.lineno < line <= held_out.end_lineno  # and it is the one of score_held_out


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
    out = evaluate_case(case("native_only"), ["s"], pres("1", "3", "4"), pres("1", "3", "4"), NONE, set(), wcvp, PARTS, REGIONS, plantae("s"))
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
    s = plantae("s")
    ok = evaluate_case(c, ["s"], pres("1", "3", "4"), pres("1"), NONE, set(), wcvp, PARTS, taxa=s)
    assert (ok["native_before"], ok["outside_before"], ok["native_after"], ok["outside_after"], ok["passed"]) == (1, 2, 1, 0, True)
    assert evaluate_case(c, ["s"], pres("1", "3"), pres("1", "3"), NONE, set(), wcvp, PARTS, taxa=s)["passed"] is False
    assert evaluate_case(c, ["s"], pres("1", "3"), pres("3"), NONE, set(), wcvp, PARTS, taxa=s)["passed"] is False
    assert evaluate_case(c, ["s"], pres("1", "3"), pres(), NONE, set(), wcvp, PARTS, taxa=s)["passed"] is False
    with pytest.raises(ValueError, match="taxa"):
        evaluate_case(c, ["s"], pres("1"), pres("1"), NONE, set(), wcvp, PARTS)


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


def made_up_ingest(where: Path, frame: pl.DataFrame, resolution: int, extra: list[str] = ()) -> Path:
    """An ingest of `frame` (country counts) and the `extra` lines of the download, with the marker of a complete ingest at `resolution`."""
    lines = ["\t".join(COLUMNS)]
    for x in frame.iter_rows(named=True):
        cell, _, kind = x["part"].rpartition(":")
        lines.append("\t".join(str(v) for v in (x["species_id"], x["countrycode"], cell, kind, *(x[c] for c in COUNT_COLUMNS))))
    where.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(where / "d.zip", "w") as archive:
        archive.writestr("d.csv", "\n".join([*lines, *extra]) + "\n")
    out = where / "out"
    mark_complete(out, where / "d.zip", resolution, ingest(where / "d.zip", out))
    return out


# The world twice over: the species of half A, and twins of them in half B with the same records, in the same groups, so that the dry run on A and the
# scoring of B each have plants with outside and native presences. One made-up download holds both halves.
TWIN = {s: f"b_{s}" for s in (*PLANTS, "heavy", "bird")}
BOTH_COUNTS = pl.concat([COUNTS, COUNTS.with_columns(pl.col("species_id").replace(TWIN))])
BOTH_REFS = refs(
    [*TAXA, *({**t, "species_id": TWIN[t["species_id"]]} for t in TAXA)],
    [*GROUPS, *((TWIN[s], g) for s, g in GROUPS)],
    [*WCVP, *((TWIN[s], r, n) for s, r, n in WCVP)],
    tdwg3=TDWG,
)
BOTH_SPLIT = pl.DataFrame(
    {"species_id": [*TWIN, *TWIN.values()], "lineage": [f"family:{s}" for s in TWIN] * 2, "half": ["A"] * len(TWIN) + ["B"] * len(TWIN)}
)
WORLD_REGIONS = pl.DataFrame([("NAT", "Native land", 10, 1), ("FAR", "Far away", 20, 2), ("ISL", "An island", 30, 3)], schema=REGION_SCHEMA, orient="row")
REFERENCE_STUBS = ("grids/ISEA3H-7.gpkg", "grids/ISEA3H-6.gpkg", "natural_earth/ne_10m_land.zip", "exports/extended/col.zip", "wgsrpd/level3.geojson")


class ScriptWorld:
    """`scripts/calibrate.py` on the twin world: a made-up ingest at `resolution`, stubs for the reference files, the loaders of the grid, the references,
    the regions, the split and the names replaced by the world's, the known cases left out, every region listed, and the fixed marker in `tmp_path`. Every call of
    `load_country_counts` (the real one, on the made-up ingest) is recorded in `loads` as (half, unsealed, species)."""

    def __init__(self, tmp_path: Path, monkeypatch, resolution: int = 7):
        self.tmp_path, self.monkeypatch = tmp_path, monkeypatch
        self.script = script = load_script()
        self.counts_dir = made_up_ingest(tmp_path / "gbif", BOTH_COUNTS, resolution)
        self.results = tmp_path / "results"
        self.reference = tmp_path / "reference"
        for name in REFERENCE_STUBS:
            (self.reference / name).parent.mkdir(parents=True, exist_ok=True)
            (self.reference / name).touch()
        self.split_path = tmp_path / "split" / "halves-315557.parquet"
        self.split_path.parent.mkdir()
        self.split_path.touch()
        self.seal = tmp_path / "split" / "held_out-315557.started"
        self.loads = []

        def spy(out_dir, split_, half, **kwargs):
            self.loads.append((half, kwargs.get("unsealed", False), sorted(split_.filter(pl.col("half") == half)["species_id"]) if half != "all" else None))
            return load_country_counts(out_dir, split_, half, **kwargs)

        grid = Grid(resolution, WORLD.cells, WORLD.adjacency)
        monkeypatch.setattr(script, "SEAL_MARKER", self.seal)
        monkeypatch.setattr(script, "load_country_counts", spy)
        monkeypatch.setattr(script, "load_grid", lambda gpkg, land: grid)
        monkeypatch.setattr(script, "load_references", lambda grid_, reference: BOTH_REFS)
        monkeypatch.setattr(script, "load_split", lambda path: BOTH_SPLIT)
        monkeypatch.setattr(script, "tdwg_regions", lambda path: WORLD_REGIONS)
        monkeypatch.setattr(script, "name_index", lambda path: pl.DataFrame())
        monkeypatch.setattr(script, "CASES", ())
        monkeypatch.setattr(script, "MIN_NATIVE", 1)  # the world's regions have a few native presences, not the 100 a region needs to be listed

    def run(self, *args):
        argv = ["--counts", str(self.counts_dir), "--results", str(self.results), "--reference", str(self.reference), "--split", str(self.split_path), *args]
        self.monkeypatch.setattr("sys.argv", ["calibrate.py", *argv])
        self.script.main()

    def report(self) -> dict:
        return json.loads((self.results / "res7" / "report.json").read_text())


@pytest.fixture
def world7(tmp_path, monkeypatch):
    return ScriptWorld(tmp_path, monkeypatch)


def test_the_dry_run_loads_half_a_only_and_writes_its_report_with_the_choice_and_the_ingest(world7):
    world7.run()
    assert world7.loads == [("A", False, sorted(TWIN))]  # half A, sealed, and half B never
    report = world7.report()
    assert (report["half"], report["resolution"]) == ("A", 7) and isinstance(report["chosen"]["setting"], int)
    assert report["ingest"] == json.loads((world7.counts_dir / "_COMPLETE").read_text())
    assert not world7.seal.exists() and not list(world7.results.glob("res*/held_out.*"))


def test_the_dry_run_reports_the_knee_walk_and_what_the_old_rule_would_have_chosen(world7):
    world7.run()
    report = world7.report()
    table = pl.read_csv(world7.results / "res7" / "settings.csv")
    assert table.height == 144
    walk = report["knee_rule"]
    assert (walk["rule"], walk["max_marginal"], walk["native_loss_limit"], walk["reading"]) == ("knee", MAX_MARGINAL, NATIVE_LOSS_LIMIT, "run")
    assert [row["setting"] for row in walk["frontier"]] == frontier(table)["setting"].to_list()
    assert report["chosen"]["setting"] == choose(table)[0] == walk["step_into"]["setting"]  # the chosen setting is the last the walk accepted, and the step into it is reported
    ids = [row["setting"] for row in walk["frontier"]]
    after = walk["step_after"]
    if after is None:
        assert ids[-1] == report["chosen"]["setting"] and walk["stopped_by"] is None
    else:  # the step after is the next row of the frontier, with its marginal cost and the test that stopped the walk
        assert after["setting"] == ids[ids.index(report["chosen"]["setting"]) + 1] and after["accepted"] is False
        assert walk["stopped_by"] and set(walk["stopped_by"]) <= {"cost", "limit"}
        assert after["marginal"] == pytest.approx(after["extra_native"] / after["extra_outside"])
    old = report["chosen_by_old_rule"]
    assert old["setting"] == choose(table, rule="most_removed")[0] and set(old["readings"]) == {"layer5", "run_all", "run_drops"}
    assert set(report["chosen"]["readings"]) == {"layer5", "run_all", "run_drops"}


def test_a_dry_run_that_fails_leaves_no_report_of_an_earlier_one(world7, monkeypatch):
    world7.run()

    def fails(*args, **kwargs):
        raise MemoryError("the dry run failed")

    monkeypatch.setattr(world7.script, "calibrate", fails)
    with pytest.raises(MemoryError):
        world7.run()
    assert not (world7.results / "res7" / "report.json").exists()  # so the scoring of half B is refused until a dry run finishes


def test_scoring_half_b_loads_the_references_and_half_bs_calibration_species_and_nothing_else(world7):
    world7.run()
    chosen = world7.report()["chosen"]["setting"]
    setting = (chosen + 1) % 144
    before = {p.name: p.stat().st_mtime_ns for p in (world7.results / "res7").iterdir()}
    world7.loads.clear()

    def not_in_this_run(*args, **kwargs):
        raise AssertionError("the scoring of half B runs no known case, no grid of settings, no region and no cross-check")

    for name in ("known_cases", "name_index", "clean", "baseline", "calibrate", "griis_cross_check", "region_loss", "tdwg_regions"):
        world7.monkeypatch.setattr(world7.script, name, not_in_this_run)
    world7.run("--score-half-b", "--setting", str(setting))
    # one load, of half B's calibration species (the bird of half B is in no plant group), unsealed; nothing of half A
    assert world7.loads == [("B", True, sorted(TWIN[s] for s in (*PLANTS, "heavy")))]
    after = {p.name: p.stat().st_mtime_ns for p in (world7.results / "res7").iterdir()}
    assert {k: v for k, v in after.items() if k in before} == before  # the dry run's files are as they were
    assert set(after) - set(before) <= {"held_out.json", "held_out.started", "held_out.log"} and {"held_out.json", "held_out.started"} <= set(after)
    held = json.loads((world7.results / "res7" / "held_out.json").read_text())
    assert set(held) == world7.script.HELD_OUT_KEYS
    assert (held["setting"], held["dry_run_chosen"], held["setting_is_dry_run_choice"]) == (setting, chosen, False)
    expected = settings_grid()[setting]
    b_counts = BOTH_COUNTS.filter(pl.col("species_id").is_in(list(TWIN.values())))
    assert {k: held[k] for k in SCORE} == score(b_counts, expected, BOTH_REFS, Grid(7, WORLD.cells, WORLD.adjacency))
    assert {k: held[k] for k in THRESHOLDS} == {k: getattr(expected, k) for k in THRESHOLDS}
    assert held["species_loaded"] == {"species": 5, "plantae": 5, "not_plantae": 0}
    # both markers, with the setting and whether it is the dry run's choice
    for marker in (world7.seal, world7.results / "res7" / "held_out.started"):
        lines = dict(line.split(" ", 1) for line in marker.read_text().splitlines())
        assert list(lines) == ["started", "commit", "setting", "dry_run_chosen", "setting_is_dry_run_choice", "counts", "zip_sha256"]
        assert (lines["setting"], lines["dry_run_chosen"], lines["setting_is_dry_run_choice"]) == (str(setting), str(chosen), "False")
        assert re.fullmatch(r"[0-9a-f]{40}|unknown", lines["commit"]) and datetime.fromisoformat(lines["started"])
    # and it is done once
    world7.loads.clear()
    with pytest.raises(SystemExit, match="was started"):
        world7.run("--score-half-b", "--setting", str(chosen))
    assert world7.loads == []


def test_scoring_half_b_records_when_the_setting_is_the_dry_runs_choice(world7):
    world7.run()
    chosen = world7.report()["chosen"]["setting"]
    world7.run("--score-half-b", "--setting", str(chosen))
    held = json.loads((world7.results / "res7" / "held_out.json").read_text())
    assert (held["setting"], held["dry_run_chosen"], held["setting_is_dry_run_choice"]) == (chosen, chosen, True)


def test_scoring_half_b_needs_a_setting_and_a_setting_needs_the_scoring(world7):
    with pytest.raises(SystemExit, match="needs --setting"):
        world7.run("--score-half-b")
    with pytest.raises(SystemExit, match="for --score-half-b only"):
        world7.run("--setting", "3")
    assert world7.loads == [] and not world7.results.exists()


def fake_dry_run(world: ScriptWorld, **changes) -> Path:
    """A `report.json` of a finished dry run at resolution 7 on the world's ingest, with `changes`, and nothing else."""
    report = {"half": "A", "resolution": 7, "chosen": {"setting": 75}, "ingest": json.loads((world.counts_dir / "_COMPLETE").read_text()), **changes}
    path = world.results / "res7" / "report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({k: v for k, v in report.items() if v is not None}))
    return path


@pytest.mark.parametrize(
    "changes, message",
    [
        (None, "no dry run"),
        ({"chosen": None}, "not a complete dry run"),
        ({"chosen": {"setting": "75"}}, "not a complete dry run"),
        ({"resolution": 6}, "not a complete dry run"),
        ({"half": "B"}, "not a complete dry run"),
        ({"ingest": {"zip_sha256": "0" * 64}}, "another ingest"),
    ],
)
def test_scoring_half_b_refuses_without_a_finished_dry_run_at_7_on_the_same_ingest(world7, changes, message):
    if changes is not None:
        fake_dry_run(world7, **changes)
    with pytest.raises(SystemExit, match=message):
        world7.run("--score-half-b", "--setting", "75")
    assert world7.loads == [] and not world7.seal.exists() and not list(world7.results.glob("res*/held_out.*"))


@pytest.mark.parametrize("setting", ["-1", "144"])
def test_scoring_half_b_refuses_a_setting_that_is_not_one_of_the_144(world7, setting):
    fake_dry_run(world7)
    with pytest.raises(SystemExit, match="144 settings"):
        world7.run("--score-half-b", "--setting", setting)
    assert world7.loads == [] and not world7.seal.exists()


@pytest.mark.parametrize("marker, message", [("held_out.json", "already scored"), ("held_out.started", "was started")])
def test_score_half_b_refuses_if_any_resolution_has_a_scoring_started_or_done(world7, marker, message):
    """The marker is written before half B is loaded and the result after, so a run that failed in between is stopped by the marker alone."""
    fake_dry_run(world7)
    for resolution in (6, 7):  # the file of any resolution stops it, and so does the one of its own
        done = world7.results / f"res{resolution}"
        done.mkdir(exist_ok=True)
        (done / marker).write_text("")
        with pytest.raises(SystemExit, match=message) as stopped:
            world7.run("--score-half-b", "--setting", "75")
        assert "Markus" in str(stopped.value) and "no flag" in str(stopped.value)  # and it says how to go on
        assert f"res{resolution}" in str(stopped.value)
        (done / marker).unlink()
    assert world7.loads == []


def test_score_half_b_refuses_on_the_fixed_marker_alone(world7, tmp_path):
    """The fixed marker next to the split holds whatever --results says: a run with results elsewhere, or after results/ was cleared, is still refused."""
    fake_dry_run(world7)
    world7.seal.write_text("started 2026-11-01T00:00:00+00:00\n")
    assert not list(world7.results.glob("res*/held_out.*"))
    with pytest.raises(SystemExit, match="was started") as stopped:
        world7.run("--score-half-b", "--setting", "75")
    assert str(world7.seal) in str(stopped.value) and "Markus" in str(stopped.value)
    assert world7.loads == []


def test_the_fixed_marker_is_next_to_the_frozen_split():
    script = load_script()
    assert script.SEAL_MARKER == ROOT / "data" / "split" / "held_out-315557.started"
    assert script.SEAL_MARKER.parent == (ROOT / "data" / "split" / "halves-315557.parquet").parent


def test_score_half_b_refuses_any_resolution_but_7_before_it_loads_anything(tmp_path, monkeypatch):
    """The resolution is the ingest's, from its marker, so stating 7 for an ingest of another does not get past it."""
    for resolution in (5, 6, 8):
        world = ScriptWorld(tmp_path / str(resolution), monkeypatch, resolution)
        with pytest.raises(SystemExit, match="resolution 7"):
            world.run("--score-half-b", "--setting", "75")
        with pytest.raises(SystemExit, match=f"resolution {resolution}, not 7"):
            world.run("--resolution", "7", "--score-half-b", "--setting", "75")
        assert world.loads == [] and not world.results.exists()  # nothing was made


def failing_load(world: ScriptWorld):
    """`score_held_out` with a loader that fails after half B is asked for, as an out-of-memory error would, and the markers it leaves."""
    markers = [world.seal, world.results / "res7" / "held_out.started"]
    attempts = []

    def fails(out_dir, split_, half, **kwargs):
        attempts.append(half)
        raise MemoryError("the load of half B failed")

    world.monkeypatch.setattr(world.script, "load_country_counts", fails)
    with pytest.raises(MemoryError):
        world.script.score_held_out(world.counts_dir, BOTH_SPLIT, BOTH_REFS, WORLD, Settings(), markers, {"setting": 75})
    return markers, attempts


def test_a_run_that_fails_after_loading_half_b_leaves_both_markers_and_a_retry_is_refused(world7):
    fake_dry_run(world7)
    markers, attempts = failing_load(world7)
    assert attempts == ["B"] and all(m.exists() for m in markers) and not (markers[1].parent / "held_out.json").exists()  # no result, but the markers
    with pytest.raises(SystemExit, match="was started"):
        world7.run("--score-half-b", "--setting", "75")
    markers[1].unlink()  # with the results' marker gone, the fixed one still stops it
    with pytest.raises(SystemExit, match="was started"):
        world7.run("--score-half-b", "--setting", "75")
    with pytest.raises(FileExistsError):  # and the function itself will not write over a marker
        world7.script.mark_started([world7.seal], {})


def test_the_guards_are_for_scoring_half_b_only(tmp_path, monkeypatch):
    """A dry run at resolution 6 is not stopped by a held_out.json or by the fixed marker, and runs."""
    world = ScriptWorld(tmp_path, monkeypatch, resolution=6)
    (world.results / "res7").mkdir(parents=True)
    (world.results / "res7" / "held_out.json").write_text("{}")
    world.seal.write_text("")
    world.run()
    assert (world.results / "res6" / "report.json").exists() and world.loads == [("A", False, sorted(TWIN))]


def test_score_held_out_loads_the_calibration_species_of_half_b_group_mates_included_and_no_other_species(tmp_path, monkeypatch):
    """On a made-up download with every species in half B but one. The calibration species are the plants and the other species of their groups, so
    the alga that shares the plants' group is loaded with them, and the bird of another group is not. The result is numbers: the score, and how many
    species were loaded and how many of those are not Plantae."""
    script = load_script()
    alga = land("alga", {"5": 10})
    everything = counts(*COUNTS.iter_rows(named=True), *alga)
    r = refs([*TAXA, taxon("alga", "Chromista", "Phaeophyceae")], [*GROUPS, ("alga", "angiosperms")], WCVP, tdwg3=TDWG)
    out = made_up_ingest(tmp_path, everything, 7, ["\t".join(str(v) for v in ("elsewhere", "AA", "1", "land", *([5] * len(COUNT_COLUMNS))))])
    split = pl.DataFrame({"species_id": [*PLANTS, "heavy", "alga", "bird", "elsewhere"], "lineage": ["f"] * 8, "half": ["B"] * 7 + ["A"]})
    calls = []
    markers = [tmp_path / "split" / "held_out-315557.started", tmp_path / "res7" / "held_out.started"]

    def spy(out_dir, split_, half, **kwargs):
        assert all(m.exists() for m in markers)  # written before half B is loaded
        calls.append((sorted(split_["species_id"]), half, kwargs))
        return load_country_counts(out_dir, split_, half, **kwargs)

    monkeypatch.setattr(script, "load_country_counts", spy)
    got = script.score_held_out(out, split, r, WORLD, Settings(), markers, {"setting": 40, "dry_run_chosen": 75})
    for marker in markers:
        lines = marker.read_text().splitlines()
        assert [line.split()[0] for line in lines] == ["started", "commit", "setting", "dry_run_chosen"]  # the time, the commit that ran it, and what it scores
        assert re.fullmatch(r"[0-9a-f]{40}|unknown", lines[1].split()[1]) and datetime.fromisoformat(lines[0].split()[1])
    assert calls == [(["alga", "heavy", "p1", "p2", "p3", "p4"], "B", {"unsealed": True})]  # the alga is loaded, the bird of half B never is
    assert got == {**score(everything, Settings(), r, WORLD), "species_loaded": {"species": 6, "plantae": 5, "not_plantae": 1}}
    assert set(got) == {*SCORE, "species_loaded"}
    assert all(isinstance(v, int | float) for k, v in got.items() if k != "species_loaded")
    assert set(got["species_loaded"]) == {"species", "plantae", "not_plantae"}  # no species, no cell: numbers only
    with pytest.raises(SealedHalf):
        load_country_counts(out, split, "B")  # and the seal itself is as it was


def test_the_script_has_the_flag_and_the_resolution_and_the_counts():
    text = (ROOT / "scripts" / "calibrate.py").read_text()
    for flag in ("--counts", "--resolution", "--score-half-b", "--setting"):
        assert flag in text
