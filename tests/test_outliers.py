"""Layer 5, the rare-outlier heuristics, on hand-made worlds.

The world is a chain of 60 cells along the equator, 3 degrees (334 km) apart and neighbours of each other, so the grid's median neighbour
spacing is 334 km and two cells are of one patch when their centres are at most 500 + 334 = 834 km apart. Far from the chain, at 50 degrees
north and south, are the cells where a rare record can sit by itself. All cells are coastal, so that a cell has two parts.
"""

import math
from dataclasses import replace

import polars as pl
import pytest

from clean_helpers import DESIGN_5, counts, gone, kept, make_grid, refs, row, taxon
from realms.clean import LAYERS, REMOVED_SCHEMA, RULES, SWITCHABLE, Cleaned, clean
from realms.outliers import (
    cell_graph,
    core_parts,
    detached_satellites,
    neighbour_spacing,
    suspect_cells,
    thin_long_range,
    under_reported,
)

KM_PER_DEGREE = 6371.0088 * math.pi / 180
L5 = frozenset({5})
ALL = frozenset(LAYERS)
VAGRANT, ZOO, PATCH, PATCH_2, GHOST_FAR, GHOST_NEAR = "100", "101", "102", "103", "104", "105"
FAR = {VAGRANT: (50.0, 0.0), ZOO: (50.0, 40.0), PATCH: (50.0, 80.0), PATCH_2: (50.0, 120.0), GHOST_FAR: (-50.0, 0.0), GHOST_NEAR: (-50.0, 3.0)}


def world(step: float = 3.0, extra: dict | None = None, edges=()):
    """60 cells on the equator `step` degrees apart, in a chain, and the cells of `extra` (name: (lat, lon)), default `FAR`."""
    extra = FAR if extra is None else extra
    cells = {str(i): (500.0, 500.0) for i in range(1, 61)} | {c: (500.0, 500.0) for c in extra}
    centres = {str(i): (0.0, (i - 1) * step) for i in range(1, 61)} | extra
    chain = [(str(i), str(i + 1)) for i in range(1, 60)]
    return make_grid(cells, chain + list(edges), centres)


WORLD = world()


def rows(species, cells: dict, **given):
    """Rows of `species` on the land part of each cell, with that many records. A cell of 0 has no row."""
    return [row(species, f"{cell}:land", n=n, **given) for cell, n in cells.items() if n > 0]


def presence(*rs) -> pl.DataFrame:
    """The input of the heuristics for rows of `row`: counts, with a `records` column."""
    return counts(*rs).with_columns(pl.col("n").alias("records"))


def big_range(species, extra: dict | None = None, **given_extra):
    """A species of 1,000 records over the first 45 cells of the chain (22 each, 32 in the first), and `extra` (cell: records) with `given_extra`."""
    chain = {str(i): 22 for i in range(1, 46)} | {"1": 32}
    out = rows(species, chain)
    for cell, n in (extra or {}).items():
        out += rows(species, {cell: n}, **given_extra)
    return out


# ---- the core

def test_core_parts_are_the_top_cells_holding_90_percent_of_the_records():
    # 60 + 30 hold 90 of the 100, so the cells with 6 and 4 are not core
    core = core_parts(presence(*rows("a", {"1": 60, "2": 30, "3": 6, "4": 4})), DESIGN_5)
    assert {"species_id", "part", "cell"} <= set(core.columns)
    assert set(core["part"]) == {"1:land", "2:land"}
    assert set(core["cell"]) == {"1", "2"}


def test_a_species_with_under_20_records_is_all_core_and_with_20_it_is_not():
    assert set(core_parts(presence(*rows("a", {"1": 10, "2": 8, "3": 1})), DESIGN_5)["part"]) == {"1:land", "2:land", "3:land"}  # 19
    assert set(core_parts(presence(*rows("b", {"1": 17, "2": 2, "3": 1})), DESIGN_5)["part"]) == {"1:land", "2:land"}  # 20: 18 are 90%


def test_among_tied_cells_the_core_is_the_lowest_cell_ids_up_to_the_line():
    # half of the 100 records is reached by two of the three cells of 30: the first two by cell id, which is arbitrary but always the same
    core = core_parts(presence(*rows("a", {"3": 30, "1": 30, "2": 30, "4": 10})), replace(DESIGN_5, core_share=0.5))
    assert set(core["cell"]) == {"1", "2"}


def test_a_long_tail_of_equal_cells_that_crosses_the_line_is_core_only_up_to_the_line():
    # 1,000 records: 17 cells of 50 hold 850, and the line is at 900, so 50 of the 150 singletons are core as well: 67 cells, not all 167
    cells = {f"big{i:02d}": 50 for i in range(17)} | {f"one{i:03d}": 1 for i in range(150)}
    core = core_parts(presence(*rows("a", cells)), DESIGN_5)
    assert core["cell"].n_unique() == 67
    assert {c for c in core["cell"] if c.startswith("big")} == {f"big{i:02d}" for i in range(17)}
    assert {c for c in core["cell"] if c.startswith("one")} == {f"one{i:03d}" for i in range(50)}


def test_a_cell_is_core_with_all_its_parts_and_countries():
    r = [row("a", "1:land", "AA", n=30), row("a", "1:sea", "BB", n=30), row("a", "2:land", "AA", n=30), row("a", "3:land", "AA", n=10)]
    core = core_parts(presence(*r), replace(DESIGN_5, core_share=0.5))
    assert set(core["part"]) == {"1:land", "1:sea"}  # the cell holds 60 of the 100, summed over its parts


def test_core_share_of_one_makes_every_cell_core():
    assert core_parts(presence(*rows("a", {"1": 600, "2": 300, "3": 99, "4": 1})), replace(DESIGN_5, core_share=1.0))["cell"].n_unique() == 4


# ---- the grid

def test_the_neighbour_spacing_is_measured_from_the_grid():
    assert neighbour_spacing(world(step=3.0)) == pytest.approx(3 * KM_PER_DEGREE, rel=1e-9)
    assert neighbour_spacing(world(step=6.0)) == pytest.approx(6 * KM_PER_DEGREE, rel=1e-9)
    # one long edge does not move the median
    chord = make_grid({"1": (1.0, 0.0), "2": (1.0, 0.0), "3": (1.0, 0.0), "4": (1.0, 0.0)}, [("1", "2"), ("2", "3"), ("1", "4")], {"1": (0, 0), "2": (0, 1), "3": (0, 2), "4": (60, 0)})
    assert neighbour_spacing(chord) == pytest.approx(KM_PER_DEGREE, rel=1e-9)
    assert neighbour_spacing(make_grid({"1": (1.0, 0.0)})) == 0.0  # no neighbours, no spacing


def test_the_cell_graph_joins_cells_within_a_distance_and_no_further():
    g = make_grid({"a": (1.0, 0.0), "b": (1.0, 0.0), "c": (1.0, 0.0)}, centres={"a": (0, 0), "b": (0, 3), "c": (0, 6)})
    km = 3 * KM_PER_DEGREE
    assert (cell_graph(g, km * 1.001).toarray() != 0).astype(int).tolist() == [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    assert cell_graph(g, km * 0.999).nnz == 0
    assert (cell_graph(g, 2 * km * 1.001).toarray() != 0).astype(int).tolist() == [[0, 1, 1], [1, 0, 1], [1, 1, 0]]


def test_the_cell_graph_measures_over_the_antimeridian_and_the_poles():
    g = make_grid({"a": (1.0, 0.0), "b": (1.0, 0.0), "c": (1.0, 0.0)}, centres={"a": (0, 179), "b": (0, -179), "c": (89, 0)})
    near = (cell_graph(g, 2.1 * KM_PER_DEGREE).toarray() != 0).astype(int).tolist()
    assert near == [[0, 1, 0], [1, 0, 0], [0, 0, 0]]  # 2 degrees across the antimeridian, not 358
    assert cell_graph(g, 91 * KM_PER_DEGREE).toarray()[0, 2] != 0  # 89 north and 0 north are 89 degrees apart


# ---- the brief's checks

def tiger(core_other=2000, here=1, zoo_total=5001, tiger_in_groups=True):
    """A tiger with 400 records in two cells and `here` in the zoo's cell, where its group holds `zoo_total` records."""
    c = counts(*rows("tiger", {"1": 200, "2": 200, ZOO: here}), *rows("other", {"1": core_other, "2": core_other, ZOO: zoo_total - here}))
    groups = [("other", "Mammals")] + ([("tiger", "Mammals")] if tiger_in_groups else [])
    return c, refs(groups=groups)


def test_an_injected_zoo_record_is_removed():
    # 1 record in a cell holding 5,000 records of its group, far from a 400-record core
    c, r = tiger()
    out = clean(c, r, WORLD, L5, DESIGN_5)
    assert gone(out, rule="under_reported") == {("tiger", f"{ZOO}:land", "AA")}
    assert gone(out) == gone(out, rule="under_reported")
    assert set(out.removed["layer"]) == {5}
    assert kept(out) | {("tiger", f"{ZOO}:land"): 1} == kept(clean(c, r, WORLD, frozenset()))
    assert kept(out)[("tiger", "1:land")] == 200 and kept(out)[("other", f"{ZOO}:land")] == 5000


@pytest.mark.parametrize(
    "core_other, here, zoo_total, removed",
    [
        (2000, 1, 5001, True),
        (2000, 5, 6000, True),  # at most 5 records
        (2000, 6, 6000, False),
        (2000, 5, 5000, False),  # a rate of 0.1% is not under 1% of the 9.1% in its core cells
        (0, 1, 999, False),  # a group that is not sampled by 1,000 records is not judged, though the rate would be low enough
        (0, 1, 1000, True),
    ],
)
def test_a_record_is_under_reported_only_if_all_four_conditions_hold(core_other, here, zoo_total, removed):
    c, r = tiger(core_other, here, zoo_total)
    assert (gone(clean(c, r, WORLD, L5, DESIGN_5), rule="under_reported") == {("tiger", f"{ZOO}:land", "AA")}) is removed


def test_a_species_with_no_group_is_not_judged_by_its_rate():
    c, r = tiger(tiger_in_groups=False)
    assert not gone(clean(c, r, WORLD, L5, DESIGN_5))


def test_under_reported_needs_20_records_in_all():
    # 19 and 20 records, the one in the zoo's cell the same
    for core, removed in ((18, False), (19, True)):
        c = counts(*rows("s", {"1": core, ZOO: 1}), *rows("other", {"1": 100, ZOO: 5000}))
        out = clean(c, refs(groups=[("s", "Mammals"), ("other", "Mammals")]), WORLD, L5, DESIGN_5)
        assert (gone(out, rule="under_reported") == {("s", f"{ZOO}:land", "AA")}) is removed, core


def test_a_vagrant_patch_is_removed():
    c = counts(*big_range("s", {VAGRANT: 1}))
    out = clean(c, refs(), WORLD, L5, DESIGN_5)
    assert gone(out, rule="detached_satellite") == {("s", f"{VAGRANT}:land", "AA")}
    assert gone(out) == gone(out, rule="detached_satellite")
    assert len(kept(out)) == 45 and sum(kept(out).values()) == 1000


def test_a_narrow_endemic_with_15_records_in_2_cells_survives():
    # the cells are 5,500 km apart, and 10 + 5 records is under the 20 that the heuristics need
    c = counts(*rows("endemic", {VAGRANT: 10, ZOO: 5}))
    r = refs(groups=[("endemic", "Birds")])
    out = clean(c, r, WORLD, L5, DESIGN_5)
    assert not gone(out) and sum(kept(out).values()) == 15
    # and with the other cell holding group records, as the zoo's does
    out = clean(counts(*rows("endemic", {VAGRANT: 10, ZOO: 5}), *rows("noise", {ZOO: 5000})), refs(groups=[("endemic", "Birds"), ("noise", "Birds")]), WORLD, L5, DESIGN_5)
    assert ("endemic", f"{ZOO}:land") in kept(out)


@pytest.mark.parametrize("datasets, years, survives", [(2, 3, True), (5, 9, True), (1, 3, False), (2, 2, False), (1, 1, False)])
def test_a_disjunct_patch_survives_with_2_datasets_and_3_years(datasets, years, survives):
    c = counts(*big_range("s", {VAGRANT: 5}, n_datasets=datasets, n_years=years))
    out = clean(c, refs(), WORLD, L5, DESIGN_5)
    assert (not gone(out)) is survives
    assert gone(out) == (set() if survives else {("s", f"{VAGRANT}:land", "AA")})


def test_the_datasets_and_years_of_a_patch_are_the_largest_over_its_parts():
    # no part of the two-cell patch has 2 datasets and 3 years, but one has 2 datasets and the other 3 years. A lower bound for the patch, and kept.
    c = counts(*big_range("s"), *rows("s", {GHOST_FAR: 5}, n_datasets=2, n_years=1), *rows("s", {GHOST_NEAR: 5}, n_datasets=1, n_years=3))
    assert not gone(clean(c, refs(), WORLD, L5, DESIGN_5))
    c = counts(*big_range("s"), *rows("s", {GHOST_FAR: 5}, n_datasets=1, n_years=1), *rows("s", {GHOST_NEAR: 5}, n_datasets=1, n_years=1))
    assert gone(clean(c, refs(), WORLD, L5, DESIGN_5)) == {("s", f"{GHOST_FAR}:land", "AA"), ("s", f"{GHOST_NEAR}:land", "AA")}


def test_a_patch_of_two_cells_is_judged_as_one():
    # 2 of 47 cells is 4.3%, and its 10 records are 0.99% of 1,010. 5 and 5 are each far under, but it is the patch that counts
    c = counts(*big_range("s", {GHOST_FAR: 5, GHOST_NEAR: 5}))
    assert {part for _, part, _ in gone(clean(c, refs(), WORLD, L5, DESIGN_5))} == {f"{GHOST_FAR}:land", f"{GHOST_NEAR}:land"}
    c = counts(*big_range("s", {GHOST_FAR: 5, GHOST_NEAR: 6}))  # 11 of 1,011 is 1.09%
    assert not gone(clean(c, refs(), WORLD, L5, DESIGN_5))


@pytest.mark.parametrize("records, removed", [(10, True), (11, False)])
def test_a_patch_of_1_percent_of_the_records_is_not_under_it(records, removed):
    # 10 of 1,010 is 0.99% and 11 of 1,011 is 1.09%, of 1 cell in 46
    out = clean(counts(*big_range("s", {VAGRANT: records})), refs(), WORLD, L5, DESIGN_5)
    assert bool(gone(out)) is removed


def test_a_thin_species_is_taken_out():
    c = counts(
        *rows("far", {"1": 1, VAGRANT: 1}),  # 2 records, 5,600 km apart
        *rows("far4", {"1": 1, "2": 1, ZOO: 1, VAGRANT: 1}),  # 4 records
        *rows("near", {"1": 1, "3": 1}),  # 667 km apart
        *rows("five", {"1": 1, "2": 2, VAGRANT: 2}),  # not under 5 records
        *rows("one", {"1": 3}),
        *rows("cell", {"1": 2, "2": 1}),
    )
    out = clean(c, refs(), WORLD, L5, DESIGN_5)
    assert {s for s, _, _ in gone(out, rule="thin_long_range")} == {"far", "far4"}
    assert gone(out) == gone(out, rule="thin_long_range")
    assert {s for s, _ in kept(out)} == {"near", "five", "one", "cell"}
    # the whole species goes, core cells and all
    assert {(s, p) for s, p, _ in gone(out, rule="thin_long_range")} == {("far", "1:land"), ("far", f"{VAGRANT}:land")} | {("far4", f"{c}:land") for c in ("1", "2", ZOO, VAGRANT)}


def test_a_thin_species_needs_cells_over_1000_km_apart():
    # cells 3, 4 and 5 are 6, 9 and 12 degrees from cell 1 on the equator: 667, 1,000.8 and 1,334 km
    for far, removed in (("3", False), ("4", True), ("5", True)):
        out = clean(counts(*rows("s", {"1": 1, far: 1})), refs(), WORLD, L5, DESIGN_5)
        assert bool(gone(out)) is removed, far
    assert not gone(clean(counts(*rows("s", {"1": 1, "4": 1})), refs(), WORLD, L5, replace(DESIGN_5, thin_km=1001)))


# ---- the guardrails

def test_core_cells_are_never_removed_even_where_they_meet_the_condition():
    vagrant = counts(*big_range("s", {VAGRANT: 1}))
    c, r = tiger()
    for table, reference, rule, cell in ((vagrant, refs(), "detached_satellite", VAGRANT), (c, r, "under_reported", ZOO)):
        assert gone(clean(table, reference, WORLD, L5, DESIGN_5), rule=rule), rule  # not core: removed
        guarded = clean(table, reference, WORLD, L5, replace(DESIGN_5, core_share=1.0))  # every cell core
        assert not gone(guarded), rule
        assert any(part == f"{cell}:land" for _, part in kept(guarded)), rule


def test_a_species_under_20_records_is_never_judged_by_the_record_share_rules():
    c = counts(*rows("s", {"1": 18, ZOO: 1}), *rows("other", {"1": 2000, ZOO: 5000}), *rows("t", {"1": 10, VAGRANT: 9}))  # 19 records each
    out = clean(c, refs(groups=[("s", "Mammals"), ("other", "Mammals"), ("t", "Mammals")]), WORLD, L5, DESIGN_5)
    assert not gone(out)


def test_the_thin_rule_is_the_exception_to_the_core_guardrail():
    out = clean(counts(*rows("s", {"1": 1, VAGRANT: 1})), refs(), WORLD, L5, replace(DESIGN_5, core_share=1.0))
    assert gone(out, rule="thin_long_range") == {("s", "1:land", "AA"), ("s", f"{VAGRANT}:land", "AA")}


# ---- grain

def test_a_cell_is_judged_on_its_records_summed_over_parts_and_countries_and_goes_whole():
    other = rows("other", {"1": 2000, "2": 2000, ZOO: 7000})
    core = rows("tiger", {"1": 200, "2": 200})
    r = refs(groups=[("tiger", "Mammals"), ("other", "Mammals")])
    # 3 + 2 records in two parts and two countries of the zoo's cell: 5 at most, and both rows go
    zoo = [row("tiger", f"{ZOO}:land", "AA", n=3), row("tiger", f"{ZOO}:sea", "BB", n=2)]
    assert gone(clean(counts(*core, *other, *zoo), r, WORLD, L5, DESIGN_5)) == {("tiger", f"{ZOO}:land", "AA"), ("tiger", f"{ZOO}:sea", "BB")}
    # 3 + 3 is 6, and each part alone is under 5: kept
    zoo = [row("tiger", f"{ZOO}:land", "AA", n=3), row("tiger", f"{ZOO}:sea", "BB", n=3)]
    assert not gone(clean(counts(*core, *other, *zoo), r, WORLD, L5, DESIGN_5))


def test_a_patch_is_5_percent_of_the_cells_not_of_the_parts():
    # the vagrant's cell has two parts. 1 cell of 20 is 5%, which is not under, and 1 of 22 is 4.5%, which is, where 2 parts of 23 would be 8.7%
    chain = {str(i): 50 for i in range(1, 20)}
    c = counts(*rows("s", chain), row("s", f"{VAGRANT}:land", n=1), row("s", f"{VAGRANT}:sea", n=1))
    assert not gone(clean(c, refs(), WORLD, L5, DESIGN_5))
    chain = {str(i): 50 for i in range(1, 22)}
    c = counts(*rows("s", chain), row("s", f"{VAGRANT}:land", n=1), row("s", f"{VAGRANT}:sea", n=1))
    assert gone(clean(c, refs(), WORLD, L5, DESIGN_5)) == {("s", f"{VAGRANT}:land", "AA"), ("s", f"{VAGRANT}:sea", "AA")}


# ---- the gap

@pytest.mark.parametrize("step, removed", [(1.0, True), (3.0, False)])
def test_the_gap_is_the_setting_plus_the_median_neighbour_spacing_of_the_grid(step, removed):
    # the same 6 degrees (667 km) between a patch and the range: apart where cells are a degree wide (gap 500 + 111 km), one patch where they are 3 (500 + 334)
    grid = world(step=step, extra={VAGRANT: (6.0, 0.0)})
    out = clean(counts(*big_range("s", {VAGRANT: 1})), refs(), grid, L5, DESIGN_5)
    assert bool(gone(out, rule="detached_satellite")) is removed


def test_the_gap_setting_moves_the_gap():
    c = counts(*big_range("s", {VAGRANT: 1}))
    grid = world(extra={VAGRANT: (6.0, 0.0)})  # 667 km from the range
    assert not gone(clean(c, refs(), grid, L5, replace(DESIGN_5, satellite_gap_km=500)))
    assert gone(clean(c, refs(), grid, L5, replace(DESIGN_5, satellite_gap_km=300)))


# ---- the rules

def everything():
    """A tiger in the zoo, a vagrant and a thin species in one table."""
    c = counts(*tiger()[0].iter_rows(named=True), *big_range("vagrant", {VAGRANT: 1}), *rows("ghost", {"1": 1, GHOST_FAR: 1}))
    return c, refs(groups=[("tiger", "Mammals"), ("other", "Mammals"), ("vagrant", "Birds"), ("ghost", "Birds")])


def test_the_rules_are_of_layer_5():
    for rule in ("detached_satellite", "under_reported", "thin_long_range"):
        assert RULES[rule] == 5 and rule in SWITCHABLE


@pytest.mark.parametrize("rule", ["detached_satellite", "under_reported", "thin_long_range"])
def test_each_rule_can_be_switched_off_on_its_own(rule):
    c, r = everything()
    on = clean(c, r, WORLD, L5, DESIGN_5)
    off = clean(c, r, WORLD, L5, replace(DESIGN_5, rules_off={rule}))
    assert gone(on, rule=rule) and not gone(off, rule=rule)
    for other in {"detached_satellite", "under_reported", "thin_long_range"} - {rule}:
        assert gone(off, rule=other) == gone(on, rule=other)
    assert off.presence["records"].sum() == on.presence["records"].sum() + on.removed.filter(pl.col("rule") == rule)["records"].sum()


def test_a_cell_that_two_rules_flag_goes_with_the_first_and_the_second_takes_it_when_that_is_off():
    # the vagrant cell holds 5,000 records of the vagrant's group: a satellite, and under-reported
    c = counts(*big_range("s", {VAGRANT: 1}), *rows("filler", {VAGRANT: 5000}))
    r = refs(groups=[("s", "Birds"), ("filler", "Birds")])
    both = {("s", f"{VAGRANT}:land", "AA")}
    assert gone(clean(c, r, WORLD, L5, DESIGN_5), rule="detached_satellite") == both
    assert not gone(clean(c, r, WORLD, L5, DESIGN_5), rule="under_reported")
    assert gone(clean(c, r, WORLD, L5, replace(DESIGN_5, rules_off={"detached_satellite"})), rule="under_reported") == both


def test_the_heuristics_are_judged_on_what_the_layers_before_them_leave():
    # the zoo's group has its 5,000 records from a species that layer 2 takes out, so the zoo cell is not well sampled
    c, r = tiger()
    r = refs([taxon("tiger"), taxon("other", extinct=True)], groups=r.groups.rows())
    assert ("tiger", f"{ZOO}:land") in kept(clean(c, r, WORLD, ALL, DESIGN_5))
    assert not gone(clean(c, r, WORLD, ALL, DESIGN_5), rule="under_reported")
    assert gone(clean(c, r, WORLD, ALL - {2}, DESIGN_5), rule="under_reported") == {("tiger", f"{ZOO}:land", "AA")}


def test_layer_5_off_removes_nothing_and_on_it_changes_nothing_where_nothing_is_an_outlier():
    c = counts(*rows("s", {"1": 30, "2": 30, "3": 30}))
    assert clean(c, refs(), WORLD, ALL, DESIGN_5).presence.equals(clean(c, refs(), WORLD, ALL - {5}, DESIGN_5).presence)
    c, r = everything()
    assert not clean(c, r, WORLD, ALL - {5}, DESIGN_5).removed.filter(pl.col("layer") == 5).height


def test_the_heuristics_return_what_they_flag():
    c, r = everything()
    p = c.with_columns(pl.col("n").alias("records"))
    satellite = detached_satellites(p, WORLD, DESIGN_5)
    assert satellite.select("species_id", "cell").rows() == [("vagrant", VAGRANT)]
    assert {"records", "patch_records", "patch_cells", "patch_datasets", "patch_years"} <= set(satellite.columns)
    assert satellite["patch_records"].to_list() == [1] and satellite["patch_cells"].to_list() == [1]
    under = under_reported(p, r.groups, DESIGN_5)
    assert under.select("species_id", "cell").rows() == [("tiger", ZOO)]
    assert {"records", "group_records", "rate", "core_rate"} <= set(under.columns)
    thin = thin_long_range(p, WORLD, DESIGN_5)
    assert sorted(thin.select("species_id", "cell").rows()) == [("ghost", "1"), ("ghost", GHOST_FAR)]
    assert thin["km"].min() > 1000


def test_the_heuristics_refuse_a_presence_without_the_columns_they_need():
    with pytest.raises(ValueError, match="n_years"):
        core_parts(presence(*rows("a", {"1": 5})).drop("n_years"), DESIGN_5)


# ---- suspect cells

def cleaned_of(present: dict, removed: dict) -> Cleaned:
    """A `Cleaned` with `present[cell]` species left in each cell, and `removed[cell]` species that layer 5 took out of it (two rows each, one a sea part)."""
    left = [(f"p{cell}_{i}", f"{cell}:land", 5) for cell, n in present.items() for i in range(n)]
    gone_rows = []
    for cell, n in removed.items():
        for i in range(n):
            for part in ("land", "sea"):
                gone_rows.append((5, "under_reported", f"r{cell}_{i}", f"{cell}:{part}", "AA", 1, None, True))
    gone_rows.append((1, "human_text", "somebody", "1:land", "AA", 9, None, True))  # layers other than 5 are not counted
    pres = pl.DataFrame(left, schema={"species_id": pl.String, "part": pl.String, "records": pl.Int64}, orient="row")
    return Cleaned(pres, pl.DataFrame(gone_rows, schema=REMOVED_SCHEMA, orient="row"))


SUSPECT_EDGES = [("1", "2"), ("1", "3"), ("1", "4"), ("11", "12"), ("11", "13"), ("21", "22"), ("31", "32"), ("41", "42")]


def suspect_grid():
    cells = ["1", "2", "3", "4", "11", "12", "13", "21", "22", "31", "32", "41", "42"]
    return make_grid({c: (1000.0, 0.0) for c in cells}, SUSPECT_EDGES, {c: (float(i), float(i)) for i, c in enumerate(cells)})


def test_a_suspect_cell_has_3_times_its_neighbours_share_and_10_points_more():
    present = {"1": 6, "2": 10, "3": 9, "4": 10, "11": 8, "12": 9, "13": 9, "21": 22, "22": 20, "31": 23, "32": 10, "41": 0, "42": 0}
    removed = {"1": 4, "2": 0, "3": 1, "4": 0, "11": 2, "12": 1, "13": 1, "21": 3, "22": 0, "31": 2, "32": 0, "41": 5, "42": 0}
    cleaned = cleaned_of(present, removed)
    before = (cleaned.presence.clone(), cleaned.removed.clone())
    out = suspect_cells(cleaned, suspect_grid())
    # 1: 40% against a median of 0. 11: 20% against 10% is only twice. 21: 12 points against 0. 31: 8 points. 41: no neighbour has a species to compare with.
    assert out["cell"].to_list() == ["1", "21"]
    first = out.row(0, named=True)
    assert (first["species"], first["removed"], first["share"], first["neighbours"], first["neighbour_median"]) == (10, 4, 0.4, 3, 0.0)
    assert {"lat", "lon"} <= set(out.columns)
    assert out["share"].to_list() == sorted(out["share"].to_list(), reverse=True)
    # a list, nothing removed
    assert cleaned.presence.equals(before[0]) and cleaned.removed.equals(before[1])


def test_with_no_layer_5_removals_there_are_no_suspect_cells():
    out = suspect_cells(cleaned_of({"1": 5, "2": 5}, {}), suspect_grid())
    assert out.height == 0 and {"cell", "share"} <= set(out.columns)


def test_a_zoo_cell_stands_out_in_the_output_of_clean():
    tigers = [f"t{i}" for i in range(12)]
    c = counts(*[r for t in tigers for r in rows(t, {"1": 200, "2": 200, ZOO: 1})], *rows("other", {"1": 2000, "2": 2000, ZOO: 5000}), *rows("bystander", {PATCH: 30}))
    r = refs(groups=[(s, "Mammals") for s in (*tigers, "other")] + [("bystander", "Birds")])
    grid = world(edges=[(ZOO, PATCH)])
    out = clean(c, r, grid, L5, DESIGN_5)
    assert len(gone(out, rule="under_reported")) == 12
    zoo = suspect_cells(out, grid)
    assert zoo["cell"].to_list() == [ZOO]
    assert zoo.row(0, named=True)["removed"] == 12 and zoo.row(0, named=True)["species"] == 13
    assert (zoo["lat"][0], zoo["lon"][0]) == FAR[ZOO]
