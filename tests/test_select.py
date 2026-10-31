import polars as pl

from realms.cluster import run_infomap
from realms.candidates import candidate_maps
from realms.network import build_network
from realms.select import bootstrap_stability, select_level
from realms.synth import planted_network


def table(rows):
    return pl.DataFrame(rows, schema={"k": pl.Int64, "stability": pl.Float64, "excess": pl.Float64}, orient="row")


def test_picks_highest_excess_among_eligible():
    selection = select_level(table([(5, 0.9, 0.10), (10, 0.85, 0.30), (20, 0.95, 0.20)]))
    assert selection.k == 10 and not selection.fell_short


def test_ineligible_high_excess_is_not_picked():
    selection = select_level(table([(5, 0.9, 0.10), (30, 0.79, 0.90)]))
    assert selection.k == 5


def test_stability_of_exactly_08_is_eligible():
    assert select_level(table([(5, 0.8, 0.10), (6, 0.5, 0.90)])).k == 5


def test_tie_within_001_takes_coarser():
    assert select_level(table([(8, 0.9, 0.295), (12, 0.9, 0.300), (20, 0.9, 0.2)])).k == 8


def test_difference_above_001_is_no_tie():
    assert select_level(table([(8, 0.9, 0.28), (12, 0.9, 0.30)])).k == 12


def test_none_eligible_takes_most_stable_and_flags_it():
    selection = select_level(table([(5, 0.5, 0.9), (9, 0.7, 0.1), (20, 0.6, 0.5)]))
    assert selection.k == 9 and selection.fell_short


def test_stability_high_on_planted_realms():
    presence, taxa, truth = planted_network(n_realms=6, cells_per_realm=10, endemic_species=30, widespread_species=5)
    areas = {c: 1e6 for c in truth}
    adjacency = {c: set() for c in truth}
    result = candidate_maps(run_infomap(build_network(presence)), presence, areas, adjacency)
    stability = bootstrap_stability(presence, None, result.maps, areas, adjacency, n_boot=3)
    assert stability[6] >= 0.8


def test_stability_low_on_noise():
    # species scattered at random have no realms, so the maps shouldn't reproduce
    import numpy as np

    rng = np.random.default_rng(1)
    cells = [f"c{i}" for i in range(40)]
    rows = [(str(c), f"s{s}") for s in range(300) for c in rng.choice(cells, size=3, replace=False)]
    presence = pl.DataFrame(rows, schema={"cell": pl.String, "species_id": pl.String}, orient="row")
    areas = {c: 1e6 for c in cells}
    adjacency = {c: set() for c in cells}
    maps = {5: {c: i % 5 for i, c in enumerate(cells)}}
    assert bootstrap_stability(presence, None, maps, areas, adjacency, n_boot=3)[5] < 0.5
