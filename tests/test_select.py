import numpy as np
import polars as pl
import pytest

import reference_impl as ref
from realms.cluster import run_infomap
from realms.candidates import DegenerateClustering, candidate_maps
from realms.matrix import cell_taxa
from realms.network import network_from
from realms.score import v_measure
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
    ct = cell_taxa(presence, None)
    result = candidate_maps(run_infomap(network_from(ct)), ct, areas, adjacency)
    stability = bootstrap_stability(ct, result.maps, areas, adjacency, n_boot=3)
    assert stability[6] >= 0.8


def test_stability_low_on_noise():
    # species scattered at random have no realms, so the maps shouldn't reproduce
    rng = np.random.default_rng(1)
    cells = [f"c{i}" for i in range(40)]
    rows = [(str(c), f"s{s}") for s in range(300) for c in rng.choice(cells, size=3, replace=False)]
    presence = pl.DataFrame(rows, schema={"cell": pl.String, "species_id": pl.String}, orient="row")
    areas = {c: 1e6 for c in cells}
    adjacency = {c: set() for c in cells}
    maps = {5: {c: i % 5 for i, c in enumerate(cells)}}
    assert bootstrap_stability(cell_taxa(presence, None), maps, areas, adjacency, n_boot=3)[5] < 0.5


def family_only_world():
    """Cells share no species at all, so only family nodes can tell the two groups apart."""
    p, t, truth = [], [], {}
    for group, family in enumerate(("FA", "FB")):
        for i in range(8):
            cell = f"{family}-{i}"
            truth[cell] = group
            for s in range(4):
                p.append((cell, f"{cell}-s{s}"))
                t.append((f"{cell}-s{s}", f"G{cell}", family))
    return (
        pl.DataFrame(p, schema={"cell": pl.String, "species_id": pl.String}, orient="row"),
        pl.DataFrame(t, schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String}, orient="row"),
        truth,
    )


def test_stability_reruns_the_callers_ranks():
    presence, taxa, truth = family_only_world()
    areas = {c: 1e6 for c in truth}
    adjacency = {c: set() for c in truth}
    with_family = bootstrap_stability(cell_taxa(presence, taxa, ("species", "family")), {2: truth}, areas, adjacency, n_boot=3)
    species_only = bootstrap_stability(cell_taxa(presence, taxa, ("species",)), {2: truth}, areas, adjacency, n_boot=3)
    assert with_family[2] > 0.8 > species_only[2]


def test_stability_reruns_the_callers_size_floor():
    # with a floor above the whole world every map collapses to one realm, which can't reproduce two
    presence, taxa, truth = planted_network(n_realms=2, cells_per_realm=6, endemic_species=20, widespread_species=0)
    areas = {c: 1e6 for c in truth}
    adjacency = {c: set() for c in truth}
    ct = cell_taxa(presence, None)
    kept = bootstrap_stability(ct, {2: truth}, areas, adjacency, n_boot=2)
    collapsed = bootstrap_stability(ct, {2: truth}, areas, adjacency, n_boot=2, floor_km2=1e9)
    assert kept[2] > 0.8 and collapsed[2] == 0.0


def reference_stability(presence, taxa, maps, areas, adjacency, n_boot, ranks, seed=1, trials=3, floor_km2=1_000_000, k_max=60):
    """`bootstrap_stability` as it was: resample the string tables, rebuild the network and the rank table from them."""
    rng = np.random.default_rng(seed)
    scores = {k: [] for k in maps}
    for b in range(n_boot):
        new_presence, new_taxa = ref._resample(presence, taxa, rng)
        partition = ref.run_infomap(ref.build_network(new_presence, new_taxa, ranks=ranks), seed=seed + b, trials=trials)
        compared = ref.rank_presence(new_presence, new_taxa, ranks)
        try:
            units = ref.resolve_pieces(partition.leaf_labels(), adjacency, compared)
            floored, _ = ref.apply_size_floor(units, areas, compared, floor_km2)
            if len(set(floored.values())) < 2:
                raise DegenerateClustering
            candidates = ref.merge_sequence(floored, compared, k_max, 2)
        except DegenerateClustering:
            candidates = {}
        for k, labels in maps.items():
            if candidates:
                scores[k].append(v_measure(labels, candidates[min(candidates, key=lambda available: (abs(available - k), available))]))
            else:
                scores[k].append(0.0)
    return {k: float(np.median(values)) for k, values in scores.items()}


@pytest.mark.parametrize("ranks", [("species",), ("species", "genus", "family")])
def test_bootstrap_stability_matches_reference(ranks):
    presence, taxa, truth = planted_network(n_realms=6, cells_per_realm=10, endemic_species=30, widespread_species=5)
    areas = {c: 1e6 for c in truth}
    adjacency = {c: set() for c in truth}
    ct = cell_taxa(presence, taxa, ranks)
    maps = candidate_maps(run_infomap(network_from(ct)), ct, areas, adjacency, k_min=2).maps
    assert len(maps) > 1
    old = reference_stability(presence, taxa, maps, areas, adjacency, n_boot=3, ranks=ranks)
    assert len(set(old.values())) > 1  # the maps are told apart, so the comparison isn't of a constant
    assert bootstrap_stability(ct, maps, areas, adjacency, n_boot=3) == old


def test_stability_needs_species_to_resample():
    presence, taxa, truth = planted_network(n_realms=2, cells_per_realm=6, endemic_species=20, widespread_species=0)
    with pytest.raises(ValueError, match="species"):
        bootstrap_stability(cell_taxa(presence, taxa, ("genus",)), {2: truth}, {c: 1e6 for c in truth}, {c: set() for c in truth}, n_boot=1)
