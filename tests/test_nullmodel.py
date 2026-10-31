import numpy as np
import pytest

import reference_impl as ref
from helpers import grid
from realms.candidates import connected_pieces
from realms.endemism import area_weighted, realm_share
from realms.matrix import cell_taxa
from realms.nullmodel import endemism_excess, excess_table, null_partitions
from realms.synth import planted_network


def striped_target():
    """A 12 x 12 grid with three realms, one of them in two separate pieces."""
    cells, adjacency = grid(12, 12)
    target = {}
    for c in cells:
        x = int(c.split(",")[0])
        target[c] = 0 if x < 4 else (1 if x < 8 else 2)
    # carve a detached block of realm 0 out of the far side of realm 2
    for x in (10, 11):
        for y in (0, 1, 2):
            target[f"{x},{y}"] = 0
    return cells, adjacency, target, {c: 1.0 for c in cells}


def realm_cells(labels):
    realms = {}
    for c, r in labels.items():
        realms.setdefault(r, set()).add(c)
    return realms


def total_pieces(labels, adjacency):
    return sum(len(connected_pieces(cells, adjacency)) for cells in realm_cells(labels).values())


def test_null_matches_realm_count_and_piece_count():
    _, adjacency, target, areas = striped_target()
    for null in null_partitions(target, adjacency, areas, n=20):
        assert set(null.values()) == set(target.values())
        assert abs(total_pieces(null, adjacency) - total_pieces(target, adjacency)) <= 1


def test_null_areas_within_one_cell_of_target():
    _, adjacency, target, areas = striped_target()
    want = {r: len(cells) for r, cells in realm_cells(target).items()}
    for null in null_partitions(target, adjacency, areas, n=20):
        got = {r: len(cells) for r, cells in realm_cells(null).items()}
        assert all(abs(got[r] - want[r]) <= 2 for r in want)


def test_null_covers_all_cells_exactly_once():
    cells, adjacency, target, areas = striped_target()
    for null in null_partitions(target, adjacency, areas, n=5):
        assert sorted(null) == sorted(cells)


def test_null_differs_from_target_and_between_draws():
    _, adjacency, target, areas = striped_target()
    first, second = null_partitions(target, adjacency, areas, n=2)
    assert first != target and first != second


def test_null_is_deterministic_for_seed():
    _, adjacency, target, areas = striped_target()
    assert null_partitions(target, adjacency, areas, n=3, seed=5) == null_partitions(target, adjacency, areas, n=3, seed=5)
    assert null_partitions(target, adjacency, areas, n=3, seed=5) != null_partitions(target, adjacency, areas, n=3, seed=6)


def test_null_handles_cells_with_no_neighbours():
    cells, adjacency = grid(4, 4)
    cells = cells + ["island"]
    adjacency = adjacency | {"island": set()}
    target = {c: (0 if c == "island" or int(c.split(",")[0]) < 2 else 1) for c in cells}
    areas = {c: 1.0 for c in cells}
    for null in null_partitions(target, adjacency, areas, n=5):
        assert sorted(null) == sorted(cells)


def planted_on_grid(shuffle_space):
    """Three planted realms on a 6 x 5 grid, in spatial blocks, or scattered so that location says nothing about the biota."""
    presence, taxa, truth = planted_network(n_realms=3, cells_per_realm=10, endemic_species=30, widespread_species=5)
    positions, adjacency = grid(6, 5)
    ordered = sorted(truth, key=lambda c: (truth[c], int(c.split("-c")[1])))
    if shuffle_space:
        ordered = [ordered[i] for i in np.random.default_rng(4).permutation(len(ordered))]
    name = dict(zip(positions, ordered))
    named_adjacency = {name[p]: {name[n] for n in adjacency[p]} for p in positions}
    return presence, taxa, truth, named_adjacency, {c: 1.0 for c in ordered}


def test_planted_map_beats_null():
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=False)
    assert endemism_excess(truth, cell_taxa(presence, taxa), adjacency, areas, n_null=19) > 0.3


def test_random_maps_do_not_beat_null():
    # with no signal in location, a map drawn from the null shows no excess on average
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=True)
    ct = cell_taxa(presence, taxa)
    maps = null_partitions(truth, adjacency, areas, n=8, seed=42)
    excess = [endemism_excess(m, ct, adjacency, areas, n_null=19, seed=3) for m in maps]
    assert abs(float(np.mean(excess))) < 0.05


def reference_excess(labels, presence, taxa, adjacency, areas, n_null, seed, tolerance=0.95):
    """`endemism_excess` as it was: one string-keyed table per rank, each turned into its own `incidence`."""
    cells = sorted(labels)
    cell_area = np.array([areas[c] for c in cells])
    realms = sorted(set(labels.values()))
    observed = np.array([realms.index(labels[c]) for c in cells])
    nulls = [np.array([realms.index(null[c]) for c in cells]) for null in null_partitions(labels, adjacency, areas, n=n_null, seed=seed)]
    excess = []
    for table in ref.presence_by_rank(presence, taxa).values():
        incidence = ref.incidence(table, cells)

        def score(cell_realm):
            return area_weighted(realm_share(incidence, cell_realm, len(realms), tolerance), cell_realm, cell_area)

        excess.append(score(observed) - float(np.mean([score(null) for null in nulls])))
    return float(np.mean(excess))


def test_endemism_excess_matches_reference():
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=False)
    for tolerance in (0.95, 1.0):
        old = reference_excess(truth, presence, taxa, adjacency, areas, n_null=9, seed=2, tolerance=tolerance)
        assert endemism_excess(truth, cell_taxa(presence, taxa), adjacency, areas, n_null=9, seed=2, tolerance=tolerance) == old


def test_endemism_excess_averages_the_ranks_of_the_matrix():
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=False)
    ct = cell_taxa(presence, taxa)
    excess = {r: endemism_excess(truth, ct.select([r]), adjacency, areas, n_null=9) for r in ("species", "genus", "family")}
    assert endemism_excess(truth, ct, adjacency, areas, n_null=9) == float(np.mean(list(excess.values())))


def test_excess_table_equals_endemism_excess():
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=False)
    ct = cell_taxa(presence, taxa)
    maps = {3: truth, 2: {c: min(r, 1) for c, r in truth.items()}, 4: null_partitions(truth, adjacency, areas, n=1, seed=8)[0]}
    old = {k: endemism_excess(maps[k], ct, adjacency, areas, n_null=9, seed=2) for k in maps}
    assert len(set(old.values())) == len(maps)
    assert excess_table(maps, ct, adjacency, areas, n_null=9, seed=2, workers=2) == old
    assert excess_table(maps, ct, adjacency, areas, n_null=9, seed=2) == old
    assert list(excess_table(maps, ct, adjacency, areas, n_null=9, seed=2, workers=8)) == list(maps)  # more workers than maps


def test_excess_table_keeps_the_tolerance():
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=False)
    ct = cell_taxa(presence, taxa)
    maps = {4: null_partitions(truth, adjacency, areas, n=1, seed=8)[0]}  # realms that cut across the planted ones
    lax = excess_table(maps, ct, adjacency, areas, n_null=9, tolerance=0.5, workers=2)
    strict = excess_table(maps, ct, adjacency, areas, n_null=9, tolerance=1.0, workers=2)
    assert lax != strict
    assert lax == {4: endemism_excess(maps[4], ct, adjacency, areas, n_null=9, tolerance=0.5)}
    assert strict == {4: endemism_excess(maps[4], ct, adjacency, areas, n_null=9, tolerance=1.0)}


def test_excess_table_needs_every_map_on_the_same_cells():
    presence, taxa, truth, adjacency, areas = planted_on_grid(shuffle_space=False)
    smaller = dict(list(truth.items())[:-1])
    with pytest.raises(ValueError, match="same cells"):
        excess_table({3: truth, 2: smaller}, cell_taxa(presence, taxa), adjacency, areas, n_null=3)
