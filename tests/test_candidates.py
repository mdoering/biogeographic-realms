import numpy as np
import pytest

import reference_impl as ref
from helpers import grid, matrix
from realms.candidates import DegenerateClustering, apply_size_floor, candidate_maps, connected_pieces, merge_sequence, resolve_pieces
from realms.cluster import Partition
from realms.matrix import RANKS, cell_taxa
from realms.synth import planted_network


def units(n, shared_with_next=0, own=5):
    """n units of one cell each, with `own` species of their own and some shared with the next unit."""
    rows = []
    for u in range(n):
        rows += [(f"c{u}", f"own{u}-{i}") for i in range(own)]
        rows += [(f"c{u}", f"link{u}-{i}") for i in range(shared_with_next)]
        if u:
            rows += [(f"c{u}", f"link{u - 1}-{i}") for i in range(shared_with_next)]
    return matrix(rows), {f"c{u}": u for u in range(n)}


def test_merge_sequence_has_one_map_per_k():
    p, labels = units(8, shared_with_next=1)
    maps = merge_sequence(labels, p, k_max=6, k_min=3)
    assert sorted(maps) == [3, 4, 5, 6]
    assert {k: len(set(m.values())) for k, m in maps.items()} == {3: 3, 4: 4, 5: 5, 6: 6}


def test_merge_sequence_starts_at_the_units_available():
    p, labels = units(4, shared_with_next=1)
    assert sorted(merge_sequence(labels, p, k_max=60, k_min=2)) == [2, 3, 4]


def test_merge_sequence_is_nested():
    p, labels = units(8, shared_with_next=1)
    maps = merge_sequence(labels, p, k_max=8, k_min=2)
    for k in range(3, 9):
        for realm in set(maps[k].values()):
            cells = [c for c, r in maps[k].items() if r == realm]
            assert len({maps[k - 1][c] for c in cells}) == 1


def test_merge_takes_the_pair_sharing_most_species():
    # c0 and c1 share 4 of 5 species, every other pair shares nothing
    rows = [("c0", f"s{i}") for i in range(5)] + [("c1", f"s{i}") for i in range(1, 6)]
    rows += [("c2", f"t{i}") for i in range(5)] + [("c3", f"u{i}") for i in range(5)]
    maps = merge_sequence({"c0": 0, "c1": 1, "c2": 2, "c3": 3}, matrix(rows), k_max=4, k_min=3)
    assert maps[3]["c0"] == maps[3]["c1"]
    assert len({maps[3]["c2"], maps[3]["c3"], maps[3]["c0"]}) == 3


def test_merge_uses_replacement_not_richness():
    # c0 is a subset of the much richer c1 (beta-sim 0); c2 and c3 share 3 of 5, which Jaccard would prefer to the subset
    rows = [("c0", f"a{i}") for i in range(2)] + [("c1", f"a{i}") for i in range(40)]
    rows += [("c2", f"b{i}") for i in range(5)] + [("c3", f"b{i}") for i in range(2, 7)]
    maps = merge_sequence({"c0": 0, "c1": 1, "c2": 2, "c3": 3}, matrix(rows), k_max=4, k_min=3)
    assert maps[3]["c0"] == maps[3]["c1"]
    assert maps[3]["c2"] != maps[3]["c3"]


def test_merged_species_sets_are_unions():
    # after c0+c1 merge, c2 is a subset of the union and so joins them, though it overlaps each only partly
    rows = [("c0", "a"), ("c0", "b"), ("c1", "a"), ("c1", "b"), ("c1", "c"), ("c1", "d"), ("c0", "e")]
    rows += [("c2", "c"), ("c2", "e"), ("c3", "x"), ("c3", "y")]
    maps = merge_sequence({"c0": 0, "c1": 1, "c2": 2, "c3": 3}, matrix(rows), k_max=4, k_min=2)
    assert maps[3]["c0"] == maps[3]["c1"]
    assert maps[2]["c2"] == maps[2]["c0"]
    assert maps[2]["c3"] != maps[2]["c0"]


def test_size_floor_merges_small_unit_and_records_it():
    # unit 2 is small and shares species with unit 0 only
    rows = [("big0", "a"), ("big0", "b"), ("big1", "x"), ("big1", "y"), ("small", "a"), ("small", "b")]
    labels = {"big0": 0, "big1": 1, "small": 2}
    areas = {"big0": 2e6, "big1": 2e6, "small": 5e5}
    floored, sub_realms = apply_size_floor(labels, areas, matrix(rows))
    assert floored == {"big0": 0, "big1": 1, "small": 0}
    assert sub_realms == {2: 0}


def test_size_floor_never_leaves_a_realm_below_floor():
    rows = [(f"c{u}", f"s{u}") for u in range(6)] + [(f"c{u}", "shared") for u in range(6)]
    labels = {f"c{u}": u for u in range(6)}
    areas = {f"c{u}": 4e5 for u in range(6)}
    floored, sub_realms = apply_size_floor(labels, areas, matrix(rows))
    for realm in set(floored.values()):
        assert sum(areas[c] for c, r in floored.items() if r == realm) >= 1e6
    assert set(sub_realms.values()) <= set(floored.values())


def test_size_floor_leaves_large_units_alone():
    rows = [("a", "s"), ("b", "t")]
    labels = {"a": 0, "b": 1}
    floored, sub_realms = apply_size_floor(labels, {"a": 2e6, "b": 2e6}, matrix(rows))
    assert floored == labels and sub_realms == {}


def pieces_fixture(shared_species, extra_rows=()):
    # realm 0 has a main piece a0-a2 and a detached piece d0. Realm 1 is b0-b2, between them and touching d0.
    rows = [(c, f"main{i}") for c in ("a0", "a1", "a2") for i in range(5)]
    rows += [(c, f"other{i}") for c in ("b0", "b1", "b2") for i in range(5)]
    rows += [("d0", f"main{i}") for i in range(shared_species)]
    rows += [("d0", "wide"), ("a0", "wide"), ("b0", "wide")]
    rows += extra_rows
    labels = {"a0": 0, "a1": 0, "a2": 0, "b0": 1, "b1": 1, "b2": 1, "d0": 0}
    adjacency = {
        "a0": {"a1"}, "a1": {"a0", "a2"}, "a2": {"a1", "b0"},
        "b0": {"a2", "b1"}, "b1": {"b0", "b2"}, "b2": {"b1", "d0"}, "d0": {"b2"},
    }
    return matrix(rows), labels, adjacency


def test_connected_pieces_largest_first():
    _, _, adjacency = pieces_fixture(0)
    assert connected_pieces({"a0", "a1", "a2", "d0"}, adjacency) == [{"a0", "a1", "a2"}, {"d0"}]


def test_piece_sharing_restricted_species_stays():
    p, labels, adjacency = pieces_fixture(shared_species=3)
    assert resolve_pieces(labels, adjacency, p) == labels


def test_piece_sharing_only_widespread_species_moves():
    # d0 shares only "wide", which isn't restricted to realm 0, so it joins the realm next to it
    p, labels, adjacency = pieces_fixture(shared_species=0)
    assert resolve_pieces(labels, adjacency, p)["d0"] == 1


def test_piece_below_min_shared_moves():
    p, labels, adjacency = pieces_fixture(shared_species=2)
    assert resolve_pieces(labels, adjacency, p, min_shared=3)["d0"] == 1


def test_largest_piece_never_moves():
    p, labels, adjacency = pieces_fixture(shared_species=0)
    resolved = resolve_pieces(labels, adjacency, p)
    assert {resolved[c] for c in ("a0", "a1", "a2")} == {0}


def test_isolated_piece_goes_to_the_most_similar_realm():
    # d0 has no neighbours at all and shares species with realm 1 only
    rows = [("a0", "x"), ("a1", "x"), ("b0", "y"), ("b1", "y"), ("c0", "z"), ("c1", "z"), ("d0", "y")]
    labels = {"a0": 0, "a1": 0, "b0": 1, "b1": 1, "c0": 2, "c1": 2, "d0": 0}
    adjacency = {"a0": {"a1"}, "a1": {"a0"}, "b0": {"b1"}, "b1": {"b0"}, "c0": {"c1"}, "c1": {"c0"}, "d0": set()}
    assert resolve_pieces(labels, adjacency, matrix(rows))["d0"] == 1


def test_candidate_maps_applies_the_steps_in_order():
    # 6 big units in a chain plus one tiny unit: the tiny one becomes a sub-realm, the maps run from 6 down to 5
    rows, paths, areas = [], {}, {}
    adjacency = {}
    names = [f"u{i}" for i in range(6)] + ["tiny"]
    for i, name in enumerate(names):
        rows += [(name, f"s{i}-{j}") for j in range(4)]
        paths[name] = (i + 1,)
        areas[name] = 3e6 if name != "tiny" else 1e5
        adjacency[name] = {names[j] for j in (i - 1, i + 1) if 0 <= j < len(names)}
    rows += [("tiny", "s0-0"), ("tiny", "s0-1")]
    result = candidate_maps(Partition(paths=paths, unassigned=[]), matrix(rows), areas, adjacency)
    assert sorted(result.maps) == [5, 6]
    tiny_label = result.units["tiny"]
    assert tiny_label in result.sub_realms
    assert result.maps[6]["tiny"] == result.maps[6]["u0"]


def test_degenerate_clustering_raises():
    rows = [(f"c{u}", f"s{u}") for u in range(3)]
    paths = {f"c{u}": (u + 1,) for u in range(3)}
    areas = {f"c{u}": 2e6 for u in range(3)}
    adjacency = {f"c{u}": set() for u in range(3)}
    with pytest.raises(DegenerateClustering):
        candidate_maps(Partition(paths=paths, unassigned=[]), matrix(rows), areas, adjacency)


def test_species_also_outside_the_map_is_not_restricted():
    # the shared species occur mostly in cells the map doesn't cover, so they can't hold the detached piece in the realm
    elsewhere = [(f"far{i}", f"main{s}") for i in range(30) for s in range(3)]
    p, labels, adjacency = pieces_fixture(shared_species=3, extra_rows=elsewhere)
    assert resolve_pieces(labels, adjacency, p)["d0"] == 1


def test_piece_sharing_nothing_with_other_realms_stays():
    # d0 has its own species and shares none with realm 1, so moving it there would be arbitrary
    rows = [("a0", "x"), ("a1", "x"), ("b0", "y"), ("b1", "y"), ("d0", "own")]
    labels = {"a0": 0, "a1": 0, "b0": 1, "b1": 1, "d0": 0}
    adjacency = {"a0": {"a1"}, "a1": {"a0"}, "b0": {"b1"}, "b1": {"b0"}, "d0": set()}
    assert resolve_pieces(labels, adjacency, matrix(rows))["d0"] == 0


def planted_case(seed):
    """Planted realms cut into 40 random units, each cell's area drawn so that most units start under the floor."""
    presence, taxa, truth = planted_network(n_realms=4, cells_per_realm=15, endemic_species=60, widespread_species=10, seed=seed)
    rng = np.random.default_rng(seed)
    cells = sorted(truth)
    labels = {c: int(u) for c, u in zip(cells, rng.integers(40, size=len(cells)))}
    areas = {c: float(a) for c, a in zip(cells, rng.uniform(2e5, 6e5, size=len(cells)))}
    return presence, taxa, labels, areas


@pytest.mark.parametrize("seed", range(10))
def test_merge_sequence_matches_reference(seed):
    presence, taxa, labels, _ = planted_case(seed)
    old = ref.merge_sequence(labels, ref.rank_presence(presence, taxa, RANKS), k_max=40, k_min=5)
    new = merge_sequence(labels, cell_taxa(presence, taxa), k_max=40, k_min=5)
    assert len(new) > 10
    assert new == old


@pytest.mark.parametrize("seed", range(10))
def test_size_floor_matches_reference(seed):
    presence, taxa, labels, areas = planted_case(seed)
    old = ref.apply_size_floor(labels, areas, ref.rank_presence(presence, taxa, RANKS))
    new = apply_size_floor(labels, areas, cell_taxa(presence, taxa))
    assert len(new[1]) > 5  # many units are absorbed, so the merges are exercised
    assert new == old


def scattered_case(seed):
    """Planted realms scattered over a 10 x 10 grid, with a fifth of the cells labelled as some other realm.

    A scattered piece of a true realm shares its endemic species with the rest of the realm and stays. A mislabelled cell doesn't,
    and moves.
    """
    presence, taxa, truth = planted_network(n_realms=4, cells_per_realm=25, endemic_species=40, widespread_species=5, seed=seed)
    positions, grid_adjacency = grid(10, 10)
    rng = np.random.default_rng(seed)
    cells = sorted(truth)
    name = dict(zip(positions, (cells[i] for i in rng.permutation(len(cells)))))
    adjacency = {name[p]: {name[n] for n in grid_adjacency[p]} for p in positions}
    labels = {c: (int(rng.integers(4)) if rng.random() < 0.2 else truth[c]) for c in cells}
    return presence, taxa, labels, adjacency


@pytest.mark.parametrize("seed", range(10))
def test_resolve_pieces_matches_reference(seed):
    presence, taxa, labels, adjacency = scattered_case(seed)
    old = ref.resolve_pieces(labels, adjacency, ref.rank_presence(presence, taxa, RANKS))
    new = resolve_pieces(labels, adjacency, cell_taxa(presence, taxa))
    realms = {}
    for cell, realm in labels.items():
        realms.setdefault(realm, set()).add(cell)
    detached = sum(len(piece) for members in realms.values() for piece in connected_pieces(members, adjacency)[1:])
    assert 0 < sum(old[c] != labels[c] for c in labels) < detached  # some pieces move and some stay, or this proves little
    assert new == old
