import numpy as np
import pytest

from helpers import grid, matrix
from realms.breaks import merge_heights, split_of, split_tests
from realms.candidates import merge_sequence
from realms.matrix import cell_taxa
from realms.synth import planted_network


def chain(n, shared_with_next=0, own=5):
    rows = []
    for u in range(n):
        rows += [(f"c{u}", f"own{u}-{i}") for i in range(own)]
        rows += [(f"c{u}", f"link{u}-{i}") for i in range(shared_with_next)]
        if u:
            rows += [(f"c{u}", f"link{u - 1}-{i}") for i in range(shared_with_next)]
    return matrix(rows), {f"c{u}": u for u in range(n)}


def test_merge_heights_gives_the_maps_of_merge_sequence():
    ct, labels = chain(8, shared_with_next=2)
    table, maps = merge_heights(labels, ct, k_min=3)
    assert maps == merge_sequence(labels, ct, k_max=8, k_min=3)
    assert table["k"].to_list() == [8, 7, 6, 5, 4]


def test_merge_heights_records_the_beta_sim_of_the_pair_merged():
    # units 0 and 1 share 4 of their 5 taxa, unit 2 shares nothing: the first merge is at 1 - 4/5, the second at 1
    rows = [("a", f"s{i}") for i in range(5)] + [("b", f"s{i}") for i in range(1, 6)] + [("c", f"t{i}") for i in range(5)]
    table, maps = merge_heights({"a": 0, "b": 1, "c": 2}, matrix(rows), k_min=1)
    assert table["height"].to_list() == pytest.approx([0.2, 1.0])
    assert (table["keep"][0], table["drop"][0]) == (0, 1)
    assert maps[2] == {"a": 0, "b": 0, "c": 2}


def test_split_of_finds_the_realm_divided():
    coarse = {"a": 0, "b": 0, "c": 1, "d": 1}
    fine = {"a": 0, "b": 0, "c": 1, "d": 2}
    assert split_of(coarse, fine) == (1, ["c", "d"], {"c": 1, "d": 2})


def test_split_of_refuses_maps_that_are_not_one_split_apart():
    with pytest.raises(ValueError, match="one realm split"):
        split_of({"a": 0, "b": 0, "c": 1, "d": 1}, {"a": 0, "b": 1, "c": 2, "d": 3})
    with pytest.raises(ValueError, match="different cells"):
        split_of({"a": 0}, {"b": 0})


def planted(shuffle_space):
    """Two planted realms on a 6 x 5 grid, as blocks or scattered so that place says nothing about the biota."""
    presence, taxa, truth = planted_network(n_realms=2, cells_per_realm=15, endemic_species=30, widespread_species=5)
    positions, adjacency = grid(6, 5)
    ordered = sorted(truth, key=lambda c: (truth[c], int(c.split("-c")[1])))
    if shuffle_space:
        ordered = [ordered[i] for i in np.random.default_rng(4).permutation(len(ordered))]
    name = dict(zip(positions, ordered))
    named = {name[p]: {name[n] for n in adjacency[p]} for p in positions}
    return cell_taxa(presence, taxa), truth, named, {c: 1.0 for c in ordered}


def test_a_planted_split_beats_random_splits():
    ct, truth, adjacency, areas = planted(shuffle_space=False)
    maps = {1: {c: 0 for c in truth}, 2: truth}
    table = split_tests(maps, ct, adjacency, areas, n_null=19)
    assert table["score"].to_list() == ["area_weighted", "equal"]
    for row in table.iter_rows(named=True):
        assert row["k"] == 2 and row["cells"] == 30
        assert row["observed"] > row["null_q95"] and row["p"] == pytest.approx(1 / 20)
        assert row["gain"] == pytest.approx(row["observed"] - row["null_mean"])


def test_a_split_of_scattered_cells_does_not_beat_random_ones_on_average():
    ct, truth, adjacency, areas = planted(shuffle_space=True)
    # an arbitrary split of a map with no signal in place
    left = {c: int(c.split("-c")[1]) % 2 for c in truth}
    maps = {1: {c: 0 for c in truth}, 2: left}
    for row in split_tests(maps, ct, adjacency, areas, n_null=19).iter_rows(named=True):
        assert row["p"] > 0.05


def test_split_tests_is_the_same_for_any_number_of_workers():
    ct, truth, adjacency, areas = planted(shuffle_space=False)
    three = {c: (2 if truth[c] == 1 and int(c.split("-c")[1]) < 5 else truth[c]) for c in truth}
    maps = {1: {c: 0 for c in truth}, 2: truth, 3: three}
    assert split_tests(maps, ct, adjacency, areas, n_null=9).equals(split_tests(maps, ct, adjacency, areas, n_null=9, workers=2))


def test_equal_weight_is_the_mean_of_the_sides_and_area_weight_follows_the_areas():
    # side 0 (cells a, b) holds s0 only, side 1 (cell c) holds s1 and w, which is also in d outside the map: shares 1 and 1/2
    ct = matrix([("a", "s0"), ("b", "s0"), ("c", "s1"), ("c", "w"), ("d", "w")])
    adjacency = {"a": {"b"}, "b": {"a", "c"}, "c": {"b"}}
    maps = {1: {"a": 0, "b": 0, "c": 0}, 2: {"a": 0, "b": 0, "c": 1}}
    table = split_tests(maps, ct, adjacency, {"a": 1.0, "b": 1.0, "c": 6.0}, n_null=3)
    observed = dict(table.select("score", "observed").iter_rows())
    assert observed["equal"] == pytest.approx(0.75)
    assert observed["area_weighted"] == pytest.approx((2 * 1.0 + 6 * 0.5) / 8)
