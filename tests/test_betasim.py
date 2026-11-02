import numpy as np
import polars as pl
import pytest
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform
from sklearn.metrics import v_measure_score

from realms.betasim import Consensus, anosim_r, beta_sim, consensus, hierarchical_labels, rank_beta_sim
from realms.matrix import cell_taxa
from realms.synth import planted_network


def test_beta_sim_nested_is_zero():
    # cell 0's species are a subset of cell 1's: no replacement at all, however different the richness
    assert beta_sim(csr_matrix([[1, 1, 0, 0, 0], [1, 1, 1, 1, 1]])).tolist() == [0.0]


def test_beta_sim_disjoint_is_one():
    assert beta_sim(csr_matrix([[1, 1, 0], [0, 0, 1]])).tolist() == [1.0]


def test_beta_sim_hand_example():
    # shared a=2, only in first b=1, only in second c=3: min(b,c)/(a+min(b,c)) = 1/3
    x = csr_matrix([[1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 1, 1]])
    assert beta_sim(x)[0] == pytest.approx(1 / 3)


def test_beta_sim_is_condensed_in_pair_order():
    x = csr_matrix([[1, 1], [1, 1], [0, 1]])
    # pairs (0,1), (0,2), (1,2)
    assert beta_sim(x).tolist() == [0.0, 0.0, 0.0]
    y = csr_matrix([[1, 0], [1, 0], [0, 1]])
    assert beta_sim(y).tolist() == [0.0, 1.0, 1.0]


def test_ward_separates_two_blocks():
    dist = np.array([0.1, 0.9, 0.9, 0.9, 0.9, 0.1])  # (0,1) close, (2,3) close
    labels = hierarchical_labels(dist, "ward", 2)
    assert labels[0] == labels[1] and labels[2] == labels[3] and labels[0] != labels[2]


def test_average_separates_two_blocks():
    dist = np.array([0.1, 0.9, 0.9, 0.9, 0.9, 0.1])
    labels = hierarchical_labels(dist, "average", 2)
    assert labels[0] == labels[1] and labels[2] == labels[3] and labels[0] != labels[2]


def test_labels_give_exactly_k_groups():
    dist = np.array([0.1, 0.9, 0.9, 0.9, 0.9, 0.1])
    assert len(set(hierarchical_labels(dist, "ward", 3).tolist())) == 3


LABELS = np.array([0, 0, 1, 1])


def test_anosim_perfect_grouping_is_one():
    # pairs (0,1)(0,2)(0,3)(1,2)(1,3)(2,3); the within pairs (0,1) and (2,3) are the two smallest
    assert anosim_r(np.array([1.0, 4.0, 5.0, 3.0, 6.0, 2.0]), LABELS) == pytest.approx(1.0)


def test_anosim_no_structure_is_zero():
    # within ranks 1 and 6, mean 3.5; between ranks 2..5, mean 3.5
    assert anosim_r(np.array([1.0, 3.0, 4.0, 2.0, 5.0, 6.0]), LABELS) == pytest.approx(0.0)


def test_anosim_matches_vegan_hand_example():
    # within ranks 1 and 4, mean 2.5; between ranks 2, 5, 3, 6, mean 4; (4 - 2.5) / (6 / 2)
    assert anosim_r(np.array([1.0, 2.0, 5.0, 3.0, 6.0, 4.0]), LABELS) == pytest.approx(0.5)


def test_anosim_ties_use_average_ranks():
    # all equal: every rank is 3.5, so no difference between within and between
    assert anosim_r(np.ones(6), LABELS) == pytest.approx(0.0)


def ranked(rows, lineage):
    """The cell-taxon matrix of `(cell, species)` rows, where `lineage` gives each species' `(genus, family)`."""
    presence = pl.DataFrame(rows, schema={"cell": pl.String, "species_id": pl.String}, orient="row")
    taxa = pl.DataFrame(
        [(s, g, f) for s, (g, f) in lineage.items()],
        schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String},
        orient="row",
    )
    return cell_taxa(presence, taxa)


LINEAGE = {"s1": ("g1", "f1"), "s2": ("g1", "f1"), "s3": ("g2", "f1"), "s4": ("g3", "f2")}
ROWS = [("A", "s1"), ("A", "s2"), ("B", "s1"), ("B", "s3"), ("C", "s4")]


def test_rank_beta_sim_is_the_mean_of_the_three_ranks():
    # A holds s1 s2 / g1 / f1, B holds s1 s3 / g1 g2 / f1, C holds s4 / g3 / f2
    # A-B: species 1 - 1/2, genus 1 - 1/1 (A's g1 is within B's), family 1 - 1/1; A-C and B-C are disjoint at every rank
    d = rank_beta_sim(ranked(ROWS, LINEAGE), ["A", "B", "C"])
    assert d.dtype == np.float32 and d.shape == (3, 3)
    assert d[0, 1] == pytest.approx((0.5 + 0.0 + 0.0) / 3)
    assert d[0, 2] == d[1, 2] == 1.0
    assert np.array_equal(d, d.T) and not d.diagonal().any()


def test_rank_beta_sim_matches_beta_sim_rank_by_rank():
    presence, taxa, _ = planted_network(n_realms=3, cells_per_realm=8, endemic_species=15, widespread_species=5)
    ct = cell_taxa(presence, taxa)
    d = rank_beta_sim(ct, ct.cells)
    each = [squareform(beta_sim(ct.select([r]).x), checks=False) for r in ("species", "genus", "family")]
    assert d == pytest.approx(sum(each) / 3, abs=1e-6)


def test_rank_beta_sim_follows_the_units_asked_for():
    ct = ranked(ROWS, LINEAGE)
    d = rank_beta_sim(ct, ["A", "B", "C"])
    assert np.array_equal(rank_beta_sim(ct, ["C", "A", "B"]), d[np.ix_([2, 0, 1], [2, 0, 1])])
    assert np.array_equal(rank_beta_sim(ct, ["B", "A"]), d[np.ix_([1, 0], [1, 0])])


def test_rank_beta_sim_unknown_unit_raises():
    with pytest.raises(ValueError, match="Z"):
        rank_beta_sim(ranked(ROWS, LINEAGE), ["A", "Z"])


def test_rank_beta_sim_unit_without_taxa_at_a_rank_is_fully_dissimilar():
    # D has a species but its genus and family are unknown, so at those ranks D is a 0/0 against everything: distance 1 there
    lineage = {**LINEAGE, "s5": (None, None)}
    d = rank_beta_sim(ranked(ROWS + [("D", "s5"), ("E", "s5")], lineage), ["A", "B", "D", "E"])
    assert d[0, 2] == 1.0
    assert d[2, 3] == pytest.approx(2 / 3)  # D and E share their one species, so 0 there and 1 at the two empty ranks
    assert not d.diagonal().any()


def test_rank_beta_sim_needs_all_three_ranks():
    species_only = cell_taxa(pl.DataFrame({"cell": ["A", "B"], "species_id": ["s1", "s1"]}), None)
    with pytest.raises(ValueError, match="genus"):
        rank_beta_sim(species_only, ["A", "B"])


def planted_distances(**kwargs):
    presence, taxa, truth = planted_network(**kwargs)
    ct = cell_taxa(presence, taxa)
    return rank_beta_sim(ct, ct.cells), [truth[c] for c in ct.cells]


@pytest.mark.parametrize("method", ["ward", "average"])
def test_consensus_recovers_planted_realms(method):
    dist, truth = planted_distances(n_realms=4, cells_per_realm=30, endemic_species=60, widespread_species=10)
    result = consensus(dist, method, [4], shuffles=8)
    assert isinstance(result[4], Consensus)
    assert v_measure_score(truth, result[4].labels) >= 0.95
    assert len(set(result[4].labels.tolist())) == 4


def tie_heavy(n=240, realms=6, pool=30, per=4, wide=2, wide_pool=100, seed=0):
    """Beta-sim of cells that each hold `per` species of their realm's pool and `wide` of a shared pool, so most pairs tie at 1.

    Returns the distances (88% of the pairs at 1 by default), the planted realm of each cell and the cell ids.
    """
    rng = np.random.default_rng(seed)
    realm = np.repeat(np.arange(realms), n // realms)
    rows, cols = [], []
    for i in range(n):
        held = np.concatenate([realm[i] * pool + rng.choice(pool, per, replace=False), realms * pool + rng.integers(0, wide_pool, wide)])
        rows += [i] * len(held)
        cols += held.tolist()
    x = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, realms * pool + wide_pool))
    return squareform(beta_sim(x)).astype(np.float32), realm, [f"u{i:04d}" for i in range(n)]


def together(labels):
    return labels[:, None] == labels[None, :]


def test_tie_heavy_matrix_makes_a_single_tree_depend_on_the_order():
    # guards the order test below: on this matrix a single tree does change with the order, so a test that it doesn't would fail
    dist, _, _ = tie_heavy(per=3, pool=60)
    assert (dist[np.triu_indices(len(dist), 1)] == 1).mean() > 0.8
    trees = []
    for seed in range(5):
        order = np.random.default_rng(seed).permutation(len(dist))
        labels = hierarchical_labels(squareform(dist[np.ix_(order, order)], checks=False), "ward", 6)
        back = np.empty(len(order), dtype=int)
        back[order] = labels
        trees.append(back)
    assert min(v_measure_score(trees[0], other) for other in trees[1:]) < 0.99


@pytest.mark.parametrize("method", ["ward", "average"])
def test_consensus_does_not_depend_on_the_input_order(method):
    dist, _, ids = tie_heavy(per=3, pool=60)  # 94% of the pairs at 1, too many for the shuffles alone to settle the groups
    straight = consensus(dist, method, [4, 6], shuffles=10, ids=ids)
    for seed in range(3):
        order = np.random.default_rng(seed).permutation(len(dist))
        moved = consensus(dist[np.ix_(order, order)], method, [4, 6], shuffles=10, ids=[ids[i] for i in order])
        for k in (4, 6):
            back = np.empty(len(order), dtype=int)
            back[order] = moved[k].labels
            # exactly the same groups (the numbers differ, as they follow the first appearance in each input order), and spread
            assert np.array_equal(together(back), together(straight[k].labels))
            assert moved[k].spread == straight[k].spread


@pytest.mark.parametrize("method", ["ward", "average"])
def test_without_ids_the_input_order_reaches_the_groups(method):
    # the contrast to the test above, and what the docstring says: with no ids the input order is the canonical one
    dist, _, _ = tie_heavy(per=3, pool=60)
    straight = consensus(dist, method, [6], shuffles=10)[6].labels
    order = np.random.default_rng(0).permutation(len(dist))
    back = np.empty(len(order), dtype=int)
    back[order] = consensus(dist[np.ix_(order, order)], method, [6], shuffles=10)[6].labels
    assert not np.array_equal(together(back), together(straight))


@pytest.mark.parametrize("method", ["ward", "average"])
def test_consensus_recovers_planted_realms_among_ties(method):
    dist, realm, ids = tie_heavy()
    result = consensus(dist, method, [6], shuffles=10, ids=ids)[6]
    assert v_measure_score(realm, result.labels) >= 0.95


def test_labels_without_ids_are_numbered_by_first_appearance_in_the_input_order():
    dist, _, ids = tie_heavy(n=60, realms=3)
    order = np.random.default_rng(3).permutation(len(dist))
    labels = consensus(dist[np.ix_(order, order)], "ward", [3], shuffles=4, ids=[ids[i] for i in order])[3].labels
    assert labels[0] == 0 and list(dict.fromkeys(labels.tolist())) == [0, 1, 2]


def blocks(sizes, within=0.1, between=0.9, seed=3):
    """Distances that are small inside each block and large between them, with no ties among the within-block pairs."""
    rng = np.random.default_rng(seed)
    block = np.repeat(np.arange(len(sizes)), sizes)
    dist = np.where(block[:, None] == block[None, :], within, between) + rng.random((len(block), len(block))) * 0.01
    dist = np.triu(dist, 1)
    dist = (dist + dist.T).astype(np.float32)
    return dist, block


@pytest.mark.parametrize("method", ["ward", "average"])
def test_spread_is_one_when_every_shuffle_agrees(method):
    dist, block = blocks([10, 15, 20])
    result = consensus(dist, method, [3], shuffles=6)[3]
    assert result.spread == pytest.approx(1.0)
    assert v_measure_score(block, result.labels) == pytest.approx(1.0)


def test_spread_drops_when_the_order_decides():
    # every pair equally far: which units end up together is only the order they come in, so the shuffles disagree
    dist = (1 - np.eye(40)).astype(np.float32)
    assert consensus(dist, "ward", [4], shuffles=6)[4].spread < 0.9


def test_consensus_labels_number_groups_by_first_appearance():
    dist, _ = blocks([6, 6, 6])
    order = np.random.default_rng(1).permutation(len(dist))
    labels = consensus(dist[np.ix_(order, order)], "ward", [3], shuffles=4)[3].labels
    assert labels.dtype.kind == "i"
    assert labels[0] == 0
    assert list(dict.fromkeys(labels.tolist())) == [0, 1, 2]


def test_consensus_is_reproducible_and_seeded():
    dist = np.random.default_rng(2).integers(0, 3, (60, 60)).astype(np.float32)  # three values, so the order decides much
    dist = np.triu(dist, 1) + np.triu(dist, 1).T
    a, b = consensus(dist, "ward", [3, 5], shuffles=5, seed=4), consensus(dist, "ward", [3, 5], shuffles=5, seed=4)
    assert all(np.array_equal(a[k].labels, b[k].labels) and a[k].spread == b[k].spread for k in (3, 5))
    assert list(a) == [3, 5] and all(a[k].k == k for k in a)
    assert consensus(dist, "ward", [3], shuffles=5, seed=5)[3].spread != a[3].spread


def test_consensus_refuses_what_it_cannot_do():
    dist, _ = blocks([5, 5])
    with pytest.raises(ValueError, match="method"):
        consensus(dist, "single", [2])
    with pytest.raises(ValueError, match="k"):
        consensus(dist, "ward", [11])
    with pytest.raises(ValueError, match="shuffles"):
        consensus(dist, "ward", [2], shuffles=1)
    with pytest.raises(ValueError, match="square"):
        consensus(dist[:, :5], "ward", [2])
    with pytest.raises(ValueError, match="symmetric"):
        consensus(np.triu(dist), "ward", [2])
    with pytest.raises(ValueError, match="ids"):
        consensus(dist, "ward", [2], ids=list("abc"))
    with pytest.raises(ValueError, match="unique"):
        consensus(dist, "ward", [2], ids=list("abcdefghia"))
