import numpy as np
import pytest
from scipy.sparse import csr_matrix

from realms.betasim import anosim_r, beta_sim, hierarchical_labels


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
