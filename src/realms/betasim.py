"""Similarity-based regionalisation, as the Kreft and Jetz framework and the Costello critique do it."""

from typing import Literal

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform
from scipy.stats import rankdata


def beta_sim(cells_by_species: csr_matrix) -> np.ndarray:
    """Condensed Simpson dissimilarity, min(b, c) / (a + min(b, c)), which ignores differences in richness.

    With a shared species and b, c unique to each cell, a + min(b, c) is the richness of the poorer cell. A pair where one
    cell has no species counts as fully dissimilar.
    """
    binary = csr_matrix(cells_by_species, dtype=np.int64)
    binary.data = np.ones_like(binary.data)
    shared = (binary @ binary.T).toarray().astype(np.float64)
    richness = np.diag(shared).copy()
    poorer = np.minimum.outer(richness, richness)
    with np.errstate(divide="ignore", invalid="ignore"):
        dissimilarity = np.where(poorer > 0, 1.0 - shared / poorer, 1.0)
    np.fill_diagonal(dissimilarity, 0.0)
    return squareform(dissimilarity, checks=False)


def hierarchical_labels(dist: np.ndarray, method: Literal["ward", "average"], k: int) -> np.ndarray:
    """Cut a hierarchical clustering into exactly k groups. SciPy's ward on distances is R's ward.D2.

    cut_tree behaves like R's cutree. fcluster's maxclust would return fewer groups when merge heights tie, which beta-sim
    produces a lot of.
    """
    return cut_tree(linkage(dist, method=method), n_clusters=k).ravel()


def anosim_r(dist: np.ndarray, labels: np.ndarray) -> float:
    """ANOSIM R, (mean rank between - mean rank within) / (M / 2), with M the number of pairs."""
    ranks = rankdata(dist)
    left, right = np.triu_indices(len(labels), k=1)
    within = labels[left] == labels[right]
    return float((ranks[~within].mean() - ranks[within].mean()) / (len(dist) / 2))
