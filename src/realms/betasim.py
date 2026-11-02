"""Similarity-based regionalisation, as the Kreft and Jetz framework and the Costello critique do it."""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.sparse import csr_matrix
from scipy.spatial.distance import squareform
from scipy.stats import rankdata
from sklearn.metrics import v_measure_score

from realms.matrix import RANKS, CellTaxa


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


def rank_beta_sim(ct: CellTaxa, units: Sequence[str]) -> np.ndarray:
    """Beta-sim between `units`, averaged over the species, genus and family columns: the mean of the three distance matrices.

    `units` are row ids of `ct`, in the order of the output, so they may be a subset or a reordering of `ct.cells`; an id that
    `ct` doesn't know raises. Returns a square float32 matrix, symmetric, with a zero diagonal. A unit that holds no taxon at
    a rank (a 0/0) is at distance 1 from every other unit at that rank, the same as `beta_sim` takes it, so the unit is not
    taken for similar to another that holds none either. A rank with no column in `ct` at all raises, as the average would
    otherwise quietly count one rank as no information.

    The matrix is built in row blocks, so one block of float32 is all that is held beside the result, which at 11,292 units
    is 510 MB, where `beta_sim` would hold several float64 copies of the square.
    """
    position = {c: i for i, c in enumerate(ct.cells)}
    unknown = [u for u in units if u not in position]
    if unknown:
        raise ValueError(f"{len(unknown)} units are not rows of the matrix, such as {unknown[:5]}")
    rows = ct.x[np.array([position[u] for u in units], dtype=np.int64)]
    n = len(units)
    total = np.zeros((n, n), dtype=np.float32)
    step = max(1, 2**24 // max(n, 1))  # blocks of about 64 MB
    for rank in RANKS:
        columns = np.flatnonzero(ct.rank == RANKS.index(rank))
        if not len(columns):
            raise ValueError(f"the matrix has no {rank} columns, so there is no average over the three ranks")
        present = csr_matrix(rows[:, columns], dtype=np.float32)
        present_t = present.T.tocsr()
        richness = np.diff(present.indptr).astype(np.float32)
        for start in range(0, n, step):
            shared = (present[start : start + step] @ present_t).toarray()  # integer counts, exact in float32
            poorer = np.minimum(richness[start : start + step, None], richness[None, :])
            with np.errstate(divide="ignore", invalid="ignore"):
                total[start : start + step] += np.where(poorer > 0, 1 - shared / poorer, np.float32(1))
    total /= 3
    np.fill_diagonal(total, 0)
    return total


@dataclass(frozen=True)
class Consensus:
    """The consensus grouping of the units into `k` groups.

    `labels` number the groups 0..k-1 by first appearance in the unit order. `spread` is the mean pairwise V-measure between
    the shuffles' own groupings at this k: 1.0 when the cell order never mattered, lower the more it did.
    """

    k: int
    labels: np.ndarray
    spread: float


def _condensed(square: np.ndarray, order: np.ndarray | None = None) -> np.ndarray:
    """The upper triangle of `square[order][:, order]` as a float64 vector, row by row, without ever copying the square."""
    n = len(square)
    order = np.arange(n) if order is None else order
    out = np.empty(n * (n - 1) // 2, dtype=np.float64)
    at = 0
    for i in range(n - 1):
        length = n - 1 - i
        out[at : at + length] = square[order[i]][order[i + 1 :]]
        at += length
    return out


def _by_first_appearance(labels: np.ndarray) -> np.ndarray:
    _, first, inverse = np.unique(labels, return_index=True, return_inverse=True)
    return np.argsort(np.argsort(first))[inverse.ravel()]


def consensus(
    dist: np.ndarray, method: Literal["ward", "average"], ks: Sequence[int], shuffles: int = 20, seed: int = 1
) -> dict[int, Consensus]:
    """A grouping at each k that does not depend on the order of the units, from Ward or UPGMA trees over shuffled orders.

    Ties in the distances, which beta-sim has by the million, are broken by the order the units come in, so a single tree
    changes with it. Here the units are put in `shuffles` random orders (a generator seeded with `seed`), each is clustered by
    `method` ("ward" or "average", UPGMA) and the one tree is cut at every k. For each k, `co[i, j]` is the share of shuffles
    that put units i and j in the same group, and the groups are an average linkage of `1 - co` cut at k. The `spread` of
    a k is how much the shuffles agree with each other there, the mean V-measure over all pairs of them.

    `dist` is the square matrix of `rank_beta_sim`. The input order still reaches the answer through the shuffles that are
    drawn for it and the ties of the last linkage, so a reordering gives the same partition when the structure is clear,
    which is what the consensus is for, but not bit for bit.

    The counts of every k are held at once, one n x n uint8 matrix each (127 MB per k at 11,292 units), and the last linkage
    for a k frees its own.
    """
    if method not in ("ward", "average"):
        raise ValueError(f"method is 'ward' or 'average', not {method!r}")
    dist = np.asarray(dist)
    n = len(dist)
    if dist.ndim != 2 or dist.shape != (n, n):
        raise ValueError(f"dist must be square, not {dist.shape}")
    ks = sorted({int(k) for k in ks})
    if not ks or ks[0] < 1 or ks[-1] > n:
        raise ValueError(f"each k must be between 1 and the {n} units, not {ks[:1] + ks[-1:]}")
    if not 2 <= shuffles <= 255:
        raise ValueError(f"shuffles must be 2 to 255, as a spread needs two groupings and the counts are bytes, not {shuffles}")
    for start in range(0, n, 2048):  # blockwise, so that no second square is made
        if not np.array_equal(dist[start : start + 2048], dist[:, start : start + 2048].T):
            raise ValueError("dist must be symmetric")

    rng = np.random.default_rng(seed)
    counts = {k: np.zeros((n, n), dtype=np.uint8) for k in ks}
    shuffled = {k: [] for k in ks}
    for _ in range(shuffles):
        order = rng.permutation(n)
        tree = linkage(_condensed(dist, order), method=method)
        cut = cut_tree(tree, n_clusters=ks)  # n x len(ks), in the shuffled order
        del tree
        for column, k in enumerate(ks):
            labels = np.empty(n, dtype=np.int64)
            labels[order] = cut[:, column]  # back to the original order
            shuffled[k].append(labels)
            for start in range(0, n, 1024):
                counts[k][start : start + 1024] += labels[start : start + 1024, None] == labels[None, :]

    result = {}
    for k in ks:
        co = counts.pop(k)
        distance = np.empty(n * (n - 1) // 2, dtype=np.float64)
        at = 0
        for i in range(n - 1):
            length = n - 1 - i
            distance[at : at + length] = 1 - co[i, i + 1 :] / shuffles
            at += length
        del co
        labels = _by_first_appearance(cut_tree(linkage(distance, method="average"), n_clusters=k).ravel())
        del distance
        pairs = [v_measure_score(a, b) for i, a in enumerate(shuffled[k]) for b in shuffled[k][i + 1 :]]
        result[k] = Consensus(k=k, labels=labels, spread=float(np.mean(pairs)))
    return result
