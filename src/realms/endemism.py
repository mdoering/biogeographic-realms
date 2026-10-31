"""How much of a realm's biota lives nowhere else."""

import numpy as np
from scipy.sparse import coo_matrix

from realms.matrix import CellTaxa, Incidence

Labels = dict[str, int]


def realm_share(table: Incidence, cell_realm: np.ndarray, n_realms: int, tolerance: float) -> np.ndarray:
    """Endemic share per realm, for a realm index per row of the table.

    A taxon is endemic to a realm when at least `tolerance` of all its cells lie inside, counting cells outside the map too. Realm
    order follows the index.
    """
    x, totals = table.x, table.totals
    membership = coo_matrix(
        (np.ones(len(cell_realm), dtype=np.int32), (np.arange(len(cell_realm)), cell_realm)), shape=(x.shape[0], n_realms)
    ).tocsr()
    counts = (x.T @ membership).tocoo()
    present = np.bincount(counts.col, minlength=n_realms)
    endemic_rows = counts.data / totals[counts.row] >= tolerance
    endemic = np.bincount(counts.col[endemic_rows], minlength=n_realms)
    return np.divide(endemic, present, out=np.zeros(n_realms), where=present > 0)


def _indexed(labels: Labels) -> tuple[list[str], np.ndarray, list[int]]:
    cells = sorted(labels)
    realms = sorted(set(labels.values()))
    index = {r: i for i, r in enumerate(realms)}
    return cells, np.array([index[labels[c]] for c in cells], dtype=np.int64), realms


def realm_endemism(labels: Labels, ct: CellTaxa, tolerance: float = 0.95) -> dict[int, float]:
    cells, cell_realm, realms = _indexed(labels)
    shares = realm_share(ct.rows(cells), cell_realm, len(realms), tolerance)
    return {r: float(s) for r, s in zip(realms, shares)}


def area_weighted(shares: np.ndarray, cell_realm: np.ndarray, cell_area: np.ndarray) -> float:
    realm_area = np.bincount(cell_realm, weights=cell_area, minlength=len(shares))
    return float((shares * realm_area).sum() / realm_area.sum())


def endemism_score(labels: Labels, ct: CellTaxa, areas: dict[str, float], tolerance: float = 0.95) -> float:
    """Endemic share averaged over realms, weighted by realm area."""
    cells, cell_realm, realms = _indexed(labels)
    shares = realm_share(ct.rows(cells), cell_realm, len(realms), tolerance)
    return area_weighted(shares, cell_realm, np.array([areas[c] for c in cells]))
