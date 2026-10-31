"""How much of a realm's biota lives nowhere else."""

import numpy as np
import polars as pl
from scipy.sparse import coo_matrix, csr_matrix

Labels = dict[str, int]


def presence_by_rank(presence: pl.DataFrame, taxa: pl.DataFrame | None) -> dict[str, pl.DataFrame]:
    """One `cell, species_id` table per rank that has values, with the genus or family name in `species_id`."""
    ranks = {"species": presence.select("cell", "species_id").unique()}
    if taxa is None:
        return ranks
    for rank in ("genus", "family"):
        if rank not in taxa.columns or taxa[rank].null_count() == taxa.height:
            continue
        rolled = presence.join(taxa.select("species_id", rank), on="species_id", how="inner").filter(pl.col(rank).is_not_null())
        ranks[rank] = rolled.select("cell", pl.col(rank).alias("species_id")).unique()
    return ranks


def incidence(presence: pl.DataFrame, cells: list[str]) -> csr_matrix:
    """Binary cells x taxa matrix, rows in the order of `cells`. Rows of cells not listed are ignored."""
    cell_frame = pl.DataFrame({"cell": cells, "row": range(len(cells))}, schema={"cell": pl.String, "row": pl.Int64})
    rows = presence.join(cell_frame, on="cell", how="inner").select("row", "species_id").unique()
    species = rows["species_id"].to_numpy(allow_copy=True)
    _, columns = np.unique(species.astype(str), return_inverse=True)
    n_taxa = int(columns.max()) + 1 if len(columns) else 0
    return coo_matrix((np.ones(len(columns), dtype=np.int32), (rows["row"].to_numpy(), columns)), shape=(len(cells), n_taxa)).tocsr()


def realm_share(x: csr_matrix, cell_realm: np.ndarray, n_realms: int, tolerance: float) -> np.ndarray:
    """Endemic share per realm, for a realm index per row of `x`.

    A taxon is endemic to a realm when at least `tolerance` of its cells lie inside. Realm order follows the index.
    """
    membership = coo_matrix(
        (np.ones(len(cell_realm), dtype=np.int32), (np.arange(len(cell_realm)), cell_realm)), shape=(x.shape[0], n_realms)
    ).tocsr()
    counts = (x.T @ membership).tocoo()
    totals = np.asarray(x.sum(axis=0)).ravel()
    present = np.bincount(counts.col, minlength=n_realms)
    endemic_rows = counts.data / totals[counts.row] >= tolerance
    endemic = np.bincount(counts.col[endemic_rows], minlength=n_realms)
    return np.divide(endemic, present, out=np.zeros(n_realms), where=present > 0)


def _indexed(labels: Labels) -> tuple[list[str], np.ndarray, list[int]]:
    cells = sorted(labels)
    realms = sorted(set(labels.values()))
    index = {r: i for i, r in enumerate(realms)}
    return cells, np.array([index[labels[c]] for c in cells], dtype=np.int64), realms


def realm_endemism(labels: Labels, presence: pl.DataFrame, tolerance: float = 0.95) -> dict[int, float]:
    cells, cell_realm, realms = _indexed(labels)
    shares = realm_share(incidence(presence, cells), cell_realm, len(realms), tolerance)
    return {r: float(s) for r, s in zip(realms, shares)}


def area_weighted(shares: np.ndarray, cell_realm: np.ndarray, cell_area: np.ndarray) -> float:
    realm_area = np.bincount(cell_realm, weights=cell_area, minlength=len(shares))
    return float((shares * realm_area).sum() / realm_area.sum())


def endemism_score(labels: Labels, presence: pl.DataFrame, areas: dict[str, float], tolerance: float = 0.95) -> float:
    """Endemic share averaged over realms, weighted by realm area."""
    cells, cell_realm, realms = _indexed(labels)
    shares = realm_share(incidence(presence, cells), cell_realm, len(realms), tolerance)
    return area_weighted(shares, cell_realm, np.array([areas[c] for c in cells]))
