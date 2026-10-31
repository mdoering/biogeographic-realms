"""How much of a realm's biota lives nowhere else."""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl
from scipy.sparse import coo_matrix, csr_matrix

from realms.network import accepted_presence, rank_values

Labels = dict[str, int]


@dataclass(frozen=True)
class Incidence:
    x: csr_matrix  # binary, the cells asked for x all taxa in the table
    totals: np.ndarray  # occupied cells per taxon across the whole table, including cells not asked for


def presence_by_rank(presence: pl.DataFrame, taxa: pl.DataFrame | None) -> dict[str, pl.DataFrame]:
    """One `cell, species_id` table per rank that has values, with the genus or family id in `species_id`."""
    occupied = accepted_presence(presence, taxa)
    ranks = {"species": occupied}
    if taxa is None:
        return ranks
    for rank in ("genus", "family"):
        if taxa[f"{rank}_id"].null_count() == taxa.height:
            continue
        rolled = occupied.join(rank_values(taxa, rank), on="species_id", how="inner")
        ranks[rank] = rolled.select("cell", pl.col("value").alias("species_id")).unique()
    return ranks


def rank_presence(presence: pl.DataFrame, taxa: pl.DataFrame | None, ranks: Sequence[str]) -> pl.DataFrame:
    """All the chosen ranks in one table, each taxon named by rank, for comparing cells or realms at every depth at once."""
    tables = presence_by_rank(presence, taxa)
    return pl.concat(
        table.select("cell", (pl.lit(f"{rank}:") + pl.col("species_id")).alias("species_id")) for rank, table in tables.items() if rank in ranks
    )


def incidence(presence: pl.DataFrame, cells: list[str]) -> Incidence:
    """Binary cells x taxa for `cells`, in that order, plus each taxon's cell count over the whole table.

    The totals include cells not listed, so a taxon that also occurs outside a map isn't mistaken for an endemic of it.
    """
    pairs = presence.filter(pl.col("species_id").is_not_null()).select("cell", "species_id").unique()
    pairs = pairs.with_columns((pl.col("species_id").rank("dense").cast(pl.Int64) - 1).alias("column"))
    n_taxa = int(pairs["column"].max()) + 1 if pairs.height else 0
    totals = np.bincount(pairs["column"].to_numpy(), minlength=n_taxa)
    cell_frame = pl.DataFrame({"cell": cells, "row": range(len(cells))}, schema={"cell": pl.String, "row": pl.Int64})
    rows = pairs.join(cell_frame, on="cell", how="inner")
    x = coo_matrix((np.ones(rows.height, dtype=np.int32), (rows["row"].to_numpy(), rows["column"].to_numpy())), shape=(len(cells), n_taxa)).tocsr()
    return Incidence(x, totals)


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
