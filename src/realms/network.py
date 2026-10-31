"""Bipartite network of grid cells and taxa at several ranks.

Species link to the cells they occur in. A genus or family links to every cell holding one of its species, weighted by the
spatial information it carries, following Infomap Bioregions 2 (Edler et al.).
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import compress

import numpy as np
import polars as pl

from realms.matrix import RANKS, SPECIES, CellTaxa, cell_taxa


@dataclass(frozen=True)
class Network:
    cells: list[str]
    taxa: list[str]
    cell_idx: np.ndarray
    taxon_idx: np.ndarray
    weight: np.ndarray


def information_weight(n_cells: int, total_cells: int) -> float:
    """1 - log(n) / log(N): full weight for a taxon in one cell, none for a taxon in every cell."""
    if total_cells == 1:
        return 1.0
    return 1.0 - math.log(n_cells) / math.log(total_cells)


def accepted_presence(presence: pl.DataFrame, taxa: pl.DataFrame | None) -> pl.DataFrame:
    """Distinct `cell, species_id` rows without null species. `taxa` lists the accepted species, and anything else stops here."""
    occupied = presence.filter(pl.col("species_id").is_not_null()).select("cell", "species_id").unique()
    return occupied if taxa is None else occupied.join(taxa, on="species_id", how="semi")


def rank_values(taxa: pl.DataFrame, rank: str) -> pl.DataFrame:
    """`species_id, value` for a rank, without nulls. The value is the taxon's id, because names repeat even within a kingdom."""
    return taxa.select("species_id", pl.col(f"{rank}_id").alias("value")).filter(pl.col("value").is_not_null())


def network_from(ct: CellTaxa) -> Network:
    """Links for every cell-taxon pair of the matrix, weighted as `information_weight` says for genera and families.

    A genus or family held by every cell has weight 0 and drops out before the nodes are numbered. A cell with no links stays in
    the network, so a cell with no species is reported rather than lost.
    """
    n_cells = len(ct.cells)
    weight_of_count = np.array([0.0] + [information_weight(k, n_cells) for k in range(1, n_cells + 1)])
    weight = np.where(ct.rank == SPECIES, 1.0, weight_of_count[ct.x.getnnz(axis=0)])
    keep = weight > 0

    cell_idx = np.repeat(np.arange(n_cells, dtype=np.int64), np.diff(ct.x.indptr))
    column = ct.x.indices
    if not keep.all():
        linked = keep[column]
        cell_idx, column = cell_idx[linked], column[linked]
    node = np.cumsum(keep) - 1
    return Network(
        cells=ct.cells,
        taxa=list(compress(ct.taxa, keep.tolist())),
        cell_idx=cell_idx,
        taxon_idx=node[column],
        weight=weight[column],
    )


def build_network(presence: pl.DataFrame, taxa: pl.DataFrame | None = None, ranks: Sequence[str] = RANKS) -> Network:
    """Build the network from `cell, species_id` rows.

    A null species_id keeps its cell in the network without a link, so a cell with no species is reported rather than lost.
    With `taxa`, only the species it lists become nodes, and one with no genus or family simply gets no node at that rank.
    """
    return network_from(cell_taxa(presence, taxa, ranks))
