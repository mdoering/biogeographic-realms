"""Bipartite network of grid cells and taxa at several ranks.

Species link to the cells they occur in. A genus or family links to every cell holding one of its species, weighted by the
spatial information it carries, following Infomap Bioregions 2 (Edler et al.).
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl


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


def _rank_links(presence: pl.DataFrame, taxa: pl.DataFrame | None, rank: str, total_cells: int) -> pl.DataFrame:
    if rank == "species":
        return presence.select(
            "cell",
            (pl.lit("species:") + pl.col("species_id")).alias("taxon"),
            pl.lit(1.0).alias("weight"),
        )
    pairs = presence.join(rank_values(taxa, rank), on="species_id", how="inner").select("cell", "value").unique()
    weights = pairs.group_by("value").len().with_columns(
        pl.col("len").map_elements(lambda n: information_weight(n, total_cells), return_dtype=pl.Float64).alias("weight")
    )
    return pairs.join(weights, on="value").select(
        "cell",
        (pl.lit(f"{rank}:") + pl.col("value")).alias("taxon"),
        "weight",
    )


def build_network(
    presence: pl.DataFrame,
    taxa: pl.DataFrame | None = None,
    ranks: Sequence[str] = ("species", "genus", "family"),
) -> Network:
    """Build the network from `cell, species_id` rows.

    A null species_id keeps its cell in the network without a link, so a cell with no species is reported rather than lost.
    With `taxa`, only the species it lists become nodes, and one with no genus or family simply gets no node at that rank.
    """
    cells = sorted(presence["cell"].unique().to_list())
    occupied = accepted_presence(presence, taxa)
    use_ranks = [r for r in ranks if r == "species" or taxa is not None]

    links = pl.concat([_rank_links(occupied, taxa, rank, len(cells)) for rank in use_ranks]).filter(pl.col("weight") > 0)
    taxon_names = sorted(links["taxon"].unique().to_list())

    cell_index = {c: i for i, c in enumerate(cells)}
    taxon_index = {t: i for i, t in enumerate(taxon_names)}
    links = links.sort("cell", "taxon")
    return Network(
        cells=cells,
        taxa=taxon_names,
        cell_idx=np.array([cell_index[c] for c in links["cell"]], dtype=np.int64),
        taxon_idx=np.array([taxon_index[t] for t in links["taxon"]], dtype=np.int64),
        weight=links["weight"].to_numpy(),
    )
