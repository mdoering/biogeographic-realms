"""The cell-taxon matrix, built once from a presence table and sliced from then on.

Rows are cells, columns are the taxa of every rank, named `<rank>:<id>` and sorted by that name. The network, the endemism counts,
the merges and the bootstrap resamples are all slices or products of this one boolean matrix, where each used to rebuild
string-keyed tables from the presence rows.
"""

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
import polars as pl
from scipy.sparse import coo_matrix, csr_matrix, hstack

RANKS = ("species", "genus", "family")
SPECIES = RANKS.index("species")


@dataclass(frozen=True)
class Incidence:
    x: csr_matrix  # binary, the cells asked for x all taxa in the table
    totals: np.ndarray  # occupied cells per taxon across the whole table, including cells not asked for


def _check_ranks(ranks: Sequence[str]) -> None:
    unknown = [r for r in ranks if r not in RANKS]
    if unknown:
        raise ValueError(f"unknown rank {unknown}, the ranks are {list(RANKS)}")


@dataclass(frozen=True)
class CellTaxa:
    """Which cells hold which taxa. `cells` and `taxa` are sorted, and `x[i, j]` says that cell `i` holds taxon `j`.

    A cell stays a row when none of its species are accepted. `rank` indexes `RANKS` for each column, and a species' genus and
    family columns are `genus_col` and `family_col` at the species' own column. Those hold -1 for a species without one, and on
    every column that isn't a species. `x` has sorted indices, which `network_from` relies on for its link order.
    """

    cells: list[str]
    taxa: list[str]
    rank: np.ndarray  # int8, per column
    x: csr_matrix  # bool, cells x taxa
    genus_col: np.ndarray  # int64, per column
    family_col: np.ndarray  # int64, per column

    def select(self, ranks: Sequence[str]) -> "CellTaxa":
        """The columns of these ranks, with the genus and family columns of species renumbered, or -1 if their rank is gone."""
        _check_ranks(ranks)
        keep = np.isin(self.rank, [RANKS.index(r) for r in ranks])
        if keep.all():
            return self
        kept = np.flatnonzero(keep)
        renumbered = np.full(len(keep), -1, dtype=np.int64)
        renumbered[kept] = np.arange(len(kept))

        def moved(parent: np.ndarray) -> np.ndarray:
            parent = parent[kept]
            return np.where(parent >= 0, renumbered[np.maximum(parent, 0)], -1)

        x = self.x[:, kept]
        x.sort_indices()
        taxa = [self.taxa[j] for j in kept.tolist()]
        return CellTaxa(self.cells, taxa, self.rank[kept], x, moved(self.genus_col), moved(self.family_col))

    def rows(self, cells: Sequence[str]) -> Incidence:
        """Binary cells x taxa for `cells`, in that order, with a cell this matrix doesn't know as an empty row.

        The totals count every cell of the matrix, so a taxon that also occurs outside the rows asked for isn't mistaken for an
        endemic of them.
        """
        position = {c: i for i, c in enumerate(self.cells)}
        asked = np.array([position.get(c, -1) for c in cells], dtype=np.int64)
        known = asked >= 0
        found = self.x[asked[known]]
        length = np.zeros(len(asked), dtype=np.int64)
        length[known] = np.diff(found.indptr)
        indptr = np.concatenate([[0], np.cumsum(length)]).astype(found.indptr.dtype)  # known rows keep their relative order
        x = csr_matrix((found.data.astype(np.int32), found.indices, indptr), shape=(len(asked), self.x.shape[1]))  # int32 counts, not bool
        return Incidence(x, self.x.getnnz(axis=0))

    def species_ids(self) -> list[str]:
        """The ids of the species columns, in column order."""
        return [self.taxa[j].removeprefix("species:") for j in np.flatnonzero(self.rank == SPECIES).tolist()]

    def resample(self, draws: np.ndarray) -> "CellTaxa":
        """Species drawn with replacement, each draw its own species so that one drawn twice counts twice.

        Draw `i` of `species_ids()[d]` becomes the column `species:<id>#<i>`. The genus and family columns are rebuilt from the
        drawn species, and one with no drawn species is gone. Columns are sorted by name again. A cell that none of the drawn
        species occupy is gone too, as it dropped out of a resampled presence table, so it isn't counted in the weights.
        """
        draws = np.asarray(draws, dtype=np.int64)
        species = np.flatnonzero(self.rank == SPECIES)
        ids = pl.Series(self.species_ids(), dtype=pl.String).gather(draws)
        names = pl.select((pl.lit("species:") + ids + pl.lit("#") + pl.int_range(len(draws)).cast(pl.String)).alias("name")).to_series()
        order = names.arg_sort().to_numpy()
        columns = species[draws[order]]
        x_species = self.x[:, columns]
        x_species.sort_indices()
        occupied = np.diff(x_species.indptr) > 0
        cells = self.cells
        if not occupied.all():
            x_species, cells = x_species[occupied], [c for c, o in zip(self.cells, occupied.tolist()) if o]

        parents = {}
        for rank, own in (("genus", self.genus_col), ("family", self.family_col)):
            parent = own[columns]
            used = np.unique(parent[parent >= 0])  # a column is never empty, so this is exactly the taxa with a drawn species
            if len(used):
                parents[rank] = ([self.taxa[j] for j in used.tolist()], np.where(parent >= 0, np.searchsorted(used, parent), -1))
        return _assemble(cells, names.gather(order).to_list(), x_species, parents)


def _assemble(
    cells: list[str], species: list[str], x_species: csr_matrix, parents: dict[str, tuple[list[str], np.ndarray]]
) -> CellTaxa:
    """Put the species columns and the genus and family columns rolled up from them side by side, in name order.

    `species` are the qualified names of the columns of `x_species`. For each rank in `parents` there are the qualified names of
    its taxa, sorted and each with a species under it, and per species the position of its taxon among them, -1 for none.
    """
    blocks = {"species": (species, x_species)}
    for rank, (names, parent) in parents.items():
        belongs = parent >= 0
        roll_up = csr_matrix(
            (np.ones(belongs.sum(), dtype=bool), (np.flatnonzero(belongs), parent[belongs])), shape=(len(species), len(names))
        )
        blocks[rank] = (names, x_species @ roll_up)

    layout = sorted(blocks, key=lambda rank: f"{rank}:")  # qualified names sort by rank first
    start, width = {}, 0
    for rank in layout:
        start[rank], width = width, width + len(blocks[rank][0])

    def parent_cols(rank: str) -> np.ndarray:
        cols = np.full(width, -1, dtype=np.int64)
        if rank in parents:
            parent = parents[rank][1]
            cols[start["species"] : start["species"] + len(species)] = np.where(parent >= 0, parent + start[rank], -1)
        return cols

    x = hstack([blocks[rank][1] for rank in layout], format="csr", dtype=bool)
    x.sort_indices()
    return CellTaxa(
        cells=cells,
        taxa=[name for rank in layout for name in blocks[rank][0]],
        rank=np.concatenate([np.full(len(blocks[rank][0]), RANKS.index(rank), dtype=np.int8) for rank in layout]),
        x=x,
        genus_col=parent_cols("genus"),
        family_col=parent_cols("family"),
    )


def _parents(species: pl.Series, taxa: pl.DataFrame, rank: str) -> tuple[list[str], np.ndarray]:
    """The qualified names of the taxa at `rank` that have a species in `species`, sorted, and each species' position among them."""
    under = (
        pl.DataFrame({"species_id": species})
        .with_row_index("column")
        .join(taxa.select("species_id", pl.col(f"{rank}_id").alias("value")).drop_nulls(), on="species_id")
        .with_columns((pl.col("value").rank("dense").cast(pl.Int64) - 1).alias("position"))
    )
    names = under.unique("position").sort("position").select((pl.lit(f"{rank}:") + pl.col("value")).alias("name"))["name"]
    position = np.full(len(species), -1, dtype=np.int64)
    position[under["column"].to_numpy()] = under["position"].to_numpy()
    return names.to_list(), position


def cell_taxa(presence: pl.DataFrame, taxa: pl.DataFrame | None, ranks: Sequence[str] = RANKS) -> CellTaxa:
    """The matrix for `cell, species_id` rows.

    A null species keeps its cell without a link. With `taxa`, which lists each accepted species once, any other species stops
    here, and its genus and family become columns where `genus_id` and `family_id` have values. Without `taxa` there are only
    species. The id of a taxon, not its name, tells taxa apart, because names repeat even within a kingdom.
    """
    _check_ranks(ranks)
    if taxa is not None and taxa["species_id"].is_duplicated().any():
        twice = taxa["species_id"].filter(taxa["species_id"].is_duplicated()).unique().sort().head(5).to_list()
        raise ValueError(f"taxa lists a species more than once, such as {twice}")

    cells = presence["cell"].drop_nulls().unique().sort()
    species = presence["species_id"].drop_nulls().unique().to_frame()
    if taxa is not None:
        species = species.join(taxa.select("species_id"), on="species_id", how="semi")
    species = species["species_id"].sort()

    pairs = (
        presence.lazy()
        .join(pl.DataFrame({"cell": cells}).with_row_index("row").lazy(), on="cell")
        .join(pl.DataFrame({"species_id": species}).with_row_index("column").lazy(), on="species_id")
        .select(pl.col("row").cast(pl.Int32), pl.col("column").cast(pl.Int32))
        .collect()
    )
    x_species = coo_matrix(
        (np.ones(pairs.height, dtype=bool), (pairs["row"].to_numpy(), pairs["column"].to_numpy())), shape=(len(cells), len(species))
    ).tocsr()  # a repeated row of the table is still one link

    parents = {}
    if taxa is not None:
        for rank in ("genus", "family"):
            if rank in ranks:
                taxon_names, position = _parents(species, taxa, rank)
                if taxon_names:
                    parents[rank] = (taxon_names, position)
    species_names = pl.select((pl.lit("species:") + species).alias("name")).to_series().to_list()
    matrix = _assemble(cells.to_list(), species_names, x_species, parents)
    return matrix if "species" in ranks else matrix.select(ranks)
