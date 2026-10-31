"""Frozen copies of the pre-speed-up code, used only as test oracles.

These are the functions as they stood at 5bcde8f, copied verbatim. They are never imported by the package and never edited:
the equivalence tests compare the integer-coded matrix code against them, so a change here would hide a change of result.
`connected_pieces`, `_beta_sim`, `_unit_matrix` and `_rank_links` are the helpers the copies call. `run_infomap` is whole, and
its `add_link` loop is what the bulk link input is checked against. `Network`, `Partition`, `Incidence`, `Labels` and
`information_weight` are data containers and a formula that the speed-up leaves alone, so they are imported.
"""

from collections.abc import Sequence

import infomap
import numpy as np
import polars as pl
from scipy.sparse import coo_matrix, csr_matrix, vstack

from realms.cluster import Partition
from realms.endemism import Labels
from realms.matrix import Incidence
from realms.network import Network, information_weight


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


def connected_pieces(cells: set[str], adjacency: dict[str, set[str]]) -> list[set[str]]:
    """Connected components of `cells`, largest first, with ties in order of their first cell."""
    pieces, seen = [], set()
    for start in sorted(cells):
        if start in seen:
            continue
        piece, stack = {start}, [start]
        while stack:
            for neighbour in adjacency.get(stack.pop(), ()):
                if neighbour in cells and neighbour not in piece:
                    piece.add(neighbour)
                    stack.append(neighbour)
        seen |= piece
        pieces.append(piece)
    return sorted(pieces, key=lambda p: (-len(p), min(p)))


def _beta_sim(a: np.ndarray, b: np.ndarray) -> float:
    poorer = min(int(a.sum()), int(b.sum()))
    return 1.0 - int((a & b).sum()) / poorer if poorer else 1.0


def resolve_pieces(
    labels: Labels,
    adjacency: dict[str, set[str]],
    presence: pl.DataFrame,
    min_shared: int = 3,
    tolerance: float = 0.95,
) -> Labels:
    """Keep a separate piece in its realm only if it really shares restricted taxa with the rest of that realm.

    Cells with few records often hold only widespread taxa, which makes distant cells look alike. A piece stays if at least
    `min_shared` of its taxa are endemic to the realm and also occur in another piece. Otherwise it moves to the adjacent
    realm it shares most with, or the most similar realm of all if nothing touches it. The largest piece never moves, and neither
    does a piece that shares nothing with any option, since moving it would be arbitrary.
    """
    cells = sorted(labels)
    index = {c: i for i, c in enumerate(cells)}
    table = incidence(presence, cells)
    x, totals = table.x.astype(bool), table.totals

    realm_cells: dict[int, set[str]] = {}
    for cell, realm in labels.items():
        realm_cells.setdefault(realm, set()).add(cell)
    realm_rows = {r: [index[c] for c in members] for r, members in realm_cells.items()}

    def taxa_of(rows: list[int]) -> np.ndarray:
        return np.asarray(x[rows].sum(axis=0)).ravel() > 0

    result = dict(labels)
    for realm, members in sorted(realm_cells.items()):
        pieces = connected_pieces(members, adjacency)
        if len(pieces) == 1:
            continue
        inside = np.asarray(x[realm_rows[realm]].sum(axis=0)).ravel()
        endemic = (inside > 0) & (inside / np.maximum(totals, 1) >= tolerance)
        for piece in pieces[1:]:
            piece_rows = [index[c] for c in piece]
            rest_rows = [index[c] for c in members - piece]
            piece_taxa = taxa_of(piece_rows)
            if int((piece_taxa & taxa_of(rest_rows) & endemic).sum()) >= min_shared:
                continue
            touching = {labels[n] for c in piece for n in adjacency.get(c, ()) if labels.get(n, realm) != realm}
            options = sorted(touching or (set(realm_cells) - {realm}))
            if options:
                best = min(options, key=lambda r: (_beta_sim(piece_taxa, taxa_of(realm_rows[r])), r))
                if _beta_sim(piece_taxa, taxa_of(realm_rows[best])) < 1.0:
                    for cell in piece:
                        result[cell] = best
    return result


def _unit_matrix(labels: Labels, presence: pl.DataFrame) -> tuple[list[int], csr_matrix]:
    cells = sorted(labels)
    units = sorted(set(labels.values()))
    position = {u: i for i, u in enumerate(units)}
    x = incidence(presence, cells).x
    membership = csr_matrix(
        (np.ones(len(cells), dtype=np.int32), (np.arange(len(cells)), [position[labels[c]] for c in cells])),
        shape=(len(cells), len(units)),
    )
    grouped = (membership.T @ x).tocsr()
    grouped.data = np.ones_like(grouped.data)
    return units, grouped


class _Merger:
    """Units that merge one pair at a time, keeping pairwise beta-sim up to date."""

    def __init__(self, labels: Labels, presence: pl.DataFrame):
        self.units, grouped = _unit_matrix(labels, presence)
        self.rows = [grouped[i] for i in range(grouped.shape[0])]
        self.shared = (grouped @ grouped.T).toarray().astype(np.float64)
        self.active = np.ones(len(self.units), dtype=bool)

    def beta_sim(self) -> np.ndarray:
        richness = np.diag(self.shared)
        poorer = np.minimum.outer(richness, richness)
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(poorer > 0, 1.0 - self.shared / poorer, 1.0)
        return beta

    def merge(self, keep: int, drop: int) -> None:
        merged = self.rows[keep] + self.rows[drop]
        merged.data = np.ones_like(merged.data)
        self.rows[keep] = merged
        self.rows[drop] = csr_matrix(merged.shape, dtype=merged.dtype)
        self.active[drop] = False
        overlap = (vstack(self.rows).tocsr() @ merged.T).toarray().ravel().astype(np.float64)
        self.shared[keep, :] = overlap
        self.shared[:, keep] = overlap
        self.shared[drop, :] = 0.0
        self.shared[:, drop] = 0.0

    def closest_pair(self, among: np.ndarray | None = None) -> tuple[int, int]:
        """The lowest beta-sim pair of active units, ties going to the lowest labels. `among` limits the first of the pair."""
        beta = self.beta_sim()
        allowed = np.outer(self.active, self.active)
        if among is None:
            allowed &= np.triu(np.ones_like(allowed), k=1)
        else:
            allowed &= np.outer(among, np.ones(len(among), dtype=bool))
            np.fill_diagonal(allowed, False)
        i, j = np.unravel_index(np.argmin(np.where(allowed, beta, np.inf)), beta.shape)
        return int(i), int(j)


def merge_sequence(labels: Labels, presence: pl.DataFrame, k_max: int = 60, k_min: int = 5) -> dict[int, Labels]:
    """One nested map per realm count, merging the pair of units that share the most taxa by beta-sim each time."""
    merger = _Merger(labels, presence)
    current = dict(labels)
    maps: dict[int, Labels] = {}
    k = len(merger.units)
    while k >= k_min:
        if k <= k_max:
            maps[k] = dict(current)
        if k == k_min:
            break
        keep, drop = merger.closest_pair()
        keep_label, drop_label = merger.units[keep], merger.units[drop]
        current = {c: (keep_label if r == drop_label else r) for c, r in current.items()}
        merger.merge(keep, drop)
        k -= 1
    return maps


def apply_size_floor(
    labels: Labels, areas: dict[str, float], presence: pl.DataFrame, floor_km2: float = 1_000_000
) -> tuple[Labels, dict[int, int]]:
    """Merge units under the floor into the realm they share the most taxa with, smallest first.

    Returns the new labels and `{sub_realm: parent}`, where the parent is the realm the unit ended up in.
    """
    merger = _Merger(labels, presence)
    current = dict(labels)
    area = {u: 0.0 for u in merger.units}
    for cell, unit in labels.items():
        area[unit] += areas[cell]
    position = {u: i for i, u in enumerate(merger.units)}
    absorbed_by: dict[int, int] = {}

    while True:
        small = [u for u in merger.units if merger.active[position[u]] and area[u] < floor_km2]
        if not small or int(merger.active.sum()) == 1:
            break
        unit = min(small, key=lambda u: (area[u], u))
        only_this = np.zeros(len(merger.units), dtype=bool)
        only_this[position[unit]] = True
        _, other = merger.closest_pair(among=only_this)
        target = merger.units[other]
        current = {c: (target if r == unit else r) for c, r in current.items()}
        merger.merge(other, position[unit])
        area[target] += area.pop(unit)
        absorbed_by[unit] = target

    def final(unit: int) -> int:
        while unit in absorbed_by:
            unit = absorbed_by[unit]
        return unit

    return current, {unit: final(unit) for unit in absorbed_by}


def _resample(presence: pl.DataFrame, taxa: pl.DataFrame | None, rng: np.random.Generator) -> tuple[pl.DataFrame, pl.DataFrame | None]:
    """Draw species with replacement. Each draw becomes its own species, so a species drawn twice counts twice."""
    species = sorted(presence["species_id"].drop_nulls().unique().to_list())
    draws = rng.integers(len(species), size=len(species))
    drawn = pl.DataFrame(
        {"species_id": [species[d] for d in draws], "draw": [f"{species[d]}#{i}" for i, d in enumerate(draws)]},
        schema={"species_id": pl.String, "draw": pl.String},
    )
    new_presence = presence.join(drawn, on="species_id").select("cell", pl.col("draw").alias("species_id"))
    if taxa is None:
        return new_presence, None
    new_taxa = taxa.join(drawn, on="species_id").drop("species_id").rename({"draw": "species_id"})
    return new_presence, new_taxa


def run_infomap(net: Network, seed: int = 1, trials: int = 10) -> Partition:
    n_cells = len(net.cells)
    linked = set(net.cell_idx.tolist())
    unassigned = sorted(net.cells[i] for i in range(n_cells) if i not in linked)
    if not linked:
        return Partition(paths={}, unassigned=unassigned)

    im = infomap.Infomap(f"--flow-model undirected --seed {seed} --num-trials {trials} --silent")
    im.bipartite_start_id = n_cells
    for cell, taxon, weight in zip(net.cell_idx.tolist(), net.taxon_idx.tolist(), net.weight.tolist()):
        im.add_link(cell, n_cells + taxon, weight)
    modules = im.run().multilevel_modules()

    paths = {net.cells[node]: tuple(path) for node, path in modules.items() if node < n_cells}
    return Partition(paths=paths, unassigned=unassigned)
