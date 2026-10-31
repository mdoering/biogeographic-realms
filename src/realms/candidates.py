"""From clustering modules to one nested map per realm count."""

from dataclasses import dataclass

import numpy as np
import polars as pl
from scipy.sparse import csr_matrix, vstack

from realms.cluster import Partition
from realms.endemism import Labels, incidence


class DegenerateClustering(Exception):
    """Too few units survive to build candidate maps from."""


@dataclass(frozen=True)
class CandidateMaps:
    maps: dict[int, Labels]
    sub_realms: dict[int, int]
    units: Labels  # the units before the size floor, which is where a sub-realm's cells can be read from


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
    realm it shares most with, or the most similar realm of all if nothing touches it. The largest piece never moves.
    """
    cells = sorted(labels)
    index = {c: i for i, c in enumerate(cells)}
    x = incidence(presence, cells).astype(bool)
    totals = np.asarray(x.sum(axis=0)).ravel()

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
                for cell in piece:
                    result[cell] = best
    return result


def _unit_matrix(labels: Labels, presence: pl.DataFrame) -> tuple[list[int], csr_matrix]:
    cells = sorted(labels)
    units = sorted(set(labels.values()))
    position = {u: i for i, u in enumerate(units)}
    x = incidence(presence, cells).astype(np.int32)
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


def candidate_maps(
    partition: Partition,
    presence: pl.DataFrame,
    areas: dict[str, float],
    adjacency: dict[str, set[str]],
    k_max: int = 60,
    k_min: int = 5,
    floor_km2: float = 1_000_000,
) -> CandidateMaps:
    """Finest modules, then separate pieces, then the size floor, then the merge sequence.

    Flooring first keeps exactly one map per realm count, and a merge only coarsens, so pieces need resolving just once.
    """
    units = resolve_pieces(partition.leaf_labels(), adjacency, presence)
    floored, sub_realms = apply_size_floor(units, areas, presence, floor_km2)
    if len(set(floored.values())) < k_min:
        raise DegenerateClustering(f"{len(set(floored.values()))} units above the size floor, need at least {k_min}")
    return CandidateMaps(maps=merge_sequence(floored, presence, k_max, k_min), sub_realms=sub_realms, units=units)
