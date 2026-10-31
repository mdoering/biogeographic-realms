"""From clustering modules to one nested map per realm count."""

from dataclasses import dataclass
from heapq import heapify, heappop, heappush

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from realms.cluster import Partition
from realms.endemism import Labels
from realms.matrix import CellTaxa


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


def _membership(labels: Labels) -> tuple[list[str], list[int], csr_matrix]:
    """The sorted cells, the sorted units, and which unit each cell belongs to as a cells x units indicator."""
    cells = sorted(labels)
    units = sorted(set(labels.values()))
    position = {u: i for i, u in enumerate(units)}
    membership = coo_matrix(
        (np.ones(len(cells), dtype=np.int32), (np.arange(len(cells)), [position[labels[c]] for c in cells])),
        shape=(len(cells), len(units)),
    ).tocsr()
    return cells, units, membership


def resolve_pieces(
    labels: Labels,
    adjacency: dict[str, set[str]],
    ct: CellTaxa,
    min_shared: int = 3,
    tolerance: float = 0.95,
) -> Labels:
    """Keep a separate piece in its realm only if it really shares restricted taxa with the rest of that realm.

    Cells with few records often hold only widespread taxa, which makes distant cells look alike. A piece stays if at least
    `min_shared` of its taxa are endemic to the realm and also occur in another piece. Otherwise it moves to the adjacent
    realm it shares most with, or the most similar realm of all if nothing touches it. The largest piece never moves, and neither
    does a piece that shares nothing with any option, since moving it would be arbitrary.
    """
    cells, realms, membership = _membership(labels)
    index = {c: i for i, c in enumerate(cells)}
    position = {r: i for i, r in enumerate(realms)}
    table = ct.rows(cells)
    x, totals = table.x, table.totals

    realm_cells: dict[int, set[str]] = {}
    for cell, realm in labels.items():
        realm_cells.setdefault(realm, set()).add(cell)
    inside_counts = (membership.T @ x).tocsr()  # the cells of each realm that hold each taxon
    realm_taxa = inside_counts > 0

    def taxa_of(realm: int) -> np.ndarray:
        return realm_taxa[position[realm]].toarray().ravel()

    result = dict(labels)
    for realm, members in sorted(realm_cells.items()):
        pieces = connected_pieces(members, adjacency)
        if len(pieces) == 1:
            continue
        inside = inside_counts[position[realm]].toarray().ravel()
        endemic = (inside > 0) & (inside / np.maximum(totals, 1) >= tolerance)
        for piece in pieces[1:]:
            in_piece = np.asarray(x[[index[c] for c in piece]].sum(axis=0)).ravel()
            piece_taxa = in_piece > 0
            rest_taxa = inside > in_piece  # held by a cell of the realm outside the piece
            if int((piece_taxa & rest_taxa & endemic).sum()) >= min_shared:
                continue
            touching = {labels[n] for c in piece for n in adjacency.get(c, ()) if labels.get(n, realm) != realm}
            options = sorted(touching or (set(realm_cells) - {realm}))
            if options:
                best = min(options, key=lambda r: (_beta_sim(piece_taxa, taxa_of(r)), r))
                if _beta_sim(piece_taxa, taxa_of(best)) < 1.0:
                    for cell in piece:
                        result[cell] = best
    return result


def _unit_matrix(labels: Labels, ct: CellTaxa) -> tuple[list[int], csr_matrix]:
    cells, units, membership = _membership(labels)
    grouped = (membership.T @ ct.rows(cells).x).tocsr()
    grouped.data = np.ones_like(grouped.data)
    return units, grouped


class _Merger:
    """Units that merge one pair at a time, keeping pairwise beta-sim up to date."""

    def __init__(self, labels: Labels, ct: CellTaxa):
        self.units, self.x = _unit_matrix(labels, ct)  # units x taxa, binary, and a merged-away unit is an empty row
        self.shared = (self.x @ self.x.T).toarray().astype(np.float64)
        self.active = np.ones(len(self.units), dtype=bool)

    def beta_sim(self) -> np.ndarray:
        richness = np.diag(self.shared)
        poorer = np.minimum.outer(richness, richness)
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(poorer > 0, 1.0 - self.shared / poorer, 1.0)
        return beta

    def merge(self, keep: int, drop: int) -> None:
        k = len(self.units)
        others = np.delete(np.arange(k), drop)
        # a k x k update: every unit keeps its row, `keep` also takes the row of `drop`, and `drop` is left empty
        update = coo_matrix((np.ones(k, dtype=np.int32), (np.append(others, keep), np.append(others, drop))), shape=(k, k)).tocsr()
        self.x = (update @ self.x).tocsr()
        self.x.data[:] = 1  # a taxon held by both units is still one taxon
        self.active[drop] = False
        overlap = (self.x @ self.x[keep].T).toarray().ravel().astype(np.float64)
        self.shared[keep, :] = overlap
        self.shared[:, keep] = overlap
        self.shared[drop, :] = 0.0
        self.shared[:, drop] = 0.0

    def closest_pair(self) -> tuple[int, int]:
        """The lowest beta-sim pair of active units, ties going to the lowest labels."""
        beta = self.beta_sim()
        allowed = np.outer(self.active, self.active)
        allowed &= np.triu(np.ones_like(allowed), k=1)
        i, j = np.unravel_index(np.argmin(np.where(allowed, beta, np.inf)), beta.shape)
        return int(i), int(j)


def merge_sequence(labels: Labels, ct: CellTaxa, k_max: int = 60, k_min: int = 5) -> dict[int, Labels]:
    """One nested map per realm count, merging the pair of units that share the most taxa by beta-sim each time."""
    merger = _Merger(labels, ct)
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


class _TaxaIndex:
    """Which unit holds each taxon, for units that merge into one another one at a time.

    The overlap of one unit with every other is then one gather over the taxa of that unit, where a units x units array would
    have to be updated for every merge. The index is the transpose of the units x taxa matrix, with the unit of each entry
    changed as units merge. When two units hold a taxon and merge, one of the two entries is dropped, so a taxon has one entry per
    unit that holds it and the entries count shared taxa.
    """

    def __init__(self, x: csr_matrix):
        """`x` is units x taxa, binary."""
        by_taxon = x.T.tocsr()
        self.n_units = x.shape[0]
        self.indptr = by_taxon.indptr.astype(np.int64)  # where each taxon's entries start
        self.holder = by_taxon.indices.astype(np.int32)  # the unit of each entry, or `n_units` for one that was dropped
        self.taxa = [[x.indices[x.indptr[i] : x.indptr[i + 1]]] for i in range(self.n_units)]  # per unit, in pieces that merging adds
        self.richness = np.diff(x.indptr).astype(np.float64)
        self.active = np.ones(self.n_units, dtype=bool)

    def absorb(self, unit: int) -> int:
        """Merge `unit` into the active unit with the lowest beta-sim to it, and return that unit's position.

        Beta-sim is `1 - shared / min(richness)` in float64, and 1.0 where the poorer unit has no taxa. Ties go to the lowest position.
        """
        self.active[unit] = False  # so that neither it nor a unit merged before can be the partner
        taxa = self.taxa[unit][0] if len(self.taxa[unit]) == 1 else np.concatenate(self.taxa[unit])
        starts = self.indptr[taxa]
        lengths = self.indptr[taxa + 1] - starts
        ends = np.cumsum(lengths)
        entries = np.repeat(starts - (ends - lengths), lengths)  # every entry of every taxon of the unit, as positions in `holder`
        entries += np.arange(len(entries))
        holders = self.holder[entries]
        shared = np.bincount(holders, minlength=self.n_units + 1)[: self.n_units].astype(np.float64)
        poorer = np.minimum(self.richness[unit], self.richness)
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(poorer > 0, 1.0 - shared / poorer, 1.0)
        partner = int(np.argmin(np.where(self.active, beta, np.inf)))

        # the unit's entry for each of its taxa moves to the partner, or is dropped where the partner holds the taxon already
        own = entries[holders == unit]  # one for each taxon, in the order of `taxa`
        both = np.zeros(len(taxa), dtype=bool)
        both[np.searchsorted(ends, np.flatnonzero(holders == partner), side="right")] = True
        self.holder[own[both]] = self.n_units
        self.holder[own[~both]] = partner
        self.taxa[partner].append(taxa[~both])
        self.taxa[unit] = []
        self.richness[partner] += len(taxa) - int(both.sum())
        self.richness[unit] = 0.0
        return partner


def apply_size_floor(
    labels: Labels, areas: dict[str, float], ct: CellTaxa, floor_km2: float = 1_000_000
) -> tuple[Labels, dict[int, int]]:
    """Merge units under the floor into the realm they share the most taxa with, smallest first.

    Returns the new labels and `{sub_realm: parent}`, where the parent is the realm the unit ended up in. The smallest area goes
    first and a tie in area goes to the lowest label, and the realm to merge into is the one with the lowest beta-sim, a tie going to
    the lowest label. It stops when no unit is under the floor or one unit is left.
    """
    units, x = _unit_matrix(labels, ct)
    index = _TaxaIndex(x)
    area = {u: 0.0 for u in units}
    for cell, unit in labels.items():
        area[unit] += areas[cell]
    position = {u: i for i, u in enumerate(units)}
    small = [(a, u) for u, a in area.items() if a < floor_km2]  # a unit is listed again when it grows, and the old entry is stale
    heapify(small)
    absorbed_by: dict[int, int] = {}

    while len(units) - len(absorbed_by) > 1 and small:
        size, unit = heappop(small)
        if unit in absorbed_by or area[unit] != size:
            continue
        target = units[index.absorb(position[unit])]
        area[target] += area.pop(unit)
        absorbed_by[unit] = target
        if area[target] < floor_km2:
            heappush(small, (area[target], target))

    def final(unit: int) -> int:
        chain = []
        while unit in absorbed_by:
            chain.append(unit)
            unit = absorbed_by[unit]
        for passed in chain:
            absorbed_by[passed] = unit
        return unit

    parents = {unit: final(unit) for unit in list(absorbed_by)}
    return {c: parents.get(r, r) for c, r in labels.items()}, parents


def candidate_maps(
    partition: Partition,
    ct: CellTaxa,
    areas: dict[str, float],
    adjacency: dict[str, set[str]],
    k_max: int = 60,
    k_min: int = 5,
    floor_km2: float = 1_000_000,
) -> CandidateMaps:
    """Finest modules, then separate pieces, then the size floor, then the merge sequence.

    Flooring first keeps exactly one map per realm count, and a merge only coarsens, so pieces need resolving just once. Pass
    the matrix the network was built from, so that pieces and merges see the same evolutionary depth.
    """
    units = resolve_pieces(partition.leaf_labels(), adjacency, ct)
    floored, sub_realms = apply_size_floor(units, areas, ct, floor_km2)
    if len(set(floored.values())) < k_min:
        raise DegenerateClustering(f"{len(set(floored.values()))} units above the size floor, need at least {k_min}")
    return CandidateMaps(maps=merge_sequence(floored, ct, k_max, k_min), sub_realms=sub_realms, units=units)
