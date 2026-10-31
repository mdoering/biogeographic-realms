"""Random maps with the same shape as a real one, to judge how much endemism a map shows beyond chance."""

import heapq

import numpy as np
import polars as pl

from realms.candidates import connected_pieces
from realms.endemism import Labels, area_weighted, incidence, realm_share

MAX_SWEEPS = 200


def _one_partition(
    pieces: list[tuple[int, float]],
    cells: list[str],
    adjacency: dict[str, set[str]],
    areas: dict[str, float],
    rng: np.random.Generator,
) -> Labels:
    assigned: dict[str, int] = {}  # cell -> index into pieces
    realm_of = [realm for realm, _ in pieces]
    remaining = [area for _, area in pieces]
    frontiers: list[list[str]] = [[] for _ in pieces]

    def touches_own_realm(cell: str, piece: int) -> bool:
        """Would this cell join the piece to another piece of the same realm?"""
        return any(realm_of[assigned[n]] == realm_of[piece] and assigned[n] != piece for n in adjacency.get(cell, ()) if n in assigned)

    def place(cell: str, piece: int) -> None:
        assigned[cell] = piece
        remaining[piece] -= areas[cell]
        frontiers[piece].extend(n for n in sorted(adjacency.get(cell, ())) if n not in assigned)

    shuffled = [cells[i] for i in rng.permutation(len(cells))]
    for piece in range(len(pieces)):
        free = (c for c in shuffled if c not in assigned)
        seed = next((c for c in free if not touches_own_realm(c, piece)), None)
        if seed is None:
            seed = next((c for c in shuffled if c not in assigned), None)
        if seed is not None:
            place(seed, piece)

    # all pieces grow together, always the one furthest below its area, so regions rarely strand each other
    targets = [area for _, area in pieces]
    heap = [(-remaining[p] / targets[p], float(rng.random()), p) for p in range(len(pieces))]
    heapq.heapify(heap)
    deferred: list[tuple[str, int]] = []
    while heap:
        _, _, piece = heapq.heappop(heap)
        frontier = frontiers[piece]
        while frontier:
            cell = frontier.pop(int(rng.integers(len(frontier))))
            if cell in assigned:
                continue
            if touches_own_realm(cell, piece):
                deferred.append((cell, piece))
                continue
            place(cell, piece)
            heapq.heappush(heap, (-remaining[piece] / targets[piece], float(rng.random()), piece))
            break

    # what's left would have merged two pieces of a realm, or nothing touches it: give it to whoever needs area most
    leftovers = [c for c in cells if c not in assigned]
    while leftovers:
        progress = []
        for cell in leftovers:
            options = {assigned[n] for n in adjacency.get(cell, ()) if n in assigned}
            if options:
                place(cell, max(sorted(options), key=lambda p: remaining[p]))
                progress.append(cell)
        if not progress:
            place(leftovers[0], max(range(len(pieces)), key=lambda p: remaining[p]))
            progress.append(leftovers[0])
        leftovers = [c for c in leftovers if c not in assigned]

    _rebalance(assigned, remaining, realm_of, adjacency, areas, shuffled)
    return {cell: realm_of[piece] for cell, piece in assigned.items()}


def _stays_connected(cell: str, piece: int, assigned: dict[str, int], adjacency: dict[str, set[str]]) -> bool:
    """A cheap local test that removing `cell` won't split its piece: its neighbours in the piece must connect within two steps.

    It can refuse a safe move, but a move it allows never cuts the piece.
    """
    own = [n for n in adjacency.get(cell, ()) if assigned[n] == piece]
    if len(own) <= 1:
        return bool(own)
    reached = {own[0]}
    stack = [own[0]]
    wanted = set(own)
    near = wanted | {m for n in own for m in adjacency.get(n, ()) if assigned[m] == piece}
    near.discard(cell)
    while stack:
        for m in adjacency.get(stack.pop(), ()):
            if m in near and m not in reached:
                reached.add(m)
                stack.append(m)
    return wanted <= reached


def _rebalance(
    assigned: dict[str, int],
    remaining: list[float],
    realm_of: list[int],
    adjacency: dict[str, set[str]],
    areas: dict[str, float],
    cells: list[str],
) -> None:
    """Pass border cells from pieces above their area to neighbours below it, until no such move improves the balance.

    Area diffuses through neighbouring pieces, so a walled-in piece recovers what it lost even with no direct route.
    """
    for _ in range(MAX_SWEEPS):
        moved = False
        for cell in cells:
            donor = assigned[cell]
            options = {assigned[n] for n in adjacency.get(cell, ())} - {donor}
            if not options:
                continue
            taker = max(sorted(options), key=lambda p: remaining[p])
            if remaining[taker] - remaining[donor] <= areas[cell]:
                continue
            joins_own_realm = any(
                realm_of[assigned[n]] == realm_of[taker] and assigned[n] != taker for n in adjacency.get(cell, ())
            )
            if joins_own_realm or not _stays_connected(cell, donor, assigned, adjacency):
                continue
            assigned[cell] = taker
            remaining[donor] += areas[cell]
            remaining[taker] -= areas[cell]
            moved = True
        if not moved:
            return


def null_partitions(
    target: Labels,
    adjacency: dict[str, set[str]],
    areas: dict[str, float],
    n: int = 99,
    seed: int = 1,
) -> list[Labels]:
    """Random maps matching the target's realm count, realm areas and number of pieces.

    Each target piece becomes a region grown from a random seed cell, all at once, always extending the piece furthest below its
    area. Leftovers join a touching piece, so a realm can end a little off its area, and two of its pieces can occasionally
    touch and become one.
    """
    rng = np.random.default_rng(seed)
    cells = sorted(target)
    realms = sorted(set(target.values()))
    pieces: list[tuple[int, float]] = []
    for realm in realms:
        members = {c for c in cells if target[c] == realm}
        pieces += [(realm, sum(areas[c] for c in piece)) for piece in connected_pieces(members, adjacency)]
    pieces.sort(key=lambda p: (-p[1], p[0]))
    return [_one_partition(pieces, cells, adjacency, areas, rng) for _ in range(n)]


def endemism_excess(
    labels: Labels,
    presence_by_rank: dict[str, pl.DataFrame],
    adjacency: dict[str, set[str]],
    areas: dict[str, float],
    n_null: int = 99,
    seed: int = 1,
    tolerance: float = 0.95,
) -> float:
    """Area-weighted endemism above the mean of the null maps, averaged equally across ranks."""
    cells = sorted(labels)
    cell_area = np.array([areas[c] for c in cells])
    realms = sorted(set(labels.values()))
    position = {r: i for i, r in enumerate(realms)}
    nulls = null_partitions(labels, adjacency, areas, n=n_null, seed=seed)
    observed_realm = np.array([position[labels[c]] for c in cells])
    null_realms = [np.array([position[null[c]] for c in cells]) for null in nulls]

    def score(x, cell_realm: np.ndarray) -> float:
        return area_weighted(realm_share(x, cell_realm, len(realms), tolerance), cell_realm, cell_area)

    excess = []
    for presence in presence_by_rank.values():
        x = incidence(presence, cells)
        excess.append(score(x, observed_realm) - float(np.mean([score(x, realm) for realm in null_realms])))
    return float(np.mean(excess))
