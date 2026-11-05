"""Nested levels of realms instead of one realm count (DECISIONS, "Nested levels", 2026-11-05).

The merge heights ("Is there a natural k?") single out no k, but a floor (five blocks about equally distinct) and an order in which the smaller
realms join. So the realms are reported at nested levels, and a unit is a realm by how distinct it is, not by the 1M km2 size floor:

- `leaf_realms`: Infomap's units, merged as `candidates.merge_sequence` merges them (the pair with the lowest beta-sim, on the union of their taxa),
  for as long as some pair is below `cut`. What is left are the finest realms. Small islands join the realm they share most with.
- `nested_merges`: the leaf realms merged on in the same way, but land only with land and sea only with sea (`sea_realms`, by area), the two joined last.
  Unconstrained, the vast ocean realm, whose coasts hold taxa from everywhere, takes in Madagascar, Australia and South and East Asia before any
  land realm does. Marine and terrestrial realms are two systems in most schemes, and the data put the sea's join with land at the top anyway.
- `at_cuts`: the map of each level, cut from the nested merges at a beta-sim. A cut never joins land and sea.

`candidates.py`, the selection and the Costello gate are unchanged: this is a separate path.
"""

import numpy as np
import polars as pl

from realms.candidates import _Merger
from realms.endemism import Labels
from realms.matrix import CellTaxa

SEA_SUFFIX = ":sea"
MERGES = {"k": pl.Int64, "height": pl.Float64, "keep": pl.Int64, "drop": pl.Int64, "kind": pl.String}


def _relabel(current: Labels, keep: int, drop: int) -> Labels:
    return {c: (keep if r == drop else r) for c, r in current.items()}


def leaf_realms(units: Labels, ct: CellTaxa, cut: float) -> tuple[Labels, pl.DataFrame]:
    """`units` merged as `merge_sequence` merges them while the closest pair's beta-sim is below `cut`, and the merges made (`MERGES`, `kind` null).

    A row's `k` is the count before its merge. The merged unit keeps the label of `keep`, as in `merge_sequence`.
    """
    merger = _Merger(units, ct)
    current = dict(units)
    rows = []
    while merger.active.sum() > 1:
        keep, drop = merger.closest_pair()
        height = float(merger.beta_sim()[keep, drop])
        if height >= cut:
            break
        rows.append({"k": int(merger.active.sum()), "height": height, "keep": merger.units[keep], "drop": merger.units[drop], "kind": None})
        current = _relabel(current, merger.units[keep], merger.units[drop])
        merger.merge(keep, drop)
    return current, pl.DataFrame(rows, schema=MERGES)


def sea_realms(labels: Labels, areas: dict[str, float], share: float = 0.5) -> dict[int, bool]:
    """Whether each realm is a sea realm: at least `share` of its area in sea parts (units whose id ends in `SEA_SUFFIX`)."""
    total: dict[int, float] = {}
    sea: dict[int, float] = {}
    for unit, realm in labels.items():
        total[realm] = total.get(realm, 0.0) + areas[unit]
        if unit.endswith(SEA_SUFFIX):
            sea[realm] = sea.get(realm, 0.0) + areas[unit]
    return {r: total[r] > 0 and sea.get(r, 0.0) / total[r] >= share for r in total}


def nested_merges(labels: Labels, ct: CellTaxa, sea: dict[int, bool]) -> tuple[pl.DataFrame, dict[int, Labels]]:
    """`labels` merged down to one realm, land only with land and sea only with sea until one of each is left, and the map at every count.

    Each step merges the pair of the same kind with the lowest beta-sim on the union of their taxa (ties to the lowest positions, as `_Merger`);
    the last merge joins land and sea. A row (`MERGES`) has the count before the merge, its height, the two labels and `kind`: "land", "sea" or "both".
    """
    merger = _Merger(labels, ct)
    kind = np.array([sea[u] for u in merger.units])
    current = dict(labels)
    maps = {int(merger.active.sum()): dict(current)}
    rows = []
    upper = np.triu(np.ones((len(merger.units),) * 2, dtype=bool), k=1)
    while merger.active.sum() > 1:
        beta = merger.beta_sim()
        allowed = np.outer(merger.active, merger.active) & upper
        same = allowed & np.equal.outer(kind, kind)
        pool = same if same.any() else allowed
        keep, drop = np.unravel_index(np.argmin(np.where(pool, beta, np.inf)), beta.shape)
        keep, drop = int(keep), int(drop)
        joined = "both" if kind[keep] != kind[drop] else ("sea" if kind[keep] else "land")
        k = int(merger.active.sum())
        rows.append({"k": k, "height": float(beta[keep, drop]), "keep": merger.units[keep], "drop": merger.units[drop], "kind": joined})
        current = _relabel(current, merger.units[keep], merger.units[drop])
        merger.merge(keep, drop)
        maps[k - 1] = dict(current)
    return pl.DataFrame(rows, schema=MERGES), maps


def at_cuts(merges: pl.DataFrame, maps: dict[int, Labels], cuts: list[float]) -> dict[float, Labels]:
    """The map of each cut: the realms before the first merge of `merges` at a height of at least the cut, or joining land and sea, which a cut never
    does (that join is the construction's, and its height says nothing). The merges are in the order made."""
    found = {}
    for cut in cuts:
        first = merges.filter((pl.col("height") >= cut) | (pl.col("kind") == "both")).head(1)
        k = int(first["k"].item()) if first.height else min(maps)
        found[cut] = maps[k]
    return found


def levels_of(units: Labels, ct: CellTaxa, areas: dict[str, float], leaf_cut: float, cuts: list[float]) -> tuple[dict[float, Labels], pl.DataFrame]:
    """The map of every level: the leaf realms at `leaf_cut`, and `at_cuts` of their nested merges at each of `cuts` (each above `leaf_cut`).
    Also the merges: the leaf merges, then the nested ones (`MERGES`, with a `stage` column)."""
    if any(c <= leaf_cut for c in cuts):
        raise ValueError(f"the cuts of the levels above the leaves must be above the leaf cut {leaf_cut}: {cuts}")
    leaves, leaf_merges = leaf_realms(units, ct, leaf_cut)
    merges, maps = nested_merges(leaves, ct, sea_realms(leaves, areas))
    levels = {leaf_cut: leaves} | at_cuts(merges, maps, cuts)
    table = pl.concat([leaf_merges.with_columns(pl.lit("leaf").alias("stage")), merges.with_columns(pl.lit("nested").alias("stage"))])
    return levels, table


def realm_recovery(reference: Labels, replicate: Labels, areas: dict[str, float]) -> dict[int, float]:
    """For each realm of `reference`, the highest Jaccard index by area with a realm of `replicate`, on the cells both maps have."""
    common = [c for c in reference if c in replicate]
    inside: dict[tuple[int, int], float] = {}
    ref_area: dict[int, float] = {}
    rep_area: dict[int, float] = {}
    for c in common:
        a, b, w = reference[c], replicate[c], areas[c]
        inside[(a, b)] = inside.get((a, b), 0.0) + w
        ref_area[a] = ref_area.get(a, 0.0) + w
        rep_area[b] = rep_area.get(b, 0.0) + w
    best = {r: 0.0 for r in set(reference.values())}
    for (a, b), w in inside.items():
        best[a] = max(best[a], w / (ref_area[a] + rep_area[b] - w))
    return best


def _boot(shared: tuple, task: tuple[int, np.ndarray]) -> list[dict]:
    from realms.candidates import resolve_pieces
    from realms.cluster import run_infomap
    from realms.network import network_from

    ct, areas, adjacency, reference, leaf_cut, cuts, options, seed, trials = shared
    b, draws = task
    resampled = ct.resample(draws)
    partition = run_infomap(network_from(resampled), seed=seed + b, trials=trials, options=options)
    units = resolve_pieces(partition.leaf_labels(), adjacency, resampled)
    levels, _ = levels_of(units, resampled, areas, leaf_cut, cuts)
    return [{"boot": b, "cut": cut, "realm": realm, "jaccard": j}
            for cut, labels in levels.items() for realm, j in realm_recovery(reference[cut], labels, areas).items()]


def bootstrap_recovery(
    ct: CellTaxa,
    areas: dict[str, float],
    adjacency: dict[str, set[str]],
    reference: dict[float, Labels],
    leaf_cut: float,
    cuts: list[float],
    options,
    n_boot: int = 20,
    seed: int = 1,
    trials: int = 3,
    workers: int = 1,
) -> pl.DataFrame:
    """How well each realm of each level of `reference` comes back when the whole path (Infomap, `resolve_pieces`, `levels_of`) runs on resampled species.

    The draws are those of `select.bootstrap_stability` (species with replacement, all drawn up front from one generator at `seed`, Infomap at `seed`
    plus the bootstrap's number), so the result does not depend on `workers`. A row: `boot`, `cut`, `realm` and its `jaccard` (`realm_recovery`).
    """
    from realms.parallel import map_tasks

    rng = np.random.default_rng(seed)
    n_species = len(ct.species_ids())
    if not n_species:
        raise ValueError("bootstrap_recovery resamples species, so the matrix needs the species rank")
    draws = [rng.integers(n_species, size=n_species) for _ in range(n_boot)]
    shared = (ct, areas, adjacency, reference, leaf_cut, cuts, options, seed, trials)
    rows = [row for rows in map_tasks(_boot, list(enumerate(draws)), shared, workers) for row in rows]
    return pl.DataFrame(rows, schema={"boot": pl.Int64, "cut": pl.Float64, "realm": pl.Int64, "jaccard": pl.Float64})
