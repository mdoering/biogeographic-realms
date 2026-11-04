"""Does the data single out a realm count? Two diagnostics of a nested sequence of maps (DECISIONS, "Is there a natural k?", 2026-11-04).

They choose nothing: `select.select_level` and `seeds.choose` still do, and this module changes neither them nor `candidates`.

- `merge_heights`: the merge sequence of `candidates.merge_sequence` (the same `_Merger`, so the same maps), with the beta-sim of the pair merged
  at each step. A break shows as a jump: merging from k to k - 1 joins units far less alike than the merges before it.
- `split_tests`: for each step from k - 1 to k realms, the endemism of the split within the realm it divides, against random splits of that realm
  into contiguous regions of the same areas and number of pieces (`nullmodel.null_partitions`). The endemism is `nullmodel`'s: the share of a
  side's taxa with at least `tolerance` of all their cells there, averaged equally over the ranks, with the two sides weighted by area or equally. A boundary the clustering drew from the
  data will usually beat random ones, so a split that does not is the telling case.
"""

import numpy as np
import polars as pl

from realms.candidates import _Merger
from realms.endemism import Labels, area_weighted, realm_share
from realms.matrix import CellTaxa
from realms.nullmodel import _rank_tables, null_partitions
from realms.parallel import map_tasks


def merge_heights(labels: Labels, ct: CellTaxa, k_min: int = 2) -> tuple[pl.DataFrame, dict[int, Labels]]:
    """Merge `labels` down to `k_min` units as `candidates.merge_sequence` does, and return one row per merge with the maps at every count.

    A row has `k`, the count before the merge, so it leads to the map with `k - 1` units; `height`, the beta-sim of the pair merged; `keep` and
    `drop`, the two units' labels (`drop`'s cells take `keep`'s label); and the taxa of each (`taxa_keep`, `taxa_drop`).
    """
    merger = _Merger(labels, ct)
    current = dict(labels)
    k = len(merger.units)
    maps: dict[int, Labels] = {k: dict(current)}
    rows = []
    while k > k_min:
        beta = merger.beta_sim()
        richness = np.diag(merger.shared)
        keep, drop = merger.closest_pair()
        keep_label, drop_label = merger.units[keep], merger.units[drop]
        rows.append({"k": k, "height": float(beta[keep, drop]), "keep": keep_label, "drop": drop_label,
                     "taxa_keep": int(richness[keep]), "taxa_drop": int(richness[drop])})
        current = {c: (keep_label if r == drop_label else r) for c, r in current.items()}
        merger.merge(keep, drop)
        k -= 1
        maps[k] = dict(current)
    schema = {"k": pl.Int64, "height": pl.Float64, "keep": pl.Int64, "drop": pl.Int64, "taxa_keep": pl.Int64, "taxa_drop": pl.Int64}
    return pl.DataFrame(rows, schema=schema), maps


def split_of(coarse: Labels, fine: Labels) -> tuple[int, list[str], Labels]:
    """The realm of `coarse` that `fine` divides in two, its cells (sorted), and their labels in `fine`. `fine` must be `coarse` with one realm split."""
    if coarse.keys() != fine.keys():
        raise ValueError("the two maps cover different cells")
    halves: dict[int, set[int]] = {}
    for cell, realm in coarse.items():
        halves.setdefault(realm, set()).add(fine[cell])
    divided = [realm for realm, labels in halves.items() if len(labels) > 1]
    if len(divided) != 1 or len(halves[divided[0]]) != 2 or len(set(fine.values())) != len(set(coarse.values())) + 1:
        raise ValueError("the finer map is not the coarser one with one realm split in two")
    parent = divided[0]
    cells = sorted(c for c, r in coarse.items() if r == parent)
    return parent, cells, {c: fine[c] for c in cells}


SCORES = ("area_weighted", "equal")


def _split_scores(tables, cell_realm: np.ndarray, cell_area: np.ndarray) -> tuple[float, float]:
    """The split's endemism averaged over the ranks: the sides' shares weighted by area, and with equal weight."""
    shares = [realm_share(table, cell_realm, 2, tolerance) for table, tolerance in tables]
    return float(np.mean([area_weighted(s, cell_realm, cell_area) for s in shares])), float(np.mean([s.mean() for s in shares]))


def _one_split(shared: tuple, task: tuple[int, Labels]) -> list[dict]:
    ct, adjacency, areas, n_null, seed, tolerance = shared
    k, sides = task
    cells = sorted(sides)
    inside = set(cells)
    local_adjacency = {c: {n for n in adjacency.get(c, ()) if n in inside} for c in cells}
    local_areas = {c: areas[c] for c in cells}
    order = sorted(set(sides.values()))
    position = {r: i for i, r in enumerate(order)}
    cell_area = np.array([areas[c] for c in cells])
    tables = [(t, tolerance) for t in _rank_tables(ct, cells)]  # the totals are the taxa's cells everywhere, so a taxon also found outside the realm is endemic to neither side
    observed = _split_scores(tables, np.array([position[sides[c]] for c in cells]), cell_area)
    nulls = null_partitions(sides, local_adjacency, local_areas, n=n_null, seed=seed)
    null_scores = np.array([_split_scores(tables, np.array([position[null[c]] for c in cells]), cell_area) for null in nulls]).reshape(-1, 2)
    side_area = [float(sum(areas[c] for c in cells if sides[c] == r)) for r in order]
    rows = []
    for i, score in enumerate(SCORES):
        column = null_scores[:, i]
        rows.append({
            "k": k, "score": score, "cells": len(cells), "area_small": min(side_area), "area_large": max(side_area),
            "observed": observed[i], "null_mean": float(column.mean()), "null_q95": float(np.quantile(column, 0.95)),
            "null_sd": float(column.std(ddof=1)) if n_null > 1 else float("nan"),
            "p": float((1 + int((column >= observed[i]).sum())) / (n_null + 1)),
        })
    return rows


def split_tests(
    maps: dict[int, Labels],
    ct: CellTaxa,
    adjacency: dict[str, set[str]],
    areas: dict[str, float],
    ks: list[int] | None = None,
    n_null: int = 99,
    seed: int = 1,
    tolerance: float = 0.95,
    workers: int = 1,
) -> pl.DataFrame:
    """For each k of `ks` (default every k with a map at k - 1), the split from `maps[k - 1]` to `maps[k]` within the realm it divides.

    Each split has a row per score of `SCORES`: the sides' endemic shares weighted by their area, as `nullmodel` weights realms, or with equal
    weight, since in a lopsided split the larger side's share is high for any boundary and swamps the area-weighted score.
    A row has `k` and `score`; the `parent` realm and its `cells`; the areas of the smaller and larger side; `observed`, the endemism of the split;
    `null_mean`, `null_q95` and `null_sd` of `n_null` random splits; `gain`, observed minus the null mean; `z`, the gain in null standard deviations;
    and `p`, the share of random splits at least as endemic (counting the observed one). The null maps of every split are drawn from `seed`.
    """
    ks = sorted(k for k in maps if k - 1 in maps) if ks is None else sorted(ks)
    tasks, parents = [], {}
    for k in ks:
        parent, _, sides = split_of(maps[k - 1], maps[k])
        parents[k] = parent
        tasks.append((k, sides))
    rows = [row for split in map_tasks(_one_split, tasks, (ct, adjacency, areas, n_null, seed, tolerance), workers) for row in split]
    table = pl.DataFrame(rows).with_columns(pl.col("k").replace_strict(parents, return_dtype=pl.Int64).alias("parent"))
    return table.with_columns((pl.col("observed") - pl.col("null_mean")).alias("gain")).with_columns((pl.col("gain") / pl.col("null_sd")).alias("z")).select(
        "k", "score", "parent", "cells", "area_small", "area_large", "observed", "null_mean", "null_q95", "null_sd", "gain", "z", "p"
    ).sort("k", "score")
