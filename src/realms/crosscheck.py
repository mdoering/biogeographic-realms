"""Comparing the headline Infomap map with other groupings of the same units, for the cross-checks (`scripts/crosscheck.py`).

The beta-sim consensus groups the network's units, so it is compared unit by unit: `unit_labels` reads the realm of each unit from `parts.parquet`,
`pairwise_v` scores groupings against each other, and `area_agreement` says how much of the area two maps put in matching realms.
"""

import itertools
from collections.abc import Mapping

import numpy as np
import polars as pl
from scipy.optimize import linear_sum_assignment

from realms.score import v_measure


def unit_labels(parts: pl.DataFrame, k: int) -> dict[str, int]:
    """Unit to realm at one k of a `parts.parquet` table, from the parts that are not imputed (each holds the unit it is in).

    The parts of a unit share its realm, so a unit with parts in two realms raises, as does a part that is not imputed and has no unit.
    """
    rows = parts.filter((pl.col("k") == k) & ~pl.col("imputed"))
    if not rows.height:
        raise ValueError(f"the parts hold no map at k={k}; they have k={sorted(parts['k'].unique())}")
    if rows["unit"].null_count():
        raise ValueError(f"{rows['unit'].null_count()} parts at k={k} are not imputed and have no unit")
    realms = rows.group_by("unit").agg(pl.col("realm").unique())
    split = realms.filter(pl.col("realm").list.len() > 1)
    if split.height:
        raise ValueError(f"units with parts in more than one realm at k={k}, such as {split['unit'].head(5).to_list()}")
    return dict(realms.select("unit", pl.col("realm").list.first()).iter_rows())


def pairwise_v(maps: Mapping) -> list[dict]:
    """`a`, `b` and `v_measure` (over the keys both hold) for every pair of the named groupings in `maps`, in their order."""
    return [{"a": a, "b": b, "v_measure": v_measure(dict(maps[a]), dict(maps[b]))} for a, b in itertools.combinations(maps, 2)]


def area_agreement(a: Mapping[str, int], b: Mapping[str, int], areas: Mapping[str, float]) -> float:
    """The share of the area of the units both maps hold that lies in realms matched one to one between them.

    The realms are matched to maximise the area they share (the Hungarian method on the table of shared area), so a realm of one map that
    has no partner in the other, when they have different numbers of realms, counts as disagreement on all its area.
    """
    shared = sorted(a.keys() & b.keys())
    missing = [u for u in shared if u not in areas]
    if missing:
        raise ValueError(f"{len(missing)} units have no area, such as {missing[:5]}")
    rows = {r: i for i, r in enumerate(sorted({a[u] for u in shared}))}
    columns = {r: i for i, r in enumerate(sorted({b[u] for u in shared}))}
    table = np.zeros((len(rows), len(columns)))
    for u in shared:
        table[rows[a[u]], columns[b[u]]] += areas[u]
    left, right = linear_sum_assignment(table, maximize=True)
    return float(table[left, right].sum() / table.sum())
