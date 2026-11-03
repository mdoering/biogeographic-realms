"""Stability over seeds (DECISIONS, "Stability over five seeds", 2026-11-03).

Twenty bootstraps give a stability that moved by 0.054 between two seeds of the same configuration, more than the margins it was chosen on.
So each configuration is selected at several seeds, and `across_seeds` takes the median of their `selection.csv` tables per k: a seed with no
map at a k counts as stability 0 there, and the excess is the median over the seeds that have the map. `choose` then applies
`select.select_level`'s rule to the pooled maps of every configuration.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import polars as pl

from realms.select import MIN_STABILITY, TIE


def across_seeds(tables: Mapping[int, pl.DataFrame]) -> pl.DataFrame:
    """`k`, `stability` (the median over the seeds, 0 where a seed has no map), `excess` (the median over the seeds that have the map), `seeds`
    (how many have it), and `stability_min` and `stability_max` (with the zeros), from each seed's `selection.csv` table, sorted by k."""
    if not tables:
        raise ValueError("no seed to take the median over")
    long = pl.concat([t.select("k", "stability", "excess").with_columns(pl.lit(seed).alias("seed")) for seed, t in tables.items()])
    ks = long.select("k").unique()
    full = ks.join(pl.DataFrame({"seed": list(tables)}), how="cross").join(long, on=["k", "seed"], how="left")
    return (
        full.group_by("k")
        .agg(
            pl.col("stability").fill_null(0.0).median().alias("stability"),
            pl.col("excess").drop_nulls().median().alias("excess"),
            pl.col("excess").is_not_null().sum().cast(pl.Int64).alias("seeds"),
            pl.col("stability").fill_null(0.0).min().alias("stability_min"),
            pl.col("stability").fill_null(0.0).max().alias("stability_max"),
        )
        .sort("k")
    )


@dataclass(frozen=True)
class Choice:
    config: str
    k: int
    fell_short: bool
    reason: str


def choose(summaries: Mapping[str, pl.DataFrame]) -> Choice:
    """Of every configuration's maps (`across_seeds` tables), those with stability at least `MIN_STABILITY`: the highest excess, the coarser map
    within `TIE`, and on a tie in k the higher excess. With no stable map, the most stable one (then the coarser, then the higher excess), flagged."""
    pooled = pl.concat([s.select("k", "stability", "excess").with_columns(pl.lit(config).alias("config")) for config, s in summaries.items()])
    eligible = pooled.filter(pl.col("stability") >= MIN_STABILITY)
    if eligible.height == 0:
        best = pooled.sort(["stability", "k", "excess"], descending=[True, False, True]).row(0, named=True)
        return Choice(best["config"], int(best["k"]), True,
                      f"no map reached stability {MIN_STABILITY}; {best['config']} at k={best['k']} is the most stable at {best['stability']:.3f}")
    top = eligible["excess"].max()
    best = eligible.filter(pl.col("excess") >= top - TIE).sort(["k", "excess"], descending=[False, True]).row(0, named=True)
    return Choice(best["config"], int(best["k"]), False,
                  f"{best['config']} at k={best['k']} has excess within {TIE} of the highest, {top:.3f}, among {eligible.height} stable maps")
