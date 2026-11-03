"""Stability over seeds (DECISIONS, "Stability over five seeds", 2026-11-03).

Twenty bootstraps give a stability that moved by 0.054 between two seeds of the same configuration, more than the margins it was chosen on.
So each configuration is selected at several seeds, and `across_seeds` takes the median of their `selection.csv` tables per k: a seed with no
map at a k counts as stability 0 there, and the excess is the median over the seeds that have the map. `choose` then applies
`select.select_level`'s rule to the pooled maps of every configuration.

`scripts/seeds.py` writes that choice to `<run>/seeds/choice.json`. `read_headline` is `explore.read_selected` with the choice applied: where it names
the configuration, its k (and whether it fell short) replace those the select stage chose at seed 1, so the gate, the cross-checks and the review
draw the map chosen over seeds.
"""

import json
from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path

import polars as pl

from realms.explore import Selected, read_selected
from realms.select import MIN_STABILITY, TIE

CHOICE = Path("seeds") / "choice.json"


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


def apply_choice(found: Selected, choice: Mapping | None, config: str) -> Selected:
    """`found` with the k and `fell_short` of `choice` (a `choice.json`) if it names `config`, its parts' `selected` column marking that k; else `found`."""
    if not choice or choice["config"] != config:
        return found
    k = int(choice["k"])
    if k not in set(found.parts["k"].unique().to_list()):
        raise ValueError(f"the choice over seeds is k={k} of {config}, but {found.directory} has maps at k={sorted(found.parts['k'].unique())}")
    return replace(found, k=k, fell_short=bool(choice["fell_short"]), parts=found.parts.with_columns((pl.col("k") == k).alias("selected")))


def read_choice(run_dir: Path) -> dict | None:
    """The `choice.json` that `scripts/seeds.py` wrote for the run `run_dir`, or None if it has none."""
    path = Path(run_dir) / CHOICE
    return json.loads(path.read_text()) if path.exists() else None


def read_headline(run_dir: Path, config: str, configs_path: Path, **kwargs) -> Selected:
    """`explore.read_selected` at seed 1 with the run's choice over seeds applied (`apply_choice`)."""
    return apply_choice(read_selected(run_dir, config, configs_path, **kwargs), read_choice(run_dir), config)
