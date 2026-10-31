"""Choosing the headline map by rules fixed in advance."""

from dataclasses import dataclass

import numpy as np
import polars as pl

from realms.candidates import DegenerateClustering, candidate_maps
from realms.cluster import run_infomap
from realms.endemism import Labels
from realms.network import build_network
from realms.score import v_measure

MIN_STABILITY = 0.8
TIE = 0.01


@dataclass(frozen=True)
class Selection:
    k: int
    fell_short: bool
    reason: str


def select_level(table: pl.DataFrame) -> Selection:
    """Highest excess among stable maps, the coarser map on a tie. With no stable map, the most stable one, flagged."""
    eligible = table.filter(pl.col("stability") >= MIN_STABILITY)
    if eligible.height == 0:
        best = table.sort(["stability", "k"], descending=[True, False]).row(0, named=True)
        return Selection(int(best["k"]), True, f"no map reached stability {MIN_STABILITY}, k={best['k']} is the most stable at {best['stability']:.3f}")
    top = eligible["excess"].max()
    k = eligible.filter(pl.col("excess") >= top - TIE)["k"].min()
    return Selection(int(k), False, f"k={k} has excess within {TIE} of the highest, {top:.3f}, among {eligible.height} stable maps")


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


def _nearest(maps: dict[int, Labels], k: int) -> Labels:
    return maps[min(maps, key=lambda available: (abs(available - k), available))]


def bootstrap_stability(
    presence: pl.DataFrame,
    taxa: pl.DataFrame | None,
    maps: dict[int, Labels],
    areas: dict[str, float],
    adjacency: dict[str, set[str]],
    n_boot: int = 20,
    seed: int = 1,
    trials: int = 3,
) -> dict[int, float]:
    """Median V-measure per realm count between each map and the same pipeline rerun on resampled species.

    A bootstrap that lacks the exact count uses its nearest. One that collapses completely scores zero.
    """
    rng = np.random.default_rng(seed)
    scores: dict[int, list[float]] = {k: [] for k in maps}
    for b in range(n_boot):
        new_presence, new_taxa = _resample(presence, taxa, rng)
        partition = run_infomap(build_network(new_presence, new_taxa), seed=seed + b, trials=trials)
        try:
            candidates = candidate_maps(partition, new_presence, areas, adjacency, k_min=2).maps
        except DegenerateClustering:
            candidates = {}
        for k, labels in maps.items():
            scores[k].append(v_measure(labels, _nearest(candidates, k)) if candidates else 0.0)
    return {k: float(np.median(values)) for k, values in scores.items()}
