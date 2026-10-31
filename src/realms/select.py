"""Choosing the headline map by rules fixed in advance."""

from dataclasses import dataclass

import numpy as np
import polars as pl

from realms.candidates import DegenerateClustering, candidate_maps
from realms.cluster import run_infomap
from realms.endemism import Labels
from realms.matrix import CellTaxa
from realms.network import network_from
from realms.parallel import map_tasks
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


def _nearest(maps: dict[int, Labels], k: int) -> Labels:
    return maps[min(maps, key=lambda available: (abs(available - k), available))]


def _rerun(shared: tuple, task: tuple[int, np.ndarray]) -> dict[int, float]:
    """The V-measure of every map against the pipeline rerun on the species of one draw."""
    ct, maps, areas, adjacency, seed, trials, floor_km2, k_max, markov_time = shared
    b, draws = task
    resampled = ct.resample(draws)
    partition = run_infomap(network_from(resampled), seed=seed + b, trials=trials, markov_time=markov_time)
    try:
        candidates = candidate_maps(partition, resampled, areas, adjacency, k_max=k_max, k_min=2, floor_km2=floor_km2).maps
    except DegenerateClustering:
        candidates = {}
    return {k: v_measure(labels, _nearest(candidates, k)) if candidates else 0.0 for k, labels in maps.items()}


def bootstrap_stability(
    ct: CellTaxa,
    maps: dict[int, Labels],
    areas: dict[str, float],
    adjacency: dict[str, set[str]],
    n_boot: int = 20,
    seed: int = 1,
    trials: int = 3,
    floor_km2: float = 1_000_000,
    k_max: int = 60,
    markov_time: float = 1.0,
    workers: int = 1,
) -> dict[int, float]:
    """Median V-measure per realm count between each map and the same pipeline rerun on resampled species.

    Each bootstrap draws species with replacement from those of `ct`, each draw its own species so that one drawn twice counts
    twice. The ranks are those of `ct`, and they, `floor_km2`, `k_max` and `markov_time` have to be the ones that produced `maps`,
    or this measures some other pipeline. A bootstrap that lacks the exact count uses its nearest. One that collapses completely
    scores zero.

    All the draws are made up front, in order, from one generator, and each bootstrap's Infomap seed is `seed` plus its number.
    So the result is the same for any number of `workers`, which run the bootstraps in separate processes.
    """
    rng = np.random.default_rng(seed)
    n_species = len(ct.species_ids())
    if not n_species:
        raise ValueError("bootstrap_stability resamples species, so the matrix needs the species rank")
    draws = [rng.integers(n_species, size=n_species) for _ in range(n_boot)]
    shared = (ct, maps, areas, adjacency, seed, trials, floor_km2, k_max, markov_time)
    scores = map_tasks(_rerun, list(enumerate(draws)), shared, workers)
    return {k: float(np.median([score[k] for score in scores])) for k in maps}
