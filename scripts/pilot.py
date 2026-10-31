"""Run the pipeline once on half A at resolution 6, weigh and time every stage, and decide from the numbers whether resolution 7 fits.

The stages, in order: load the taxa, the split, the counts and the grid, then `cell_taxa` and `network_from` on the parts of the grid
(land and sea part of a cell are separate nodes), one-trial Infomap at Markov times 1.0 and 0.75, and a few more (`SURVEY`) to see what lies
between them. Then `candidate_maps` on each, the smallest number of units first, with `resolve_pieces`, `apply_size_floor` and
`merge_sequence` timed apart. The first of the Markov times that is not degenerate goes on to a 10-trial Infomap, a bootstrap and
`excess_table`. Half A only: `load_presence` raises for half B, and nothing here asks it to.

Every stage's wall time and the peak memory (RSS, GB of 10^9 bytes) go to `results/pilot/res<r>/measure.json` and `pilot.log`, written as
the run goes. A process only reports a peak for its whole life, so Infomap, `candidate_maps`, the taxa load (which peaks at 12 to 15 GB),
and the parents of the bootstrap's and `excess_table`'s workers each run in a child process of their own (Infomap loads its links from a
`.npy`), and the workers record their own peaks. `null_partitions` divides by zero for a piece of parts with no area, which are kept: the
stage is then run again without them, only to time it and weigh its workers. `--project-only` projects again from a saved `measure.json`.

`project` carries the measures over to resolution 7 with the conservative growth factor of 2.5 from resolution 6 to 7, and with the 1.38
measured on Denmark, which is a lower bound because Denmark spans only 8 cells.

The counts are a download, so run it as `uv run python -I scripts/pilot.py --counts data/gbif/<key> --workers N`.
"""

import argparse
import json
import math
import multiprocessing
import os
import resource
import shutil
import sys
import time
import traceback
from collections import Counter
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack, contextmanager
from functools import partial
from pathlib import Path
from unittest import mock

import numpy as np
import polars as pl

from realms import candidates, cluster, nullmodel, parallel, select
from realms.candidates import DegenerateClustering, candidate_maps
from realms.cluster import Partition, run_infomap
from realms.grid import ZERO_KM2, load_grid, part_adjacency, parts
from realms.matrix import cell_taxa
from realms.network import Network, network_from
from realms.nullmodel import excess_table
from realms.occurrences import load_presence
from realms.reference import load_taxa
from realms.select import bootstrap_stability
from realms.split import load_split

SEED = 1
MARKOV_TIMES = (1.0, 0.75)  # the first that gives candidate maps is used from the 10-trial Infomap onwards
SURVEY = (0.97, 0.95, 0.9, 0.85, 0.8)  # between them: the module counts and what `candidate_maps` costs for them, and nothing after
K_MAX, K_MIN, FLOOR_KM2 = 60, 5, 1_000_000
BOOT_TRIALS, BOOT_PROBE = 3, 2  # a bootstrap of 3 trials, and 2 of them in 2 workers, to see a worker's peak
N_NULL = 99
HEADLINE_TRIALS = 10
TOP_MODULES = 40  # of the top-level modules, the largest ones are listed

# The projection to resolution 7.
HEADLINE = 2.0  # the headline holds all species, which is about twice the links of half A
GROWTHS = {"2.5": 2.5, "1.38": 1.38}  # links, resolution 6 to 7: the conservative ruling R3, and Denmark's 117,960 -> 163,091 pairs
CELLS_RATIO = 21_872 / 7_292  # ISEA3H cells at resolution 7 over resolution 6, `scripts/check_grid.py`
BUDGET_GB, FREE_GB, HOURS = 40, 10, 24  # the most that one process or all the workers may hold, the disk to leave, the time limit
BOOTSTRAPS = 20  # of a selection run, as `bootstrap_stability` has it by default
CORES = os.cpu_count() or 1


def peak_gb() -> float:
    """The most memory this process has had resident so far. `ru_maxrss` is in bytes on macOS and kilobytes on Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)


class Run:
    """The measurements and the log, written as the run goes, so that a failure late on keeps what was measured."""

    def __init__(self, directory: Path):
        directory.mkdir(parents=True, exist_ok=True)
        self.path = directory / "measure.json"
        self.log = (directory / "pilot.log").open("w")
        self.data: dict = {"stages": {}}

    def say(self, text: str) -> None:
        line = f"{time.strftime('%H:%M:%S')} {text}"
        print(line, flush=True)
        self.log.write(line + "\n")
        self.log.flush()

    def save(self) -> None:
        self.path.write_text(json.dumps(self.data, indent=2) + "\n")

    @contextmanager
    def stage(self, name: str):
        """Yields the stage's entry to add figures to. `peak_rss_gb` is this process's peak so far, at the end of the stage."""
        entry: dict = {}
        started = time.perf_counter()
        yield entry
        entry = {"seconds": round(time.perf_counter() - started, 1), **entry, "peak_rss_gb": round(peak_gb(), 2)}
        self.data["stages"][name] = entry
        self.say(f"{name}: " + ", ".join(f"{key} {value}" for key, value in entry.items() if not isinstance(value, (dict, list))))
        self.save()


def in_child(function: Callable, *args):
    """`function(*args)` in a child process of its own, which is spawned, so nothing is inherited and its peak memory is its own."""
    with ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("spawn")) as pool:
        return pool.submit(function, *args).result()


def _load_taxa(path: Path) -> tuple[pl.DataFrame, float]:
    return load_taxa(path), peak_gb()


def _infomap(links_path: Path, cells: list[str], seed: int, trials: int, markov_time: float) -> dict:
    """In a child: `run_infomap` on the links of a `.npy` of `cell_idx, taxon_idx, weight`, as the parent would run it on its own network."""
    links = np.load(links_path)
    net = Network(cells, [], links[:, 0].astype(np.int64), links[:, 1].astype(np.int64), links[:, 2].copy())
    del links
    loaded = {"peak_before_gb": peak_gb()}
    real_load = cluster.load_infomap

    def timed_load(*args, **kwargs):
        started = time.perf_counter()
        infomap = real_load(*args, **kwargs)
        loaded["load_seconds"], loaded["peak_after_load_gb"] = time.perf_counter() - started, peak_gb()
        return infomap

    started = time.perf_counter()
    with mock.patch.object(cluster, "load_infomap", timed_load):
        partition = run_infomap(net, seed=seed, trials=trials, markov_time=markov_time)
    return {**loaded, "seconds": time.perf_counter() - started, "peak_gb": peak_gb(), "paths": partition.paths, "unassigned": partition.unassigned}


@contextmanager
def timed_candidates(timings: dict[str, dict]):
    """Within it, the three stages of `candidate_maps` time themselves and note how many units each is given."""

    def timed(name: str, function: Callable) -> Callable:
        def wrapper(labels, *args, **kwargs):
            started = time.perf_counter()
            result = function(labels, *args, **kwargs)
            timings[name] = {"seconds": round(time.perf_counter() - started, 2), "units_in": len(set(labels.values()))}
            return result

        return wrapper

    with ExitStack() as stack:
        for name in ("resolve_pieces", "apply_size_floor", "merge_sequence"):
            stack.enter_context(mock.patch.object(candidates, name, timed(name, getattr(candidates, name))))
        yield


def _candidates(partition: Partition, ct, areas: dict[str, float], adjacency: dict[str, set[str]]) -> dict:
    """In a child: `candidate_maps`, with its three stages timed apart. A degenerate clustering comes back as its message."""
    timings: dict[str, dict] = {}
    started = time.perf_counter()
    try:
        with timed_candidates(timings):
            found = candidate_maps(partition, ct, areas, adjacency, k_max=K_MAX, k_min=K_MIN, floor_km2=FLOOR_KM2)
    except DegenerateClustering as error:
        return {"degenerate": str(error), "seconds": time.perf_counter() - started, "peak_gb": peak_gb(), "timings": timings}
    return {"found": found, "seconds": time.perf_counter() - started, "peak_gb": peak_gb(), "timings": timings}


def _bootstrap(ct, maps: dict, areas: dict, adjacency: dict, markov_time: float) -> dict:
    """In a child, which is clean: `bootstrap_stability` in `BOOT_PROBE` workers, with the peak and time of each, and this process's peak."""
    started = time.perf_counter()
    with worker_peaks(select) as seen:
        scores = bootstrap_stability(
            ct, maps, areas, adjacency, n_boot=BOOT_PROBE, seed=SEED, trials=BOOT_TRIALS, floor_km2=FLOOR_KM2, k_max=K_MAX,
            markov_time=markov_time, workers=BOOT_PROBE,
        )
    return {"scores": scores, "workers": seen, "seconds": time.perf_counter() - started, "peak_gb": peak_gb()}


def _excess(maps: dict, ct, adjacency: dict, areas: dict, workers: int) -> dict:
    """In a child, which is clean: `excess_table` in `workers` workers, with the peak and time of each, and this process's peak.

    A division by zero comes back as its message instead: `null_partitions` divides by the area of each piece of a realm, and a piece
    of parts that have no area has none.
    """
    started = time.perf_counter()
    try:
        with worker_peaks(nullmodel) as seen:
            excess = excess_table(maps, ct, adjacency, areas, n_null=N_NULL, seed=SEED, workers=workers)
    except ZeroDivisionError as error:
        return {"error": f"ZeroDivisionError: {error}", "seconds": time.perf_counter() - started, "peak_gb": peak_gb()}
    return {"excess": excess, "workers": seen, "seconds": time.perf_counter() - started, "peak_gb": peak_gb()}


def structure(partition: Partition, part_info: pl.DataFrame) -> dict:
    """What Infomap found: how many modules at each level, the share of the largest, and the top-level modules one by one.

    `part_info` is `part, kind, area_km2, lat, lon` for every part. The latitude and longitude of a module are the medians of its parts'
    centres, which is only a pointer to where it lies, as a module can span the antimeridian. The `TOP_MODULES` largest are listed.
    """
    paths = partition.paths
    top = Counter(path[0] for path in paths.values())
    leaf = Counter(paths.values())
    table = (
        pl.DataFrame({"part": list(paths), "module": [path[0] for path in paths.values()], "leaf": [str(path) for path in paths.values()]})
        .join(part_info, on="part")
        .with_columns((pl.col("kind") == "sea").alias("sea"))
    )
    by_module = (
        table.group_by("module")
        .agg(
            pl.len().alias("parts"),
            pl.col("area_km2").sum().round(0).alias("area_km2"),
            pl.col("sea").mean().round(3).alias("sea_share"),
            pl.col("lat").median().round(1).alias("median_lat"),
            pl.col("lon").median().round(1).alias("median_lon"),
            pl.col("leaf").n_unique().alias("leaf_modules"),
        )
        .sort("parts", "module", descending=[True, False])
    )
    leaf_area = table.group_by("leaf").agg(pl.col("area_km2").sum())["area_km2"]
    n, area = len(paths), table["area_km2"].sum()
    return {
        "parts": n,
        "unassigned": len(partition.unassigned),
        "levels": max(len(path) for path in paths.values()),
        "top_modules": len(top),
        "leaf_modules": len(leaf),
        "largest_top_share_of_parts": round(max(top.values()) / n, 4),
        "largest_leaf_share_of_parts": round(max(leaf.values()) / n, 4),
        "largest_top_share_of_area": round(by_module["area_km2"].max() / area, 4),
        "leaf_modules_of_one_part": sum(1 for size in leaf.values() if size == 1),
        "leaf_modules_under_floor": int((leaf_area < FLOOR_KM2).sum()),
        "top_level_modules": by_module.head(TOP_MODULES).to_dicts(),
    }


def _with_peak(function: Callable, shared, task):
    """In a worker: the result of one task, with the seconds it took and the peak of this worker so far."""
    started = time.perf_counter()
    result = function(shared, task)
    return result, {"seconds": round(time.perf_counter() - started, 1), "peak_gb": round(peak_gb(), 2), "pid": os.getpid()}


@contextmanager
def worker_peaks(module):
    """Replaces `module.map_tasks` with one that also collects, per task, what `_with_peak` measured. Yields that list."""
    seen: list[dict] = []

    def traced(function, tasks, shared, workers=1):
        out = parallel.map_tasks(partial(_with_peak, function), tasks, shared, workers)
        seen.extend(info for _, info in out)
        return [result for result, _ in out]

    with mock.patch.object(module, "map_tasks", traced):
        yield seen


def project(data: dict, growth: float) -> dict:
    """What a selection run at resolution 7 would take, from the measurements of this one, for links that grow by `growth`.

    Every stage is scaled from its measured figure, and wherever the scaling is uncertain it is the larger of the options, so the
    result is a conservative one. The headline has `HEADLINE` times the links of half A, so its links are `f` times those measured,
    and `c` is how many times more cells the grid has.

    - Infomap and the resampling are proportional to links, and so is their memory, with no fixed share. A 10-trial Infomap costs
      as much as 10 times the 1-trial run, or the 10-trial run measured if that is more.
    - `candidate_maps` merges once per unit, and units follow cells, so there are `c` times the merges. A merge goes over a dense
      units x units array (`c` squared) and over the links (`f`), and is as slow as the slower. Those dense arrays need `c` squared
      the memory too, where the rest needs `f`.
    - A bootstrap is a 3-trial Infomap, the resample and the network, which grow with `f`, and a `candidate_maps`, which grows as above.
    - `excess_table` is a Python loop over cells and a product over links, so it grows as the larger of the two, and as the number of
      maps, of which a larger grid can have all `K_MAX - K_MIN + 1`.
    - Workers: as many as fit in `BUDGET_GB` next to the parent, and at most one fewer than the cores.
    """
    stages, infomap = data["stages"], data["infomap"]
    f, c = HEADLINE * growth, CELLS_RATIO
    chosen = data["chosen_markov_time"]
    merge_scale, dense_scale = c * max(c**2, f), max(c**2, f)

    # (a) the headline Infomap run, the worst of the Infomap runs measured, alone and next to the matrix
    measured_ten = stages["infomap_10_trials"]
    child_peak = max(measured_ten["child_peak_rss_gb"], *(run["peak_gb"] for run in infomap.values()))
    alone = f * child_peak
    inside = f * (child_peak + data["ct_gb"])

    # (b) one selection run: the 10-trial Infomap, the candidate maps, 20 bootstraps and the excess
    trial = infomap[f"{chosen}"]["seconds"]
    ten_trials = f * max(HEADLINE_TRIALS * trial, measured_ten["child_seconds"])
    cand = stages[f"candidate_maps_{chosen}"]
    candidate = cand["child_seconds"] * merge_scale
    candidate_peak = cand["child_peak_rss_gb"] * dense_scale

    boot_workers = stages["bootstrap"]["workers"]
    boot = max(worker["seconds"] for worker in boot_workers)
    linear = f * max(boot - cand["child_seconds"], 0.0)  # what a bootstrap does besides its candidate maps
    per_boot = linear + candidate
    parent = f * stages["bootstrap"]["child_peak_rss_gb"]  # the process that runs the workers, which holds the matrix and the maps
    worker = max(f * max(w["peak_gb"] for w in boot_workers), candidate_peak)
    boot_room = (BUDGET_GB - parent) / worker  # workers that fit next to the parent
    boot_fit = max(1, min(CORES - 1, math.floor(boot_room)))
    waves = math.ceil(BOOTSTRAPS / boot_fit)

    which = "excess_table" if "workers" in stages["excess_table"] else "excess_table_without_zero_area_parts"
    excess = stages[which]
    excess_room = (BUDGET_GB - f * excess["child_peak_rss_gb"]) / (f * max(w["peak_gb"] for w in excess["workers"]))
    excess_fit = max(1, min(CORES - 1, math.floor(excess_room)))
    more_maps = (K_MAX - K_MIN + 1) / excess["maps"]
    excess_time = excess["tasks_seconds_sum"] * max(f, c) * more_maps / excess_fit  # the work of every map, over as many workers as fit

    total = ten_trials + candidate + waves * per_boot + excess_time
    without_candidates = ten_trials + waves * linear + excess_time

    # (c) the zip and Parquet of the download, which grow with its rows, and so with the pairs
    disk = data["disk"]
    needed = growth * (disk["zip_bytes"] + disk["parquet_bytes"]) / 1e9
    left = disk["free_gb"] - needed

    a = alone <= BUDGET_GB and inside <= BUDGET_GB
    b = total / 3600 <= HOURS and boot_room >= 1 and excess_room >= 1  # not even one worker is a failure, whatever the time
    cc = left >= FREE_GB
    return {
        "growth": growth,
        "link_factor": f,
        "links": round(f * data["links"]),
        "a_infomap_peak_gb": {"alone": round(alone, 1), "with_matrix": round(inside, 1), "limit": BUDGET_GB, "ok": a},
        "b_selection_hours": {
            "total": round(total / 3600, 2),
            "infomap_10_trials": round(ten_trials / 3600, 2),
            "candidate_maps": round(candidate / 3600, 2),
            "bootstraps": round(waves * per_boot / 3600, 2),
            "one_bootstrap_hours": round(per_boot / 3600, 2),
            "bootstrap_workers": boot_fit,
            "excess_table": round(excess_time / 3600, 2),
            "excess_workers": excess_fit,
            "without_candidate_maps": round(without_candidates / 3600, 2),
            "limit": HOURS,
            "ok": b,
        },
        "memory_gb": {
            "candidate_maps": round(candidate_peak, 1),
            "bootstrap_worker": round(worker, 1),
            "bootstrap_parent": round(parent, 1),
            "excess_worker": round(f * max(w["peak_gb"] for w in excess["workers"]), 1),
            "limit": BUDGET_GB,
        },
        "c_disk_left_gb": {"after_download": round(left, 1), "download_gb": round(needed, 2), "minimum": FREE_GB, "ok": cc},
        "excess_measured_on": which,
        "go": bool(a and b and cc),
    }


def add_projection(data: dict) -> None:
    """The projection for every growth, and the decision on the conservative one, which a measure file can be given again."""
    data["projection"] = {name: project(data, growth) for name, growth in GROWTHS.items()}
    data["decision"] = "go to resolution 7" if data["projection"]["2.5"]["go"] else "stay at resolution 6"


def disk_numbers(counts: Path) -> dict:
    """The size of the download as a zip and as Parquet, and the disk that is free."""
    return {
        "zip_bytes": counts.with_suffix(".zip").stat().st_size,
        "parquet_bytes": sum(p.stat().st_size for p in counts.glob("part-*.parquet")),
        "free_gb": round(shutil.disk_usage(counts).free / 1e9, 1),
    }


def infomap_stage(run: Run, links_path: Path, cells: list[str], trials: int, markov_time: float, part_info: pl.DataFrame) -> Partition:
    """Cluster in a child process and describe what it found. Returns the partition."""
    name = f"infomap_{markov_time}"
    with run.stage(name) as entry:
        result = in_child(_infomap, links_path, cells, SEED, trials, markov_time)
        partition = Partition(result["paths"], result["unassigned"])
        found = structure(partition, part_info)
        run.data["infomap"][str(markov_time)] = {
            "trials": trials,
            "markov_time": markov_time,
            "seconds": round(result["seconds"], 1),
            "load_seconds": round(result["load_seconds"], 1),
            "peak_gb": round(result["peak_gb"], 2),
            "peak_after_load_gb": round(result["peak_after_load_gb"], 2),
            "peak_before_gb": round(result["peak_before_gb"], 2),
            "structure": found,
        }
        entry.update(
            child_seconds=round(result["seconds"], 1),
            load_seconds=round(result["load_seconds"], 1),
            child_peak_rss_gb=round(result["peak_gb"], 2),
            child_peak_after_load_gb=round(result["peak_after_load_gb"], 2),
            top_modules=found["top_modules"],
            leaf_modules=found["leaf_modules"],
            levels=found["levels"],
            largest_top_share_of_parts=found["largest_top_share_of_parts"],
            largest_leaf_share_of_parts=found["largest_leaf_share_of_parts"],
        )
    return partition


def candidate_stage(run: Run, markov_time: float, partition: Partition, ct, areas: dict, adjacency: dict):
    """`candidate_maps` in a child process, with its three stages timed. Returns the maps, or None for a degenerate clustering."""
    with run.stage(f"candidate_maps_{markov_time}") as entry:
        result = in_child(_candidates, partition, ct, areas, adjacency)
        timings = result["timings"]
        entry.update(child_seconds=round(result["seconds"], 1), child_peak_rss_gb=round(result["peak_gb"], 2))
        if "apply_size_floor" in timings:
            entry["units_after_pieces"] = timings["apply_size_floor"]["units_in"]
        for name, call in timings.items():
            entry[f"{name}_seconds"] = call["seconds"]
        if "degenerate" in result:
            entry["degenerate"] = result["degenerate"]
            return None
        found = result["found"]
        entry.update(
            units_after_floor=timings["merge_sequence"]["units_in"],
            merges=timings["merge_sequence"]["units_in"] - K_MIN,
            floor_merges=len(found.sub_realms),
            maps=len(found.maps),
            k_max=max(found.maps),
            k_min=min(found.maps),
        )
    return found


def excess_stage(run: Run, name: str, maps: dict, ct, adjacency: dict, areas: dict, workers: int) -> bool:
    """`excess_table` in a child, with the peak of each worker. A division by zero is recorded and returns False."""
    with run.stage(name) as entry:
        result = in_child(_excess, maps, ct, adjacency, areas, workers)
        entry.update(maps=len(maps), process_workers=workers, parts=len(next(iter(maps.values()))), child_seconds=round(result["seconds"], 1))
        entry["child_peak_rss_gb"] = round(result["peak_gb"], 2)
        if "error" in result:
            entry["error"] = result["error"] + ", in the null maps, from a piece of a realm with no area"
            run.say(f"{name} failed: {result['error']}")
            return False
        excess = result["excess"]
        entry.update(n_null=N_NULL, excess_min=round(min(excess.values()), 4), excess_max=round(max(excess.values()), 4))
        entry["workers"] = result["workers"]
        entry["worker_peak_rss_gb"] = max(w["peak_gb"] for w in result["workers"])
        entry["tasks_seconds_sum"] = round(sum(w["seconds"] for w in result["workers"]), 1)
    return True


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument("--workers", type=int, required=True, help="processes for excess_table")
    parser.add_argument("--resolution", type=int, default=6, help="the ISEA3H resolution of the counts, default 6")
    parser.add_argument("--markov-times", type=float, nargs="+", default=MARKOV_TIMES, help="the first with candidate maps is used for the rest, default 1.0 0.75")
    parser.add_argument("--survey", type=float, nargs="*", default=SURVEY, help="more Markov times, for their module counts and candidate maps only")
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="default data/reference")
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--results", type=Path, default=Path("results/pilot"), help="where res<r>/ goes, default results/pilot")
    parser.add_argument("--project-only", action="store_true", help="measure nothing: project again from the measure.json of an earlier run")
    args = parser.parse_args(argv)

    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{args.resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    release = reference / "exports" / "extended" / "col.zip"
    directory = args.results / f"res{args.resolution}"
    missing = [str(p) for p in (args.counts, args.split, gpkg, land, release, args.counts.with_suffix(".zip")) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")

    if args.project_only:
        measured = json.loads((directory / "measure.json").read_text())
        add_projection(measured)
        (directory / "measure.json").write_text(json.dumps(measured, indent=2) + "\n")
        print(json.dumps(measured["projection"], indent=2))
        print(f"decision at growth 2.5: {measured['decision']}")
        return

    run = Run(directory)
    run.data.update(
        half="A",
        resolution=args.resolution,
        counts=str(args.counts),
        workers=args.workers,
        machine={"cores": CORES, "ram_gb": round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9, 1)},
        infomap={},
    )
    try:
        measure(run, args, directory, gpkg, land, release)
    except BaseException:
        run.say("FAILED\n" + traceback.format_exc())
        raise


def measure(run: Run, args: argparse.Namespace, directory: Path, gpkg: Path, land: Path, release: Path) -> None:
    run.say(f"pilot on half A, resolution {args.resolution}, {args.counts}, {args.workers} workers, {CORES} cores")

    with run.stage("load_taxa") as entry:  # in a child: the load peaks far above anything after it
        taxa, child_peak = in_child(_load_taxa, release)
        entry.update(species=taxa.height, child_peak_rss_gb=round(child_peak, 2))
    with run.stage("load_split") as entry:
        split = load_split(args.split)
        taxa_a = taxa.join(split.filter(pl.col("half") == "A").select("species_id"), on="species_id", how="semi")
        entry["species_in_half"] = taxa_a.height
    del taxa
    with run.stage("load_grid") as entry:
        grid = load_grid(gpkg, land)
        entry["cells"] = grid.cells.height
    with run.stage("load_presence") as entry:
        presence = load_presence(args.counts, split, "A")
        present = presence.filter(pl.col("n") > 0)  # the presence of layer 1: a clean record
        entry.update(rows_with_records=presence.height, rows=present.height, species=present["species_id"].n_unique(), parts=present["part"].n_unique())
        del presence, split

    with run.stage("parts") as entry:
        names = present["part"].unique().sort()
        part_table = parts(grid, names)
        areas = dict(zip(part_table["part"], part_table["area_km2"], strict=True))
        adjacency = part_adjacency(grid, names)
        zero = part_table.filter(pl.col("area_km2") < ZERO_KM2)
        zero_rows = present.join(zero.select("part"), on="part", how="semi")
        part_info = part_table.join(grid.cells.select("cell", "lat", "lon"), on="cell").select("part", "kind", "area_km2", "lat", "lon")
        entry.update(
            parts=part_table.height,
            land_parts=int((part_table["kind"] == "land").sum()),
            sea_parts=int((part_table["kind"] == "sea").sum()),
            zero_area_parts=zero.height,
            zero_area_land_parts=int((zero["kind"] == "land").sum()),
            zero_area_sea_parts=int((zero["kind"] == "sea").sum()),
            zero_area_rows=zero_rows.height,
            zero_area_share_of_rows=round(zero_rows.height / present.height, 5),
            zero_area_share_of_records=round(zero_rows["n"].sum() / present["n"].sum(), 5),
            parts_without_neighbours=sum(1 for near in adjacency.values() if not near),
            adjacent_pairs=sum(len(near) for near in adjacency.values()) // 2,
        )
    with run.stage("cell_taxa") as entry:
        ct = cell_taxa(present.select("species_id", pl.col("part").alias("cell")), taxa_a)
        run.data["ct_gb"] = round((ct.x.data.nbytes + ct.x.indices.nbytes + ct.x.indptr.nbytes) / 1e9, 3)
        entry.update(
            cells=len(ct.cells),
            taxa=len(ct.taxa),
            links=int(ct.x.nnz),
            species=int((ct.rank == 0).sum()),
            genera=int((ct.rank == 1).sum()),
            families=int((ct.rank == 2).sum()),
            matrix_gb=run.data["ct_gb"],
        )
    del present, taxa_a
    with run.stage("network_from") as entry:
        net = network_from(ct)
        entry.update(nodes=len(net.cells) + len(net.taxa), cells=len(net.cells), taxon_nodes=len(net.taxa), links=len(net.weight))
        run.data["links"], run.data["nodes"] = len(net.weight), len(net.cells) + len(net.taxa)
    cells = net.cells
    links_path = directory / "links.npy"
    with run.stage("save_links") as entry:
        np.save(links_path, np.column_stack([net.cell_idx, net.taxon_idx, net.weight]))
        entry["file_gb"] = round(links_path.stat().st_size / 1e9, 3)
    del net

    times = [*args.markov_times, *(m for m in args.survey if m not in args.markov_times)]
    partitions = {markov: infomap_stage(run, links_path, cells, 1, markov, part_info) for markov in times}
    # the fewest units first, so that the slow ones come last and everything before them is measured and saved
    found = {}
    for markov in sorted(times, key=lambda m: len(set(partitions[m].paths.values()))):
        found[markov] = candidate_stage(run, markov, partitions[markov], ct, areas, adjacency)
    chosen = next((m for m in args.markov_times if found[m] is not None), None)
    if chosen is None:
        run.say("every Markov time is degenerate, so the bootstrap and excess_table cannot be measured")
        return
    run.data["chosen_markov_time"] = chosen
    maps = found[chosen].maps
    run.say(f"the later stages use Markov time {chosen}, with {len(maps)} maps")

    with run.stage("infomap_10_trials") as entry:  # the headline's setting, to see what ten trials add
        result = in_child(_infomap, links_path, cells, SEED, HEADLINE_TRIALS, chosen)
        entry.update(
            markov_time=chosen,
            child_seconds=round(result["seconds"], 1),
            child_peak_rss_gb=round(result["peak_gb"], 2),
            child_peak_after_load_gb=round(result["peak_after_load_gb"], 2),
            leaf_modules=len(set(result["paths"].values())),
        )

    with run.stage("bootstrap") as entry:
        result = in_child(_bootstrap, ct, maps, areas, adjacency, chosen)
        scores = result["scores"]
        entry.update(n_boot=BOOT_PROBE, trials=BOOT_TRIALS, markov_time=chosen, maps_scored=len(scores), child_seconds=round(result["seconds"], 1))
        entry.update(child_peak_rss_gb=round(result["peak_gb"], 2), median_stability=round(float(np.median(list(scores.values()))), 3))
        entry["workers"] = result["workers"]
        entry["worker_peak_rss_gb"] = max(w["peak_gb"] for w in result["workers"])
        entry["one_bootstrap_seconds"] = max(w["seconds"] for w in result["workers"])

    if not excess_stage(run, "excess_table", maps, ct, adjacency, areas, args.workers):
        # `null_partitions` divides by the area of each piece of a realm, and a piece of parts without area has none. The measure goes on
        # without those parts, only to time the stage and weigh its workers: the maps it scores are not the pilot's
        kept = {part for part, area in areas.items() if area >= ZERO_KM2}
        reduced = {k: {part: realm for part, realm in labels.items() if part in kept} for k, labels in maps.items()}
        near = {part: adjacency[part] & kept for part in kept}
        excess_stage(run, "excess_table_without_zero_area_parts", reduced, ct, near, {part: areas[part] for part in kept}, args.workers)

    run.data["disk"] = disk_numbers(args.counts)
    add_projection(run.data)
    run.save()
    run.say(json.dumps(run.data["projection"], indent=2))
    run.say(f"decision at growth 2.5: {run.data['decision']}")


if __name__ == "__main__":
    main()
