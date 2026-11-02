"""The exploration pipeline: cleaned counts to network units, to Infomap, to candidate maps, to a screen, to a selected map. Every stage is cached.

A run is one WCVP mode and the cleaning rules switched off, at the resolution of the counts: `<results>/res<r>/<run>/`.

- `cleaned/`: `clean` on half A's counts (`presence.parquet`, `summary.csv` and `summary_by_group.csv`, the removals by layer and rule), and what the later
  stages need of the references, so that they never load them: `major.parquet` (the major group of each species, computed once here) and `taxa.parquet`
  (genus and family ids).
- `units/`: `network_units` (`members`, `presence`, `areas`, `adjacency`, `excluded` and `stats`), `removed_by_region.csv`, and `cells.parquet`, the centres
  and areas of the fine grid, which is all that imputing the excluded parts needs.
- `<configuration>/`: Infomap with the options of one named configuration, then `candidate_maps`.
  - **screen** (`trials_screen` trials, `n_null_screen` null maps): `screen.csv`, the endemism excess of every k.
  - **select** (10 trials, 20 bootstraps of 3 trials, 99 null maps): `selection.csv` with the stability and the excess of every k, and `selected.json`, from `select_level`.

  Both write `partition.parquet` (unit, leaf, top), `structure.json`, `candidates.parquet` (every k, with the realm of each unit of the network), `parts.parquet` (the
  realm of each part for every k, with the selected one marked) and `log.json`. They write the same files, so the files are those of the stage that ran last, and
  `structure.json` says which (`stage`) and with how many trials. A screen that runs again removes what select wrote, which it overwrote.

**The units are the network's cells.** The areas and adjacency that `candidate_maps`, the bootstrap and the null maps are given are the units', and the parts that were
left out (`Units.excluded`) are in none of them. They get their realm from `units.impute`, in `parts.parquet` only, flagged as imputed.

**Caching.** A stage has a stamp (`stage.json`, `stage-screen.json`, `stage-select.json`) holding the hash of its inputs, and it is skipped when its outputs exist with
the same hash. The hash of cleaning is that of the cleaning settings, the layers, the mode and the rules switched off (all in `Settings`), the frozen GRIIS code table,
the `_COMPLETE` marker of the counts and the manifest of the references. Each later stage hashes the stage before it with its own settings (the unit settings, then the
configuration and the parameters), so a change anywhere reruns what follows it.

**A degenerate clustering is a result.** `DegenerateClustering` (too few units above the size floor) and a network with no link at all are written to `structure.json`
as `degenerate` with the reason, with empty maps, and never raised.

**Half A only.** The one loader of counts is `load_country_counts(..., "A")`, and nothing here can ask for another half.
"""

import hashlib
import json
import logging
import multiprocessing
import os
import resource
import sys
import time
import tomllib
from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import asdict, dataclass, fields, replace
from datetime import UTC, datetime
from functools import partial
from pathlib import Path
from typing import Any

import polars as pl

from realms import nullmodel, parallel
from realms import select as selection
from realms.candidates import CandidateMaps, DegenerateClustering, candidate_maps
from realms.clean import LAYERS, SWITCHABLE, WCVP_MODES, Settings, clean, frozen_list_codes
from realms.clean_report import removal_summary
from realms.cluster import InfomapOptions, Partition, run_infomap
from realms.grid import Grid, load_grid
from realms.matrix import CellTaxa, cell_taxa
from realms.natives import load_references
from realms.network import network_from
from realms.nullmodel import excess_table
from realms.occurrences import load_country_counts, read_complete, resolution_of
from realms.select import bootstrap_stability, select_level
from realms.split import load_split
from realms.units import UnitSettings, Units, impute, is_coarse, major_groups, network_units, part_regions, removed_by_region

log = logging.getLogger(__name__)

STAGES = ("clean", "units", "screen", "select")
CLEANED, UNITS = "cleaned", "units"
CLEAN_FILES = ("presence.parquet", "summary.csv", "summary_by_group.csv", "major.parquet", "taxa.parquet")
UNIT_FILES = ("members.parquet", "presence.parquet", "areas.parquet", "adjacency.parquet", "excluded.parquet", "stats.parquet", "removed_by_region.csv", "cells.parquet")
MAP_FILES = ("partition.parquet", "structure.json", "candidates.parquet", "parts.parquet")
SCREEN_FILES = (*MAP_FILES, "screen.csv")
SELECT_FILES = (*MAP_FILES, "selection.csv", "selected.json")
BUDGET_GB = 40.0  # the most that the parent and all its workers may hold together
MARGIN = 1.25  # a worker's later tasks may need more than the first one measured

CANDIDATE_SCHEMA = {"unit": pl.String, "k": pl.Int64, "realm": pl.Int64, "sub_realm": pl.Int64, "sub_realm_of": pl.Int64}
PARTS_SCHEMA = {"part": pl.String, "k": pl.Int64, "realm": pl.Int64, "unit": pl.String, "imputed": pl.Boolean, "selected": pl.Boolean}
UNLINKED = "unlinked"  # why a unit that Infomap gave no module to is left out of the maps


@dataclass(frozen=True)
class Params:
    """What the stages take besides the configuration. The defaults are the design's: 10 trials for the headline run, 20 bootstraps of 3 trials, 99 null maps, and
    the screen's `trials_screen` and `n_null_screen`, which `configs.toml` sets. The candidates are 60 down to 5 realms with a floor of 1M km2."""

    seed: int = 1
    k_max: int = 60
    k_min: int = 5
    floor_km2: float = 1_000_000
    trials_screen: int = 3
    n_null_screen: int = 19
    trials_select: int = 10
    n_boot: int = 20
    boot_trials: int = 3
    n_null_select: int = 99


def load_configs(path: Path) -> tuple[dict[str, InfomapOptions], Params]:
    """The named `InfomapOptions` of a `configs.toml` (`[configs.<name>]`, each key a field of `InfomapOptions`), and its `trials_screen` and `n_null_screen`."""
    data = tomllib.loads(Path(path).read_text())
    named = data.pop("configs", {})
    screen = {k: data.pop(k) for k in ("trials_screen", "n_null_screen") if k in data}
    if data:
        raise ValueError(f"{path} has keys that are not trials_screen, n_null_screen or [configs.<name>]: {sorted(data)}")
    if not named:
        raise ValueError(f"{path} has no configuration: add a [configs.<name>] table, empty for Infomap's defaults")
    known = {f.name for f in fields(InfomapOptions)}
    configs = {}
    for name, table in named.items():
        unknown = sorted(set(table) - known)
        if unknown:
            raise ValueError(f"configuration {name!r} of {path} has keys that are not InfomapOptions fields: {unknown}")
        configs[name] = InfomapOptions(**table)
    return configs, Params(**screen)


def run_name(mode: str, rules_off: Iterable[str] = ()) -> str:
    """`native_only`, or `native_only-no-griis` with the rules switched off in sorted order: the directory of a run."""
    if mode not in WCVP_MODES:
        raise ValueError(f"the WCVP mode is one of {WCVP_MODES}, not {mode!r}")
    unknown = sorted(set(rules_off) - set(SWITCHABLE))
    if unknown:
        raise ValueError(f"rules to switch off are of {SWITCHABLE}, not {unknown}")
    return "-".join([mode, *(f"no-{rule}" for rule in sorted(set(rules_off)))])


def peak_gb() -> float:
    """The most memory this process has had resident so far. `ru_maxrss` is in bytes on macOS and kilobytes on Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)


# ------------------------------------------------------------------ stamps: what makes a stage current


def digest(value: Any) -> str:
    """The sha256 of the JSON of `value` with sorted keys. Sets and paths are written as sorted lists and text."""

    def plain(item: Any) -> Any:
        if isinstance(item, (set, frozenset)):
            return sorted(item)
        if isinstance(item, Path):
            return str(item)
        raise TypeError(f"{type(item)} is not hashable here")

    return hashlib.sha256(json.dumps(value, sort_keys=True, default=plain).encode()).hexdigest()


def read_stamp(directory: Path, name: str = "stage") -> dict:
    path = Path(directory) / f"{name}.json"
    return json.loads(path.read_text()) if path.exists() else {}


def is_current(directory: Path, name: str, hash_: str, files: Sequence[str]) -> bool:
    """Whether the stage's stamp has this hash and all its files are there."""
    return read_stamp(directory, name).get("hash") == hash_ and all((Path(directory) / f).exists() for f in files)


def write_stamp(directory: Path, name: str, hash_: str, inputs: dict, log_entry: dict) -> None:
    stamp = {"hash": hash_, "written": datetime.now(UTC).isoformat(timespec="seconds"), "inputs": inputs, "log": log_entry}
    (Path(directory) / f"{name}.json").write_text(json.dumps(stamp, indent=2, default=str) + "\n")


def clean_inputs(settings: Settings, layers: Iterable[int], list_codes: pl.DataFrame, marker: dict, manifest: str | None = None) -> dict:
    """What the cleaning stage's hash is of: the settings (the mode and the rules switched off are in them), the layers, the content of the frozen code table, the
    counts' `_COMPLETE` marker and the references' manifest (as a digest)."""
    return {
        "settings": asdict(settings) | {"rules_off": sorted(settings.rules_off)},
        "layers": sorted(layers),
        "list_codes": digest(list_codes.sort(list_codes.columns).rows()),
        "counts": marker,
        "manifest": manifest,
    }


def _elapsed(started: float) -> dict:
    return {"seconds": round(time.perf_counter() - started, 1), "peak_rss_gb": round(peak_gb(), 2)}


# ------------------------------------------------------------------ cleaned counts


def clean_stage(
    out: Path, counts: pl.DataFrame, refs, grid: Grid, settings: Settings, layers: Iterable[int], list_codes: pl.DataFrame, inputs: dict, *, force: bool = False
) -> dict:
    """`clean` on the counts, and the parts of the references that the later stages use, into `<out>/cleaned/`. Skipped if the hash of `inputs` is the stamp's."""
    directory = Path(out) / CLEANED
    hash_ = digest(inputs)
    if not force and is_current(directory, "stage", hash_, CLEAN_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "stage.json").unlink(missing_ok=True)
    cleaned = clean(counts, refs, grid, frozenset(layers), settings, list_codes)
    species = cleaned.presence.select("species_id").unique()
    major = major_groups(refs.taxa, refs.groups).join(species, on="species_id", how="semi").sort("species_id")
    taxa = refs.taxa.select("species_id", "genus_id", "family_id").join(species, on="species_id", how="semi").sort("species_id")
    cleaned.presence.write_parquet(directory / "presence.parquet")
    removal_summary(cleaned, refs, by=()).write_csv(directory / "summary.csv")
    removal_summary(cleaned, refs, by=("tax_group",)).write_csv(directory / "summary_by_group.csv")
    major.write_parquet(directory / "major.parquet")
    taxa.write_parquet(directory / "taxa.parquet")
    entry = {
        **_elapsed(started),
        "presence_rows": cleaned.presence.height,
        "species": species.height,
        "records": int(cleaned.presence["records"].sum()),
        "removed_rows": cleaned.removed.height,
        "major_group_species": major.height,
    }
    write_stamp(directory, "stage", hash_, inputs, entry)
    _note(Path(out), "clean", entry)
    return {"skipped": False, "hash": hash_, **entry}


def _note(out: Path, stage: str, entry: dict) -> None:
    """Add a stage's figures to the run's `log.json`."""
    path = Path(out) / "log.json"
    data = json.loads(path.read_text()) if path.exists() else {}
    data[stage] = entry
    path.write_text(json.dumps(data, indent=2, default=str) + "\n")


# ------------------------------------------------------------------ units


def units_stage(
    out: Path, grid: Grid, coarse: Grid, coarse_gpkg: Path, unit_settings: UnitSettings, regions: pl.DataFrame, parent: str, *, force: bool = False
) -> dict:
    """`network_units` on the cleaned presence, and its files, into `<out>/units/`. `parent` is the hash of the cleaning stage. Skipped if the hash is the stamp's."""
    directory = Path(out) / UNITS
    inputs = {"cleaned": parent, "unit_settings": asdict(unit_settings)}
    hash_ = digest(inputs)
    if not force and is_current(directory, "stage", hash_, UNIT_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "stage.json").unlink(missing_ok=True)
    presence = pl.read_parquet(Path(out) / CLEANED / "presence.parquet")
    major = pl.read_parquet(Path(out) / CLEANED / "major.parquet")
    units = network_units(presence, grid, coarse, major, unit_settings, coarse_gpkg=coarse_gpkg)
    write_units(directory, units)
    by_region = removed_by_region(units, regions)
    by_region.write_csv(directory / "removed_by_region.csv")
    grid.cells.write_parquet(directory / "cells.parquet")
    network = units.stats.filter(pl.col("status") == "network")
    entry = {
        **_elapsed(started),
        "units": len(units.areas),
        "coarse_units": sum(1 for u in units.areas if is_coarse(u)),
        "parts_in_network": len(units.members),
        "parts_excluded": units.excluded.height,
        "network_records": int(network["records"].sum() or 0),
    }
    write_stamp(directory, "stage", hash_, inputs, entry)
    _note(Path(out), "units", entry)
    return {"skipped": False, "hash": hash_, **entry}


def write_units(directory: Path, units: Units) -> None:
    pl.DataFrame({"part": list(units.members), "unit": list(units.members.values())}, schema={"part": pl.String, "unit": pl.String}).write_parquet(directory / "members.parquet")
    units.presence.write_parquet(directory / "presence.parquet")
    pl.DataFrame({"unit": list(units.areas), "area_km2": list(units.areas.values())}, schema={"unit": pl.String, "area_km2": pl.Float64}).write_parquet(directory / "areas.parquet")
    pairs = [(u, n) for u, near in units.adjacency.items() for n in sorted(near)]
    pl.DataFrame(pairs, schema={"unit": pl.String, "neighbour": pl.String}, orient="row").write_parquet(directory / "adjacency.parquet")
    units.excluded.write_parquet(directory / "excluded.parquet")
    units.stats.write_parquet(directory / "stats.parquet")


def read_units(out: Path) -> Units:
    """The `Units` that `units_stage` wrote, with the unit settings of its stamp."""
    directory = Path(out) / UNITS
    stamp = read_stamp(directory)
    if not stamp:
        raise FileNotFoundError(f"{directory} has no stage.json: run the units stage first")
    areas = pl.read_parquet(directory / "areas.parquet")
    adjacency: dict[str, set[str]] = {u: set() for u in areas["unit"]}
    for unit, neighbour in pl.read_parquet(directory / "adjacency.parquet").iter_rows():
        adjacency[unit].add(neighbour)
    return Units(
        members=dict(pl.read_parquet(directory / "members.parquet").iter_rows()),
        presence=pl.read_parquet(directory / "presence.parquet"),
        areas=dict(areas.iter_rows()),
        adjacency=adjacency,
        excluded=pl.read_parquet(directory / "excluded.parquet"),
        stats=pl.read_parquet(directory / "stats.parquet"),
        settings=UnitSettings(**stamp["inputs"]["unit_settings"]),
    )


def read_light_grid(out: Path) -> Grid:
    """The fine grid's cells (centres and areas) without its adjacency, which is all `units.impute` reads."""
    cells = pl.read_parquet(Path(out) / UNITS / "cells.parquet")
    return Grid(0, cells, {c: set() for c in cells["cell"]})


def network_matrix(out: Path, units: Units | None = None) -> CellTaxa:
    """The cell-taxon matrix of the units of the network (they are its cells), with the genus and family of each species from the cleaned stage."""
    units = units or read_units(out)
    taxa = pl.read_parquet(Path(out) / CLEANED / "taxa.parquet")
    return cell_taxa(units.presence.select(pl.col("unit").alias("cell"), "species_id"), taxa)


def labelled_units(units: Units, labelled: set[str]) -> Units:
    """`units` without the ones Infomap gave no module to (a unit with no link), which are left out like the thin ones: in `excluded`, as `unlinked`, to be imputed.

    The `records` and `species` of such a part are those of its unit. Returns `units` itself if every unit has a module.
    """
    lost = set(units.areas) - labelled
    if not lost:
        return units
    parts = [p for p, u in units.members.items() if u in lost]
    stats = units.stats.filter(pl.col("unit").is_in(lost)).select("unit", "records", "species")
    unlinked = (
        pl.DataFrame({"part": parts, "tried_unit": [units.members[p] for p in parts]}, schema={"part": pl.String, "tried_unit": pl.String})
        .join(stats, left_on="tried_unit", right_on="unit")
        .with_columns(pl.lit(UNLINKED).alias("reason"))
        .select(*units.excluded.columns)
        .cast(units.excluded.schema)
    )
    return replace(
        units,
        members={p: u for p, u in units.members.items() if u in labelled},
        areas={u: a for u, a in units.areas.items() if u in labelled},
        adjacency={u: near & labelled for u, near in units.adjacency.items() if u in labelled},
        excluded=pl.concat([units.excluded, unlinked]).sort("part"),
    )


# ------------------------------------------------------------------ workers


def fit_workers(requested: int, parent_gb: float, worker_gb: float, budget_gb: float = BUDGET_GB) -> int:
    """How many workers fit in `budget_gb` next to the parent, each at its measured peak plus `MARGIN`. At least one, at most `requested`."""
    return max(1, min(requested, int((budget_gb - parent_gb) // (worker_gb * MARGIN))))


def _with_peak(function: Callable, shared: Any, task: Any) -> tuple[Any, dict]:
    """In a worker: the result of one task, with the seconds it took and the peak of this process so far."""
    started = time.perf_counter()
    result = function(shared, task)
    return result, {"seconds": round(time.perf_counter() - started, 2), "peak_gb": round(peak_gb(), 3), "pid": os.getpid()}


@contextmanager
def measured_workers(requested: int, seen: list[dict], used: list[int], budget_gb: float = BUDGET_GB):
    """Within it, the bootstrap's and the null maps' tasks record the seconds and peak of each (in `seen`), and run in as many workers as fit.

    With more than one worker requested and more than one task, the first two tasks run in two workers to measure a worker's peak, and the rest in `fit_workers`
    of them next to the parent (`used` gets the number). A task is a pure function of its arguments, so the results are the same for any number of workers.
    """

    def traced(function, tasks, shared, workers=1):
        if requested <= 1:
            used.append(1)
            return [function(shared, task) for task in tasks]  # in this process, whose peak is the stage's
        wrapped = partial(_with_peak, function)
        if len(tasks) < 2:
            out = [wrapped(shared, task) for task in tasks]
            used.append(1)
        else:
            out = parallel.map_tasks(wrapped, tasks[:2], shared, 2)
            fitted = fit_workers(requested, peak_gb(), max(info["peak_gb"] for _, info in out), budget_gb)
            used.append(fitted)
            if len(tasks) > 2:
                out += parallel.map_tasks(wrapped, tasks[2:], shared, fitted)
        seen.extend(info for _, info in out)
        return [result for result, _ in out]

    original = selection.map_tasks, nullmodel.map_tasks
    selection.map_tasks = nullmodel.map_tasks = traced
    try:
        yield
    finally:
        selection.map_tasks, nullmodel.map_tasks = original


# ------------------------------------------------------------------ one configuration


def structure(partition: Partition, areas: dict[str, float], n_units: int) -> dict:
    """What Infomap found: modules (the finest), the top-level ones, the levels, the largest module's share of units and area, and single-unit modules."""
    paths = partition.paths
    if not paths:
        return {
            "units": n_units, "unassigned": len(partition.unassigned), "modules": 0, "top_modules": 0, "levels": 0, "single_unit_modules": 0,
            "largest_module_share_of_units": None, "largest_module_share_of_area": None,
            "largest_top_module_share_of_units": None, "largest_top_module_share_of_area": None,
        }
    total = sum(areas[u] for u in paths)

    def shares(key: Callable[[tuple[int, ...]], Any]) -> tuple[int, Counter, float]:
        members = Counter(key(path) for path in paths.values())
        size: Counter = Counter()
        for unit, path in paths.items():
            size[key(path)] += areas[unit]
        return len(members), members, max(size.values()) / total

    modules, leaf, leaf_area = shares(lambda path: path)
    top_modules, top, top_area = shares(lambda path: path[0])
    return {
        "units": n_units,
        "unassigned": len(partition.unassigned),
        "modules": modules,
        "top_modules": top_modules,
        "levels": max(len(path) for path in paths.values()),
        "single_unit_modules": sum(1 for n in leaf.values() if n == 1),
        "largest_module_share_of_units": round(max(leaf.values()) / len(paths), 4),
        "largest_module_share_of_area": round(leaf_area, 4),
        "largest_top_module_share_of_units": round(max(top.values()) / len(paths), 4),
        "largest_top_module_share_of_area": round(top_area, 4),
    }


def _cluster(units: Units, ct: CellTaxa, options: InfomapOptions, trials: int, params: Params, timings: dict) -> tuple[Partition, CandidateMaps | None, str | None]:
    """Infomap, then the candidate maps. The reason is the message of a `DegenerateClustering`, or of a network that has no link at all, and the maps are then `None`."""
    started = time.perf_counter()
    net = network_from(ct)
    partition = run_infomap(net, seed=params.seed, trials=trials, options=options)
    timings["infomap"] = round(time.perf_counter() - started, 2)
    if not partition.paths:
        return partition, None, "no unit has a link to cluster on"
    started = time.perf_counter()
    try:
        found = candidate_maps(partition, ct, units.areas, units.adjacency, k_max=params.k_max, k_min=params.k_min, floor_km2=params.floor_km2)
    except DegenerateClustering as error:
        return partition, None, str(error)
    finally:
        timings["candidate_maps"] = round(time.perf_counter() - started, 2)
    return partition, found, None


def _write_maps(
    directory: Path, stage: str, trials: int, partition: Partition, shape: dict, found: CandidateMaps | None, reason: str | None, units: Units, grid: Grid, selected: int | None
) -> None:
    """The files that screen and select share. The maps of a degenerate clustering are empty files."""
    leaf, top = partition.leaf_labels(), partition.top_labels()
    ids = sorted(units.areas)
    pl.DataFrame(
        {"unit": ids, "leaf": [leaf.get(u) for u in ids], "top": [top.get(u) for u in ids]}, schema={"unit": pl.String, "leaf": pl.Int64, "top": pl.Int64}
    ).write_parquet(directory / "partition.parquet")
    shape = {"stage": stage, "trials": trials, **shape, "degenerate": reason is not None, "degenerate_reason": reason}
    if found is None:
        (directory / "structure.json").write_text(json.dumps(shape, indent=2) + "\n")
        pl.DataFrame(schema=CANDIDATE_SCHEMA).write_parquet(directory / "candidates.parquet")
        pl.DataFrame(schema=PARTS_SCHEMA).write_parquet(directory / "parts.parquet")
        return
    sub_realms = found.sub_realms
    cells = sorted(next(iter(found.maps.values())))
    sub = [found.units[c] if found.units[c] in sub_realms else None for c in cells]
    rows = []
    for k, labels in sorted(found.maps.items()):
        rows.append(
            pl.DataFrame(
                {"unit": cells, "k": k, "realm": [labels[c] for c in cells], "sub_realm": sub, "sub_realm_of": [labels[c] if s is not None else None for c, s in zip(cells, sub, strict=True)]},
                schema=CANDIDATE_SCHEMA,
            )
        )
    pl.concat(rows).write_parquet(directory / "candidates.parquet")
    shape = {**shape, "candidate_ks": [min(found.maps), max(found.maps)], "sub_realms": len(sub_realms)}
    (directory / "structure.json").write_text(json.dumps(shape, indent=2) + "\n")
    pl.concat(_part_rows(found.maps, units, grid, selected)).write_parquet(directory / "parts.parquet")


def _part_rows(maps: dict[int, dict[str, int]], units: Units, grid: Grid, selected: int | None) -> list[pl.DataFrame]:
    """For each k, the realm of every part of the network and of every excluded part, imputed (`impute`), with the unit of a part of the network."""
    network = labelled_units(units, set(next(iter(maps.values()))))
    frames = []
    for k, labels in sorted(maps.items()):
        realms, imputed = impute(labels, network, grid)
        parts = sorted(realms)
        frames.append(
            pl.DataFrame(
                {
                    "part": parts,
                    "k": k,
                    "realm": [realms[p] for p in parts],
                    "unit": [network.members.get(p) for p in parts],
                    "imputed": [p in imputed for p in parts],
                    "selected": k == selected,
                },
                schema=PARTS_SCHEMA,
            )
        )
    return frames


def _config_dir(out: Path, name: str, stages_to_clear: Sequence[str]) -> Path:
    directory = Path(out) / name
    directory.mkdir(parents=True, exist_ok=True)
    for stage in stages_to_clear:
        (directory / f"stage-{stage}.json").unlink(missing_ok=True)
    return directory


def _log(directory: Path, stage: str, entry: dict) -> None:
    path = directory / "log.json"
    data = json.loads(path.read_text()) if path.exists() else {}
    data[stage] = entry
    path.write_text(json.dumps(data, indent=2, default=str) + "\n")


def _table_rows(table: dict[int, float], name: str) -> pl.DataFrame:
    return pl.DataFrame({"k": list(table), name: list(table.values())}, schema={"k": pl.Int64, name: pl.Float64}).sort("k")


def config_inputs(parent: str, options: InfomapOptions, params: Params, stage: str) -> dict:
    """What a configuration's stamp is a hash of: the units stage, the options, the parameters and which stage it is."""
    return {"units": parent, "options": asdict(options), "params": asdict(params), "stage": stage}


def screen_stage(
    out: Path, name: str, units: Units, ct: CellTaxa, grid: Grid, options: InfomapOptions, params: Params, parent: str, *, workers: int = 1, force: bool = False
) -> dict:
    """Infomap with `trials_screen` trials, `candidate_maps`, and the endemism excess of every k over `n_null_screen` null maps, into `<out>/<name>/`.

    `parent` is the hash of the units stage. A degenerate clustering is written to `structure.json` and never raised.
    """
    inputs = config_inputs(parent, options, params, "screen")
    hash_ = digest(inputs)
    directory = Path(out) / name
    if not force and is_current(directory, "stage-screen", hash_, SCREEN_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    directory = _config_dir(out, name, ("screen", "select"))  # it writes the files that a select result has, so that result is out of date
    for stale in ("selected.json", "selection.csv"):
        (directory / stale).unlink(missing_ok=True)
    timings: dict = {}
    seen: list[dict] = []
    used: list[int] = []
    partition, found, reason = _cluster(units, ct, options, params.trials_screen, params, timings)
    if found is None:
        pl.DataFrame(schema={"k": pl.Int64, "excess": pl.Float64}).write_csv(directory / "screen.csv")
    else:
        begun = time.perf_counter()
        with measured_workers(workers, seen, used):
            excess = excess_table(found.maps, ct, units.adjacency, units.areas, n_null=params.n_null_screen, seed=params.seed, workers=workers)
        timings["excess_table"] = round(time.perf_counter() - begun, 2)
        _table_rows(excess, "excess").write_csv(directory / "screen.csv")
    shape = structure(partition, units.areas, len(units.areas))
    _write_maps(directory, "screen", params.trials_screen, partition, shape, found, reason, units, grid, None)
    entry = {**_elapsed(started), "stages": timings, "workers_requested": workers, "workers_used": max(used, default=1), "workers": seen, "degenerate": reason is not None}
    _log(directory, "screen", entry)
    write_stamp(directory, "stage-screen", hash_, inputs, {k: v for k, v in entry.items() if k != "workers"})
    return {"skipped": False, "hash": hash_, "degenerate": reason is not None, **entry}


def select_stage(
    out: Path, name: str, units: Units, ct: CellTaxa, grid: Grid, options: InfomapOptions, params: Params, parent: str, *, workers: int = 1, force: bool = False
) -> dict:
    """Infomap with `trials_select` trials, `candidate_maps`, `bootstrap_stability`, `excess_table` over `n_null_select` null maps, and `select_level`, into `<out>/<name>/`.

    `parent` is the hash of the units stage. A degenerate clustering is written to `structure.json` and `selected.json` and never raised.
    """
    inputs = config_inputs(parent, options, params, "select")
    hash_ = digest(inputs)
    directory = Path(out) / name
    if not force and is_current(directory, "stage-select", hash_, SELECT_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    directory = _config_dir(out, name, ("select",))
    timings: dict = {}
    seen: list[dict] = []
    used: list[int] = []
    partition, found, reason = _cluster(units, ct, options, params.trials_select, params, timings)
    chosen: int | None = None
    if found is None:
        pl.DataFrame(schema={"k": pl.Int64, "stability": pl.Float64, "excess": pl.Float64, "selected": pl.Boolean}).write_csv(directory / "selection.csv")
        selected = {"degenerate": True, "reason": reason}
    else:
        with measured_workers(workers, seen, used):
            begun = time.perf_counter()
            stability = bootstrap_stability(
                ct, found.maps, units.areas, units.adjacency, n_boot=params.n_boot, seed=params.seed, trials=params.boot_trials,
                floor_km2=params.floor_km2, k_max=params.k_max, options=options, workers=workers,
            )
            timings["bootstrap_stability"] = round(time.perf_counter() - begun, 2)
            begun = time.perf_counter()
            excess = excess_table(found.maps, ct, units.adjacency, units.areas, n_null=params.n_null_select, seed=params.seed, workers=workers)
            timings["excess_table"] = round(time.perf_counter() - begun, 2)
        table = _table_rows(stability, "stability").join(_table_rows(excess, "excess"), on="k")
        choice = select_level(table)
        chosen = choice.k
        table.with_columns((pl.col("k") == chosen).alias("selected")).write_csv(directory / "selection.csv")
        selected = {"degenerate": False, "k": choice.k, "fell_short": choice.fell_short, "reason": choice.reason, "realms": choice.k, "units": len(units.areas)}
    (directory / "selected.json").write_text(json.dumps(selected, indent=2) + "\n")
    shape = structure(partition, units.areas, len(units.areas))
    _write_maps(directory, "select", params.trials_select, partition, shape, found, reason, units, grid, chosen)
    entry = {**_elapsed(started), "stages": timings, "workers_requested": workers, "workers_used": max(used, default=1), "workers": seen, "degenerate": reason is not None}
    _log(directory, "select", entry)
    write_stamp(directory, "stage-select", hash_, inputs, {k: v for k, v in entry.items() if k != "workers"})
    return {"skipped": False, "hash": hash_, "degenerate": reason is not None, **({"k": chosen} if chosen is not None else {}), **entry}


# ------------------------------------------------------------------ jobs: load what a stage needs, and only if it will run


@dataclass(frozen=True)
class CleanSpec:
    counts: Path
    reference: Path
    split: Path
    frozen_dir: Path
    out: Path
    settings: Settings
    layers: tuple[int, ...]
    inputs: dict
    force: bool = False


@dataclass(frozen=True)
class UnitsSpec:
    reference: Path
    out: Path
    resolution: int
    unit_settings: UnitSettings
    parent: str
    force: bool = False


def clean_job(spec: CleanSpec) -> dict:
    """Load the references, the grid and half A's counts, and run `clean_stage`. Nothing is loaded if the stage is current."""
    hash_ = digest(spec.inputs)
    if not spec.force and is_current(Path(spec.out) / CLEANED, "stage", hash_, CLEAN_FILES):
        return {"skipped": True, "hash": hash_}
    resolution = resolution_of(spec.counts)
    list_codes = frozen_list_codes(resolution, spec.frozen_dir)
    grid = load_grid(spec.reference / "grids" / f"ISEA3H-{resolution}.gpkg", spec.reference / "natural_earth" / "ne_10m_land.zip")
    if grid.resolution != resolution:
        raise ValueError(f"the grid file holds resolution {grid.resolution}, not {resolution}")
    refs = load_references(grid, spec.reference)
    split = load_split(spec.split)
    counts = load_country_counts(spec.counts, split, "A")
    del split
    return clean_stage(spec.out, counts, refs, grid, spec.settings, spec.layers, list_codes, spec.inputs, force=spec.force)


def _manifest_path(reference: Path, item: str) -> Path:
    manifest = json.loads((reference / "manifest.json").read_text())
    entry = manifest.get(item)
    if not entry or not entry.get("path"):
        raise ValueError(f"{reference / 'manifest.json'} has no {item}: run scripts/fetch_reference.py")
    return Path(entry["path"])


def units_job(spec: UnitsSpec) -> dict:
    """Load the grids and the WGSRPD regions (the res-5 file and the level 3 polygons from the manifest), and run `units_stage`. Nothing is loaded if it is current."""
    directory = Path(spec.out) / UNITS
    hash_ = digest({"cleaned": spec.parent, "unit_settings": asdict(spec.unit_settings)})
    if not spec.force and is_current(directory, "stage", hash_, UNIT_FILES):
        return {"skipped": True, "hash": hash_}
    land = spec.reference / "natural_earth" / "ne_10m_land.zip"
    grid = load_grid(_manifest_path(spec.reference, f"isea3h_{spec.resolution}"), land)
    coarse_gpkg = _manifest_path(spec.reference, f"isea3h_{spec.unit_settings.coarse_resolution}")
    coarse = load_grid(coarse_gpkg, land)
    regions = part_regions(grid, _manifest_path(spec.reference, "wgsrpd_level3"))
    return units_stage(spec.out, grid, coarse, coarse_gpkg, spec.unit_settings, regions, spec.parent, force=spec.force)


def _measured(function: Callable, spec: Any) -> tuple[Any, float]:
    return function(spec), peak_gb()


def in_child(function: Callable, spec: Any) -> tuple[Any, float]:
    """`function(spec)` in a child process of its own, which is spawned, so that its peak memory is its own. Returns the result and the child's peak."""
    with ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("spawn")) as pool:
        return pool.submit(_measured, function, spec).result()


# ------------------------------------------------------------------ the whole run


@dataclass(frozen=True)
class Plan:
    counts: Path
    reference: Path
    split: Path
    frozen_dir: Path
    results: Path
    configs: Path
    mode: str
    rules_off: tuple[str, ...]
    stage: str
    only: tuple[str, ...] = ()
    workers: int = 1
    force: bool = False
    isolate: bool = True  # clean and units each in a child process of their own, so that the parent never holds their peak and each has one to report


def run_pipeline(plan: Plan, unit_settings: UnitSettings = UnitSettings()) -> dict:
    """Every stage up to `plan.stage`, each skipped when its inputs are those of its stamp. `plan.force` reruns the stage that was asked for.

    Returns the directory of the run and each stage's result.
    """
    if plan.stage not in STAGES:
        raise ValueError(f"the stage is one of {STAGES}, not {plan.stage!r}")
    resolution = resolution_of(plan.counts)  # before anything is loaded
    list_codes = frozen_list_codes(resolution, plan.frozen_dir)
    marker = read_complete(plan.counts)
    configs, params = load_configs(plan.configs)
    missing = sorted(set(plan.only) - set(configs))
    if missing:
        raise ValueError(f"--only names {missing}, and {plan.configs} has {sorted(configs)}")
    names = list(plan.only) or list(configs)
    settings = Settings(wcvp_mode=plan.mode, rules_off=frozenset(plan.rules_off))
    out = Path(plan.results) / f"res{resolution}" / run_name(plan.mode, plan.rules_off)
    out.mkdir(parents=True, exist_ok=True)
    manifest = plan.reference / "manifest.json"
    inputs = clean_inputs(settings, LAYERS, list_codes, marker, hashlib.sha256(manifest.read_bytes()).hexdigest() if manifest.exists() else None)
    wanted = STAGES.index(plan.stage)
    force = {stage: plan.force and stage == plan.stage for stage in STAGES}
    call = in_child if plan.isolate else (lambda function, spec: (function(spec), peak_gb()))
    done: dict = {"directory": out}

    begun = time.perf_counter()
    result, peak = call(clean_job, CleanSpec(plan.counts, plan.reference, plan.split, plan.frozen_dir, out, settings, tuple(LAYERS), inputs, force["clean"]))
    done["clean"] = result
    if not result["skipped"]:
        # `seconds` is the stage itself and `wall_seconds` the job with the loading of the references, the grid and the counts, whose peak is the child's
        done["clean"] = result | {"wall_seconds": round(time.perf_counter() - begun, 1), "child_peak_rss_gb": round(peak, 2)}
        _note(out, "clean", {**read_stamp(out / CLEANED)["log"], "wall_seconds": done["clean"]["wall_seconds"], "child_peak_rss_gb": round(peak, 2)})
    if wanted >= 1:
        begun = time.perf_counter()
        result, peak = call(units_job, UnitsSpec(plan.reference, out, resolution, unit_settings, done["clean"]["hash"], force["units"]))
        done["units"] = result
        if not result["skipped"]:
            done["units"] = result | {"wall_seconds": round(time.perf_counter() - begun, 1), "child_peak_rss_gb": round(peak, 2)}
            _note(out, "units", {**read_stamp(out / UNITS)["log"], "wall_seconds": done["units"]["wall_seconds"], "child_peak_rss_gb": round(peak, 2)})
    if wanted >= 2:
        parent = done["units"]["hash"]
        stages = [stage for stage in ("screen", "select") if STAGES.index(stage) <= wanted]
        files = {"screen": SCREEN_FILES, "select": SELECT_FILES}

        def current(name: str, stage: str) -> bool:
            return is_current(out / name, f"stage-{stage}", digest(config_inputs(parent, configs[name], params, stage)), files[stage])

        needed = [n for n in names if plan.force or not all(current(n, stage) for stage in stages)]
        units = read_units(out)
        grid = read_light_grid(out)
        ct = network_matrix(out, units) if needed else None
        for name in names:
            done[name] = {}
            if name not in needed:
                done[name] = {stage: {"skipped": True} for stage in stages}
                continue
            done[name]["screen"] = screen_stage(out, name, units, ct, grid, configs[name], params, parent, workers=plan.workers, force=force["screen"])
            if wanted >= 3:
                done[name]["select"] = select_stage(out, name, units, ct, grid, configs[name], params, parent, workers=plan.workers, force=force["select"])
    return done
