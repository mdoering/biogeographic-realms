"""The exploration pipeline: cleaned counts to network units, to Infomap, to candidate maps, to a screen, to a selected map. Every stage is cached.

A run is one WCVP mode and the cleaning rules switched off, at the resolution of the counts: `<results>/res<r>/<run>/`. A run may also merge duplicate genus and family
names (`matrix.duplicate_names`) in the cell-taxon matrix: it is a run of its own, named with `-merged-names`, whose cleaned taxa carry the names that merging needs.

- `cleaned/`: `clean` on half A's counts (`presence.parquet`, `summary.csv` and `summary_by_group.csv`, the removals by layer and rule), and what the later
  stages need of the references, so that they never load them: `major.parquet` (the major group of each species, computed once here) and `taxa.parquet`
  (genus and family ids, and in a run that merges duplicate names also `kingdom`, `order`, `family` and `genus`).
- `units/`: `network_units` (`members`, `presence`, `areas`, `adjacency`, `excluded` and `stats`), `removed_by_region.csv`, and `cells.parquet`, the centres
  and areas of the fine grid, which is all that imputing the excluded parts needs.
- `<configuration>/` (`<configuration>-s<N>/` for a seed N other than 1): Infomap with the options of one named configuration, then `candidate_maps`.
  - **screen** (`trials_screen` trials, `n_null_screen` null maps), in `<configuration>/`: `screen.csv`, the endemism excess of every k.
  - **select** (10 trials, 20 bootstraps of 3 trials, 99 null maps), in `<configuration>/select/`: `selection.csv` with the stability and the excess of every k, and
    `selected.json`, from `select_level`.

  Each writes `partition.parquet` (unit, leaf, top), `structure.json`, `candidates.parquet` (every k, with the realm of each unit of the network), `parts.parquet` (the
  realm of each part for every k, with the selected one marked) and `log.json` into its own directory, so neither overwrites the other's. `read_selected` reads a select
  result back, checked against the stamps. Runs made before select had a directory of its own are moved by `migrate_select`.

**The units are the network's cells.** The areas and adjacency that `candidate_maps`, the bootstrap and the null maps are given are the units', and the parts that were
left out (`Units.excluded`) are in none of them. They get their realm from `units.impute`, in `parts.parquet` only, flagged as imputed.

**Caching.** A stage has a stamp (`stage.json`, `stage-screen.json`, `stage-select.json`) holding the hash of its inputs, and it is skipped when its outputs exist with
the same hash. The hash of cleaning is that of the cleaning settings, the layers, the mode and the rules switched off (all in `Settings`), the frozen GRIIS code table,
the `_COMPLETE` marker of the counts and the manifest of the references. Each later stage hashes the stage before it with its own settings (the unit settings, then the
configuration and the parameters, the seed among them), so a change of an input anywhere reruns what follows it. Merging duplicate names enters the hash of cleaning and
of each configuration as `merge_duplicate_names` only when it is on, so that the hashes of the runs on ids are what they were before the option existed; so do a stage's
`METHOD_VERSION` when it is not 1, and the half with the split file's sha256 when it is not half A.

The hashes are of the inputs, not of what the stage before wrote. So a stage that writes its outputs again (`--force`, or a changed input) first deletes the stamps of
every stage of the run that depends on it (`invalidate_after`): cleaning those of the units and of every configuration, the units those of every configuration. A screen
and a select depend on the units only, and invalidate nothing. A change of what a stage does, with its inputs the same, raises its `METHOD_VERSION`.

**Units that Infomap puts in no module** (a unit with no link) leave the scores: `assigned_only` takes them out of the matrix, areas and adjacency that `candidate_maps`, the bootstrap and the null maps are given, so
the scored cells are the mapped ones. `structure.json` lists them (`unassigned`, `unassigned_units`, and a `warning` above 1% of the units), and `parts.parquet` imputes them as it does the excluded parts.

**A degenerate clustering is a result.** `DegenerateClustering` (too few units above the size floor) and a network with no link at all are written to `structure.json`
as `degenerate` with the reason, with empty maps, and never raised.

**Half A only.** The one loader of counts is `load_country_counts(..., "A")`, and nothing here can ask for another half: `Plan` refuses any other (sealed until M3).
The half enters the hashes and the run's name only when it is not A, so that a later run on another half can never take a stamp or a directory of half A's.
"""

import hashlib
import json
import logging
import multiprocessing
import os
import re
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

import numpy as np
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
from realms.occurrences import VISITS, load_country_counts, read_complete, resolution_of
from realms.select import bootstrap_stability, select_level
from realms.split import load_split
from realms.units import UnitSettings, Units, impute, is_coarse, major_groups, network_units, part_regions, removed_by_region

log = logging.getLogger(__name__)

STAGES = ("clean", "units", "screen", "select")
CLEANED, UNITS = "cleaned", "units"
CLEAN_FILES = ("presence.parquet", "summary.csv", "summary_by_group.csv", "major.parquet", "taxa.parquet")
UNIT_FILES = ("members.parquet", "presence.parquet", "areas.parquet", "adjacency.parquet", "excluded.parquet", "stats.parquet", "removed_by_region.csv", "cells.parquet")
NAME_COLUMNS = ("kingdom", "order", "family", "genus")  # what `matrix.duplicate_names` needs of the taxa
MERGED_SUFFIX = "merged-names"
# what counts as a row's evidence (DECISIONS, "Species-days", 2026-11-06): its records, or its visits, which only a download with visits has
EVIDENCE = ("records", "visits")
VISITS_SUFFIX = "visits"
# the counts that are records among `n`, scaled to visits by a row's visits per record when visits are the evidence
RECORD_SUBSETS = ("n_nc", "n_ebird", "n_inat", "n_synonym", "n_human_text")
MAP_FILES = ("partition.parquet", "structure.json", "candidates.parquet", "parts.parquet")
SCREEN_FILES = (*MAP_FILES, "screen.csv")
SELECT_FILES = (*MAP_FILES, "selection.csv", "selected.json")
SELECT_DIR = "select"  # a select result is in `<configuration>/select/`, apart from the screen's files in `<configuration>/`
HALF = "A"  # the one half this module runs on; `Plan` refuses any other until M3
BUDGET_GB = 40.0  # the most that the parent and all its workers may hold together, the default of `Plan.budget_gb` (`--budget-gb`)
MARGIN = 1.25  # a worker's later tasks may need more than the first one measured

CANDIDATE_SCHEMA = {"unit": pl.String, "k": pl.Int64, "realm": pl.Int64, "sub_realm": pl.Int64, "sub_realm_of": pl.Int64}
PARTS_SCHEMA = {"part": pl.String, "k": pl.Int64, "realm": pl.Int64, "unit": pl.String, "imputed": pl.Boolean, "selected": pl.Boolean}
UNLINKED = "unlinked"  # why a unit that Infomap gave no module to is left out of the maps

# The version of what each stage does. Raise a stage's number when its behaviour changes with its inputs the same (a fix to `clean`, to `network_units`, to the
# clustering or the scores): the new number enters that stage's hash, so it and every stage after it run again. It is in the hash only when it is not 1, so that
# the hashes of every run made before it existed stay what they were.
METHOD_VERSION = {"clean": 1, "units": 1, "screen": 1, "select": 1}


def _versioned(inputs: dict, stage: str) -> dict:
    """`inputs` with the stage's `METHOD_VERSION`, if that is not 1."""
    version = METHOD_VERSION[stage]
    return inputs | ({"method_version": version} if version != 1 else {})


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
    """The named `InfomapOptions` of a `configs.toml` (`[configs.<name>]`, each key a field of `InfomapOptions`), and its `trials_screen` and `n_null_screen`.

    A name may not end in `-s<N>` or be `select`, which are the directories of a seed (`config_dir`) and of a select result."""
    data = tomllib.loads(Path(path).read_text())
    named = data.pop("configs", {})
    clashing = sorted(n for n in named if n == SELECT_DIR or re.fullmatch(r".*-s\d+", n))
    if clashing:
        raise ValueError(f"{path} has configurations named like a seed's or the select directory: {clashing}")
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


def run_name(mode: str, rules_off: Iterable[str] = (), merge_duplicate_names: bool = False, half: str = HALF, evidence: str = "records") -> str:
    """`native_only`, or `native_only-no-griis` with the rules switched off in sorted order, and `-merged-names` after them for a run that merges duplicate
    genus and family names: the directory of a run. `-visits` follows for a run whose evidence is visits. A half other than A ends it with `-half-<half>`,
    so that such a run never shares half A's directory."""
    if evidence not in EVIDENCE:
        raise ValueError(f"the evidence is one of {EVIDENCE}, not {evidence!r}")
    if mode not in WCVP_MODES:
        raise ValueError(f"the WCVP mode is one of {WCVP_MODES}, not {mode!r}")
    unknown = sorted(set(rules_off) - set(SWITCHABLE))
    if unknown:
        raise ValueError(f"rules to switch off are of {SWITCHABLE}, not {unknown}")
    other_half = [f"half-{half}"] if half != HALF else []
    visits = [VISITS_SUFFIX] if evidence == "visits" else []
    return "-".join([mode, *(f"no-{rule}" for rule in sorted(set(rules_off))), *([MERGED_SUFFIX] if merge_duplicate_names else []), *visits, *other_half])


def config_dir(name: str, seed: int = 1) -> str:
    """The directory of a configuration in a run: its name, and `-s<seed>` after it for a seed other than 1, so that another seed never overwrites seed 1's."""
    return name if seed == 1 else f"{name}-s{seed}"


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


def clean_inputs(
    settings: Settings, layers: Iterable[int], list_codes: pl.DataFrame, marker: dict, manifest: str | None = None, merge_duplicate_names: bool = False,
    half: str = HALF, split_sha256: str | None = None, evidence: str = "records",
) -> dict:
    """What the cleaning stage's hash is of: the settings (the mode and the rules switched off are in them), the layers, the content of the frozen code table, the
    counts' `_COMPLETE` marker and the references' manifest (as a digest). `merge_duplicate_names` is in it only when it is on, as the taxa then carry names.

    A half other than A is in it with the sha256 of the split file (`split_sha256`, then required), so that its stamps can never be half A's. Half A adds nothing,
    which keeps the hashes of M1. `METHOD_VERSION["clean"]` is in it when it is not 1."""
    if half != HALF and not split_sha256:
        raise ValueError(f"half {half!r} needs the sha256 of the split file in the hash")
    inputs = {
        "settings": asdict(settings) | {k: sorted(getattr(settings, k)) for k in ("rules_off", "land_classes", "marine_exceptions")},
        "layers": sorted(layers),
        "list_codes": digest(list_codes.sort(list_codes.columns).rows()),
        "counts": marker,
        "manifest": manifest,
    }
    inputs |= {"merge_duplicate_names": True} if merge_duplicate_names else {}
    inputs |= {"evidence": evidence} if evidence != "records" else {}  # records add nothing, which keeps the hashes of M1
    inputs |= {"half": half, "split": split_sha256} if half != HALF else {}
    return _versioned(inputs, "clean")


def units_inputs(parent: str, unit_settings: UnitSettings) -> dict:
    """What the units stage's hash is of: the cleaning stage's hash and the unit settings, and `METHOD_VERSION["units"]` when it is not 1."""
    return _versioned({"cleaned": parent, "unit_settings": asdict(unit_settings)}, "units")


def invalidate_after(out: Path, stage: str) -> list[Path]:
    """Delete the stamps of every stage of the run `out` that depends on `stage`, which is about to write its outputs again, so that they report not current.

    Cleaning takes the units and every configuration's screen and select with it, the units every configuration's; a screen or a select takes nothing, as the
    two have directories of their own. Returns the stamps deleted."""
    out = Path(out)
    stamps: list[Path] = []
    if stage == "clean":
        stamps.append(out / UNITS / "stage.json")
    if stage in ("clean", "units"):
        stamps += sorted(out.glob("*/stage-screen.json")) + sorted(out.glob(f"*/{SELECT_DIR}/stage-select.json"))
    deleted = [path for path in stamps if path.exists()]
    for path in deleted:
        path.unlink()
        log.info(f"{path} deleted: {stage} is written again")
    return deleted


def _elapsed(started: float) -> dict:
    return {"seconds": round(time.perf_counter() - started, 1), "peak_rss_gb": round(peak_gb(), 2)}


# ------------------------------------------------------------------ cleaned counts


def clean_stage(
    out: Path, counts: pl.DataFrame, refs, grid: Grid, settings: Settings, layers: Iterable[int], list_codes: pl.DataFrame, inputs: dict, *, force: bool = False
) -> dict:
    """`clean` on the counts, and the parts of the references that the later stages use, into `<out>/cleaned/`. Skipped if the hash of `inputs` is the stamp's.

    When `inputs` merge duplicate names (`clean_inputs`), the taxa keep the `NAME_COLUMNS` as well, which `matrix.duplicate_names` needs.
    """
    directory = Path(out) / CLEANED
    hash_ = digest(inputs)
    if not force and is_current(directory, "stage", hash_, CLEAN_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "stage.json").unlink(missing_ok=True)
    invalidate_after(Path(out), "clean")
    cleaned = clean(counts, refs, grid, frozenset(layers), settings, list_codes)
    species = cleaned.presence.select("species_id").unique()
    major = major_groups(refs.taxa, refs.groups).join(species, on="species_id", how="semi").sort("species_id")
    kept = ("species_id", "genus_id", "family_id", *(NAME_COLUMNS if inputs.get("merge_duplicate_names") else ()))
    taxa = refs.taxa.select(*kept).join(species, on="species_id", how="semi").sort("species_id")
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
    inputs = units_inputs(parent, unit_settings)
    hash_ = digest(inputs)
    if not force and is_current(directory, "stage", hash_, UNIT_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "stage.json").unlink(missing_ok=True)
    invalidate_after(Path(out), "units")
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


def network_matrix(out: Path, units: Units | None = None, merge_duplicate_names: bool = False) -> CellTaxa:
    """The cell-taxon matrix of the units of the network (they are its cells), with the genus and family of each species from the cleaned stage.

    With `merge_duplicate_names`, the ids of a duplicate genus or family name are one column (`matrix.cell_taxa`), which needs the names that only the cleaning of
    a run that merges keeps.
    """
    units = units or read_units(out)
    taxa = pl.read_parquet(Path(out) / CLEANED / "taxa.parquet")
    missing = [c for c in NAME_COLUMNS if c not in taxa.columns]
    if merge_duplicate_names and missing:
        raise ValueError(f"{Path(out) / CLEANED / 'taxa.parquet'} has no {missing}, which merging duplicate names needs: clean with --merge-duplicate-names")
    return cell_taxa(units.presence.select(pl.col("unit").alias("cell"), "species_id"), taxa, merge_duplicate_names=merge_duplicate_names)


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


@dataclass(frozen=True)
class Scored:
    """What the maps are scored on: the matrix, areas and adjacency of the units that have a module, and the ones Infomap gave none (sorted)."""

    ct: CellTaxa
    areas: dict[str, float]
    adjacency: dict[str, set[str]]
    unassigned: list[str]


def assigned_only(partition: Partition, ct: CellTaxa, units: Units) -> Scored:
    """The matrix, areas and adjacency without the units that Infomap put in no module (`Partition.unassigned`), so that the scored cells are the mapped ones.

    Such a unit has no label in the maps, so it must not count in the endemism scores, the null maps or the bootstraps either. Its row is gone from the
    matrix, and it is in no neighbour set; a taxon column that no remaining unit holds goes too, with the genus and family columns renumbered. The units
    themselves are left as they are: the parts of an unassigned unit get their realm from `labelled_units`.
    """
    lost = set(partition.unassigned)
    if not lost:
        return Scored(ct, units.areas, units.adjacency, [])
    keep = np.array([c not in lost for c in ct.cells], dtype=bool)
    x = ct.x[np.flatnonzero(keep)]
    used = np.flatnonzero(x.getnnz(axis=0) > 0)
    renumbered = np.full(x.shape[1], -1, dtype=np.int64)
    renumbered[used] = np.arange(len(used))

    def moved(parent: np.ndarray) -> np.ndarray:
        parent = parent[used]
        return np.where(parent >= 0, renumbered[np.maximum(parent, 0)], -1)

    x = x[:, used]
    x.sort_indices()
    restricted = CellTaxa([c for c in ct.cells if c not in lost], [ct.taxa[j] for j in used.tolist()], ct.rank[used], x, moved(ct.genus_col), moved(ct.family_col))
    return Scored(
        restricted,
        {u: a for u, a in units.areas.items() if u not in lost},
        {u: near - lost for u, near in units.adjacency.items() if u not in lost},
        sorted(lost),
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


UNASSIGNED_WARNING_SHARE = 0.01  # more than this share of the units without a module is worth a warning


def structure(partition: Partition, areas: dict[str, float], n_units: int) -> dict:
    """What Infomap found: modules (the finest), the top-level ones, the levels, the largest module's share of units and area, and single-unit modules.

    `unassigned` counts the units that are in no module and `unassigned_units` lists them: they are out of the scores (`assigned_only`). If they are more than
    1% of `units` there is also a `warning`.
    """
    paths = partition.paths
    lost = {"unassigned": len(partition.unassigned), "unassigned_units": list(partition.unassigned)}
    if len(partition.unassigned) > UNASSIGNED_WARNING_SHARE * n_units:
        lost["warning"] = f"{len(partition.unassigned)} of {n_units} units are in no module, more than 1%: they are out of the scores and the bootstraps, and imputed in parts.parquet"
    if not paths:
        return {
            "units": n_units, **lost, "modules": 0, "top_modules": 0, "levels": 0, "single_unit_modules": 0,
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
        **lost,
        "modules": modules,
        "top_modules": top_modules,
        "levels": max(len(path) for path in paths.values()),
        "single_unit_modules": sum(1 for n in leaf.values() if n == 1),
        "largest_module_share_of_units": round(max(leaf.values()) / len(paths), 4),
        "largest_module_share_of_area": round(leaf_area, 4),
        "largest_top_module_share_of_units": round(max(top.values()) / len(paths), 4),
        "largest_top_module_share_of_area": round(top_area, 4),
    }


def _cluster(units: Units, ct: CellTaxa, options: InfomapOptions, trials: int, params: Params, timings: dict) -> tuple[Partition, CandidateMaps | None, str | None, Scored]:
    """Infomap, then the candidate maps on the units that have a module (`assigned_only`), which are returned for the scores that follow.

    The reason is the message of a `DegenerateClustering`, or of a network that has no link at all, and the maps are then `None`.
    """
    started = time.perf_counter()
    net = network_from(ct)
    partition = run_infomap(net, seed=params.seed, trials=trials, options=options)
    timings["infomap"] = round(time.perf_counter() - started, 2)
    scored = assigned_only(partition, ct, units)
    if not partition.paths:
        return partition, None, "no unit has a link to cluster on", scored
    started = time.perf_counter()
    try:
        found = candidate_maps(partition, scored.ct, scored.areas, scored.adjacency, k_max=params.k_max, k_min=params.k_min, floor_km2=params.floor_km2)
    except DegenerateClustering as error:
        return partition, None, str(error), scored
    finally:
        timings["candidate_maps"] = round(time.perf_counter() - started, 2)
    return partition, found, None, scored


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
    pl.concat(part_rows(found.maps, units, grid, selected)).write_parquet(directory / "parts.parquet")


def part_rows(maps: dict[int, dict[str, int]], units: Units, grid: Grid, selected: int | None) -> list[pl.DataFrame]:
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


def stage_dir(out: Path, name: str, stage: str, seed: int = 1) -> Path:
    """Where a stage of a configuration writes: `<out>/<config_dir>/` for the screen, `<out>/<config_dir>/select/` for the select stage."""
    directory = Path(out) / config_dir(name, seed)
    return directory / SELECT_DIR if stage == "select" else directory


def _fresh_dir(directory: Path, stage: str) -> Path:
    """The stage's directory, made if need be, without the stage's stamp: a stage that fails half way leaves no stamp over the files it wrote."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"stage-{stage}.json").unlink(missing_ok=True)
    return directory


def _log(directory: Path, stage: str, entry: dict) -> None:
    path = directory / "log.json"
    data = json.loads(path.read_text()) if path.exists() else {}
    data[stage] = entry
    path.write_text(json.dumps(data, indent=2, default=str) + "\n")


def _table_rows(table: dict[int, float], name: str) -> pl.DataFrame:
    return pl.DataFrame({"k": list(table), name: list(table.values())}, schema={"k": pl.Int64, name: pl.Float64}).sort("k")


def config_inputs(parent: str, options: InfomapOptions, params: Params, stage: str, merge_duplicate_names: bool = False) -> dict:
    """What a configuration's stamp is a hash of: the units stage, the options, the parameters (the seed among them) and which stage it is, `merge_duplicate_names`
    only when it is on, as the matrix then merges them, and the stage's `METHOD_VERSION` only when it is not 1."""
    inputs = {"units": parent, "options": asdict(options), "params": asdict(params), "stage": stage}
    return _versioned(inputs | ({"merge_duplicate_names": True} if merge_duplicate_names else {}), stage)


def screen_stage(
    out: Path, name: str, units: Units, ct: CellTaxa, grid: Grid, options: InfomapOptions, params: Params, parent: str, *, workers: int = 1, force: bool = False,
    merge_duplicate_names: bool = False, budget_gb: float = BUDGET_GB,
) -> dict:
    """Infomap with `trials_screen` trials, `candidate_maps`, and the endemism excess of every k over `n_null_screen` null maps, into `<out>/<name>/` (`config_dir`).

    `parent` is the hash of the units stage, and `merge_duplicate_names` says that `ct` merges them (`network_matrix`). A degenerate clustering is written to
    `structure.json` and never raised. The select stage writes elsewhere (`stage_dir`), so a screen leaves a select result alone.
    """
    inputs = config_inputs(parent, options, params, "screen", merge_duplicate_names)
    hash_ = digest(inputs)
    directory = stage_dir(out, name, "screen", params.seed)
    if not force and is_current(directory, "stage-screen", hash_, SCREEN_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    _fresh_dir(directory, "screen")
    timings: dict = {}
    seen: list[dict] = []
    used: list[int] = []
    partition, found, reason, scored = _cluster(units, ct, options, params.trials_screen, params, timings)
    if found is None:
        pl.DataFrame(schema={"k": pl.Int64, "excess": pl.Float64}).write_csv(directory / "screen.csv")
    else:
        begun = time.perf_counter()
        with measured_workers(workers, seen, used, budget_gb):
            excess = excess_table(found.maps, scored.ct, scored.adjacency, scored.areas, n_null=params.n_null_screen, seed=params.seed, workers=workers)
        timings["excess_table"] = round(time.perf_counter() - begun, 2)
        _table_rows(excess, "excess").write_csv(directory / "screen.csv")
    shape = structure(partition, units.areas, len(units.areas))
    _write_maps(directory, "screen", params.trials_screen, partition, shape, found, reason, units, grid, None)
    entry = {**_elapsed(started), "stages": timings, "workers_requested": workers, "workers_used": max(used, default=1), "workers": seen, "degenerate": reason is not None}
    _log(directory, "screen", entry)
    write_stamp(directory, "stage-screen", hash_, inputs, {k: v for k, v in entry.items() if k != "workers"})
    return {"skipped": False, "hash": hash_, "degenerate": reason is not None, **entry}


def select_stage(
    out: Path, name: str, units: Units, ct: CellTaxa, grid: Grid, options: InfomapOptions, params: Params, parent: str, *, workers: int = 1, force: bool = False,
    merge_duplicate_names: bool = False, budget_gb: float = BUDGET_GB,
) -> dict:
    """Infomap with `trials_select` trials, `candidate_maps`, `bootstrap_stability`, `excess_table` over `n_null_select` null maps, and `select_level`, into
    `<out>/<name>/select/` (`stage_dir`).

    `parent` is the hash of the units stage, and `merge_duplicate_names` says that `ct` merges them (`network_matrix`). A degenerate clustering is written to
    `structure.json` and `selected.json` and never raised.
    """
    inputs = config_inputs(parent, options, params, "select", merge_duplicate_names)
    hash_ = digest(inputs)
    directory = stage_dir(out, name, "select", params.seed)
    if not force and is_current(directory, "stage-select", hash_, SELECT_FILES):
        return {"skipped": True, "hash": hash_}
    started = time.perf_counter()
    _fresh_dir(directory, "select")
    timings: dict = {}
    seen: list[dict] = []
    used: list[int] = []
    partition, found, reason, scored = _cluster(units, ct, options, params.trials_select, params, timings)
    chosen: int | None = None
    if found is None:
        pl.DataFrame(schema={"k": pl.Int64, "stability": pl.Float64, "excess": pl.Float64, "selected": pl.Boolean}).write_csv(directory / "selection.csv")
        selected = {"degenerate": True, "reason": reason}
    else:
        with measured_workers(workers, seen, used, budget_gb):
            begun = time.perf_counter()
            stability = bootstrap_stability(
                scored.ct, found.maps, scored.areas, scored.adjacency, n_boot=params.n_boot, seed=params.seed, trials=params.boot_trials,
                floor_km2=params.floor_km2, k_max=params.k_max, options=options, workers=workers,
            )
            timings["bootstrap_stability"] = round(time.perf_counter() - begun, 2)
            begun = time.perf_counter()
            excess = excess_table(found.maps, scored.ct, scored.adjacency, scored.areas, n_null=params.n_null_select, seed=params.seed, workers=workers)
            timings["excess_table"] = round(time.perf_counter() - begun, 2)
        table = _table_rows(stability, "stability").join(_table_rows(excess, "excess"), on="k")
        choice = select_level(table)
        chosen = choice.k
        table.with_columns((pl.col("k") == chosen).alias("selected")).write_csv(directory / "selection.csv")
        # `units` counts the units that were scored, those Infomap gave a module (`assigned_only`); an unassigned unit is in `structure.json`
        selected = {"degenerate": False, "k": choice.k, "fell_short": choice.fell_short, "reason": choice.reason, "realms": choice.k, "units": len(scored.areas)}
    (directory / "selected.json").write_text(json.dumps(selected, indent=2) + "\n")
    shape = structure(partition, units.areas, len(units.areas))
    _write_maps(directory, "select", params.trials_select, partition, shape, found, reason, units, grid, chosen)
    entry = {**_elapsed(started), "stages": timings, "workers_requested": workers, "workers_used": max(used, default=1), "workers": seen, "degenerate": reason is not None}
    _log(directory, "select", entry)
    write_stamp(directory, "stage-select", hash_, inputs, {k: v for k, v in entry.items() if k != "workers"})
    return {"skipped": False, "hash": hash_, "degenerate": reason is not None, **({"k": chosen} if chosen is not None else {}), **entry}


# ------------------------------------------------------------------ reading a select result back


@dataclass(frozen=True)
class Selected:
    """A select result, read back by `read_selected`: its directory and stamp, the selected k and whether it fell short of the stability bar, `parts.parquet`,
    `candidates.parquet`, `selection.csv`, the partition of `partition.parquet` and whether the run merges duplicate names."""

    directory: Path
    stamp: dict
    k: int
    fell_short: bool
    parts: pl.DataFrame
    candidates: pl.DataFrame
    selection: pl.DataFrame
    partition: Partition
    merged: bool


def read_partition(directory: Path) -> Partition:
    """The `Partition` of a stage's `partition.parquet`: the path of a unit is (top, leaf), and a unit with no leaf is unassigned."""
    table = pl.read_parquet(Path(directory) / "partition.parquet")
    return Partition(
        paths={u: (top, leaf) for u, leaf, top in table.filter(pl.col("leaf").is_not_null()).iter_rows()},
        unassigned=sorted(table.filter(pl.col("leaf").is_null())["unit"].to_list()),
    )


def read_selected(run_dir: Path, config: str, configs_path: Path, *, seed: int = 1, merge_duplicate_names: bool | None = None) -> Selected:
    """The select result of configuration `config` (of `configs_path`) in the run `run_dir`, if it is current, else a `ValueError` that says why.

    Current means what `run_pipeline` would find: the units stage is of the cleaning stage's hash, and the select stamp has the hash of `config_inputs` over the units
    stamp's hash, the configuration's options, the parameters of `configs_path` at `seed`, and `merge_duplicate_names` (by default, whether the run's name ends in
    `-merged-names`), with all its files there. Also checked: the result is not degenerate, and `parts.parquet` marks the k of `selected.json`.
    """
    run_dir = Path(run_dir)
    configs, params = load_configs(configs_path)
    if config not in configs:
        raise ValueError(f"{configs_path} has no configuration {config!r}")
    merge = run_dir.name.endswith(f"-{MERGED_SUFFIX}") if merge_duplicate_names is None else merge_duplicate_names
    cleaned, units = read_stamp(run_dir / CLEANED), read_stamp(run_dir / UNITS)
    if not cleaned or not units:
        raise ValueError(f"{run_dir} has no stamp of the cleaning or the units stage")
    if units["inputs"].get("cleaned") != cleaned["hash"]:
        raise ValueError(f"{run_dir}: the units stage is of cleaning {units['inputs'].get('cleaned')}, the cleaning stage is {cleaned['hash']}")
    directory = stage_dir(run_dir, config, "select", seed)
    hash_ = digest(config_inputs(units["hash"], configs[config], replace(params, seed=seed), "select", merge))
    if not is_current(directory, "stage-select", hash_, SELECT_FILES):
        found = read_stamp(directory, "stage-select")
        why = "has no select stamp" if not found else f"has a select stamp of {found.get('hash')}, not {hash_}, or misses one of {list(SELECT_FILES)}"
        raise ValueError(f"{directory} {why}: it is not the current select result of {config} at seed {seed}")
    stamp = read_stamp(directory, "stage-select")
    chosen = json.loads((directory / "selected.json").read_text())
    if chosen.get("degenerate"):
        raise ValueError(f"{directory} is degenerate: {chosen}")
    parts = pl.read_parquet(directory / "parts.parquet")
    flagged = parts.filter(pl.col("selected"))["k"].unique().to_list()
    if flagged != [chosen["k"]]:
        raise ValueError(f"{directory}: selected.json says k={chosen['k']}, parts.parquet marks {flagged}")
    return Selected(
        directory=directory, stamp=stamp, k=chosen["k"], fell_short=chosen["fell_short"], parts=parts, candidates=pl.read_parquet(directory / "candidates.parquet"),
        selection=pl.read_csv(directory / "selection.csv"), partition=read_partition(directory), merged=merge,
    )


# ------------------------------------------------------------------ the move of select results into select/ (once, for the runs of M1)

MIGRATED_NOTE = """This configuration was run before the select stage had a directory of its own.

Screen and select wrote the same map files here (partition.parquet, structure.json, candidates.parquet, parts.parquet), and the select stage, which ran
last, overwrote the screen's; log.json held an entry of each. explore.migrate_select moved the select result into select/ ({moved}) and deleted
stage-screen.json, so that a screen runs again and writes its files here, and the pipeline no longer takes select's files for the screen's.

screen.csv was left in place: it is still the screen's record (the endemism excess of every k at the screen's trials and null maps). The partition and maps it was
scored on were overwritten and are gone.
"""


def migrate_select(root: Path) -> list[dict]:
    """Move the select result of every configuration under `root` (`<root>/res<r>/<run>/<config>/` with a `stage-select.json`) into `<config>/select/`.

    The select files and its stamp move (`log.json` is split: its `select` entry goes with them), `stage-screen.json` is deleted, since select overwrote the screen's
    map files, and a `MIGRATED.txt` says so; `screen.csv` stays, as the screen's record. Idempotent: a configuration without a select stamp at its top is left alone,
    and one with a select stamp both at its top and in `select/` is an error. Returns what was done, per configuration.
    """
    done = []
    for stamp in sorted(Path(root).glob("res*/*/*/stage-select.json")):
        directory = stamp.parent
        target = directory / SELECT_DIR
        if (target / "stage-select.json").exists():
            raise ValueError(f"{directory} has a select stamp both here and in {SELECT_DIR}/: sort it out by hand")
        target.mkdir(exist_ok=True)
        moved = []
        for name in (*SELECT_FILES, "stage-select.json"):
            if (directory / name).exists():
                (directory / name).rename(target / name)
                moved.append(name)
        log_path = directory / "log.json"
        if log_path.exists():
            entries = json.loads(log_path.read_text())
            if "select" in entries:
                (target / "log.json").write_text(json.dumps({"select": entries.pop("select")}, indent=2, default=str) + "\n")
                moved.append("log.json (its select entry)")
                log_path.write_text(json.dumps(entries, indent=2, default=str) + "\n")
        screen = directory / "stage-screen.json"
        removed = screen.exists()
        screen.unlink(missing_ok=True)
        (directory / "MIGRATED.txt").write_text(MIGRATED_NOTE.format(moved=", ".join(moved)))
        done.append({"directory": str(directory), "moved": moved, "screen_stamp_removed": removed})
        log.info(f"{directory}: moved {moved} into {SELECT_DIR}/" + (", deleted stage-screen.json" if removed else ""))
    return done


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
    evidence: str = "records"


def evidence_counts(counts: pl.DataFrame, evidence: str) -> pl.DataFrame:
    """The counts the cleaning reads, with `COUNTS` and nothing else. With `records`, they are as downloaded (visits dropped, if the download has them).

    With `visits` (DECISIONS, "Species-days", 2026-11-06) a row's visits stand in for its records: `n` becomes `n_visits` and `n_open` `n_visits_open`, and
    the counts that are records among `n` (`RECORD_SUBSETS`, `n_human_text` among them, which layer 1 subtracts) are scaled by the row's visits per record
    and rounded, as the download has no visits for them. `n_raw` and the distinct counts stay."""
    if evidence not in EVIDENCE:
        raise ValueError(f"the evidence is one of {EVIDENCE}, not {evidence!r}")
    keep = [c for c in counts.columns if c not in VISITS]
    if evidence == "records":
        return counts.select(keep)
    missing = [c for c in VISITS if c not in counts.columns]
    if missing:
        raise ValueError(f"visits are the evidence, and the counts have no {missing}: the download was made without visits")
    ratio = pl.when(pl.col("n") > 0).then(pl.col("n_visits") / pl.col("n")).otherwise(0.0)
    scaled = [(pl.col(c) * ratio).round(0).cast(pl.Int64).alias(c) for c in RECORD_SUBSETS]
    return counts.with_columns(*scaled, pl.col("n_visits").alias("n"), pl.col("n_visits_open").alias("n_open")).select(keep)


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
    counts = evidence_counts(load_country_counts(spec.counts, split, "A"), spec.evidence)
    del split
    return clean_stage(spec.out, counts, refs, grid, spec.settings, spec.layers, list_codes, spec.inputs, force=spec.force)


def manifest_path(reference: Path, item: str) -> Path:
    """The path of `item` in the references' `manifest.json` (the grids, the WGSRPD regions, the land), or a `ValueError` that says to fetch it.

    Unlike `natives._manifest_path`, which takes the manifest's content and returns None for a missing item, this reads the file and raises."""
    manifest = json.loads((reference / "manifest.json").read_text())
    entry = manifest.get(item)
    if not entry or not entry.get("path"):
        raise ValueError(f"{reference / 'manifest.json'} has no {item}: run scripts/fetch_reference.py")
    return Path(entry["path"])


def units_job(spec: UnitsSpec) -> dict:
    """Load the grids and the WGSRPD regions (the res-5 file and the level 3 polygons from the manifest), and run `units_stage`. Nothing is loaded if it is current."""
    directory = Path(spec.out) / UNITS
    hash_ = digest(units_inputs(spec.parent, spec.unit_settings))
    if not spec.force and is_current(directory, "stage", hash_, UNIT_FILES):
        return {"skipped": True, "hash": hash_}
    land = spec.reference / "natural_earth" / "ne_10m_land.zip"
    grid = load_grid(manifest_path(spec.reference, f"isea3h_{spec.resolution}"), land)
    coarse_gpkg = manifest_path(spec.reference, f"isea3h_{spec.unit_settings.coarse_resolution}")
    coarse = load_grid(coarse_gpkg, land)
    regions = part_regions(grid, manifest_path(spec.reference, "wgsrpd_level3"))
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
    merge_duplicate_names: bool = False  # a run of its own: `matrix.cell_taxa` merges the ids of a duplicate genus or family name
    seed: int = 1  # Infomap's, the bootstraps' and the null maps' seed (`Params.seed`); another seed writes to `<configuration>-s<seed>/`
    budget_gb: float = BUDGET_GB  # the most the parent and the workers of the bootstraps and the null maps may hold together
    half: str = HALF  # half A only: any other is refused here, so nothing can run on it by accident
    evidence: str = "records"  # or "visits", for a download made with them: a run of its own, `<run>-visits`

    def __post_init__(self) -> None:
        if self.half != HALF:
            raise ValueError(f"half {self.half!r} is sealed until M3: the exploration runs on half {HALF} only")
        if self.seed < 1:
            raise ValueError(f"the seed is 1 or more, not {self.seed}")
        if self.budget_gb <= 0:
            raise ValueError(f"the memory budget is more than 0 GB, not {self.budget_gb}")


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
    params = replace(params, seed=plan.seed)
    missing = sorted(set(plan.only) - set(configs))
    if missing:
        raise ValueError(f"--only names {missing}, and {plan.configs} has {sorted(configs)}")
    names = list(plan.only) or list(configs)
    settings = Settings(wcvp_mode=plan.mode, rules_off=frozenset(plan.rules_off))
    merge = plan.merge_duplicate_names
    out = Path(plan.results) / f"res{resolution}" / run_name(plan.mode, plan.rules_off, merge, plan.half, plan.evidence)
    out.mkdir(parents=True, exist_ok=True)
    manifest = plan.reference / "manifest.json"
    split_sha256 = hashlib.sha256(Path(plan.split).read_bytes()).hexdigest() if plan.half != HALF else None  # the split file's bytes, never a half of the counts
    inputs = clean_inputs(
        settings, LAYERS, list_codes, marker, hashlib.sha256(manifest.read_bytes()).hexdigest() if manifest.exists() else None, merge, plan.half, split_sha256,
        plan.evidence,
    )
    wanted = STAGES.index(plan.stage)
    force = {stage: plan.force and stage == plan.stage for stage in STAGES}
    call = in_child if plan.isolate else (lambda function, spec: (function(spec), peak_gb()))
    done: dict = {"directory": out}

    begun = time.perf_counter()
    result, peak = call(clean_job, CleanSpec(plan.counts, plan.reference, plan.split, plan.frozen_dir, out, settings, tuple(LAYERS), inputs, force["clean"], plan.evidence))
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
            directory = stage_dir(out, name, stage, params.seed)
            return is_current(directory, f"stage-{stage}", digest(config_inputs(parent, configs[name], params, stage, merge)), files[stage])

        needed = [n for n in names if plan.force or not all(current(n, stage) for stage in stages)]
        units = read_units(out)
        grid = read_light_grid(out)
        ct = network_matrix(out, units, merge) if needed else None
        for name in names:
            done[name] = {}
            if name not in needed:
                done[name] = {stage: {"skipped": True} for stage in stages}
                continue
            done[name]["screen"] = screen_stage(
                out, name, units, ct, grid, configs[name], params, parent, workers=plan.workers, force=force["screen"], merge_duplicate_names=merge,
                budget_gb=plan.budget_gb,
            )
            if wanted >= 3:
                done[name]["select"] = select_stage(
                    out, name, units, ct, grid, configs[name], params, parent, workers=plan.workers, force=force["select"], merge_duplicate_names=merge,
                    budget_gb=plan.budget_gb,
                )
    return done
