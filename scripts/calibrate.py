"""Calibrate the layer 5 thresholds on half A's WCVP plants, and check the known cases on half A's fully cleaned data; or, in a run of its own, score one setting on half B.

**The dry run** (without `--score-half-b`) reads the ingested download (`data/gbif/<key>`), the frozen split and the references, and does, for half A only (it never
loads half B):

1. the known cases of the design (`realms.calibrate.CASES`) on all five layers at the design's thresholds with `wcvp_mode` `native_only`:
   those whose family is in half A are run and pass or fail, and the others are listed as waiting for step 9, with their family and half;
2. the calibration: layers 1, 2 and 4 once and then layer 5 at each of the 81 settings (satellite record share 0.5, 1 and 2%, satellite
   cell share 2, 5 and 10%, under-reported rate factor 0.5, 1 and 2%, group minimum 300, 1,000 and 3,000 records), with layer 3 off, on the plants
   WCVP gives a native region and the other species of their taxonomic groups (the group totals need them). Each setting is measured by the share of
   outside presences it removes and the share of native presences it loses. The limit of 5% native loss is read over the whole run of layers 1, 2, 4
   and 5 counting drops only (a presence a fold moved is moved, not lost), as the controller ruled, provisionally; the chosen setting removes the
   most outside presences within it. The three readings (layer 5 alone, the whole run counting everything, the whole run counting drops only) are
   reported for the chosen setting and the design's;
3. the cross-check of layer 5 against GRIIS on half A's animals (`realms.calibrate.griis_cross_check`, `griis_cross_check.json`): layers 1, 2, 4 and 5 with layer 3
   off, at the chosen setting and at the design's, reporting the share of the animal presences that GRIIS lists as introduced (by its own establishment means,
   from the lists with distributions) that layer 5 takes, against the share of all the other animal presences. Reports only;
4. the native loss of the chosen setting per TDWG level 3 region (`regions.csv`, layer 5's reading), whose top 20 regions of at least 100 native presences
   show whether islands lose more.

It writes `results/calibrate/res<r>/`: `settings.csv` (the 81 settings), `regions.csv`, `known_cases.json`, `griis_cross_check.json`, `calibrate.log`, and `report.json` last,
with the `chosen` setting and the ingest's marker. References are loaded once and need about 16 GB at the peak.

**The scoring of half B** (`--score-half-b --setting N`) is the only call before the freeze that opens half B, listed in `DECISIONS.md`, and a run of its own, so that the
setting is Markus's choice after he has seen the dry run, and half A's data is not in memory when half B is loaded. It scores setting `N` of the 81 on half B's calibration
species: the plants and the other species of their groups, which are not all plants (some are of other kingdoms or of none, as in the group `eukaryotes`). It loads the
references and those species of half B, and nothing else: no known cases, no grid of settings, no regions, no cross-check, and nothing of half A. It prints and saves rates
and counts only, and nothing per taxon. Before anything is loaded it refuses:

- any resolution but 7, the ingest's, from its marker;
- unless the dry run at resolution 7 is complete: `res7/report.json` there, at resolution 7, with its `chosen` setting, made from the same ingest (the zip's sha256);
- if a scoring was started or done: any `res*/held_out.started` or `res*/held_out.json` under `--results`, or the fixed marker `data/split/held_out-315557.started` next to the
  frozen split, which holds whatever `--results` says. A run that fails after the load has to be decided on by Markus, and there is no flag to override it.

It writes both markers (exclusive creation, with the time, the commit, `N`, and whether `N` is the dry run's choice) before half B is loaded, and `res7/held_out.json`
(`N`, the dry run's choice, whether they are the same, the thresholds, the score and how many species were loaded) and `res7/held_out.log` after. The content of the fixed
marker is to be copied into `DECISIONS.md`, which is tracked, as soon as it is written.

`--griis-without-status skip|apply` is what layer 3 does with the GRIIS lists that have no distribution table in the known cases (see `scripts/clean_counts.py`),
`skip` by default. It changes nothing in the calibration, where layer 3 is off. The binding calibration is at the analysis resolution (7), where the counts of that download are
the input: `uv run python -I scripts/calibrate.py --counts data/gbif/<key> --resolution 7` (the dry run), and then, once, after Markus has read it,
`uv run python -I scripts/calibrate.py --counts data/gbif/<key> --resolution 7 --score-half-b --setting N`. At resolution 6 the dry run is a rehearsal.
"""

import argparse
import gc
import json
import logging
import resource
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

from realms.calibrate import (
    CASES,
    GRID,
    NATIVE_LOSS_LIMIT,
    SCORE,
    baseline,
    calibrate,
    calibration_counts,
    calibration_species,
    choose,
    griis_cross_check,
    kingdom_counts,
    known_cases,
    readings,
    region_loss,
    score,
    settings_grid,
    tdwg_regions,
    top_regions,
    within_limit,
)
from realms.clean import GRIIS_WITHOUT_STATUS, LAYERS, Settings, clean, outlier_rows
from realms.grid import load_grid
from realms.natives import load_references, name_index
from realms.occurrences import IncompleteIngest, load_country_counts, read_complete, resolution_of
from realms.split import load_split

SCORE_RESOLUTION = 7  # the resolution at which half B is scored, once
# the marker of a scoring of half B that does not depend on --results, next to the frozen split of release 315557, which it is about
SEAL_MARKER = Path(__file__).resolve().parent.parent / "data" / "split" / "held_out-315557.started"
HELD_OUT_KEYS = {"setting", "dry_run_chosen", "setting_is_dry_run_choice", *GRID, *SCORE, "species_loaded"}
MIN_NATIVE = 100  # presences a region needs to be listed in the top of the native loss
SMALL_PARTS = 5  # a region touching this many land parts or fewer is small: islands and island groups, but not only
TOP = 20
LEADERS = 10


def peak_gb() -> float:
    """The most memory this process has had resident so far. `ru_maxrss` is in bytes on macOS and kilobytes on Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)


def commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=Path(__file__).parent, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def mark_started(markers: list[Path], lines: dict[str, object]) -> None:
    """Write each of `markers`, with the time, the commit and `lines`, to say that half B is about to be loaded. They stay whether the run goes on or fails.

    `--score-half-b` refuses if any of them, or a `held_out.json`, exists, so a run that fails after the load cannot be retried by running the script again.
    Made with exclusive creation: it raises if a marker is already there, before half B is loaded.
    """
    text = "".join(f"{k} {v}\n" for k, v in {"started": datetime.now(UTC).isoformat(timespec="seconds"), "commit": commit(), **lines}.items())
    for marker in markers:
        marker.parent.mkdir(parents=True, exist_ok=True)
        with marker.open("x") as f:
            f.write(text)


def score_held_out(counts_dir: Path, split: pl.DataFrame, refs, grid, setting: Settings, markers: list[Path], lines: dict[str, object]) -> dict:
    """Score `setting` on half B's calibration species. The only call before the freeze that opens half B, listed in DECISIONS.md.

    The calibration species are the plants WCVP gives a native region **and every other species of their taxonomic groups**, which are not all plants:
    `under_reported` needs the group's records per cell, so the score is like for like with half A's. The split is cut to those species before it goes
    to `load_country_counts`, so no other species of half B is loaded. `species_loaded` says how many of them there are and how many are not Plantae,
    counted from the split and the taxonomy alone. What comes back is that and the counts and rates of `score`: nothing per taxon, and nothing is
    printed or saved by this function but the markers: they are written (`mark_started`, with `lines`) before half B is loaded, so that a run that fails after
    the load leaves a trace that `refuse_a_second_scoring` stops a retry on.
    """
    species = split.filter(pl.col("half") == "B", pl.col("species_id").is_in(calibration_species(refs).implode()))
    loaded = kingdom_counts(species["species_id"], refs.taxa)
    mark_started(markers, lines)
    counts_b = load_country_counts(counts_dir, species, "B", unsealed=True)
    return {**score(counts_b, setting, refs, grid), "species_loaded": loaded}


def described(index: int, setting: Settings, start, refs, grid, met: bool | None = None) -> dict:
    """A setting of the grid with its thresholds and its three readings (`readings`), and whether it met the limit if that is known."""
    out = {"setting": index, **{k: getattr(setting, k) for k in GRID}, "readings": readings(start, refs, grid, setting)}
    return out if met is None else {**out, "limit_met": met}


def refuse_a_second_scoring(args: argparse.Namespace, resolution: int) -> None:
    """`--score-half-b` is made once, at resolution 7: stop before anything is loaded if it is asked otherwise, or if a scoring was started or done.

    A scoring was started or done if any resolution under `--results` has `held_out.started` or `held_out.json`, or the fixed `SEAL_MARKER` exists, whatever
    `--results` says. The markers are written before half B is loaded and `held_out.json` after, so a run that fails in between is stopped too. There is no
    flag to override: a retry after a failed run is for Markus to decide.
    """
    if resolution != SCORE_RESOLUTION:
        raise SystemExit(f"FAILED --score-half-b runs at resolution {SCORE_RESOLUTION}, the analysis resolution, and not at {resolution}: the call that opens half B is made once")
    how = "The call that opens half B is made once before the freeze and there is no flag to run it again: a retry after a failed run needs Markus's explicit decision."
    if SEAL_MARKER.exists():
        raise SystemExit(f"FAILED a scoring of half B was started, see {SEAL_MARKER}: it may have loaded half B, and a second run would be a second load. {how}")
    if started := sorted(args.results.glob("res*/held_out.started")):
        raise SystemExit(f"FAILED a scoring of half B was started, see {started[0]}: it may have loaded half B and then failed, and a retry would be a second load. {how}")
    if done := sorted(args.results.glob("res*/held_out.json")):
        raise SystemExit(f"FAILED half B was already scored, see {done[0]}. {how}")


def dry_run_choice(args: argparse.Namespace, marker: dict) -> int:
    """The setting the dry run at resolution 7 chose, from its `report.json`, which it writes last. Refuses unless that run is complete and was made from the ingest of `marker`."""
    report_path = args.results / f"res{SCORE_RESOLUTION}" / "report.json"
    first = f"Run the dry run first, `scripts/calibrate.py --counts {args.counts}` without --score-half-b, and choose the setting from it"
    if not report_path.exists():
        raise SystemExit(f"FAILED there is no dry run at resolution {SCORE_RESOLUTION}: {report_path} is missing. {first}")
    report = json.loads(report_path.read_text())
    chosen = report.get("chosen", {}).get("setting")
    if report.get("resolution") != SCORE_RESOLUTION or report.get("half") != "A" or not isinstance(chosen, int):
        raise SystemExit(f"FAILED {report_path} is not a complete dry run on half A at resolution {SCORE_RESOLUTION} with a `chosen` setting. {first}")
    if report.get("ingest", {}).get("zip_sha256") != marker["zip_sha256"]:
        raise SystemExit(f"FAILED {report_path} was made from another ingest than {args.counts} (the zip's sha256 differs). {first}")
    return chosen


def score_half_b(args: argparse.Namespace, resolution: int, marker: dict) -> None:
    """The run of `--score-half-b --setting N`: the guards, then the references and half B's calibration species, and nothing else."""
    refuse_a_second_scoring(args, resolution)
    chosen = dry_run_choice(args, marker)
    grid_of_settings = settings_grid()
    if not 0 <= args.setting < len(grid_of_settings):
        raise SystemExit(f"FAILED --setting is one of the {len(grid_of_settings)} settings of the grid, 0 to {len(grid_of_settings) - 1}, not {args.setting}")
    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    missing = [str(p) for p in (args.split, gpkg, land) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    directory = args.results / f"res{resolution}"
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(), logging.FileHandler(directory / "held_out.log", mode="a")])
    say = logging.getLogger("calibrate_script").info
    setting = grid_of_settings[args.setting]
    lines = {"setting": args.setting, "dry_run_chosen": chosen, "setting_is_dry_run_choice": args.setting == chosen, "counts": args.counts, "zip_sha256": marker["zip_sha256"]}
    say(f"scoring setting {args.setting} ({', '.join(f'{k} {getattr(setting, k)}' for k in GRID)}) on half B's calibration species; the dry run chose {chosen}")

    begun = time.perf_counter()
    grid = load_grid(gpkg, land)
    if grid.resolution != resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {resolution}")
    refs = load_references(grid, reference)
    gc.collect()
    split = load_split(args.split)
    say(f"references {refs.taxa.height} species, missing {list(refs.missing)} (peak {peak_gb():.1f} GB); the one unsealed call before the freeze follows")
    score_ = score_held_out(args.counts, split, refs, grid, setting, [SEAL_MARKER, directory / "held_out.started"], lines)
    held = {
        "setting": args.setting,
        "dry_run_chosen": chosen,
        "setting_is_dry_run_choice": args.setting == chosen,
        **{k: getattr(setting, k) for k in GRID},
        **score_,
    }
    assert set(held) == HELD_OUT_KEYS  # numbers only: nothing per taxon is printed or saved
    (directory / "held_out.json").write_text(json.dumps(held, indent=2) + "\n")
    say(f"half B's calibration species, rates and counts: {held}")
    say(f"wrote {directory / 'held_out.json'}, {time.perf_counter() - begun:.0f} s, peak {peak_gb():.1f} GB. Copy {SEAL_MARKER} into DECISIONS.md now")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument(
        "--resolution", type=int, help="the ISEA3H resolution of the counts, if stated: it must be that of the ingest's marker, which is used. The binding calibration is at 7"
    )
    parser.add_argument(
        "--griis-without-status", choices=GRIIS_WITHOUT_STATUS, default="skip", help="the known cases' layer 3 on a GRIIS list with no distribution table: `skip` (default) or `apply` it name-only"
    )
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="default data/reference")
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--results", type=Path, default=Path("results/calibrate"), help="where res<r>/ goes, default results/calibrate")
    parser.add_argument(
        "--score-half-b",
        action="store_true",
        help="instead of the dry run, score setting --setting on half B's calibration species (WCVP plants and the other species of their groups): the one call before the freeze that opens half B. Resolution 7 only, after the dry run at 7, run once",
    )
    parser.add_argument("--setting", type=int, help="with --score-half-b, and only with it: the setting of the 81 to score, by its number in the dry run's settings.csv")
    args = parser.parse_args()
    if args.score_half_b and args.setting is None:
        raise SystemExit("FAILED --score-half-b needs --setting N, the setting of the dry run at resolution 7 that Markus chose")
    if args.setting is not None and not args.score_half_b:
        raise SystemExit("FAILED --setting is for --score-half-b only: the dry run chooses its setting itself")

    try:
        marker = read_complete(args.counts)
        resolution = resolution_of(args.counts, args.resolution)  # before the references, which take a minute and 16 GB
    except (IncompleteIngest, ValueError) as error:
        raise SystemExit(f"FAILED {error}") from error
    if args.score_half_b:
        score_half_b(args, resolution, marker)
        return
    dry_run(args, resolution, marker)


def dry_run(args: argparse.Namespace, resolution: int, marker: dict) -> None:
    """The known cases, the 81 settings, the regions and the cross-check, on half A. Never loads half B, and writes `report.json` last."""
    directory = args.results / f"res{resolution}"
    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    col = reference / "exports" / "extended" / "col.zip"
    level3 = reference / "wgsrpd" / "level3.geojson"
    missing = [str(p) for p in (args.counts, args.split, gpkg, land, col, level3) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "report.json").unlink(missing_ok=True)  # written last, so that it is there only for a run that finished: the scoring of half B reads it
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(), logging.FileHandler(directory / "calibrate.log", mode="w")])
    say = logging.getLogger("calibrate_script").info

    started = time.perf_counter()
    seconds: dict[str, float] = {}
    report: dict = {
        "half": "A",
        "resolution": resolution,
        "counts": str(args.counts),
        "calibration_layers": [1, 2, 4, 5],
        "native_loss_limit": NATIVE_LOSS_LIMIT,
        "known_cases_griis_without_status": args.griis_without_status,
        "ingest": marker,
    }

    def lap(name: str, since: float) -> None:
        seconds[name] = round(time.perf_counter() - since, 1)
        say(f"{name}: {seconds[name]:.0f} s (peak {peak_gb():.1f} GB)")

    begun = time.perf_counter()
    grid = load_grid(gpkg, land)
    if grid.resolution != resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {resolution}")
    refs = load_references(grid, reference)
    gc.collect()
    regions = tdwg_regions(level3)
    split = load_split(args.split)
    counts = load_country_counts(args.counts, split, "A")  # half A, not unsealed: half B is opened only by score_held_out, in a run of its own
    say(f"grid {grid.cells.height} cells at resolution {grid.resolution}; references {refs.taxa.height} species, missing {list(refs.missing)}; half A {counts.height} rows")
    lap("load", begun)

    # ---- the known cases, on everything half A has, all five layers at the design's thresholds
    begun = time.perf_counter()
    index = name_index(col)
    cleaned = clean(counts, refs, grid, frozenset(LAYERS), Settings(wcvp_mode="native_only", griis_without_status=args.griis_without_status))
    say(f"cleaned half A with all layers: {cleaned.presence.height} presences of {cleaned.presence['species_id'].n_unique()} species")
    cases = known_cases(CASES, index, split, counts, cleaned, refs, regions)
    del index, cleaned
    gc.collect()
    (directory / "known_cases.json").write_text(json.dumps(cases, indent=2) + "\n")
    report["known_cases"] = cases
    for c in cases:
        if c["status"] == "run":
            numbers = {k: v for k, v in c.items() if k.startswith(("area_", "outside_", "native_")) or k in ("species", "presences_before", "presences_after")}
            say(f"known case {c['case']} ({c['claim']}): {'PASS' if c['passed'] else 'FAIL'} {numbers}, removed by {c['removed_by_rule']}, left in {c['leftover_regions']}")
        else:
            say(f"known case {c['case']}: family {c['family']} is in half {c['half']}, {c['status']}")
    lap("known_cases", begun)

    # ---- the calibration
    begun = time.perf_counter()
    grid_of_settings = settings_grid()
    start = baseline(counts, refs, grid, grid_of_settings[0])
    plants = calibration_counts(counts, refs)
    report["population"] = {
        "species": plants["species_id"].n_unique(),
        "presences_after_layers_1_2_4": start.presence.height,
        "status_presences": start.status.height,
    }
    say(f"baseline (layers 1, 2 and 4) on {report['population']['species']} species: {start.presence.height} presences, {start.status.height} with a status")
    lap("baseline", begun)
    begun = time.perf_counter()
    table = calibrate(counts, grid_of_settings, refs, grid, start)
    lap("layer_5_per_setting", begun)
    table = table.with_columns(within_limit(reading="run").alias("within_limit"), within_limit(reading="layer5").alias("within_limit_layer5"))
    table.write_csv(directory / "settings.csv")
    index_chosen, met = choose(table)  # the whole run counting drops only
    index_layer5, met_layer5 = choose(table, reading="layer5")
    index_design = grid_of_settings.index(Settings())
    chosen = grid_of_settings[index_chosen]
    report["chosen"] = described(index_chosen, chosen, start, refs, grid, met)
    report["design_setting"] = described(index_design, grid_of_settings[index_design], start, refs, grid)
    report["chosen_by_layer5_reading"] = described(index_layer5, grid_of_settings[index_layer5], start, refs, grid, met_layer5)
    report["within_limit"] = {"run": int(table["within_limit"].sum()), "layer5": int(table["within_limit_layer5"].sum())}
    ranked = table.sort(["within_limit", "run_outside_removed", "run_native_lost", "setting"], descending=[True, True, False, False]).head(LEADERS)
    report["top_settings"] = ranked.to_dicts()
    with pl.Config(tbl_rows=LEADERS + 5, tbl_cols=24, fmt_float="mixed"):
        say(f"top {LEADERS} of the 81 settings (within the limit first, then most outside removed, whole run counting drops only):\n{ranked}")
    say(f"settings that lose at most {NATIVE_LOSS_LIMIT:.0%} of native presences: {report['within_limit']['run']} of 81 counting the whole run's drops, {report['within_limit']['layer5']} by layer 5 alone")
    for name in ("chosen", "design_setting", "chosen_by_layer5_reading"):
        say(f"{name}: {report[name]}")

    # ---- the chosen setting in detail: native loss per region, by layer 5's reading
    begun = time.perf_counter()
    removed = outlier_rows(start.presence, refs, grid, chosen)
    losses = region_loss(start.status, removed, refs.wcvp, refs.tdwg3_parts, regions)
    losses.write_csv(directory / "regions.csv")
    eligible = losses.filter(pl.col("native") >= MIN_NATIVE)
    top = top_regions(losses, MIN_NATIVE, TOP)
    by_size = {}
    for name, part in (("small", eligible.filter(pl.col("parts") <= SMALL_PARTS)), ("large", eligible.filter(pl.col("parts") > SMALL_PARTS))):
        by_size[name] = {"regions": part.height, "native": int(part["native"].sum()), "lost": int(part["lost"].sum()), "loss": part["lost"].sum() / part["native"].sum()}
    report["regions"] = {
        "with_native_presences": losses.height,
        "eligible": eligible.height,
        "min_native": MIN_NATIVE,
        "median_loss": eligible["loss"].median(),
        "small_parts": SMALL_PARTS,
        "by_size": by_size,
        "top": top.to_dicts(),
    }
    say(f"regions of at least {MIN_NATIVE} native presences by size (small: up to {SMALL_PARTS} land parts), summed over regions: {by_size}")
    with pl.Config(tbl_rows=TOP + 5, tbl_cols=8, fmt_str_lengths=30):
        say(f"the {TOP} regions with the highest native loss among {eligible.height} with at least {MIN_NATIVE} native presences:\n{top}")
    lap("chosen_detail", begun)

    # ---- layer 5 against GRIIS, on half A's animals: reports only
    begun = time.perf_counter()
    cross = {}
    for name, setting in (("chosen", chosen), ("design", grid_of_settings[index_design])):
        found = griis_cross_check(counts, refs, grid, setting)
        cross[name] = {"setting": index_chosen if name == "chosen" else index_design, **{k: getattr(setting, k) for k in GRID}, **found}
        say(
            f"GRIIS cross-check at the {name} setting: animal presences listed as introduced {found['listed']['presences_removed']} of {found['listed']['presences']} taken by layer 5 "
            f"({found['listed']['presences_share']:.2%}), the other animal presences {found['other']['presences_removed']} of {found['other']['presences']} ({found['other']['presences_share']:.2%}); "
            f"records {found['listed']['records_share']:.2%} against {found['other']['records_share']:.2%}"
        )
    (directory / "griis_cross_check.json").write_text(json.dumps(cross, indent=2) + "\n")
    report["griis_cross_check"] = cross
    lap("griis_cross_check", begun)

    report["seconds"] = {**seconds, "total": round(time.perf_counter() - started, 1)}
    report["peak_gb"] = round(peak_gb(), 1)
    (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    say(f"wrote {directory / 'report.json'}, {report['seconds']['total']:.0f} s, peak {report['peak_gb']} GB")


if __name__ == "__main__":
    main()
