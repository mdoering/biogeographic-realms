"""Calibrate the layer 5 thresholds on half A's WCVP plants, and check the known cases on half A's fully cleaned data.

Reads the ingested download (`data/gbif/<key>`), the frozen split and the references, and does, for half A only:

1. the known cases of the design (`realms.calibrate.CASES`) on all five layers at the design's thresholds with `wcvp_mode` `native_only`:
   those whose family is in half A are run and pass or fail, and the others are listed as waiting for step 9, with their family and half;
2. the calibration: layers 1, 2 and 4 once and then layer 5 at each of the 81 settings (satellite record share 0.5, 1 and 2%, satellite
   cell share 2, 5 and 10%, under-reported rate factor 0.5, 1 and 2%, group minimum 300, 1,000 and 3,000 records), with layer 3 off, on the plants
   WCVP gives a native region and the groups of those. Each setting is measured by the share of outside presences it removes and the share of native
   presences it loses; the chosen one removes the most outside presences while losing at most 5% of the native ones;
3. the native loss of the chosen setting per TDWG level 3 region (`regions.csv`), whose top 20 regions of at least 100 native presences show whether
   islands lose more;
4. with `--score-half-b`, the one setting chosen on A scored on half B's plants. That is the only call before the freeze that opens half B. It is
   listed in `DECISIONS.md`, prints global numbers only, and refuses to run if `held_out.json` is already there.

Everything goes to `results/calibrate/res<r>/`: `settings.csv` (the 81 settings), `regions.csv`, `known_cases.json`, `report.json`, `calibrate.log`,
and `held_out.json` after a scoring of half B. The binding calibration is at the analysis resolution (7), where the counts of that download are the
input: `uv run python -I scripts/calibrate.py --counts data/gbif/<key> --resolution 7 --score-half-b`. At resolution 6 it is a dry run.
References are loaded once and need about 16 GB at the peak.
"""

import argparse
import gc
import json
import logging
import resource
import sys
import time
from pathlib import Path

import polars as pl

from realms.calibrate import (
    CASES,
    EPS,
    NATIVE_LOSS_LIMIT,
    GRID,
    baseline,
    calibrate,
    calibration_counts,
    calibration_species,
    choose,
    known_cases,
    pipeline_tally,
    region_loss,
    score,
    settings_grid,
    tdwg_regions,
    top_regions,
)
from realms.clean import LAYERS, Settings, clean, outlier_rows
from realms.grid import load_grid
from realms.natives import load_references, name_index
from realms.occurrences import load_country_counts
from realms.split import load_split

MIN_NATIVE = 100  # presences a region needs to be listed in the top of the native loss
SMALL_PARTS = 5  # a region touching this many land parts or fewer is small: islands and island groups, but not only
TOP = 20
LEADERS = 10


def peak_gb() -> float:
    """The most memory this process has had resident so far. `ru_maxrss` is in bytes on macOS and kilobytes on Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)


def score_held_out(counts_dir: Path, split: pl.DataFrame, refs, grid, setting: Settings) -> dict:
    """Score `setting` on half B's plants. The only call before the freeze that opens half B, listed in DECISIONS.md.

    Only the plants are read: the split is cut to the calibration species before it goes to `load_country_counts`, so no other species of half B is
    ever loaded. What comes back is counts and rates, and nothing per taxon.
    """
    plants = split.filter(pl.col("species_id").is_in(calibration_species(refs).implode()))
    counts_b = load_country_counts(counts_dir, plants, "B", unsealed=True)
    return score(counts_b, setting, refs, grid)


def setting_row(row: dict) -> dict:
    return {k: row[k] for k in ("setting", *GRID, "outside", "outside_removed", "native", "native_lost", "outside_rate", "native_loss")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument("--resolution", type=int, default=6, help="the ISEA3H resolution of the counts, default 6. The binding calibration is at 7")
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="default data/reference")
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--results", type=Path, default=Path("results/calibrate"), help="where res<r>/ goes, default results/calibrate")
    parser.add_argument("--score-half-b", action="store_true", help="also score the chosen setting on half B's plants: the one call before the freeze that opens half B, run once")
    args = parser.parse_args()

    directory = args.results / f"res{args.resolution}"
    if args.score_half_b and (directory / "held_out.json").exists():
        raise SystemExit(f"FAILED half B's plants were already scored, see {directory / 'held_out.json'}: that call is made once before the freeze")
    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{args.resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    col = reference / "exports" / "extended" / "col.zip"
    level3 = reference / "wgsrpd" / "level3.geojson"
    missing = [str(p) for p in (args.counts, args.split, gpkg, land, col, level3) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    directory.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(), logging.FileHandler(directory / "calibrate.log", mode="w")])
    say = logging.getLogger("calibrate_script").info

    started = time.perf_counter()
    seconds: dict[str, float] = {}
    report: dict = {"half": "A", "resolution": args.resolution, "counts": str(args.counts), "calibration_layers": [1, 2, 4, 5], "native_loss_limit": NATIVE_LOSS_LIMIT}

    def lap(name: str, since: float) -> None:
        seconds[name] = round(time.perf_counter() - since, 1)
        say(f"{name}: {seconds[name]:.0f} s (peak {peak_gb():.1f} GB)")

    begun = time.perf_counter()
    grid = load_grid(gpkg, land)
    refs = load_references(grid, reference)
    gc.collect()
    regions = tdwg_regions(level3)
    split = load_split(args.split)
    counts = load_country_counts(args.counts, split, "A")  # half A, not unsealed: half B is opened only by score_held_out
    say(f"grid {grid.cells.height} cells at resolution {grid.resolution}; references {refs.taxa.height} species, missing {list(refs.missing)}; half A {counts.height} rows")
    lap("load", begun)

    # ---- the known cases, on everything half A has, all five layers at the design's thresholds
    begun = time.perf_counter()
    index = name_index(col)
    cleaned = clean(counts, refs, grid, frozenset(LAYERS), Settings(wcvp_mode="native_only"))
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
    table = table.with_columns((pl.col("native_lost") <= NATIVE_LOSS_LIMIT * pl.col("native") + EPS).alias("within_limit"))
    table.write_csv(directory / "settings.csv")
    index_chosen, met = choose(table)
    chosen = grid_of_settings[index_chosen]
    row = setting_row(table.filter(pl.col("setting") == index_chosen).row(0, named=True))
    design = setting_row(table.filter(pl.col("setting") == grid_of_settings.index(Settings())).row(0, named=True))
    report["chosen"] = {**row, "limit_met": met}
    report["design_setting"] = design
    report["within_limit"] = int(table["within_limit"].sum())
    ranked = table.sort(["within_limit", "outside_removed", "native_lost", "setting"], descending=[True, True, False, False]).head(LEADERS)
    report["top_settings"] = ranked.to_dicts()
    with pl.Config(tbl_rows=LEADERS + 5, tbl_cols=14, fmt_float="mixed"):
        say(f"top {LEADERS} of the 81 settings (within the limit first, then most outside removed):\n{ranked}")
    say(f"{report['within_limit']} of 81 settings lose at most {NATIVE_LOSS_LIMIT:.0%} of native presences; limit met by the chosen one: {met}")
    say(f"chosen setting {index_chosen}: {row}")
    say(f"design setting {grid_of_settings.index(Settings())}: {design}")

    # ---- the chosen setting in detail
    begun = time.perf_counter()
    removed = outlier_rows(start.presence, refs, grid, chosen)
    final = start.presence.join(removed.select("species_id", "part"), on=["species_id", "part"], how="anti")
    raw = plants.filter(pl.col("n_raw") > 0).select("species_id", "part")
    report["whole_run"] = pipeline_tally(raw, final, refs.wcvp, refs.tdwg3_parts)
    whole = report["whole_run"]
    say(f"the whole run (layers 1, 2, 4 and 5) on the raw presences: {whole}, native lost {whole['native_lost'] / whole['native']:.3%}, outside removed {whole['outside_removed'] / whole['outside']:.3%}")
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

    # ---- half B, once
    if args.score_half_b:
        begun = time.perf_counter()
        say("scoring the chosen setting on half B's plants: the one unsealed call before the freeze")
        held = {"setting": index_chosen, **{k: getattr(chosen, k) for k in GRID}, **score_held_out(args.counts, split, refs, grid, chosen)}
        (directory / "held_out.json").write_text(json.dumps(held, indent=2) + "\n")
        report["held_out"] = held
        say(f"half B's plants: {held}")
        lap("held_out", begun)

    report["seconds"] = {**seconds, "total": round(time.perf_counter() - started, 1)}
    report["peak_gb"] = round(peak_gb(), 1)
    (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    say(f"wrote {directory / 'report.json'}, {report['seconds']['total']:.0f} s, peak {report['peak_gb']} GB")


if __name__ == "__main__":
    main()
