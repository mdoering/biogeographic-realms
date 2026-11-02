"""Clean half A's counts through layers 1 to 5, once per WCVP mode, and report what each layer, rule and mode removed.

Reads the ingested download (`data/gbif/<key>`), the frozen split and the references, and for each of `--modes` runs `clean` with every
layer on. Half A only: `load_country_counts` raises for half B and nothing here asks it to. Everything goes to
`results/clean/res<r>/<mode>/`: `presence.parquet`, `removed.parquet` (every row taken out or moved, with its layer and rule),
`unplaced.parquet`, `summary.csv` (per layer, rule, group and country), `summary_by_group.csv` and `summary_by_rule.csv`, and
`suspect_cells.csv`, the cells where layer 5 took far more of the species than it did next door (listed for inspection, nothing is removed
for them). The report of the whole run is `results/clean/res<r>/report.json`, with the removals per layer and rule in records and in
presences, layer 5 in records, presences, species and cells per rule with its suspect cells, the effect of each WCVP mode, the fold of
zero-area parts, the Caspian check and the provisionally accepted species whose name also has an accepted species.

`--griis-without-status` says what layer 3 does with a GRIIS list that has no distribution table, which is 250 of the 304 lists until they have an
extended export: `skip`, the default, does not apply them and `apply` applies them name-only, as every list was before GRIIS was read by its own
establishment means. The report has what the lists say (`griis_lists`: entries, introduced, unstated, native and uncertain, for the lists with a table).

**The layer 5 thresholds are the defaults of `Settings`, which are setting 132 of the calibration** (Markus, 2026-11-01): record share 2%, cell share 20%, rate factor 0.5% and group
minimum 300 (the report has them as `layer_5_thresholds`). **Which country code each sub-national GRIIS list applies to is the frozen table of the resolution**, `reference/griis_uses_code-res<r>.csv`
(`scripts/griis_codes.py`: the code with most rows in the list's cells, derived once from half A's counts), so that this run, the calibration and the later all-species run apply the same lists; the script
refuses to start without it. The report has every place with that code and its share of the rows in the list's cells (`list_codes`), and flags those under 80%.

References are loaded once and need about 16 GB at the peak (`load_taxa`), so the modes share them. The counts are a download:
`uv run python -I scripts/clean_counts.py --counts data/gbif/<key> [--modes native_only introduced off] [--griis-without-status skip|apply]`.
Run it once for each setting with its own `--results` to have both. The resolution is the ingest's, read from its marker (`_COMPLETE`) before anything is
loaded; `--resolution`, if given, must agree with it.
"""

import argparse
import gc
import json
import logging
import resource
import sys
import time
from dataclasses import replace
from pathlib import Path

import polars as pl

from realms.clean import GRIIS_WITHOUT_STATUS, LAYERS, LIST_CODES_FILE, WCVP_MODES, Settings, clean, code_rows, frozen_list_codes, list_code_shares
from realms.clean_report import provisional_twin_records, provisional_twins, removal_summary
from realms.grid import KINDS, ZERO_KM2, cell_of, load_grid, parts
from realms.natives import load_references, name_index
from realms.occurrences import IncompleteIngest, load_country_counts, resolution_of
from realms.outliers import suspect_cells
from realms.split import load_split

FROZEN_DIR = Path(__file__).resolve().parent.parent / "reference"  # the repository's `reference/`, which is tracked: the frozen decisions are part of the method
# the layer 5 thresholds that the calibration varies, which the report states
THRESHOLDS = ("satellite_record_share", "satellite_cell_share", "under_rate_factor", "group_min_records")
# The Caspian Sea is no land in Natural Earth and no country in GADM. The box is on cell centres, wide enough for the cells of 265 km that
# overlap it, and "Caspian" below means the cells of it with 1,000 km2 of Natural Earth sea or more: anything smaller is some other water.
CASPIAN_BOX = {"lon": (44.5, 56.5), "lat": (35.0, 49.5)}
CASPIAN_SEA_KM2 = 1_000


def peak_gb() -> float:
    """The most memory this process has had resident so far. `ru_maxrss` is in bytes on macOS and kilobytes on Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)


def twins_of(col: Path) -> pl.DataFrame:
    index = name_index(col)
    twins = provisional_twins(index)
    del index
    gc.collect()
    return twins


def caspian(grid, counts: pl.DataFrame, cleaned) -> dict:
    """Whether the Caspian's cells have Natural Earth sea area, and whether their records stay in the sea part.

    Before and after are the same measure, the clean records `n - n_human_text` that layer 1 counts (`records` in the result).
    """
    box = grid.cells.filter(
        pl.col("lon").is_between(*CASPIAN_BOX["lon"]), pl.col("lat").is_between(*CASPIAN_BOX["lat"]), pl.col("sea_km2") >= ZERO_KM2
    ).sort("lat", "lon")
    big = box.filter(pl.col("sea_km2") >= CASPIAN_SEA_KM2)
    names = pl.concat([counts["part"], cleaned.presence["part"], cleaned.removed["part"], cleaned.removed["to_part"].drop_nulls()]).unique()
    where = parts(grid, names).select("part", "cell", "kind")

    def per_cell(frame: pl.DataFrame, column: str, name: str) -> pl.DataFrame:
        return (
            frame.join(where, on="part")
            .join(big.select("cell"), on="cell", how="semi")
            .group_by("cell", "kind")
            .agg(pl.col(column).sum().alias(name), pl.len().alias(f"{name}_rows"))
        )

    raw = per_cell(counts.select("part", (pl.col("n") - pl.col("n_human_text")).alias("records")), "records", "records_in")
    after = per_cell(cleaned.presence, "records", "records_out")
    folded = per_cell(cleaned.removed.filter(pl.col("rule") == "zero_area_part"), "records", "zero_area_records")
    table = (
        big.join(pl.DataFrame({"kind": list(KINDS)}), how="cross")
        .join(raw, on=["cell", "kind"], how="left")
        .join(after, on=["cell", "kind"], how="left")
        .join(folded, on=["cell", "kind"], how="left")
        .fill_null(0)
        .sort("lat", "lon", "kind")
    )
    sea = table.filter(pl.col("kind") == "sea")
    land = table.filter(pl.col("kind") == "land")
    return {
        "cells_in_box_with_sea": box.height,
        "caspian_cells": big.height,
        "sea_km2": [round(v, 1) for v in big["sea_km2"].to_list()],
        "land_km2": [round(v, 1) for v in big["land_km2"].to_list()],
        "box_cells_below_threshold_km2": [round(v, 1) for v in box.filter(pl.col("sea_km2") < CASPIAN_SEA_KM2)["sea_km2"].to_list()],
        "records_in_on_sea_parts": int(sea["records_in"].sum()),
        "records_in_on_land_parts": int(land["records_in"].sum()),
        "records_out_on_sea_parts": int(sea["records_out"].sum()),
        "records_out_on_land_parts": int(land["records_out"].sum()),
        "records_folded_by_zero_area_part": int(table["zero_area_records"].sum()),
        "cells": table.to_dicts(),
    }


def layer_5(cleaned, suspects: pl.DataFrame) -> dict:
    """What layer 5 took per rule, and the suspect cells.

    `rows` are rows of the country counts, `presences` the (species, part) pairs it left with nothing, `records` the clean records it took, `species` the
    species it took anything from and `cells` the (species, cell) pairs it judged out. `species_lost` are the species left with no presence at all.
    """
    taken = cleaned.removed.filter(pl.col("layer") == 5).with_columns(cell_of(pl.col("part")).alias("cell"))
    by_rule = taken.group_by("rule").agg(
        pl.len().alias("rows"),
        pl.struct("species_id", "part").filter(pl.col("presence_lost")).n_unique().alias("presences"),
        pl.col("records").sum().alias("records"),
        pl.col("species_id").n_unique().alias("species"),
        pl.struct("species_id", "cell").n_unique().alias("cells"),
    ).sort("rule")
    left = cleaned.presence.select("species_id").unique()
    lost = taken.select("species_id").unique().join(left, on="species_id", how="anti")
    thin = taken.filter(pl.col("rule") == "thin_long_range")
    return {
        "rows": taken.height,
        "presences": int(taken.filter(pl.col("presence_lost")).select("species_id", "part").unique().height),
        "records": int(taken["records"].sum() or 0),
        "species": taken["species_id"].n_unique(),
        "species_lost": lost.height,
        "by_rule": by_rule.to_dicts(),
        "thin_long_range_species": thin["species_id"].n_unique(),
        "suspect_cells": suspects.height,
        "top_suspect_cells": suspects.head(10).to_dicts(),
    }


def zero_area(grid, cleaned) -> dict:
    fold = cleaned.removed.filter(pl.col("rule") == "zero_area_part")
    kinds = parts(grid, fold["part"].unique()).select("part", "kind")
    by_kind = fold.join(kinds, on="part").group_by("kind").agg(
        pl.col("part").n_unique().alias("parts"), pl.len().alias("rows"), pl.col("records").sum().alias("records"), pl.col("species_id").n_unique().alias("species")
    )
    return {
        "parts": fold["part"].n_unique(),
        "rows": fold.height,
        "records": int(fold["records"].sum() or 0),
        "presences_lost": int(fold.filter(pl.col("presence_lost")).select("species_id", "part").unique().height),
        **{row["kind"]: {k: row[k] for k in ("parts", "rows", "records", "species")} for row in by_kind.iter_rows(named=True)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument("--resolution", type=int, help="the ISEA3H resolution of the counts, if stated: it must be that of the ingest's marker, which is used")
    parser.add_argument("--modes", nargs="+", choices=WCVP_MODES, default=list(WCVP_MODES), help="the WCVP modes to run, default all three")
    parser.add_argument(
        "--griis-without-status", choices=GRIIS_WITHOUT_STATUS, default="skip", help="a GRIIS list with no distribution table: `skip` it (default) or `apply` it name-only"
    )
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="default data/reference")
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--results", type=Path, default=Path("results/clean"), help="where res<r>/ goes, default results/clean")
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts, args.resolution)  # before the references, which take a minute and 16 GB
    except (IncompleteIngest, ValueError) as error:
        raise SystemExit(f"FAILED {error}") from error
    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    col = reference / "exports" / "extended" / "col.zip"
    missing = [str(p) for p in (args.counts, args.split, gpkg, land, col) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    try:
        list_codes = frozen_list_codes(resolution, FROZEN_DIR)  # the code each sub-national list applies to, frozen once from half A: before anything is loaded or written
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(f"FAILED {error}") from error
    directory = args.results / f"res{resolution}"
    directory.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(), logging.FileHandler(directory / "clean.log", mode="w")])
    say = logging.getLogger("clean_counts").info

    settings = Settings(griis_without_status=args.griis_without_status)  # the defaults of `Settings` for everything but the modes: the thresholds are setting 132
    report: dict = {
        "half": "A",
        "resolution": resolution,
        "counts": str(args.counts),
        "layers": list(LAYERS),
        "griis_without_status": args.griis_without_status,
        "layer_5_thresholds": {k: getattr(settings, k) for k in THRESHOLDS},
        "list_codes": {"file": LIST_CODES_FILE.format(resolution=resolution), "places": list_codes.height, "other_than_own": int((list_codes["code"] != list_codes["country"]).sum())},
        "modes": {},
    }
    started = time.perf_counter()

    twins = twins_of(col)
    say(f"provisional twins: {twins.height} provisionally accepted species share a name with an accepted one (peak {peak_gb():.1f} GB)")
    grid = load_grid(gpkg, land)
    if grid.resolution != resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {resolution}")
    say(f"grid: {grid.cells.height} cells at resolution {grid.resolution}")
    refs = load_references(grid, reference)
    gc.collect()
    say(f"references loaded, {refs.taxa.height} species, missing {list(refs.missing)} (peak {peak_gb():.1f} GB)")
    split = load_split(args.split)
    counts = load_country_counts(args.counts, split, "A")  # never unsealed: half B is not looked at before the freeze
    del split
    clean_rows = counts.filter((pl.col("n") - pl.col("n_human_text")) > 0)
    report["input"] = {
        "rows": counts.height,
        "species": counts["species_id"].n_unique(),
        "n_raw": int(counts["n_raw"].sum()),
        "n": int(counts["n"].sum()),
        "clean_records": int((counts["n"] - counts["n_human_text"]).sum()),
        "clean_rows": clean_rows.height,
        "clean_presences": clean_rows.select("species_id", "part").unique().height,
        "clean_species": clean_rows["species_id"].n_unique(),
    }
    say(f"counts of half A: {report['input']} (peak {peak_gb():.1f} GB)")
    # every sub-national list with the code it applies to and how much of the rows in its cells carry it, from the counts of this run (half A, which is where the table was frozen)
    shares = list_code_shares(refs.admin1_cells, code_rows(counts), list_codes)
    report["list_codes"]["places_detail"] = shares.to_dicts()
    report["list_codes"]["low_share"] = shares.filter(pl.col("low")).select("country", "subdivision").to_dicts()
    for r in shares.iter_rows(named=True):
        say(f"list code {r['country']}/{r['subdivision']} applies to {r['code']}: {r['with_code']:,} of {r['rows']:,} rows in its cells ({'no rows' if r['share'] is None else format(r['share'], '.1%')}){'  LOW' if r['low'] else ''}")
    report["provisional_twins"] = {
        "ids": twins.height,
        "with_one_accepted_twin": twins.filter(pl.col("accepted_ids") == 1).height,
        "with_several_accepted_twins": twins.filter(pl.col("accepted_ids") > 1).height,
        "in_half_a_counts": provisional_twin_records(counts, twins),
        "in_half_a_clean_rows": provisional_twin_records(clean_rows, twins),
    }
    say(f"provisional twins in half A: {report['provisional_twins']}")
    del clean_rows
    unplaced_griis = refs.griis_no_cells.height + (refs.griis_name_only_no_cells.height if args.griis_without_status == "apply" else 0)
    report["unplaced_in_references"] = {"griis_pairs": unplaced_griis, "wrims_pairs": refs.wrims_no_country.height}
    with_table = refs.griis_lists.filter(pl.col("table"))
    report["griis_lists"] = {
        "with_table": with_table.height,
        "without_table": refs.griis_lists.height - with_table.height,
        "entries_without_table": int(refs.griis_lists.filter(~pl.col("table"))["entries"].sum() or 0),
        **{k: int(with_table[k].sum() or 0) for k in ("entries", "introduced", "unstated", "native", "uncertain")},
        "pairs_applied": refs.griis.height,
        "pairs_name_only": refs.griis_name_only.height,
        "admin1_cells": refs.admin1_cells.height,
        "admin1_cells_by_via": dict(refs.admin1_cells.group_by("via").len().iter_rows()),
    }
    say(f"GRIIS: {report['griis_lists']}; without a table: {args.griis_without_status}")

    for mode in args.modes:
        say(f"--- wcvp_mode {mode}")
        begun = time.perf_counter()
        out = directory / mode
        out.mkdir(exist_ok=True)
        cleaned = clean(counts, refs, grid, frozenset(LAYERS), replace(settings, wcvp_mode=mode), list_codes)
        seconds = time.perf_counter() - begun
        cleaned.presence.write_parquet(out / "presence.parquet")
        cleaned.removed.write_parquet(out / "removed.parquet")
        cleaned.unplaced.write_parquet(out / "unplaced.parquet")
        summary = removal_summary(cleaned, refs)
        summary.write_csv(out / "summary.csv")
        by_group = removal_summary(cleaned, refs, by=("tax_group",))
        by_group.write_csv(out / "summary_by_group.csv")
        by_rule = removal_summary(cleaned, refs, by=())
        by_rule.write_csv(out / "summary_by_rule.csv")
        suspects = suspect_cells(cleaned, grid)
        suspects.write_csv(out / "suspect_cells.csv")
        survivors = cleaned.presence["species_id"].n_unique()
        entry = {
            "seconds": round(seconds, 1),
            "peak_gb": round(peak_gb(), 1),
            "presence_rows": cleaned.presence.height,
            "presence_species": survivors,
            "presence_records": int(cleaned.presence["records"].sum()),
            "removed_rows": cleaned.removed.height,
            "by_rule": by_rule.to_dicts(),
            "by_group": by_group.to_dicts(),
            "layer_4_species": cleaned.sides,
            "layer_5": layer_5(cleaned, suspects),
            "zero_area": zero_area(grid, cleaned),
            "caspian": caspian(grid, counts, cleaned),
        }
        report["modes"][mode] = entry
        say(f"{mode}: {seconds:.0f} s, peak {peak_gb():.1f} GB, presence {entry['presence_rows']} rows of {survivors} species, {entry['presence_records']} records")
        with pl.Config(tbl_rows=60, tbl_cols=12, fmt_str_lengths=40):
            say(str(by_rule.filter(~pl.col("unplaced")).drop("unplaced")))
            say(str(by_rule.filter(pl.col("unplaced")).drop("unplaced", "presences")))
        say(f"layer 4 species: {entry['layer_4_species']}")
        say(f"layer 5: {entry['layer_5']['rows']} rows, {entry['layer_5']['records']} records, {entry['layer_5']['species']} species, {entry['layer_5']['suspect_cells']} suspect cells")
        say(f"zero-area fold: {entry['zero_area']}")
        del cleaned, summary, by_group, by_rule, suspects
        gc.collect()

    report["seconds"] = round(time.perf_counter() - started, 1)
    report["peak_gb"] = round(peak_gb(), 1)
    (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    say(f"wrote {directory / 'report.json'}, {report['seconds']:.0f} s, peak {report['peak_gb']} GB")


if __name__ == "__main__":
    main()
