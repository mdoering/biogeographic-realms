"""Freeze which country code each sub-national GRIIS list applies to, from half A's counts, so that every run applies the same lists.

A sub-national list (the US states, Svalbard, the French overseas places) is for the rows in the cells that overlap its subdivision, and the cells are bigger than the
places, so a border cell holds a neighbour country's rows too. The list applies to the rows that carry its **effective country code**: the code that has most rows among
the rows in its cells (counted by rows, the rows with no code not counting), with a tie going to the list's own country code and then to the lower code. For the US states
that is `US`; for Svalbard, a list of `NO` whose records GBIF codes `SJ`, it is `SJ`; for Reunion, a list of `FR` whose records are `RE`, it is `RE`. Without a frozen table
`clean` reads that off the counts it is given, so that half A, a calibration subset and the later all-species run could apply a list to different codes, by what each happens to
hold. This derives the table once, from half A's counts at the ingest's resolution, with the rule of `realms.clean.derive_list_codes`, and writes it to
`reference/griis_uses_code-res<r>.csv` (`country`, `subdivision`, `code`, sorted by country and subdivision) in the repository, where it is tracked. It prints every place with
its code and the share of the rows in its cells that carry it, and warns of each place under 80%. `scripts/clean_counts.py` and `scripts/calibrate.py` read the table and refuse
to run without it.

Half A only (`load_country_counts` raises for half B, and nothing here asks it to). It needs the grid and the admin-1 cells of the lists, not the references. A table that is
already there is left alone if it is the same, and **refused if it differs or cannot be read**, unless `--overwrite` says it is meant to: the table is a frozen decision, and a
different one changes what every later run removes.

`uv run python -I scripts/griis_codes.py --counts data/gbif/<key> [--out reference] [--overwrite]`. The resolution is the ingest's, read from its marker (`_COMPLETE`);
`--resolution`, if given, must agree with it.
"""

import argparse
import logging
import sys
from pathlib import Path

import polars as pl

from realms.clean import LIST_CODES_FILE, LOW_SHARE, code_rows, derive_list_codes, list_code_shares, read_list_codes
from realms.grid import load_grid
from realms.natives import load_admin1_cells
from realms.occurrences import IncompleteIngest, load_country_counts, resolution_of
from realms.split import load_split

FROZEN_DIR = Path(__file__).resolve().parent.parent / "reference"  # the repository's `reference/`, which is tracked: the decision is part of the method
log = logging.getLogger("griis_codes")


def inputs(counts_dir: Path, resolution: int, reference: Path, split_path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The admin-1 cells of the lists and the rows of half A's counts by country code and cell (`realms.clean.code_rows`), for `derive_list_codes`."""
    gpkg = reference / "grids" / f"ISEA3H-{resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    grid = load_grid(gpkg, land)
    if grid.resolution != resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {resolution}")
    cells = load_admin1_cells(grid, reference)
    counts = load_country_counts(counts_dir, load_split(split_path), "A")  # never unsealed: half B is not looked at before the freeze
    return cells, code_rows(counts)


def derive(counts_dir: Path, resolution: int, reference: Path, split_path: Path) -> pl.DataFrame:
    """The table of `resolution` from half A's counts in `counts_dir`: `realms.clean.derive_list_codes` of the admin-1 cells and the codes half A uses there."""
    cells, seen = inputs(counts_dir, resolution, reference, split_path)
    return derive_list_codes(cells, seen)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument("--resolution", type=int, help="the ISEA3H resolution of the counts, if stated: it must be that of the ingest's marker, which is used")
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="default data/reference")
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--out", type=Path, default=FROZEN_DIR, help="where griis_uses_code-res<r>.csv goes, default reference/ of the repository")
    parser.add_argument("--overwrite", action="store_true", help="replace a frozen table that differs from the one derived: a different decision changes what every later run removes")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stdout)

    try:
        resolution = resolution_of(args.counts, args.resolution)  # before anything is loaded
    except (IncompleteIngest, ValueError) as error:
        raise SystemExit(f"FAILED {error}") from error
    gpkg = args.reference / "grids" / f"ISEA3H-{resolution}.gpkg"
    land = args.reference / "natural_earth" / "ne_10m_land.zip"
    missing = [str(p) for p in (args.counts, args.split, gpkg, land, args.reference / "manifest.json") if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    path = args.out / LIST_CODES_FILE.format(resolution=resolution)
    cells, seen = inputs(args.counts, resolution, args.reference, args.split)
    table = derive_list_codes(cells, seen)
    shares = list_code_shares(cells, seen, table)
    log.info("resolution %d, half A: %d places of the sub-national lists, each with its effective code and the share of the rows in its cells (those with no code included) that carry it:", resolution, table.height)
    for r in shares.iter_rows(named=True):
        other = "" if r["next_code"] is None else f", next {r['next_code']} {r['next_rows']:,}"
        share = "no rows" if r["share"] is None else f"{r['share']:.1%}"
        log.info("  %s/%s %s: %s of %s rows%s%s", r["country"], r["subdivision"], r["code"], share, f"{r['rows']:,}", other, "  LOW" if r["low"] else "")
    low = shares.filter(pl.col("low"))
    if low.height:
        log.warning("%d places have an effective code with under %.0f%% of the rows in their cells: %s", low.height, LOW_SHARE * 100, ", ".join(f"{c}/{d}" for c, d in low.select("country", "subdivision").iter_rows()))
    own = table.filter(pl.col("code") != pl.col("country"))
    log.info("lists applied to another code than their own: %s", ", ".join(f"{c}/{d} as {k}" for c, d, k in own.iter_rows()) or "none")

    if path.exists():
        try:
            before = read_list_codes(path).sort("country", "subdivision")
        except ValueError as error:
            before, unreadable = None, str(error)
        if before is not None and before.equals(table.sort("country", "subdivision")):
            log.info("%s is there and is the same: nothing written", path)
            return
        if before is None:
            changed = f"it cannot be read as a table ({unreadable})"
        else:
            diff = before.join(table, on=["country", "subdivision"], how="full", coalesce=True, suffix="_new").filter(
                pl.col("code").is_null() | pl.col("code_new").is_null() | (pl.col("code") != pl.col("code_new"))
            )
            changed = f"it differs in {diff.height} places, such as {', '.join(f'{c}/{d}' for c, d in diff.select('country', 'subdivision').head(5).iter_rows())}"
        if not args.overwrite:
            raise SystemExit(f"FAILED {path} is a frozen decision, and {changed}. Keep it, or pass --overwrite to replace it and log why in DECISIONS.md")
        log.warning("%s replaced: %s", path, changed)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write_csv(path)
    log.info("wrote %s", path)


if __name__ == "__main__":
    main()
