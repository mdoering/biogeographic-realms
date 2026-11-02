"""Freeze where GBIF uses the country code of a sub-national GRIIS list, from half A's counts, so that every run applies the same lists.

Layer 3 decides each sub-national list (the US states, the French overseas places, and so on) either **by code**, where the list's country code occurs among
the rows in its cells, so that only rows with that code count, or **by cells alone**, where it does not (an overseas place carries the sovereign's code in its
list, `FR`, and its own in the records, `RE` or `GP`). Without a frozen table `clean` reads that from the counts it is given, so that half A, a calibration
subset and the later all-species run could decide a list differently, by what each happens to hold. This derives the decision once, from half A's counts at the
ingest's resolution, with the rule of `realms.clean.derive_uses_code`, and writes it to `reference/griis_uses_code-res<r>.csv` (`country`, `subdivision`,
`uses_code`, sorted by country and subdivision) in the repository, where it is tracked. `scripts/clean_counts.py` and `scripts/calibrate.py` read it and refuse to run
without it.

Half A only (`load_country_counts` raises for half B, and nothing here asks it to). It needs the grid and the admin-1 cells of the lists, not the references. A table that is
already there is left alone if it is the same, and **refused if it differs**, unless `--overwrite` says it is meant to: the table is a frozen decision, and a different one
changes what every later run removes.

`uv run python -I scripts/griis_codes.py --counts data/gbif/<key> [--out reference] [--overwrite]`. The resolution is the ingest's, read from its marker (`_COMPLETE`);
`--resolution`, if given, must agree with it.
"""

import argparse
import logging
import sys
from pathlib import Path

import polars as pl

from realms.clean import FROZEN_DIR, USES_CODE_FILE, codes_seen, derive_uses_code, read_uses_code
from realms.grid import load_grid
from realms.natives import load_admin1_cells
from realms.occurrences import IncompleteIngest, load_country_counts, resolution_of
from realms.split import load_split

log = logging.getLogger("griis_codes")


def derive(counts_dir: Path, resolution: int, reference: Path, split_path: Path) -> pl.DataFrame:
    """The table of `resolution` from half A's counts in `counts_dir`: `realms.clean.derive_uses_code` of the admin-1 cells and the codes half A uses there."""
    gpkg = reference / "grids" / f"ISEA3H-{resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    grid = load_grid(gpkg, land)
    if grid.resolution != resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {resolution}")
    cells = load_admin1_cells(grid, reference)
    counts = load_country_counts(counts_dir, load_split(split_path), "A")  # never unsealed: half B is not looked at before the freeze
    return derive_uses_code(cells, codes_seen(counts))


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
    path = args.out / USES_CODE_FILE.format(resolution=resolution)
    table = derive(args.counts, resolution, args.reference, args.split)
    by_code = table.filter(pl.col("uses_code")).height
    log.info("resolution %d, half A: %d places of the sub-national lists, %d decided by code and %d by cells alone", resolution, table.height, by_code, table.height - by_code)
    alone = table.filter(~pl.col("uses_code"))
    log.info("by cells alone: %s", ", ".join(f"{c}/{d}" for c, d in alone.select("country", "subdivision").iter_rows()) or "none")

    if path.exists():
        before = read_uses_code(path)
        if before.equals(table):
            log.info("%s is there and is the same: nothing written", path)
            return
        changed = before.join(table, on=["country", "subdivision"], how="full", coalesce=True, suffix="_new").filter(
            pl.col("uses_code").is_null() | pl.col("uses_code_new").is_null() | (pl.col("uses_code") != pl.col("uses_code_new"))
        )
        if not args.overwrite:
            raise SystemExit(
                f"FAILED {path} is a frozen decision and differs from the one derived ({changed.height} places, such as "
                f"{', '.join(f'{c}/{d}' for c, d in changed.select('country', 'subdivision').head(5).iter_rows())}). "
                "Keep it, or pass --overwrite to replace it and log why in DECISIONS.md"
            )
        log.warning("%s replaced: %d places differ", path, changed.height)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write_csv(path)
    log.info("wrote %s", path)


if __name__ == "__main__":
    main()
