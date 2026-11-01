"""What GRIIS took out of half A, and what to look at in it: totals, what changed between two runs, and the likely feral forms of native species.

Reads the removal log of a cleaning run (`removed.parquet` of `scripts/clean_counts.py`, any WCVP mode: the GRIIS rules take the same rows in all three) and half
A's counts, and writes to `--out` (default `results/clean/res<r>/griis_review/`):

- `griis_totals.json`: the rows, presences, records and species the two GRIIS rules took, and each rule alone (`realms.clean_report.griis_totals`);
- `griis_change.csv` with `--before`: the 20 species whose GRIIS removals changed most in records between that run's log and this one's;
- `griis_neighbours.csv`: the 20 largest GRIIS removals, by records, of species that also have clean, unlisted records in at least three neighbouring countries
  (`realms.clean_report.griis_feral_candidates`). Two countries are neighbours where GBIF has rows of both in the same or adjacent cells of the grid. These are the
  likely feral or released forms of native species, such as the Mallard in Denmark and Sweden: GRIIS calls them introduced, so **this is for reading, and no rule changes**.
  With `--min-neighbour-records` a neighbour must hold that many of the species' records (default 1, which is any), which leaves out the neighbours that are a few
  stray records (`griis_neighbours_min<N>.csv` then).

Half A only: `load_country_counts` raises for half B and nothing here asks it to. The names are COL's, from the release's export.
`uv run python -I scripts/griis_review.py --counts data/gbif/<key> --removed results/clean/res6/native_only/removed.parquet [--before results/clean-before/res6/native_only/removed.parquet]`.
"""

import argparse
import json
from pathlib import Path

import polars as pl

from realms.clean_report import griis_change, griis_feral_candidates, griis_totals
from realms.grid import load_grid
from realms.natives import SYNONYM, name_index
from realms.occurrences import load_country_counts
from realms.split import load_split

TOP = 20


def names_of(col: Path, species: pl.Series) -> pl.DataFrame:
    """`species_id` and `name` of the accepted species `species` names, from the COL export."""
    index = name_index(col).filter(pl.col("status") != SYNONYM).select("species_id", "name").unique("species_id")
    return pl.DataFrame({"species_id": species.unique()}).join(index, on="species_id", how="left")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument("--removed", type=Path, required=True, help="the removed.parquet of a run of scripts/clean_counts.py")
    parser.add_argument("--before", type=Path, help="another run's removed.parquet to compare it with")
    parser.add_argument("--min-neighbour-records", type=int, default=1, help="the records of the species a neighbouring country must hold to count, default 1 (any)")
    parser.add_argument("--resolution", type=int, default=6)
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--out", type=Path, help="default results/clean/res<r>/griis_review")
    args = parser.parse_args()

    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{args.resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    col = reference / "exports" / "extended" / "col.zip"
    missing = [str(p) for p in (args.counts, args.removed, args.split, gpkg, land, col, *([args.before] if args.before else [])) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    out = args.out or Path("results/clean") / f"res{args.resolution}" / "griis_review"
    out.mkdir(parents=True, exist_ok=True)

    removed = pl.read_parquet(args.removed)
    totals = {"removed": str(args.removed), **griis_totals(removed)}
    before = pl.read_parquet(args.before) if args.before else None
    if before is not None:
        totals["before"] = {"removed": str(args.before), **griis_totals(before)}
    (out / "griis_totals.json").write_text(json.dumps(totals, indent=2) + "\n")
    print(json.dumps(totals, indent=2))

    if before is not None:
        change = griis_change(before, removed, TOP)
        change = change.join(names_of(col, change["species_id"]), on="species_id", how="left").select("species_id", "name", "before", "after", "change")
        change.write_csv(out / "griis_change.csv")
        with pl.Config(tbl_rows=TOP + 2, tbl_cols=6, fmt_str_lengths=40):
            print(f"\nthe {TOP} species whose GRIIS removals changed most (records):\n{change}")

    grid = load_grid(gpkg, land)
    split = load_split(args.split)
    counts = load_country_counts(args.counts, split, "A")  # never unsealed
    found = griis_feral_candidates(counts, removed, grid, TOP, min_records=args.min_neighbour_records)
    found = found.join(names_of(col, found["species_id"]), on="species_id", how="left")
    found = found.with_columns(pl.col("neighbour_countries").list.join(" ")).select("species_id", "name", "country", "records", "rows", "neighbours", "neighbour_records", "neighbour_countries")
    suffix = "" if args.min_neighbour_records <= 1 else f"_min{args.min_neighbour_records}"
    found.write_csv(out / f"griis_neighbours{suffix}.csv")
    with pl.Config(tbl_rows=TOP + 2, tbl_cols=8, fmt_str_lengths=70):
        print(f"\nthe {TOP} largest GRIIS removals of species clean and unlisted in at least 3 neighbouring countries (records):\n{found}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
