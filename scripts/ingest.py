"""Stream a GBIF SQL download of occurrence counts into Parquet and check it against the reference data.

Writes the Parquet parts to `data/gbif/<key>/` and the QA numbers to `results/pilot/res<r>/qa.json`, and exits non-zero unless all three
pass marks hold:

1. no row has a species outside the accepted species of the COL release
2. no row has a cell code that is not in the grid
3. synonyms are between 2% and 5% of the clean records (the API's COL facets put them at 3.6% of all records on 2026-10-31)

When they hold, it writes `data/gbif/<key>/_COMPLETE` (the rows, the resolution and the zip's sha256), without which `load_presence` and
`load_country_counts` refuse the parts: a run that is killed or fails its QA leaves none.

The numbers are global: nothing lists a species, a lineage or a per-taxon figure, so half B stays unseen. The pass marks are not to be
loosened to make a run pass.

`--mark-complete` writes the marker for an ingest made before the marker existed, without ingesting again: it checks that the parts hold the rows
of `qa.json`, that the QA passed, and that every cell of the parts is in the grid of `--resolution`, and only then writes it.

The download is untrusted input, so run it as `uv run python -I scripts/ingest.py data/gbif/<key>.zip --resolution 6`.
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import polars as pl

from realms.grid import load_grid
from realms.occurrences import COMPLETE, IncompleteIngest, ingest, mark_complete, part_rows, qa, scan_counts
from realms.reference import load_taxa

SYNONYM_SHARE = (0.02, 0.05)


def problems(report: dict) -> list[str]:
    """The pass marks the QA report fails, as messages. An empty list is a pass."""
    found = []
    unaccepted = report["unaccepted"]
    if unaccepted["rows"]:
        found.append(
            f"unaccepted species: {unaccepted['species']:,} species in {unaccepted['rows']:,} rows ({unaccepted['n']:,} clean records) "
            "are not accepted species of the release"
        )
    unknown = report["unknown_cells"]
    if unknown["rows"]:
        found.append(f"unknown cells: {unknown['cells']:,} codes in {unknown['rows']:,} rows ({unknown['n']:,} clean records) are not in the grid")
    share = report["synonym_share"]
    low, high = SYNONYM_SHARE
    if share is None:
        found.append("synonym share: there are no clean records")
    elif not low <= share <= high:
        found.append(f"synonym share: {share:.2%} of the clean records is outside {low:.0%} to {high:.0%}")
    return found


def peak_rss_mb() -> float:
    """The most memory this process has used so far. `ru_maxrss` is in bytes on macOS and kilobytes on Linux."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1 << 20 if sys.platform == "darwin" else 1 << 10)


def mark_earlier_ingest(download: Path, out_dir: Path, qa_path: Path, gpkg: Path, land: Path, resolution: int) -> dict:
    """Write the marker for the parts of an ingest made before the marker existed, once they are shown to be that ingest. Ingests nothing.

    The parts must hold the rows of `qa.json` (its `rows` and its `ingest.rows`), `qa.json` must pass the marks, and every cell of the parts must be in
    the grid of `resolution`. Global numbers only.
    """
    if (out_dir / COMPLETE).exists():
        raise SystemExit(f"FAILED {out_dir / COMPLETE} is there already")
    missing = [str(p) for p in (download, qa_path, gpkg, land) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    report = json.loads(qa_path.read_text())
    try:
        rows = part_rows(out_dir)
    except FileNotFoundError as error:
        raise SystemExit(f"FAILED {error}") from error
    if not rows == report["rows"] == report["ingest"]["rows"]:
        raise SystemExit(f"FAILED the parts hold {rows:,} rows, and {qa_path} says {report['rows']:,} (ingest {report['ingest']['rows']:,})")
    failed = problems(report)
    if failed:
        raise SystemExit(f"FAILED {qa_path} fails the pass marks: " + "; ".join(failed))
    grid = load_grid(gpkg, land)
    if grid.resolution != resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {resolution}")
    cells = scan_counts(out_dir).select(pl.col("cell").unique()).collect()
    unknown = cells.join(grid.cells.select("cell"), on="cell", how="anti").height
    if unknown:
        raise SystemExit(f"FAILED {unknown:,} cell codes of the parts are not in the grid of resolution {resolution}")
    print(f"checked: {rows:,} rows as in {qa_path}, which passes; all {cells.height:,} cells of the parts are in the grid of resolution {resolution}")
    try:
        return mark_complete(out_dir, download, resolution, rows, how="mark-complete")
    except IncompleteIngest as error:
        raise SystemExit(f"FAILED {error}") from error


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("download", type=Path, help="the SQL download, data/gbif/<key>.zip")
    parser.add_argument("--resolution", type=int, required=True, help="the ISEA3H resolution of its cells")
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="the reference data, default data/reference")
    parser.add_argument("--results", type=Path, default=Path("results/pilot"), help="where res<r>/qa.json goes, default results/pilot")
    parser.add_argument(
        "--mark-complete", action="store_true", help="write the completion marker for an earlier ingest of this download, checked against its qa.json and the grid, without ingesting again"
    )
    args = parser.parse_args(argv)
    reference = args.reference
    gpkg = reference / "grids" / f"ISEA3H-{args.resolution}.gpkg"
    land = reference / "natural_earth" / "ne_10m_land.zip"
    release = reference / "exports" / "extended" / "col.zip"
    out_dir = args.download.parent / args.download.stem
    qa_path = args.results / f"res{args.resolution}" / "qa.json"

    if args.mark_complete:
        marker = mark_earlier_ingest(args.download, out_dir, qa_path, gpkg, land, args.resolution)
        print(f"wrote {out_dir / COMPLETE}: {json.dumps(marker)}")
        return

    # before the ingest, which takes minutes
    missing = [str(p) for p in (args.download, gpkg, land, release) if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")

    started = time.perf_counter()
    try:
        rows = ingest(args.download, out_dir)
    except ValueError as error:
        raise SystemExit(f"FAILED {error}") from error
    seconds = time.perf_counter() - started
    ingest_rss = peak_rss_mb()  # the high-water mark, read before the grid and the taxa are loaded
    parquet_bytes = sum(p.stat().st_size for p in out_dir.glob("part-*.parquet"))
    print(f"ingest: {rows:,} rows in {seconds:.0f} s, peak {ingest_rss:,.0f} MB, {parquet_bytes / 1e6:,.0f} MB of Parquet in {out_dir}")

    grid = load_grid(gpkg, land)
    if grid.resolution != args.resolution:
        raise SystemExit(f"FAILED {gpkg.name} holds resolution {grid.resolution}, not {args.resolution}")
    report = qa(scan_counts(out_dir), load_taxa(release), grid)
    report["ingest"] = {
        "rows": rows,
        "seconds": round(seconds, 1),
        "peak_rss_mb": round(ingest_rss),
        "parquet_bytes": parquet_bytes,
        "peak_rss_mb_with_qa": round(peak_rss_mb()),
    }
    qa_path.parent.mkdir(parents=True, exist_ok=True)
    qa_path.write_text(json.dumps(report, indent=2) + "\n")

    share = f"{report['synonym_share']:.2%}" if report["synonym_share"] is not None else "-"
    print(f"qa: {report['rows']:,} rows, {report['n_raw']:,} records, {report['n']:,} clean, synonyms {report['n_synonym']:,} ({share})")
    print(f"    unaccepted species: {report['unaccepted']['species']:,} species, {report['unaccepted']['rows']:,} rows, {report['unaccepted']['n']:,} clean records")
    print(f"    unknown cells: {report['unknown_cells']['cells']:,} codes, {report['unknown_cells']['rows']:,} rows, {report['unknown_cells']['n']:,} clean records")
    for kind, zero in report["zero_area"].items():
        print(f"    {kind} part of cells with no {kind}: {zero['cells']:,} cells, {zero['rows']:,} rows, {zero['n']:,} clean records")
    print(f"    written to {qa_path}")

    failed = problems(report)
    if failed:
        raise SystemExit("FAILED the pass marks: " + "; ".join(failed))
    marker = mark_complete(out_dir, args.download, args.resolution, rows)
    print(f"    wrote {out_dir / COMPLETE}: {marker['rows']:,} rows at resolution {marker['resolution']}, zip sha256 {marker['zip_sha256']}")
    print("ok")


if __name__ == "__main__":
    main()
