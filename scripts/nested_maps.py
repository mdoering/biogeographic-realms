"""Maps and a table of the nested levels written by `scripts/nested.py` (DECISIONS, "Nested levels: the cuts", 2026-11-06).

`uv run python -I scripts/nested_maps.py --counts data/gbif/<key>`

Reads `res<r>/<run>/nested/` (`parts.parquet`, `realms.csv`, `levels.json`) and the run's `cells.parquet`; the grid and the land come from the references'
manifest. Writes to `res<r>/<run>/nested/maps/`, with `run.log`:

- `level<n>.png` and `.gpkg` (`hexmaps.plot_parts`), one per level, 1 the coarsest: each realm numbered, the realms of the finest level outlined
  inside it, imputed parts hatched. The coarsest level has a colour of its own per realm (`review.distinct_colours`); each finer level takes its
  colours from the next coarser one (`review.match_colours`), so a realm keeps its parent's colour and a realm that splits off gets another;
  once the palette runs out, a realm that splits off takes a tint of its parent's colour (`review.tint`).
- `tables.md`: per level, each realm's area (parts, imputed ones included), sea share, mean centre, parent and bootstrap recovery, from `realms.csv`.

Half A only; it reads the explore run's outputs and loads no counts.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import polars as pl

from realms.boundaries import edge_adjacency, edge_lengths
from realms.explore import UNITS, manifest_path, peak_gb
from realms.hexmaps import part_geometries, plot_parts
from realms.occurrences import IncompleteIngest, resolution_of
from realms.review import distinct_colours, markdown, match_colours, part_areas, part_neighbours, realm_neighbours, tint

log = logging.getLogger("nested_maps")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only")
    parser.add_argument("--nested", default="nested", help="the directory of the run that scripts/nested.py wrote, default nested")
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    args = parser.parse_args()
    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    run_dir = args.results / f"res{resolution}" / args.run
    source = run_dir / args.nested
    out = source / "maps"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
                        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="w")])
    started = time.perf_counter()

    meta = json.loads((source / "levels.json").read_text())
    cut_of = {int(n): c for n, c in meta["levels"].items()}
    parts = pl.read_parquet(source / "parts.parquet")
    realms = pl.read_csv(source / "realms.csv")
    levels = sorted(cut_of)
    log.info(f"---- nested maps: resolution {resolution}, run {args.run}, configuration {meta['config']}, levels {cut_of}")

    gpkg = manifest_path(args.reference, f"isea3h_{resolution}")
    land = manifest_path(args.reference, "natural_earth_land")
    areas = part_areas(pl.read_parquet(run_dir / UNITS / "cells.parquet"))
    wanted = set(parts["part"])
    geoms = part_geometries(gpkg, land, wanted)
    adjacency = part_neighbours(wanted, edge_adjacency(edge_lengths(gpkg)))
    log.info(f"{len(geoms)} part geometries: {time.perf_counter() - started:.0f} s")

    labels = {n: dict(parts.filter(pl.col("k") == n).select("part", "realm").iter_rows()) for n in levels}
    finest = labels[levels[-1]]
    colours: dict[int, dict[int, str]] = {levels[0]: distinct_colours(labels[levels[0]], areas)}
    parent_of = {(n, r): p for n, r, p in realms.select("level", "realm", "parent").iter_rows() if p is not None}
    for coarser, n in zip(levels, levels[1:]):
        found = match_colours(labels[coarser], colours[coarser], labels[n], areas, realm_neighbours(labels[n], adjacency))
        # the palette runs out: a realm whose colour another realm of this level already has, and which did not have it at the coarser level,
        # takes a tint of its parent's colour instead, so an island realm reads as split off from its parent
        kept = {r: c for r, c in found.items() if colours[coarser].get(r) == c}
        for r in sorted(found):
            if r not in kept and any(c == found[r] for q, c in found.items() if q != r):
                found[r] = tint(colours[coarser][parent_of[(n, r)]], 0.55)
        colours[n] = found
    for n in levels:
        rows = parts.filter(pl.col("k") == n)
        count = len(set(labels[n].values()))
        title = (f"{args.run}, {meta['config']}, level {n} of {len(levels)}: {count} realms, cut {cut_of[n]} (beta-sim), resolution {resolution}; "
                 f"the {len(set(finest.values()))} realms of the finest level outlined, imputed parts hatched")
        sub = finest if n != levels[-1] else None
        plot_parts(geoms, labels[n], out / f"level{n}", title, sub_realms=sub, imputed=set(rows.filter(pl.col("imputed"))["part"]),
                   adjacency=adjacency, colours=colours[n])
        log.info(f"level {n} ({count} realms): colours " + ", ".join(f"{r} {c}" for r, c in sorted(colours[n].items())))

    shown = realms.with_columns((pl.col("area_km2") / 1e6).alias("area Mkm2"), (pl.col("sea_share") * 100).alias("sea %"),
                                (pl.col("recovered") * 100).alias("recovered %"))
    lines = [f"## Nested levels, {args.run} {meta['config']} (resolution {resolution})\n",
             f"Leaf cut {meta['leaf_cut']}, cuts {meta['cuts']}. Area is the realm's units (`scripts/nested.py`); the centre is the mean of its units' cell "
             "centres weighted by area; recovery is over the resamples of `nested.bootstrap_recovery` (median best Jaccard by area, and the share "
             "of resamples where it is at least 0.5).\n"]
    for n in levels:
        lines.append(f"\n### Level {n}: cut {cut_of[n]}, {realms.filter(pl.col('level') == n).height} realms\n")
        lines.append(markdown(shown.filter(pl.col("level") == n).sort("area_km2", descending=True)
                              .select("realm", "units", "area Mkm2", "sea %", "lat", "lon", "parent", "jaccard_median", "recovered %"),
                              {"area Mkm2": ".2f", "sea %": ".1f", "lat": ".1f", "lon": ".1f", "jaccard_median": ".3f", "recovered %": ".0f"}))
    (out / "tables.md").write_text("\n".join(lines) + "\n")
    log.info(f"wrote {out}: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


if __name__ == "__main__":
    main()
