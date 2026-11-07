"""Which classes carry the realm signal (DECISIONS, "Which groups carry the realm signal", 2026-11-07).

`uv run python -I scripts/group_signal.py --counts data/gbif/<key> [--run native_only] [--nested nested] [--level finest] [--n-null 99]`

Reads the realms of one level of `scripts/nested.py`'s output (default the finest) and the matrix the select stage scored, as `scripts/nested.py` does,
and the classification of the release (`reference.load_taxa`). For each class, `nested.class_signal`: the share of its species (those in at least two
cells) endemic to a realm, against 99 random maps of the same realm areas and pieces. A class is dropped only if it has at least `--min-species` such
species and is not above the null (p > `--alpha`); a smaller class is never judged, so rarely sampled groups and their endemics stay.

Writes `res<r>/<run>/<nested>/group_signal.csv` and logs the classes to drop. Half A only. It removes nothing: a run without them is a decision.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import polars as pl

from realms.candidates import resolve_pieces
from realms.explore import assigned_only, network_matrix, peak_gb, read_selected, read_units
from realms.nested import class_signal
from realms.occurrences import IncompleteIngest, resolution_of
from realms.reference import load_taxa

log = logging.getLogger("group_signal")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only")
    parser.add_argument("--nested", default="nested", help="the directory scripts/nested.py wrote in the run, default nested")
    parser.add_argument("--level", default="finest", help="the level to test against, a number (1 the coarsest) or finest")
    parser.add_argument("--n-null", type=int, default=99)
    parser.add_argument("--min-species", type=int, default=30)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    run_dir = args.results / f"res{resolution}" / args.run
    source = run_dir / args.nested
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
                        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(source / "group_signal.log", mode="w")])
    meta = json.loads((source / "levels.json").read_text())
    levels = {int(n): c for n, c in meta["levels"].items()}
    level = max(levels) if args.level == "finest" else int(args.level)
    config = meta["config"]
    log.info(f"---- group signal: resolution {resolution}, run {args.run}, configuration {config}, level {level} (cut {levels[level]}), "
             f"{args.n_null} null maps, at least {args.min_species} species in two cells or more, alpha {args.alpha}")
    started = time.perf_counter()

    found = read_selected(run_dir, config, args.configs)
    units = read_units(run_dir)
    ct = network_matrix(run_dir, units, found.merged)
    scored = assigned_only(found.partition, ct, units)
    del ct
    leaves = resolve_pieces(found.partition.leaf_labels(), scored.adjacency, scored.ct)
    # the units of the level: a unit's realm is that of its parts in parts.parquet (all of a unit's parts share it)
    parts = pl.read_parquet(source / "parts.parquet").filter((pl.col("k") == level) & pl.col("unit").is_not_null())
    unit_realm = dict(parts.group_by("unit").agg(pl.col("realm").first()).iter_rows())
    labels = {u: unit_realm[u] for u in leaves if u in unit_realm}
    if len(labels) != len(leaves):
        raise SystemExit(f"FAILED {len(leaves) - len(labels)} units have no realm at level {level}")
    manifest = json.loads((args.reference / "manifest.json").read_text())
    taxa = load_taxa(Path(manifest["col"]["path"]))
    species_class = {s: f"{k} {c}" for s, k, c in taxa.select("species_id", "kingdom", "class").iter_rows() if c}
    log.info(f"{len(labels)} units in {len(set(labels.values()))} realms: {time.perf_counter() - started:.0f} s")

    table = class_signal(labels, scored.ct, scored.adjacency, scored.areas, species_class, n_null=args.n_null, min_species=args.min_species, alpha=args.alpha)
    table = table.sort("species", descending=True)
    table.write_csv(source / "group_signal.csv")
    with pl.Config(tbl_rows=60, float_precision=3, tbl_width_chars=200):
        log.info(f"the 40 largest classes:\n{table.head(40)}")
        log.info(f"to drop ({table['drop'].sum()} of {table.height} classes, {table.filter(pl.col('drop'))['species'].sum():,} species):\n{table.filter(pl.col('drop'))}")
    testable = table.filter(pl.col("testable"))
    log.info(f"{testable.height} classes testable, {table.height - testable.height} too small to judge "
             f"({table.filter(~pl.col('testable'))['species'].sum():,} species); {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


if __name__ == "__main__":
    main()
