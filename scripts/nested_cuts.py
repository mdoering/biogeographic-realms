"""The cuts of the nested levels by the two rules (DECISIONS, "Nested levels: the cut rules", 2026-11-06).

`uv run python -I scripts/nested_cuts.py --counts data/gbif/<key> [--base-cut 0.3] [--n-boot 199] [--n-null 99] [--workers 4]`

Reads the select stage at seed 1 of one configuration (default the choice over seeds) of one run (default native_only), as `scripts/nested.py` does.
Infomap's finest modules, their pieces resolved, are merged to the realms at `--base-cut` (`nested.leaf_realms`), which are then merged land with land
and sea with sea (`nested.nested_merges`). Every place in that merge order where a cut could fall is a candidate level:

- **rule 2** (`nested.robust_levels`): the merges are replayed on `--n-boot` resamples of the species (`nested.height_replicates`); a level is robust
  where one cut separates the merges made from the merges to come in at least 95% of them. Its cut is the middle of the reference's gap.
- **rule 6** (`nested.realm_excess`): at each robust level, each realm's endemic share against the same realm in `--n-null` random maps of the same
  realm areas and pieces; a realm passes at p <= 0.05, and a level passes when all its realms do.

Writes to `res<r>/<run>/nested-cuts/` (or `--out`), with `run.log`: `merges.csv`, `heights.csv` (resamples x merges), `levels.csv` (every candidate,
rule 2), `excess.csv` (rule 6, per robust level and realm) and `cuts.json`: the robust levels that pass, their cuts, and the command of
`scripts/nested.py` that draws them, and `reported`, the levels a report gives (`nested.reported_levels`). Half A only. It chooses by the rules and nothing else.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

from realms.candidates import resolve_pieces
from realms.explore import assigned_only, network_matrix, peak_gb, read_selected, read_units
from realms.nested import height_replicates, leaf_realms, nested_merges, realm_excess, reported_levels, robust_levels, sea_realms
from realms.occurrences import IncompleteIngest, resolution_of
from realms.seeds import read_choice

log = logging.getLogger("nested_cuts")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only")
    parser.add_argument("--config", default=None, help="the Infomap configuration, default the one of <run>/seeds/choice.json")
    parser.add_argument("--base-cut", type=float, default=0.3, help="the units merge freely below this, and the candidate levels lie above it")
    parser.add_argument("--n-boot", type=int, default=199)
    parser.add_argument("--n-null", type=int, default=99)
    parser.add_argument("--share", type=float, default=0.95, help="rule 2: the share of resamples in which one cut must fit")
    parser.add_argument("--alpha", type=float, default=0.05, help="rule 6: a realm passes at p <= alpha")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    run_dir = args.results / f"res{resolution}" / args.run
    out = args.out or run_dir / "nested-cuts"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
                        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="a")])
    config = args.config or (read_choice(run_dir) or {}).get("config")
    if not config:
        raise SystemExit(f"FAILED no --config and no choice over seeds in {run_dir}")
    log.info(f"---- nested cuts: resolution {resolution}, run {args.run}, configuration {config}, base cut {args.base_cut}, "
             f"{args.n_boot} resamples (share {args.share}), {args.n_null} null maps (alpha {args.alpha})")
    started = time.perf_counter()

    found = read_selected(run_dir, config, args.configs)
    units = read_units(run_dir)
    ct = network_matrix(run_dir, units, found.merged)
    scored = assigned_only(found.partition, ct, units)
    del ct
    leaves = resolve_pieces(found.partition.leaf_labels(), scored.adjacency, scored.ct)
    base, _ = leaf_realms(leaves, scored.ct, args.base_cut)
    merges, maps = nested_merges(base, scored.ct, sea_realms(base, scored.areas))
    merges.write_csv(out / "merges.csv")
    log.info(f"{len(set(leaves.values()))} units, {len(set(base.values()))} realms at the base cut, {merges.height} nested merges")

    within = merges.filter(pl.col("kind") != "both")
    replicates = height_replicates(scored.ct, base, within, n_boot=args.n_boot, workers=args.workers)
    pl.DataFrame(replicates, schema=[f"m{i}" for i in range(within.height)], orient="row").write_csv(out / "heights.csv")
    levels = robust_levels(within, replicates, args.share)
    lo, hi = np.quantile(replicates, [0.025, 0.975], axis=0)
    log.info(f"rule 2 on {args.n_boot} resamples: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")
    with pl.Config(tbl_rows=80, float_precision=3, tbl_width_chars=200):
        log.info(f"merge heights, reference and 95% of resamples:\n{within.select('k', 'height', 'keep', 'drop').with_columns(pl.Series('q025', lo), pl.Series('q975', hi))}")
        log.info(f"candidate levels:\n{levels}")

    start_k = int(within["k"][0])
    rows = []
    for level in levels.filter(pl.col("robust")).iter_rows(named=True):
        labels = maps[start_k - level["done"]]
        table = realm_excess(labels, scored.ct, scored.adjacency, scored.areas, n_null=args.n_null, seed=1)
        rows.append(table.with_columns(pl.lit(level["realms"]).alias("realms"), pl.lit(level["cut"]).alias("cut")))
        log.info(f"rule 6 at {level['realms']} realms (cut {level['cut']:.3f}): {(table['p'] <= args.alpha).sum()} of {table.height} realms pass; "
                 f"{time.perf_counter() - started:.0f} s")
    excess = pl.concat(rows) if rows else pl.DataFrame()
    excess.write_csv(out / "excess.csv")
    passing = excess.group_by("realms", "cut").agg((pl.col("p") <= args.alpha).all().alias("passes")) if rows else pl.DataFrame({"realms": [], "cut": [], "passes": []})
    levels = levels.join(passing.select("realms", "passes"), on="realms", how="left")
    levels.write_csv(out / "levels.csv")
    with pl.Config(tbl_rows=80, float_precision=3, tbl_width_chars=200):
        log.info(f"rule 6, per realm:\n{excess}")

    chosen = levels.filter(pl.col("robust") & pl.col("passes").fill_null(False)).sort("cut", descending=True)
    cuts = [round(c, 3) for c in chosen["cut"].to_list()]
    reported = reported_levels(chosen["realms"].to_list())
    command = (f"uv run python -I scripts/nested.py --counts {args.counts} --config {config} --leaf-cut {cuts[-1]} --cuts {' '.join(map(str, cuts[:-1]))}"
               if cuts else None)
    (out / "cuts.json").write_text(json.dumps({"config": config, "base_cut": args.base_cut, "n_boot": args.n_boot, "share": args.share, "n_null": args.n_null,
                                               "alpha": args.alpha, "levels": chosen.select("realms", "cut").to_dicts(), "reported": reported, "nested_command": command}, indent=2) + "\n")
    log.info(f"levels by both rules: {chosen.select('realms', 'cut').rows()}; reported (the coarsest, the end of its run, the finest): {reported}; {command}")
    log.info(f"nested cuts: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


if __name__ == "__main__":
    main()
