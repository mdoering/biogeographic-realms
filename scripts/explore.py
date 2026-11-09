"""Explore on half A: clean the counts, make the network units, cluster, screen the candidate maps and select one, each stage cached.

`uv run python -I scripts/explore.py --counts data/gbif/<key> --wcvp-mode native_only --configs explore/configs.toml --stage units --workers 4`

The stages run in order up to `--stage`, and a stage whose outputs are there with the same input hash is skipped (`--force` reruns the stage given). The outputs are
`results/explore/res<r>/<run>/`, where `<run>` is the mode and the rules switched off, such as `native_only` or `native_only-no-griis`, with `-merged-names` after them
under `--merge-duplicate-names`; `realms.explore` describes the files. Each configuration has a directory, `<config>/` (`<config>-s<N>/` under `--seed N` for an N
other than 1), with the screen's files, and the select stage writes into its `select/`.

- `clean`: layers 1 to 5 of `realms.clean` on half A's counts, the references' major groups and taxa. Loads the references: the child's peak was 17.9 to 24.1 GB on
  the runs of M1 at resolution 7 (`log.json` of each run), the most for `introduced` and `native_only-no-extinct-no-no_wild_range`.
- `units`: the network units (coarse units for sparse parts, the composition gate, the parts left out) and where the left-out parts are.
- `screen`: Infomap with the configuration's options, `candidate_maps` and the endemism excess of every k, at `trials_screen` trials and `n_null_screen` null maps.
- `select`: the headline run, 10 trials, 20 bootstraps of 3 trials and 99 null maps, then `select_level`.

Cleaning and the units run in a child process of their own, so that each one's peak memory is reported (`log.json` of the run, and the stamp of the stage) and the parent
holds none of it. `--workers N` is how many processes the bootstraps and the null maps may use at most: the first two tasks of each are run in two workers to measure a
worker's peak, and no more workers start than fit in `--budget-gb` (default 40) next to the parent. Each worker's peak is in the stage's `log.json`.

`--rank R` (species, genus or family) makes a run of that rank on its own, `<run>-rank-<R>`: it copies the cleaning and the units of `<run>`, which must be there,
and clusters, merges, resamples and scores on that rank's taxa alone, as a run of every rank does on species (`explore.rank_alone`; DECISIONS, "Realms from one
rank on its own", 2026-11-09).

`--seed N` runs Infomap, the bootstraps and the null maps at seed N (default 1). The seed is in the hash of the screen and the select, and a seed other than 1
writes to `<config>-s<N>/`, so it never overwrites seed 1's result. A stage that writes its outputs again deletes the stamps of the stages that depend on it.

Half A only: the one loader of counts asks for half A, and nothing here can unseal half B. The resolution is the ingest's, read from its marker (`_COMPLETE`) before anything is
loaded; the frozen table of the sub-national GRIIS lists is the one of that resolution in `reference/`.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

from realms.clean import SWITCHABLE, WCVP_MODES
from realms.explore import BUDGET_GB, STAGES, Plan, run_pipeline
from realms.occurrences import IncompleteIngest

FROZEN_DIR = Path(__file__).resolve().parent.parent / "reference"  # the repository's `reference/`, which is tracked: the frozen decisions are part of the method


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the Parquet parts of an ingested download, data/gbif/<key>")
    parser.add_argument("--wcvp-mode", choices=WCVP_MODES, required=True, help="how layer 3 uses WCVP")
    parser.add_argument("--rules-off", nargs="*", choices=SWITCHABLE, default=[], help="cleaning rules to switch off; they are named in the run's directory")
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"), help="the named Infomap configurations, default explore/configs.toml")
    parser.add_argument("--stage", choices=STAGES, required=True, help="run every stage up to this one")
    parser.add_argument("--only", nargs="*", default=[], help="the configurations to run, default all")
    parser.add_argument("--workers", type=int, default=1, help="processes for the bootstraps and the null maps, at most: fewer start if they would not fit")
    parser.add_argument("--results", type=Path, default=Path("results/explore"), help="where res<r>/<run>/ goes, default results/explore")
    parser.add_argument("--force", action="store_true", help="run the stage given again even if its inputs are unchanged")
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="default data/reference")
    parser.add_argument("--split", type=Path, default=Path("data/split/halves-315557.parquet"))
    parser.add_argument("--in-process", action="store_true", help="run cleaning and the units in this process, not in children (their peaks are then this process's)")
    parser.add_argument(
        "--merge-duplicate-names", action="store_true",
        help="merge the ids of a duplicate genus or family name in the matrix (matrix.duplicate_names): a run of its own, <run>-merged-names",
    )
    parser.add_argument("--evidence", choices=("records", "visits"), default="records",
                        help="what counts as a row's evidence: its records, or its visits for a download made with them (a run of its own, <run>-visits)")
    parser.add_argument("--eukaryotes-only", action="store_true",
                        help="leave out COL's bacteria, archaea and viruses (and species with no kingdom): a run of its own, <run>-eukaryotes")
    parser.add_argument("--rank", choices=("species", "genus", "family"), default=None,
                        help="cluster on one rank on its own (explore.rank_alone), on the cleaning and units of the run without it: a run of its own, <run>-rank-<rank>")
    parser.add_argument("--seed", type=int, default=1, help="the seed of Infomap, the bootstraps and the null maps, default 1; another writes to <config>-s<N>/")
    parser.add_argument(
        "--budget-gb", type=float, default=BUDGET_GB,
        help=f"the memory the parent and the workers may hold together, default {BUDGET_GB:g}; cleaning, in its own child before any worker, peaked at 18-24 GB at resolution 7",
    )
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("FAILED --workers is at least 1")
    if args.seed < 1:
        raise SystemExit("FAILED --seed is at least 1")
    if args.budget_gb <= 0:
        raise SystemExit("FAILED --budget-gb is more than 0")

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S", stream=sys.stdout)
    missing = [str(p) for p in (args.counts, args.split, args.configs, args.reference / "manifest.json") if not p.exists()]
    if missing:
        raise SystemExit(f"FAILED these files are missing: {', '.join(missing)}")
    plan = Plan(
        counts=args.counts, reference=args.reference, split=args.split, frozen_dir=FROZEN_DIR, results=args.results, configs=args.configs, mode=args.wcvp_mode,
        rules_off=tuple(args.rules_off), stage=args.stage, only=tuple(args.only), workers=args.workers, force=args.force, isolate=not args.in_process,
        merge_duplicate_names=args.merge_duplicate_names, seed=args.seed, budget_gb=args.budget_gb, evidence=args.evidence, eukaryotes_only=args.eukaryotes_only,
        rank=args.rank,
    )
    try:
        done = run_pipeline(plan)
    except (IncompleteIngest, FileNotFoundError, ValueError) as error:
        raise SystemExit(f"FAILED {error}") from error
    out = done.pop("directory")
    print(f"run {out}")
    for stage, result in done.items():
        print(f"{stage}: {json.dumps(result, default=str)[:600]}")


if __name__ == "__main__":
    main()
