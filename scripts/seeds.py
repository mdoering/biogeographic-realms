"""Stability over seeds: the configuration and k chosen on the median over seeds (DECISIONS, "Stability over five seeds", 2026-11-03).

`uv run python -I scripts/seeds.py --counts data/gbif/<key>`

Reads the select stage of each configuration (`--configs-to-judge`, default the seven distinct ones of "Infomap's scale") at each seed (`--seeds`,
default 1 to 5) of one run (`--run`, default native_only) in `results/explore/res<r>/<run>/<config>[-s<N>]/select/`, at the resolution of the counts'
`_COMPLETE` marker, with `explore.read_selected` (which refuses a result that is not current). Every pair must be there: a missing one fails the
script, naming all that are missing. Writes to `res<r>/<run>/seeds/` (or `--out`), and appends to its `run.log`:

- `per_seed.csv`: per configuration, seed and k, the stability, the excess and whether that seed's select stage chose the k.
- `across.csv`: per configuration and k, `seeds.across_seeds`: the median stability (0 where a seed has no map), the median excess over the seeds
  that have the map, how many have it, and the lowest and highest stability.
- `choice.json`: `seeds.choose` over every configuration's `across.csv` rows: the configuration, k, whether it fell short, and why.
- `agreement.csv`: the V-measure by part (`circularity.compare`) of seed 1's map of the chosen configuration at the chosen k against each other
  seed's map at that k, where it has one. Seed 1's map is the headline map.

Half A only: it reads the explore runs' outputs, which are half A's, and loads no counts.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import polars as pl

from realms.circularity import compare, labels_at
from realms.explore import Selected, read_selected
from realms.occurrences import IncompleteIngest, resolution_of
from realms.seeds import across_seeds, choose

CONFIGS = ("default", "mt075", "mt080", "mt085", "mt090", "mt095", "regularized")
log = logging.getLogger("seeds")


def read_all(run_dir: Path, configs: list[str], seeds: list[int], configs_path: Path) -> dict[str, dict[int, Selected]]:
    found: dict[str, dict[int, Selected]] = {}
    missing = []
    for config in configs:
        for seed in seeds:
            try:
                found.setdefault(config, {})[seed] = read_selected(run_dir, config, configs_path, seed=seed)
            except (ValueError, FileNotFoundError) as error:
                missing.append(f"{config} seed {seed}: {error}")
    if missing:
        raise ValueError(f"{len(missing)} of {len(configs) * len(seeds)} select results are missing or not current:\n  " + "\n  ".join(missing))
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only", help="the headline run, default native_only (DECISIONS, the circularity gate)")
    parser.add_argument("--configs-to-judge", nargs="+", default=list(CONFIGS))
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--out", type=Path, default=None, help="where the outputs go, default res<r>/<run>/seeds/")
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    run_dir = args.results / f"res{resolution}" / args.run
    out = args.out or run_dir / "seeds"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="a")],
    )
    log.info(f"---- seeds: resolution {resolution} (from {args.counts}), run {args.run}, configurations {args.configs_to_judge}, seeds {args.seeds}")
    try:
        found = read_all(run_dir, args.configs_to_judge, args.seeds, args.configs)
    except ValueError as error:
        raise SystemExit(f"FAILED {error}") from error

    per_seed = pl.concat([
        s.selection.select("k", "stability", "excess", "selected").with_columns(pl.lit(config).alias("config"), pl.lit(seed).alias("seed"))
        for config, by_seed in found.items() for seed, s in by_seed.items()
    ]).select("config", "seed", "k", "stability", "excess", "selected").sort("config", "seed", "k")
    per_seed.write_csv(out / "per_seed.csv")

    summaries = {config: across_seeds({seed: s.selection for seed, s in by_seed.items()}) for config, by_seed in found.items()}
    across = pl.concat([s.with_columns(pl.lit(config).alias("config")) for config, s in summaries.items()])
    across = across.select("config", *[c for c in across.columns if c != "config"]).sort("config", "k")
    across.write_csv(out / "across.csv")
    with pl.Config(tbl_rows=100, float_precision=4):
        log.info(f"across seeds:\n{across}")

    choice = choose(summaries)
    (out / "choice.json").write_text(json.dumps({"config": choice.config, "k": choice.k, "fell_short": choice.fell_short, "reason": choice.reason,
                                                  "seeds": args.seeds, "configs": args.configs_to_judge}, indent=2) + "\n")
    log.info(f"choice: {choice}")

    headline = found[choice.config][args.seeds[0]]
    reference = labels_at(headline.parts, choice.k)
    rows = []
    for seed, s in found[choice.config].items():
        if seed == args.seeds[0]:
            continue
        if choice.k in set(s.parts["k"].unique()):
            rows.append({"seed": seed, **compare(reference, labels_at(s.parts, choice.k))})
        else:
            rows.append({"seed": seed, "v_measure": None, "shared": None, "only_a": None, "only_b": None})
    agreement = pl.DataFrame(rows, schema={"seed": pl.Int64, "v_measure": pl.Float64, "shared": pl.Int64, "only_a": pl.Int64, "only_b": pl.Int64})
    agreement.write_csv(out / "agreement.csv")
    log.info(f"seed {args.seeds[0]}'s map of {choice.config} at k={choice.k} against the other seeds:\n{agreement}")


if __name__ == "__main__":
    main()
