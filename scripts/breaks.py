"""Does the data single out a realm count? Merge heights and split tests of each configuration's maps (DECISIONS, "Is there a natural k?", 2026-11-04).

`uv run python -I scripts/breaks.py --counts data/gbif/<key>`

Reads the select stage at seed 1 of each configuration (`--configs-to-judge`, default the seven of the five-seed rule) of one run (`--run`, default
native_only) in `results/explore/res<r>/<run>/<config>/select/`, at the resolution of the counts' `_COMPLETE` marker, with `explore.read_selected`.
From its partition and the matrix it scored (`explore.network_matrix` restricted by `explore.assigned_only`), it rebuilds the steps of
`candidates.candidate_maps`: the units (`resolve_pieces` of Infomap's finest modules) and the units above the size floor (`apply_size_floor`).
It checks that merging the floored units gives the maps of `candidates.parquet` at every k there, and fails if not. Then, for both hierarchies,
"units" (no floor) and "floored":

- `heights.csv`: `breaks.merge_heights` down to one realm, with the area of the two units merged.
- `splits.csv`: `breaks.split_tests` (both scores) with `--n-null` null maps, at every k of the floored hierarchy and up to `--split-k-max` of the units'.

Writes to `res<r>/<run>/breaks/` (or `--out`) and appends to its `run.log`. It chooses nothing, and the selection does not read it.
Half A only: it reads the explore runs' outputs, which are half A's.
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import polars as pl

from realms.breaks import merge_heights, split_tests
from realms.candidates import apply_size_floor, resolve_pieces
from realms.explore import assigned_only, load_configs, network_matrix, peak_gb, read_selected, read_units
from realms.occurrences import IncompleteIngest, resolution_of

CONFIGS = ("default", "mt075", "mt080", "mt085", "mt090", "mt095", "regularized")
log = logging.getLogger("breaks")


def same_grouping(a: dict[str, int], b: dict[str, int]) -> bool:
    """Whether two maps of the same cells group them the same way, whatever their labels."""
    if a.keys() != b.keys():
        return False
    forward: dict[int, int] = {}
    backward: dict[int, int] = {}
    return all(forward.setdefault(a[c], b[c]) == b[c] and backward.setdefault(b[c], a[c]) == a[c] for c in a)


def with_areas(table: pl.DataFrame, maps: dict[int, dict[str, int]], areas: dict[str, float]) -> pl.DataFrame:
    rows = []
    for r in table.iter_rows(named=True):
        area: dict[int, float] = {}
        for cell, unit in maps[r["k"]].items():
            area[unit] = area.get(unit, 0.0) + areas[cell]
        rows.append({"area_keep": area[r["keep"]], "area_drop": area[r["drop"]]})
    return pl.concat([table, pl.DataFrame(rows, schema={"area_keep": pl.Float64, "area_drop": pl.Float64})], how="horizontal")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only", help="the headline run, default native_only (DECISIONS, the circularity gate)")
    parser.add_argument("--configs-to-judge", nargs="+", default=list(CONFIGS))
    parser.add_argument("--split-k-max", type=int, default=30, help="the units' hierarchy is split-tested up to this k, default 30")
    parser.add_argument("--n-null", type=int, default=99)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--out", type=Path, default=None, help="where the outputs go, default res<r>/<run>/breaks/")
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    run_dir = args.results / f"res{resolution}" / args.run
    out = args.out or run_dir / "breaks"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="a")],
    )
    _, params = load_configs(args.configs)
    log.info(f"---- breaks: resolution {resolution} (from {args.counts}), run {args.run}, configurations {args.configs_to_judge}, "
             f"floor {params.floor_km2:,.0f} km2, {args.n_null} null maps, units split-tested to k={args.split_k_max}")
    selected = {}
    for config in args.configs_to_judge:
        try:
            selected[config] = read_selected(run_dir, config, args.configs)
        except (ValueError, FileNotFoundError) as error:
            raise SystemExit(f"FAILED {error}") from error
    units = read_units(run_dir)
    ct = network_matrix(run_dir, units, next(iter(selected.values())).merged)
    log.info(f"matrix: {len(ct.cells)} units, {len(ct.taxa)} taxa, peak {peak_gb():.1f} GB")

    heights, splits = [], []
    for config, found in selected.items():
        started = time.perf_counter()
        scored = assigned_only(found.partition, ct, units)
        leaves = resolve_pieces(found.partition.leaf_labels(), scored.adjacency, scored.ct)
        floored, _ = apply_size_floor(leaves, scored.areas, scored.ct, params.floor_km2)
        log.info(f"{config}: {len(set(found.partition.leaf_labels().values()))} modules, {len(set(leaves.values()))} units, "
                 f"{len(set(floored.values()))} above the floor; selected k={found.k} at seed 1")
        for hierarchy, labels in (("units", leaves), ("floored", floored)):
            table, maps = merge_heights(labels, scored.ct, k_min=1)
            if hierarchy == "floored":
                for k in sorted(found.candidates["k"].unique().to_list()):
                    direct = dict(found.candidates.filter(pl.col("k") == k).select("unit", "realm").iter_rows())
                    if not same_grouping(maps[k], direct):
                        raise SystemExit(f"FAILED {config}: the floored merge at k={k} is not the map of candidates.parquet")
                log.info(f"{config}: the floored merge gives the maps of candidates.parquet at k={sorted(found.candidates['k'].unique().to_list())}")
            heights.append(with_areas(table, maps, scored.areas).with_columns(pl.lit(config).alias("config"), pl.lit(hierarchy).alias("hierarchy")))
            top = max(maps) if hierarchy == "floored" else min(max(maps), args.split_k_max)
            ks = [k for k in range(2, top + 1) if k - 1 in maps and k in maps]
            tested = split_tests(maps, scored.ct, scored.adjacency, scored.areas, ks=ks, n_null=args.n_null, seed=params.seed, workers=args.workers)
            splits.append(tested.with_columns(pl.lit(config).alias("config"), pl.lit(hierarchy).alias("hierarchy")))
            with pl.Config(tbl_rows=80, float_precision=3, tbl_width_chars=200):
                log.info(f"{config} {hierarchy}: heights\n{heights[-1].filter(pl.col('k') <= 16).select('k', 'height', 'taxa_keep', 'taxa_drop', 'area_keep', 'area_drop')}")
                log.info(f"{config} {hierarchy}: splits\n{tested.select('k', 'score', 'parent', 'cells', 'area_small', 'observed', 'null_mean', 'null_q95', 'gain', 'z', 'p')}")
        log.info(f"{config}: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")

    order = ["config", "hierarchy"]
    pl.concat(heights).select(*order, pl.exclude(order)).write_csv(out / "heights.csv")
    pl.concat(splits).select(*order, pl.exclude(order)).write_csv(out / "splits.csv")
    log.info(f"wrote {out / 'heights.csv'} and {out / 'splits.csv'}")


if __name__ == "__main__":
    main()
