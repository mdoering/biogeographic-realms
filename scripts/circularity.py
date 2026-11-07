"""The circularity gate on half A: how the cleaning runs' maps differ, where their boundaries move, and whether a WCVP mode follows TDWG edges.

`uv run python -I scripts/circularity.py --counts data/gbif/<key>`

Reads the select stage of one configuration (`--config`, default `mt090`) of each run, `results/explore/res<r>/<run><suffix>/<config>/select/`, at the resolution
of the counts' `_COMPLETE` marker, with `explore.read_selected` (which refuses a result that is not current) and `seeds.read_headline`, which puts the k chosen over five seeds (`<run>/seeds/choice.json`) in place of seed 1's where the choice names the configuration; the grid and the TDWG level 3 regions come from the references' manifest. Writes to `res<r>/circularity/`:

- `runs.csv`: each run's selected k, its stability and whether it fell short of 0.8, its parts, and its land and sea boundary km at that k.
- `pairs.csv`: for every pair of runs, the V-measure at each run's selected k and at the common k (`circularity.common_k`), the parts both hold
  and those only one holds, and the km of land and sea boundary that moves (`circularity.moved`) at both.
- `moves.csv`: the moved land boundary of every pair by TDWG level 1 region (`circularity.moves_by_region`), at the selected and the common k.
- `near_by_region.csv`: by TDWG level 1 region, the km of land edges and the share near TDWG edges (`circularity.near_by_region`): every edge
  between two land parts of the base (`land edges`, and `land edges, both centres` for those with both centres in the region), each run's
  boundary at its selected k, and what each mode's boundary adds to and removes from the base's.
- `gate.csv`: `boundaries.gate` with the first run as the base, every other run a mode, TDWG level 3 regions, each map at its selected k, read
  twice: rotations kept within `land_tolerance=0.5` of the regions' boundary over land, and every rotation (`land_tolerance=None`).
- `run.log`: what was read, the headline mode (`circularity.headline`, on the 0.5 reading) and the timings.

Half A only: it reads the explore runs' outputs, which are half A's, and loads no counts.
"""

import argparse
import itertools
import json
import logging
import sys
import time
from pathlib import Path

import geopandas as gpd
import polars as pl

from realms.boundaries import boundary, edge_adjacency, edge_lengths, edge_share, gate, near_cells, rotations
from realms.circularity import (
    BASE,
    centre_regions,
    common_k,
    compare,
    headline,
    labels_at,
    land_edges,
    moved,
    moves_by_region,
    near_by_region,
)
from realms.explore import manifest_path
from realms.seeds import read_headline
from realms.grid import read_centres
from realms.natives import TDWG_LEVEL3_CODE
from realms.occurrences import IncompleteIngest, resolution_of

RUNS = ("off", "introduced", "native_only", "native_only-no-extinct-no-no_wild_range", "native_only-no-griis-no-griis_subnational")
TOLERANCES = (0.5, None)  # the first is the reading the headline rests on


def tolerance_name(tolerance: float | None) -> str:
    return "none" if tolerance is None else f"{tolerance:g}"

log = logging.getLogger("circularity")


def read_run(run_dir: Path, config: str, configs: Path) -> dict:
    """The selected map of one run and its candidate ks: the current select result (`seeds.read_headline`), whose parts must be the same at every k."""
    found = read_headline(run_dir, config, configs)
    parts = found.parts
    sets = {frozenset(group["part"]) for _, group in parts.group_by("k")}
    if len(sets) != 1:
        raise ValueError(f"{found.directory}: parts.parquet holds other parts at different k")
    stability = found.selection.filter(pl.col("k") == found.k)["stability"].item()
    return {"parts": parts, "k": found.k, "fell_short": found.fell_short, "stability": stability, "candidates": sorted(parts["k"].unique().to_list())}


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--config", default="mt090", help="the Infomap configuration of the runs, default mt090 (DECISIONS, over five seeds, 2026-11-03)")
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"), help="the named Infomap configurations, default explore/configs.toml")
    parser.add_argument("--runs", nargs="+", default=list(RUNS), help="the runs, the base (WCVP off) first")
    parser.add_argument("--suffix", default="", help="appended to each run's directory, -eukaryotes for the eukaryotes-only runs (DECISIONS, \"The scope is eukaryotes\", 2026-11-07); the runs keep their names in the outputs")
    parser.add_argument("--rotations", type=int, default=99)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    args = parser.parse_args()
    if args.runs[0] != BASE:
        raise SystemExit(f"FAILED the first run is the base, {BASE}, not {args.runs[0]}")

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    root = args.results / f"res{resolution}"
    out = root / "circularity"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="w")],
    )
    started = time.perf_counter()
    timings: dict[str, float] = {}

    def lap(name: str, since: float) -> float:
        timings[name] = round(time.perf_counter() - since, 1)
        log.info(f"{name}: {timings[name]} s")
        return time.perf_counter()

    gpkg = manifest_path(args.reference, f"isea3h_{resolution}")
    level3 = manifest_path(args.reference, "wgsrpd_level3")
    land = manifest_path(args.reference, "natural_earth_land")
    log.info(f"resolution {resolution} (from {args.counts}), configuration {args.config}, grid {gpkg}, regions {level3}, land {land}")
    log.info(f"runs {args.runs} (directories <run>{args.suffix}), base {BASE}, {args.rotations} rotations, seed {args.seed}")

    t = time.perf_counter()
    try:
        runs = {run: read_run(root / f"{run}{args.suffix}", args.config, args.configs) for run in args.runs}
    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    k = common_k({run: r["candidates"] for run, r in runs.items()}, {run: r["k"] for run, r in runs.items()})
    for run, r in runs.items():
        log.info(f"{run}: selected k={r['k']}, stability {r['stability']:.3f}, fell short {r['fell_short']}, candidates {r['candidates']}, {r['parts']['part'].n_unique()} parts")
    log.info(f"common k={k}")
    selected = {run: labels_at(r["parts"], r["k"]) for run, r in runs.items()}
    common = {run: labels_at(r["parts"], k) for run, r in runs.items()}
    t = lap("read the runs", t)

    edges = edge_lengths(gpkg)
    adjacency = edge_adjacency(edges)
    regions = gpd.read_file(level3).sort_values(TDWG_LEVEL3_CODE, kind="stable")
    region_of = centre_regions(read_centres(gpkg), regions)
    log.info(f"{len(edges)} edges; {sum(v is not None for v in region_of.values())} of {len(region_of)} cell centres in a TDWG level 1 region")
    t = lap("edges and regions", t)

    lines = {
        (at, kind, run): boundary(labels, kind, edges)
        for at, maps in (("selected", selected), ("common", common))
        for run, labels in maps.items()
        for kind in ("land", "sea")
    }
    pl.DataFrame(
        [
            (run, r["k"], r["stability"], r["fell_short"], r["parts"]["part"].n_unique(), r["parts"].filter(pl.col("k") == r["k"])["imputed"].sum(),
             lines[("selected", "land", run)]["km"].sum(), lines[("selected", "sea", run)]["km"].sum())
            for run, r in runs.items()
        ],
        schema=["run", "k", "stability", "fell_short", "parts", "imputed_parts", "land_boundary_km", "sea_boundary_km"], orient="row",
    ).write_csv(out / "runs.csv")

    pairs, moves = [], []
    for a, b in itertools.combinations(args.runs, 2):
        at_selected, at_common = compare(selected[a], selected[b]), compare(common[a], common[b])
        row = {
            "run_a": a, "run_b": b, "k_a": runs[a]["k"], "k_b": runs[b]["k"], "common_k": k,
            "v_selected": at_selected["v_measure"], "v_common": at_common["v_measure"],
            "shared_parts": at_selected["shared"], "only_a": at_selected["only_a"], "only_b": at_selected["only_b"],
        }
        for at in ("selected", "common"):
            for kind in ("land", "sea"):
                shift = moved(lines[(at, kind, a)], lines[(at, kind, b)])
                row[f"{kind}_moved_km_{at}"] = shift["km"].sum()
                if kind == "land":
                    by_region = moves_by_region(shift, region_of)
                    total = shift["km"].sum()
                    moves.append(by_region.with_columns(
                        pl.lit(a).alias("run_a"), pl.lit(b).alias("run_b"), pl.lit(at).alias("at"),
                        (pl.col("km") / total if total > 0 else pl.lit(None, dtype=pl.Float64)).alias("share"),
                    ))
        pairs.append(row)
        log.info(f"{a} / {b}: V {row['v_selected']:.3f} at k={row['k_a']}/{row['k_b']}, {row['v_common']:.3f} at k={k}; land moved {row['land_moved_km_selected']:,.0f} km, sea {row['sea_moved_km_selected']:,.0f} km")
    pl.DataFrame(pairs).write_csv(out / "pairs.csv")
    pl.concat(moves).select("run_a", "run_b", "at", "region", "km", "edges", "share").write_csv(out / "moves.csv")
    t = lap("pairs and moves", t)

    # how near the TDWG edges the land boundaries are, by region: the base rate of every edge between two land parts, each run's boundary at
    # its selected k, and what each mode adds to and removes from the base's boundary
    near = near_cells(regions, gpkg, adjacency)
    key = ["cell_a", "cell_b"]
    base_line = lines[("selected", "land", BASE)]
    possible = land_edges(edges, selected[BASE])
    tables = [("land edges", BASE, near_by_region(possible, near, region_of)), ("land edges, both centres", BASE, near_by_region(possible, near, region_of, both=True))]
    tables += [("boundary", run, near_by_region(lines[("selected", "land", run)], near, region_of)) for run in args.runs]
    for run in args.runs[1:]:
        line = lines[("selected", "land", run)]
        tables.append(("added", run, near_by_region(line.join(base_line.select(key), on=key, how="anti"), near, region_of)))
        tables.append(("removed", run, near_by_region(base_line.join(line.select(key), on=key, how="anti"), near, region_of)))
    pl.concat([table.select(pl.lit(what).alias("what"), pl.lit(run).alias("run"), pl.all()) for what, run, table in tables]).write_csv(out / "near_by_region.csv")
    t = lap("near by region", t)

    base_share = edge_share(base_line, near)
    readings = []
    for tolerance in TOLERANCES:
        name = tolerance_name(tolerance)
        turns = rotations(regions, args.rotations, args.seed, land_tolerance=tolerance, land=land)
        t = lap(f"rotations, land_tolerance={name}", t)
        rows = gate(selected, BASE, regions, turns, edges=edges, gpkg=gpkg, adjacency=adjacency)
        t = lap(f"gate, land_tolerance={name}", t)
        readings.append(rows.with_columns(pl.lit(name).alias("land_tolerance"), pl.lit(len(turns)).alias("rotations"), pl.lit(base_share).alias("base_share")))
        for mode, share, lift, q95, passes in rows.iter_rows():
            log.info(f"gate {name}: {mode} share {share:.4f} (off {base_share:.4f}), raise {lift:+.4f}, null q95 {q95:+.4f}, {'passes' if passes else 'FAILS'}")
    gates = pl.concat(readings).select("land_tolerance", "rotations", "mode", "share", "base_share", "raise", "null_q95", "passes")
    gates.write_csv(out / "gate.csv")

    names = [tolerance_name(tolerance) for tolerance in TOLERANCES]
    verdicts = {name: dict(gates.filter(pl.col("land_tolerance") == name).select("mode", "passes").iter_rows()) for name in names}
    first = names[0]
    for other in names[1:]:
        differ = [mode for mode in verdicts[first] if verdicts[first][mode] != verdicts[other][mode]]
        if differ:
            log.info(f"THE READINGS DIFFER for {differ}: {first} {[verdicts[first][m] for m in differ]}, {other} {[verdicts[other][m] for m in differ]}")
    log.info(f"headline mode (land_tolerance={first}): {headline(verdicts[first])}; " + ", ".join(f"at land_tolerance={n} it would be {headline(verdicts[n])}" for n in names[1:]))
    timings["total"] = round(time.perf_counter() - started, 1)
    log.info(f"timings {json.dumps(timings)}")


if __name__ == "__main__":
    main()
