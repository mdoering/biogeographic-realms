"""Nested levels of realms for one configuration, with how well each realm comes back on resampled species (DECISIONS, "Nested levels", 2026-11-05).

`uv run python -I scripts/nested.py --counts data/gbif/<key> --leaf-cut 0.5 --cuts 0.65 0.72`

Reads the select stage at seed 1 of one configuration (`--config`, default the one of `<run>/seeds/choice.json`) of one run (`--run`, default native_only)
with `explore.read_selected`, and the matrix it scored (`explore.network_matrix` restricted by `explore.assigned_only`). Infomap's finest modules,
their pieces resolved (`candidates.resolve_pieces`), are the units of `nested.levels_of`: the leaf realms at `--leaf-cut`, then the land and sea realms
merged apart and cut at each of `--cuts`. Writes to `res<r>/<run>/nested/` (or `--out`), and appends to its `run.log`:

- `merges.csv`: every merge, the leaf ones then the nested ones (`nested.MERGES` and `stage`), with the area and sea share of the two realms merged.
- `realms.csv`: per level (`level`, 1 the coarsest, and its `cut`), each realm's area, sea share, mean centre (latitude and longitude weighted by area)
  and its realm at the next coarser level (`parent`), and, from `--n-boot` resamples (`nested.bootstrap_recovery`, Infomap at `--boot-trials` trials),
  the median Jaccard index of its best match and the share of resamples where that is at least one half.
- `parts.parquet`: the realm of every part at every level (`explore.part_rows`, `k` holding the level), imputed parts flagged, for maps.
- `recovery.csv`: the bootstrap's rows.

Half A only. It chooses nothing: the cuts are the caller's, recorded in DECISIONS.
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
from realms.explore import UNITS, assigned_only, load_configs, network_matrix, part_rows, peak_gb, read_light_grid, read_selected, read_units
from realms.nested import SEA_SUFFIX, bootstrap_recovery, levels_of
from realms.occurrences import IncompleteIngest, resolution_of
from realms.seeds import read_choice

log = logging.getLogger("nested")


def realm_table(levels: dict[float, dict[str, int]], areas: dict[str, float], centre: dict[str, tuple[float, float]]) -> pl.DataFrame:
    cuts = sorted(levels, reverse=True)  # the coarsest level first
    rows = []
    for number, cut in enumerate(cuts, start=1):
        labels = levels[cut]
        coarser = levels[cuts[number - 2]] if number > 1 else None
        by: dict[int, list[str]] = {}
        for unit, realm in labels.items():
            by.setdefault(realm, []).append(unit)
        for realm, units in sorted(by.items()):
            w = np.array([areas[u] for u in units])
            lat = np.array([centre[u][0] for u in units])
            lon = np.array([centre[u][1] for u in units])
            sea = sum(areas[u] for u in units if u.endswith(SEA_SUFFIX))
            parents = {coarser[u] for u in units} if coarser else set()
            rows.append({
                "level": number, "cut": cut, "realm": realm, "units": len(units), "area_km2": float(w.sum()), "sea_share": sea / float(w.sum()),
                "lat": float(np.average(lat, weights=w)), "lon": float(np.average(lon, weights=w)),
                "parent": parents.pop() if len(parents) == 1 else None,
            })
    return pl.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only")
    parser.add_argument("--config", default=None, help="the Infomap configuration, default the one of <run>/seeds/choice.json")
    parser.add_argument("--leaf-cut", type=float, required=True)
    parser.add_argument("--cuts", type=float, nargs="*", default=[], help="the cuts of the levels above the leaves")
    parser.add_argument("--n-boot", type=int, default=20)
    parser.add_argument("--boot-trials", type=int, default=3)
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
    out = args.out or run_dir / "nested"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
                        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="a")])
    choice = read_choice(run_dir)
    config = args.config or (choice or {}).get("config")
    if not config:
        raise SystemExit(f"FAILED no --config and no choice over seeds in {run_dir}")
    configs, params = load_configs(args.configs)
    log.info(f"---- nested: resolution {resolution}, run {args.run}, configuration {config}, leaf cut {args.leaf_cut}, cuts {args.cuts}, "
             f"{args.n_boot} resamples of {args.boot_trials} trials")
    started = time.perf_counter()
    found = read_selected(run_dir, config, args.configs)
    units = read_units(run_dir)
    ct = network_matrix(run_dir, units, found.merged)
    scored = assigned_only(found.partition, ct, units)
    del ct
    leaves = resolve_pieces(found.partition.leaf_labels(), scored.adjacency, scored.ct)
    levels, merges = levels_of(leaves, scored.ct, scored.areas, args.leaf_cut, args.cuts)
    log.info(f"{len(set(leaves.values()))} units; levels: " + ", ".join(f"cut {c}: {len(set(m.values()))} realms" for c, m in sorted(levels.items(), reverse=True)))

    current = dict(leaves)  # the leaf merges lead from the units to the leaf realms, and the nested ones go on from there
    extra = []
    for r in merges.iter_rows(named=True):
        area: dict[int, float] = {}
        sea: dict[int, float] = {}
        for u, x in current.items():
            if x in (r["keep"], r["drop"]):
                area[x] = area.get(x, 0.0) + scored.areas[u]
                sea[x] = sea.get(x, 0.0) + (scored.areas[u] if u.endswith(SEA_SUFFIX) else 0.0)
        extra.append({"area_keep": area[r["keep"]], "area_drop": area[r["drop"]], "sea_keep": sea[r["keep"]] / area[r["keep"]], "sea_drop": sea[r["drop"]] / area[r["drop"]]})
        current = {u: (r["keep"] if x == r["drop"] else x) for u, x in current.items()}
    schema = {"area_keep": pl.Float64, "area_drop": pl.Float64, "sea_keep": pl.Float64, "sea_drop": pl.Float64}
    merges = pl.concat([merges, pl.DataFrame(extra, schema=schema)], how="horizontal")
    merges.write_csv(out / "merges.csv")
    with pl.Config(tbl_rows=60, float_precision=3, tbl_width_chars=220):
        log.info(f"merges above 0.4:\n{merges.filter(pl.col('height') >= 0.4)}")

    cells = pl.read_parquet(run_dir / UNITS / "cells.parquet").select("cell", "lat", "lon")
    position = {c: (la, lo) for c, la, lo in cells.iter_rows()}
    members: dict[str, list[str]] = {}
    for part, unit in units.members.items():
        members.setdefault(unit, []).append(part)
    centre = {}
    for unit in leaves:
        pts = [position[p.split(":")[0]] for p in members.get(unit, [unit]) if p.split(":")[0] in position]
        centre[unit] = (float(np.mean([p[0] for p in pts])), float(np.mean([p[1] for p in pts]))) if pts else (float("nan"), float("nan"))
    realms = realm_table(levels, scored.areas, centre)

    if args.n_boot:
        recovery = bootstrap_recovery(scored.ct, scored.areas, scored.adjacency, levels, args.leaf_cut, args.cuts, configs[config],
                                      n_boot=args.n_boot, seed=params.seed, trials=args.boot_trials, workers=args.workers)
        recovery.write_csv(out / "recovery.csv")
        summary = recovery.group_by("cut", "realm").agg(pl.col("jaccard").median().alias("jaccard_median"), (pl.col("jaccard") >= 0.5).mean().alias("recovered"))
        realms = realms.join(summary, on=["cut", "realm"], how="left")
    realms.write_csv(out / "realms.csv")
    with pl.Config(tbl_rows=80, float_precision=3, tbl_width_chars=220):
        log.info(f"realms:\n{realms}")

    cuts = sorted(levels, reverse=True)
    frames = part_rows({n: levels[c] for n, c in enumerate(cuts, start=1)}, units, read_light_grid(run_dir), None)
    pl.concat(frames).write_parquet(out / "parts.parquet")
    (out / "levels.json").write_text(json.dumps({"config": config, "leaf_cut": args.leaf_cut, "cuts": args.cuts,
                                                 "levels": {str(n): c for n, c in enumerate(cuts, start=1)}}, indent=2) + "\n")
    log.info(f"nested: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


if __name__ == "__main__":
    main()
