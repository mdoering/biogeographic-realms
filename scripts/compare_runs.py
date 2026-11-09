"""A run's merge order cut at another run's realm counts, and how far the maps agree: as C2 compares half B with half A (PREREGISTRATION.md, section 9).

`uv run python -I scripts/compare_runs.py --counts data/gbif/<key> --runs native_only-eukaryotes-rank-family ... [--against native_only-eukaryotes]`

Each run's levels are rebuilt as `scripts/nested_cuts.py` builds them, from its `nested-cuts/cuts.json` (configuration and base cut): the finest modules of
its select result, their pieces resolved, merged freely to the base cut (`nested.leaf_realms`) and then land with land and sea with sea down to one realm
(`nested.nested_merges`), which gives a map at every realm count. Each run's map is compared with the reference run's (`--against`) at the same count, over
the units both map:

- at the reference's reported levels, and at each of the run's own reported levels;
- by the V-measure over units (`score.v_measure`, as the five-seed agreement) and by area agreement (`crosscheck.area_agreement`: the share of the area in
  realms matched one to one).

Writes `res<r>/<against>/compare-runs.csv` (or `--out`) and logs it. Half A only; it chooses nothing.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import polars as pl

from realms.candidates import resolve_pieces
from realms.crosscheck import area_agreement
from realms.explore import assigned_only, network_matrix, peak_gb, read_selected, read_units
from realms.nested import leaf_realms, nested_merges, sea_realms
from realms.occurrences import IncompleteIngest, resolution_of
from realms.score import v_measure

log = logging.getLogger("compare_runs")


def merge_order(run_dir: Path, configs: Path) -> tuple[dict[int, dict[str, int]], dict[str, float], dict]:
    """The run's map at every realm count of its nested merges, the areas of its scored units, and its `cuts.json`."""
    cuts = json.loads((run_dir / "nested-cuts/cuts.json").read_text())
    found = read_selected(run_dir, cuts["config"], configs)
    units = read_units(run_dir)
    scored = assigned_only(found.partition, network_matrix(run_dir, units, found.merged), units)
    leaves = resolve_pieces(found.partition.leaf_labels(), scored.adjacency, scored.ct)
    base, _ = leaf_realms(leaves, scored.ct, cuts["base_cut"])
    _, maps = nested_merges(base, scored.ct, sea_realms(base, scored.areas))
    return maps, scored.areas, cuts


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--against", default="native_only-eukaryotes")
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    root = args.results / f"res{resolution}"
    out = args.out or root / args.against / "compare-runs.csv"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S", stream=sys.stdout)
    started = time.perf_counter()
    reference, ref_areas, ref_cuts = merge_order(root / args.against, args.configs)
    log.info(f"---- compare runs: resolution {resolution}, against {args.against} ({ref_cuts['config']}, reported {ref_cuts['reported']}), runs {args.runs}")
    rows = []
    for run in args.runs:
        maps, areas, cuts = merge_order(root / run, args.configs)
        own = {lv["realms"] for lv in cuts["levels"]}
        counts = sorted(set(ref_cuts["reported"]) | set(cuts["reported"]))
        for k in counts:
            if k not in maps or k not in reference:
                rows.append({"run": run, "realms": k, "why": "reference" if k in ref_cuts["reported"] else "own", "v_measure": None, "area_agreement": None,
                             "units": 0, "run_level": k in own, "run_config": cuts["config"]})
                continue
            a, b = reference[k], maps[k]
            shared = a.keys() & b.keys()
            rows.append({"run": run, "realms": k, "why": "both" if k in ref_cuts["reported"] and k in cuts["reported"] else ("reference" if k in ref_cuts["reported"] else "own"),
                         "v_measure": v_measure({u: a[u] for u in shared}, {u: b[u] for u in shared}),
                         "area_agreement": area_agreement(a, b, {u: ref_areas.get(u, areas.get(u)) for u in shared}),
                         "units": len(shared), "run_level": k in own, "run_config": cuts["config"]})
        log.info(f"{run}: {cuts['config']}, levels {[lv['realms'] for lv in cuts['levels']]}, reported {cuts['reported']}; {time.perf_counter() - started:.0f} s")
    table = pl.DataFrame(rows)
    table.write_csv(out)
    with pl.Config(tbl_rows=100, float_precision=3, tbl_width_chars=200):
        log.info(f"against {args.against}:\n{table}")
    log.info(f"compare runs: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB, wrote {out}")


if __name__ == "__main__":
    main()
