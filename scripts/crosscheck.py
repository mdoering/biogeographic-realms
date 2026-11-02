"""The cross-checks of the headline map on half A: a beta-sim consensus of its units, and the run that merges duplicate genus and family names.

`uv run python -I scripts/crosscheck.py --counts data/gbif/<key> --parts betasim names`

Reads the select stage of one configuration (`--config`, default `mt080`) of the headline run (`--run`, default `native_only`) in
`results/explore/res<r>/`, at the resolution of the counts' `_COMPLETE` marker. Writes to `res<r>/crosscheck/`, and appends to its `run.log`:

- **betasim**: the matrix the select stage scored (`explore.network_matrix`, restricted by `explore.assigned_only` to the units Infomap gave a
  module), `betasim.rank_beta_sim` over its units and `betasim.consensus` with Ward and with UPGMA (`average`) at `--ks`, for each of `--seeds`,
  always with the unit ids. Infomap's map of a unit is that of its parts that are not imputed (`crosscheck.unit_labels`), checked against
  `candidates.parquet`.
  - `betasim_v.csv`: per method, seed, shuffles and k, the consensus's `spread` and its V-measure against Infomap's map at that k, if it has one.
  - `betasim_seeds.csv`: per method, shuffles and k, the V-measure between the consensus of each pair of seeds.
  - `betasim_area.csv`: per method, seed and shuffles, the share of area on which the consensus and Infomap's selected map agree at its k,
    realms matched one to one (`crosscheck.area_agreement`).
  - `betasim_labels.parquet`: the consensus group of every unit, per method, seed, shuffles and k.
  A pair of seeds below `--seed-floor` at the selected k runs that method again at that k only, with `--more-shuffles` shuffles and the same seeds.
- **names**: the same configuration of the run that merges duplicate names (`--merged`, default `<run>-merged-names`) against the run on ids.
  - `names_v.csv`: the V-measure by part (`circularity.compare`) at each run's selected k and at the headline's selected k, with the parts both
    hold and those only one holds.
  - `names_selection.csv`: stability and excess per k of both runs, the differences (merged minus ids), and which k each selected.
  - `names_duplicates.csv`: the duplicate-name groups (`matrix.duplicate_names`) among the merged run's cleaned taxa, per rank.

Half A only: it reads the explore runs' outputs, which are half A's, and loads no counts.
"""

import argparse
import json
import logging
import sys
import time
from typing import Literal
from dataclasses import asdict
from pathlib import Path

import polars as pl

from realms.betasim import consensus, rank_beta_sim
from realms.circularity import check_select_stamp, compare, labels_at
from realms.cluster import Partition
from realms.crosscheck import area_agreement, pairwise_v, unit_labels
from realms.explore import CLEANED, MERGED_SUFFIX, UNITS, assigned_only, load_configs, network_matrix, peak_gb, read_stamp, read_units
from realms.matrix import duplicate_names
from realms.occurrences import IncompleteIngest, resolution_of

METHODS: tuple[Literal["ward", "average"], ...] = ("ward", "average")  # average is UPGMA
log = logging.getLogger("crosscheck")


def read_select(directory: Path, options: dict) -> dict:
    """The select stage of a configuration: its stamp (checked to be select's, with these options), selected k, parts and selection table."""
    stamp = read_stamp(directory, "stage-select")
    check_select_stamp(stamp, options, str(directory))
    chosen = json.loads((directory / "selected.json").read_text())
    if chosen.get("degenerate"):
        raise ValueError(f"{directory} is degenerate: {chosen}")
    parts = pl.read_parquet(directory / "parts.parquet")
    flagged = parts.filter(pl.col("selected"))["k"].unique().to_list()
    if flagged != [chosen["k"]]:
        raise ValueError(f"{directory}: selected.json says k={chosen['k']}, parts.parquet marks {flagged}")
    units_hash = read_stamp(directory.parent / UNITS).get("hash")
    if stamp["inputs"]["units"] != units_hash:
        raise ValueError(f"{directory}: the select stamp is of units {stamp['inputs']['units']}, the units stage is {units_hash}")
    return {
        "stamp": stamp, "k": chosen["k"], "parts": parts, "selection": pl.read_csv(directory / "selection.csv"),
        "merged": bool(stamp["inputs"].get("merge_duplicate_names")),
    }


def betasim(run_dir: Path, config: str, headline: dict, ks: list[int], seeds: list[int], shuffles: int, more: int, floor: float, out: Path) -> None:
    started = time.perf_counter()
    k_sel = headline["k"]
    units = read_units(run_dir)
    ct = network_matrix(run_dir, units, headline["merged"])
    table = pl.read_parquet(run_dir / config / "partition.parquet")
    partition = Partition(
        paths={u: (top, leaf) for u, leaf, top in table.filter(pl.col("leaf").is_not_null()).iter_rows()},
        unassigned=sorted(table.filter(pl.col("leaf").is_null())["unit"].to_list()),
    )
    scored = assigned_only(partition, ct, units)  # the select stage's own restriction to the units Infomap gave a module
    ids = list(scored.ct.cells)
    log.info(f"matrix: {len(ct.cells)} units, {len(ct.taxa)} taxa; {len(partition.unassigned)} unassigned; scored {len(ids)} units, {len(scored.ct.taxa)} taxa")

    # Infomap's map of each unit, from the parts that are not imputed, checked against the candidates
    candidates = pl.read_parquet(run_dir / config / "candidates.parquet")
    infomap_ks = sorted(headline["parts"]["k"].unique().to_list())
    infomap = {k: unit_labels(headline["parts"], k) for k in infomap_ks}
    for k, labels in infomap.items():
        direct = dict(candidates.filter(pl.col("k") == k).select("unit", "realm").iter_rows())
        if labels != direct:
            raise ValueError(f"at k={k} the unit labels from parts.parquet differ from candidates.parquet")
        if set(labels) != set(ids):
            raise ValueError(f"at k={k} Infomap labels {len(labels)} units, the scored matrix has {len(ids)}")
    log.info(f"Infomap's unit labels at k={infomap_ks} from parts.parquet (not imputed), equal to candidates.parquet; selected k={k_sel}")
    areas = scored.areas

    t = time.perf_counter()
    dist = rank_beta_sim(scored.ct, ids)
    del ct, scored, units
    log.info(f"rank_beta_sim: {dist.shape}, {time.perf_counter() - t:.0f} s, peak {peak_gb():.1f} GB")

    v_rows, seed_rows, area_rows, label_frames = [], [], [], []

    def run(method: str, seed: int, n: int, at: list[int]) -> dict[int, dict[str, int]]:
        t = time.perf_counter()
        found = consensus(dist, method, at, shuffles=n, seed=seed, ids=ids)
        log.info(f"consensus {method} seed {seed}, {n} shuffles, k={at[0]}-{at[-1]}: {time.perf_counter() - t:.0f} s, peak {peak_gb():.1f} GB")
        maps = {}
        for k, c in found.items():
            maps[k] = dict(zip(ids, c.labels.tolist(), strict=True))
            v = compare(maps[k], infomap[k])["v_measure"] if k in infomap else None
            v_rows.append({"method": method, "seed": seed, "shuffles": n, "k": k, "groups": len(set(c.labels.tolist())), "spread": c.spread, "v_infomap": v})
            label_frames.append(pl.DataFrame({"unit": ids, "method": method, "seed": seed, "shuffles": n, "k": k, "group": c.labels}))
        area_rows.append({"method": method, "seed": seed, "shuffles": n, "k": k_sel, "area_agreement": area_agreement(maps[k_sel], infomap[k_sel], areas),
                          "v_infomap": compare(maps[k_sel], infomap[k_sel])["v_measure"]})
        return maps

    for method in METHODS:
        by_seed = {seed: run(method, seed, shuffles, ks) for seed in seeds}
        for k in ks:
            for row in pairwise_v({seed: maps[k] for seed, maps in by_seed.items()}):
                seed_rows.append({"method": method, "shuffles": shuffles, "k": k, "seed_a": row["a"], "seed_b": row["b"], "v_measure": row["v_measure"]})
        at_sel = [r for r in seed_rows if r["method"] == method and r["shuffles"] == shuffles and r["k"] == k_sel]
        for r in v_rows:
            if r["method"] == method and r["k"] == k_sel:
                log.info(f"{method} seed {r['seed']} k={k_sel}: V against Infomap {r['v_infomap']:.3f}, spread {r['spread']:.3f}")
        log.info(f"{method} seeds at k={k_sel}: " + ", ".join(f"{r['seed_a']}/{r['seed_b']} {r['v_measure']:.3f}" for r in at_sel))
        if min(r["v_measure"] for r in at_sel) < floor:
            log.info(f"{method}: a pair of seeds is below {floor} at k={k_sel}; again at k={k_sel} only with {more} shuffles")
            again = {seed: run(method, seed, more, [k_sel]) for seed in seeds}
            more_rows = [
                {"method": method, "shuffles": more, "k": k_sel, "seed_a": row["a"], "seed_b": row["b"], "v_measure": row["v_measure"]}
                for row in pairwise_v({seed: maps[k_sel] for seed, maps in again.items()})
            ]
            seed_rows += more_rows
            log.info(f"{method} seeds at k={k_sel}, {more} shuffles: " + ", ".join(f"{r['seed_a']}/{r['seed_b']} {r['v_measure']:.3f}" for r in more_rows))

    pl.DataFrame(v_rows).write_csv(out / "betasim_v.csv")
    pl.DataFrame(seed_rows).write_csv(out / "betasim_seeds.csv")
    pl.DataFrame(area_rows).write_csv(out / "betasim_area.csv")
    pl.concat(label_frames).write_parquet(out / "betasim_labels.parquet")
    for r in area_rows:
        log.info(f"area agreement {r['method']} seed {r['seed']} ({r['shuffles']} shuffles) at k={k_sel}: {r['area_agreement']:.3f}")
    log.info(f"betasim: {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


def names(root: Path, run: str, merged_run: str, config: str, options: dict, headline: dict, out: Path) -> None:
    merged = read_select(root / merged_run / config, options)
    if not merged["merged"] or headline["merged"]:
        raise ValueError(f"{merged_run} must merge duplicate names and {run} must not: their stamps say {merged['merged']} and {headline['merged']}")
    k_ids, k_merged, k_sel = headline["k"], merged["k"], headline["k"]
    rows = []
    for at, ka, kb in (("selected", k_ids, k_merged), (f"k={k_sel}", k_sel, k_sel)):
        found = compare(labels_at(headline["parts"], ka), labels_at(merged["parts"], kb))
        rows.append({"at": at, "k_ids": ka, "k_merged": kb, **found})
        log.info(f"names {at}: V {found['v_measure']:.3f} at k={ka}/{kb}, {found['shared']} parts shared, {found['only_a']} only on ids, {found['only_b']} only merged")
    pl.DataFrame(rows).rename({"shared": "shared_parts", "only_a": "only_ids", "only_b": "only_merged"}).write_csv(out / "names_v.csv")

    table = headline["selection"].rename({"stability": "stability_ids", "excess": "excess_ids", "selected": "selected_ids"}).join(
        merged["selection"].rename({"stability": "stability_merged", "excess": "excess_merged", "selected": "selected_merged"}), on="k", how="full", coalesce=True
    ).sort("k").with_columns(
        (pl.col("stability_merged") - pl.col("stability_ids")).alias("stability_diff"), (pl.col("excess_merged") - pl.col("excess_ids")).alias("excess_diff")
    ).select("k", "stability_ids", "stability_merged", "stability_diff", "excess_ids", "excess_merged", "excess_diff", "selected_ids", "selected_merged")
    table.write_csv(out / "names_selection.csv")
    for r in table.iter_rows(named=True):
        log.info(f"names k={r['k']}: stability {r['stability_ids']} / {r['stability_merged']}, excess {r['excess_ids']} / {r['excess_merged']} (diff {r['excess_diff']})")
    log.info(f"names: selected k={k_ids} on ids, k={k_merged} merged")

    taxa = pl.read_parquet(root / merged_run / CLEANED / "taxa.parquet")
    groups = duplicate_names(taxa)
    summary = groups.group_by("rank").agg(pl.len().alias("groups"), pl.col("ids").list.len().sum().alias("ids"), (pl.col("ids").list.len() - 1).sum().alias("ids_merged_away")).sort("rank")
    totals = pl.DataFrame({"rank": ["genus", "family"], "distinct_ids": [taxa["genus_id"].n_unique(), taxa["family_id"].n_unique()]})
    summary = totals.join(summary, on="rank", how="left").fill_null(0)
    summary.write_csv(out / "names_duplicates.csv")
    for r in summary.iter_rows(named=True):
        log.info(f"duplicate {r['rank']} names: {r['groups']} groups of {r['ids']} ids, {r['ids_merged_away']} merged away of {r['distinct_ids']} ids")
    same = pl.read_parquet(root / run / CLEANED / "presence.parquet").equals(pl.read_parquet(root / merged_run / CLEANED / "presence.parquet"))
    same_units = all(
        pl.read_parquet(root / run / UNITS / f).equals(pl.read_parquet(root / merged_run / UNITS / f)) for f in ("members.parquet", "areas.parquet", "excluded.parquet")
    )
    log.info(f"the cleaned presence of both runs is {'the same' if same else 'DIFFERENT'}; their units are {'the same' if same_units else 'DIFFERENT'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--parts", nargs="+", choices=("betasim", "names"), default=["betasim", "names"])
    parser.add_argument("--run", default="native_only", help="the headline run, default native_only (DECISIONS, the circularity gate)")
    parser.add_argument("--merged", default=None, help="the run that merges duplicate names, default <run>-merged-names")
    parser.add_argument("--config", default="mt080", help="the Infomap configuration, default mt080 (DECISIONS, Infomap's scale)")
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--ks", type=int, nargs=2, default=(5, 12), metavar=("FIRST", "LAST"), help="the consensus's k, default 5 to 12")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--shuffles", type=int, default=20)
    parser.add_argument("--seed-floor", type=float, default=0.9, help="a pair of seeds below this at the selected k runs the method again with more shuffles")
    parser.add_argument("--more-shuffles", type=int, default=100)
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    root = args.results / f"res{resolution}"
    out = root / "crosscheck"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="a")],
    )
    log.info(f"---- crosscheck {args.parts}: resolution {resolution} (from {args.counts}), run {args.run}, configuration {args.config}")
    configs, _ = load_configs(args.configs)
    if args.config not in configs:
        raise SystemExit(f"FAILED {args.configs} has no configuration {args.config}")
    options = asdict(configs[args.config])
    try:
        headline = read_select(root / args.run / args.config, options)
        if "betasim" in args.parts:
            ks = list(range(args.ks[0], args.ks[1] + 1))
            log.info(f"betasim: methods {METHODS}, k={ks[0]}-{ks[-1]}, seeds {args.seeds}, {args.shuffles} shuffles, ids passed")
            betasim(root / args.run, args.config, headline, ks, args.seeds, args.shuffles, args.more_shuffles, args.seed_floor, out)
        if "names" in args.parts:
            names(root, args.run, args.merged or f"{args.run}-{MERGED_SUFFIX}", args.config, options, headline, out)
    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error


if __name__ == "__main__":
    main()
