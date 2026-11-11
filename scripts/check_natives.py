"""Join the reference sources to COL ids and report how well each one matched, then load everything as `References`.

WCVP, each GRIIS list and WRiMS are matched to COL species by name. A source with more than 10% of its names unmatched (ambiguous names
are counted apart, as they are never guessed) is sent to GBIF's match API v2 for those names, at most 5 requests a second, and what GBIF
places on a COL species is saved to `data/reference/gbif_fallback.parquet` for `load_references` to use. This only reads from GBIF.

Writes results/natives.json, with what each GRIIS list says its entries are (`griis_means`) and the references' counts. Usage: uv run python -I scripts/check_natives.py [--resolution 6] [--no-references]
"""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import polars as pl

from realms.grid import load_grid
from realms.lists import NO_WILD_RANGE, resolve_no_wild_range
from realms.natives import (
    FALLBACK_FILE,
    UNMATCHED,
    gbif_fallback,
    griis_lists,
    griis_names,
    load_references,
    match_names,
    match_stats,
    name_index,
    read_fallback,
    wcvp_distributions,
    wcvp_names,
    write_fallback,
    wrims_distributions,
    wrims_names,
)

ROOT = Path("data/reference")
OUT = Path("results/natives.json")
LAND = ROOT / "natural_earth" / "ne_10m_land.zip"
FAIL_SHARE = 0.10


def peak_gb() -> float:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 1e9 if sys.platform == "darwin" else peak / 1e6


def line(label: str, stats: dict, extra: str = "") -> str:
    return f"  {label:28} {stats['names']:>7} names  matched {stats['matched_share']:6.1%}  ambiguous {stats['ambiguous']:>5}  unmatched {stats[UNMATCHED]:>5} ({stats[UNMATCHED] / max(stats['names'], 1):.1%}) {extra}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=6)
    parser.add_argument("--no-references", action="store_true", help="skip loading References, which needs about 15 GB for load_taxa")
    args = parser.parse_args()

    manifest = json.loads((ROOT / "manifest.json").read_text())
    started = time.perf_counter()
    index = name_index(Path(manifest["col"]["path"]))
    print(f"COL {manifest['col']['release']}: {index.height:,} names, {index.filter(pl.col('status') != 'synonym')['species_id'].n_unique():,} species ({time.perf_counter() - started:.0f} s)")

    sources: dict[str, pl.DataFrame] = {"WCVP": wcvp_names(Path(manifest["wcvp"]["path"]))}
    if manifest["wrims"].get("path"):
        sources["WRiMS"] = wrims_names(Path(manifest["wrims"]["path"]))
    griis_paths = {k.removeprefix("griis/"): Path(v["path"]) for k, v in manifest.items() if k.startswith("griis/") and v.get("path")}
    listed = griis_names(list(griis_paths.values()))
    print(f"{len(griis_paths)} GRIIS lists, {listed.height:,} species entries; {sorted(k for k, v in manifest.items() if v.get('status') == 'missing')} missing")
    everything = griis_names(list(griis_paths.values()), species_only=False)
    other_ranks = everything.filter(pl.col("rank") != "species")
    print("  accepted entries of other ranks dropped:", dict(sorted(other_ranks.group_by("rank").len().iter_rows(), key=lambda r: -r[1])[:8]), f"({other_ranks.height:,} in all)")

    # what each list says its entries are: introduced (an empty value counts), native or uncertain, or nothing for a list with no distribution table
    held = griis_lists(listed)
    kinds = ("entries", "introduced", "unstated", "native", "uncertain")
    with_table, without_table = held.filter(pl.col("table")), held.filter(~pl.col("table"))
    totals = {k: int(with_table[k].sum()) for k in kinds}
    print(
        f"  GRIIS entries by their own establishment means: {len(with_table)} lists with a distribution table hold {totals['entries']:,} species entries, "
        f"{totals['introduced']:,} introduced, {totals['unstated']:,} with no value, {totals['native']:,} native and {totals['uncertain']:,} uncertain, "
        f"the last three left out; {len(without_table)} lists with none hold {int(without_table['entries'].sum()):,} entries that say nothing"
    )
    out = pl.col("native") + pl.col("uncertain") + pl.col("unstated")
    left_out = with_table.filter(out > 0).sort(out, descending=True)
    print(
        "  lists with entries left out:",
        ", ".join(f"{r['alias']} {r['native']}n/{r['uncertain']}u/{r['unstated']}e of {r['entries']}" for r in left_out.head(12).iter_rows(named=True)),
    )

    report: dict = {
        "release": manifest["col"]["release"],
        "sources": {},
        "griis_lists": {},
        "griis_means": {"lists_with_table": with_table.height, "lists_without_table": without_table.height, "entries_without_table": int(without_table["entries"].sum()), **totals},
        "griis_means_by_list": with_table.to_dicts(),
    }
    matched = {name: match_names(frame, index) for name, frame in sources.items()}
    matched["GRIIS"] = match_names(listed, index)
    print("\nMatch rates (distinct names, then the names of one list each)")
    for name, frame in matched.items():
        report["sources"][name] = stats = match_stats(frame)
        print(line(name, stats))
    per_list = {alias: match_stats(rows) for (alias,), rows in matched["GRIIS"].group_by(["alias"])}
    report["griis_lists"] = per_list
    weak = {alias: s for alias, s in per_list.items() if s["names"] and s[UNMATCHED] / s["names"] > FAIL_SHARE}
    print(f"  GRIIS lists with more than {FAIL_SHARE:.0%} unmatched: {len(weak)} of {len(per_list)}")
    for alias, s in sorted(weak.items()):
        print(line(alias, s))

    # a source whose unmatched names exceed the share goes to GBIF, those names alone
    to_ask = [pl.DataFrame(schema={"name": pl.String, "kingdom": pl.String})]
    for name in ("WCVP", "WRiMS"):
        if name in matched and report["sources"][name][UNMATCHED] / max(report["sources"][name]["names"], 1) > FAIL_SHARE:
            to_ask.append(matched[name].filter(pl.col("status") == UNMATCHED).select("name", "kingdom"))
    for alias in weak:
        to_ask.append(matched["GRIIS"].filter(pl.col("alias") == alias, pl.col("status") == UNMATCHED).select("name", "kingdom"))
    names = pl.concat(to_ask).unique()
    if names.height:
        print(f"\nAsking GBIF's match API about {names.height:,} names of the sources above ({names.height / 5 / 60:.1f} min at 5 a second)")
        found = gbif_fallback(names, index)
        write_fallback(found, ROOT / FALLBACK_FILE, manifest["col"]["release"])
        report["gbif_fallback"] = {"asked": names.height, "placed": found.height}
        print(f"  GBIF placed {found.height:,} of them on a COL species, saved to {ROOT / FALLBACK_FILE}")
        for name in [*sources, "GRIIS"]:
            frame = matched[name]
            after = match_stats(match_names(frame, index, found))
            report["sources"][name]["after_fallback"] = after
            print(line(f"{name} after GBIF", after))
    else:
        print("\nNo source exceeds the share, so GBIF is not asked and no fallback file is written")
        (ROOT / FALLBACK_FILE).unlink(missing_ok=True)

    # what the joins hold, counting rows as well as names: a name that fails may carry many rows
    fallback = read_fallback(ROOT / FALLBACK_FILE, manifest["col"]["release"], index) if (ROOT / FALLBACK_FILE).exists() else None
    wcvp = wcvp_distributions(Path(manifest["wcvp"]["path"]))
    by_row = {"WCVP": wcvp.join(match_names(sources["WCVP"], index, fallback).select("wcvp_id", "status"), on="wcvp_id", how="left")}
    print(f"\nWCVP: {wcvp.height:,} native or introduced level 3 rows of {wcvp['wcvp_id'].n_unique():,} taxa; native {wcvp['native'].sum():,}")
    if "WRiMS" in sources:
        dist = wrims_distributions(Path(manifest["wrims"]["path"]))
        by_row["WRiMS"] = dist.join(match_names(sources["WRiMS"], index, fallback).select("wrims_id", "status"), on="wrims_id", how="left")
        print(f"WRiMS: {dist.height:,} introduced area rows of {dist['wrims_id'].n_unique():,} taxa in {dist['area_id'].n_unique()} areas")
    by_row["GRIIS"] = match_names(listed, index, fallback)
    for name, rows in by_row.items():
        counts = dict(rows.group_by("status").len().iter_rows())
        total = rows.height
        report["sources"][name]["rows"] = {"total": total, **{k: counts.get(k, 0) for k in ("matched", "ambiguous", "unmatched")}, "no_species_row": counts.get(None, 0)}
        print(f"  {name:6} rows: {total:,}  matched {counts.get('matched', 0) / total:.1%}  ambiguous {counts.get('ambiguous', 0) / total:.1%}  unmatched {counts.get('unmatched', 0) / total:.1%}  of a taxon that is not a species {counts.get(None, 0) / total:.1%}")
    ids, left_out = resolve_no_wild_range(index)
    report["no_wild_range"] = {"entries": len(NO_WILD_RANGE), "species_ids": len(ids), "left_out": left_out}
    print(f"NO_WILD_RANGE: {len(NO_WILD_RANGE)} entries, {len(ids)} COL ids, left out at run time: {left_out or 'none'}")

    if not args.no_references:
        grid = load_grid(ROOT / "grids" / f"ISEA3H-{args.resolution}.gpkg", LAND)
        refs = load_references(grid, ROOT)
        report["references"] = {
            "resolution": args.resolution,
            "taxa": refs.taxa.height,
            "groups": refs.groups.height,
            "wcvp": {"rows": refs.wcvp.height, "species": refs.wcvp["species_id"].n_unique(), "native": int(refs.wcvp["native"].sum())},
            "griis": {"rows": refs.griis.height, "species": refs.griis["species_id"].n_unique(), "lists": refs.griis.select("country", "subdivision").unique().height},
            "griis_name_only": {
                "rows": refs.griis_name_only.height,
                "species": refs.griis_name_only["species_id"].n_unique(),
                "lists": refs.griis_name_only.select("country", "subdivision").unique().height,
            },
            "wrims": {"rows": refs.wrims.height, "species": refs.wrims["species_id"].n_unique(), "with_country": refs.wrims["country"].drop_nulls().len()},
            "land_parts": refs.tdwg3_parts["part"].n_unique(),
            "land_parts_without_region": refs.tdwg3_parts.filter(pl.col("tdwg3").is_null()).height,
            "admin1": {"cells": refs.admin1_cells.height, "subdivisions": refs.admin1_cells.select("country", "subdivision").unique().height, "by": dict(refs.admin1_cells.group_by("via").len().iter_rows())},
            "no_wild_range": len(refs.no_wild_range),
            "missing": list(refs.missing),
        }
        unplaced = refs.griis_no_cells.select("country", "subdivision").unique().sort("country", "subdivision")
        report["references"]["admin1"]["subdivisions_without_cells"] = [f"{c}-{s}" for c, s in unplaced.iter_rows()]
        report["references"]["admin1"]["pairs_without_cells"] = refs.griis_no_cells.height
        name_only_unplaced = refs.griis_name_only_no_cells.select("country", "subdivision").unique().sort("country", "subdivision")
        report["references"]["admin1"]["name_only_subdivisions_without_cells"] = [f"{c}-{s}" for c, s in name_only_unplaced.iter_rows()]
        report["references"]["admin1"]["name_only_pairs_without_cells"] = refs.griis_name_only_no_cells.height
        report["references"]["wrims"]["pairs_without_country"] = refs.wrims_no_country.height
        print("\nReferences:")
        for key, value in report["references"].items():
            print(f"  {key}: {value}")
        print(f"peak RSS of this process: {peak_gb():.1f} GB, {time.perf_counter() - started:.0f} s in all")
        report["peak_rss_gb"] = round(peak_gb(), 1)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, default=str))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
