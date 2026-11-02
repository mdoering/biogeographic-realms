"""The review package of the headline map on half A: maps, and a table per realm.

`uv run python -I scripts/review_maps.py --counts data/gbif/<key>`

Reads the select stage of one configuration (`--config`, default `mt080`) of the headline run (`--run`, default `native_only`) and of the base
run (`--base`, default `off`) in `results/explore/res<r>/<run>/<config>/select/`, with `explore.read_selected` (which refuses a result that is not
current), at the resolution of the counts' `_COMPLETE` marker; the grid, the land and the TDWG
level 3 regions come from the references' manifest. Writes to `res<r>/review/`, with `run.log` and its timings:

- `k<k>.png` and `.gpkg` (`hexmaps.plot_parts`) for the selected k and the candidates either side, with sub-realm outlines, imputed parts
  hatched, a colour of its own for each realm of the selected map (`review.distinct_colours`) and the same colour for a realm across the
  three (`review.match_colours`, by largest shared area).
- `k<k>-no-floor.png` and `.gpkg` for each of `--no-floor-ks` (default 19 and 30): `candidates.merge_sequence` on Infomap's modules before the
  size floor (`review.unfloored_maps`), imputed as the select stage imputes, kept for M4's comparison and never used to choose; their parts
  and units are in `no_floor_parts.parquet` and `no_floor_units.parquet`, and each realm's area, parts and mean cell centre (longitude as a mean direction) in
  `no_floor_realms.csv` (a realm of a few island cells is hidden under its number on the map).
- `realm-<r>-k<k>.png`, one per realm of the selected map (`review.plot_realm`): its sub-realms, imputed parts and the units above the floor it merges.
- `boundary-<region>-k<k>.png`: the land boundary of the run against the base at the selected k in one TDWG level 1 region (default Northern
  America), edges kept, added and removed (`review.edge_changes`), with TDWG level 3 outlines; the edges are in `boundary-<region>-k<k>.csv`.
- `realms-k<k>.csv` (`review.realm_table`, endemic species at `--tolerance`), `sub_realms-k<k>.csv`, `summary.csv` (k, realms, stability and
  excess of the maps drawn) and `tables.md`, the same as pipe tables.

Half A only: it reads the explore runs' outputs, which are half A's, and loads no counts.
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import geopandas as gpd
import polars as pl
import shapely

from realms.boundaries import boundary, edge_adjacency, edge_lengths
from realms.candidates import apply_size_floor, merge_sequence
from realms.circularity import attribute, centre_regions
from realms.explore import UNITS, assigned_only, manifest_path, network_matrix, part_rows, peak_gb, read_light_grid, read_selected, read_units
from realms.grid import read_centres
from realms.hexmaps import part_geometries, plot_parts
from realms.natives import TDWG_LEVEL3_CODE
from realms.occurrences import IncompleteIngest, resolution_of
from realms.review import (
    distinct_colours,
    edge_changes,
    edge_lines,
    endemic_species,
    markdown,
    match_colours,
    part_areas,
    part_neighbours,
    plot_boundary_change,
    plot_realm,
    precise,
    realm_neighbours,
    realm_table,
    sub_realm_table,
    unfloored_maps,
)
from realms.units import TDWG_LEVEL1, TDWG_LEVEL1_CODE

log = logging.getLogger("review")
REGION_VIEW = {"Northern America": ((-180.0, 12.0, -10.0, 84.0), "+proj=laea +lat_0=50 +lon_0=-95 +datum=WGS84")}


def at(table: pl.DataFrame, k: int) -> pl.DataFrame:
    rows = table.filter(pl.col("k") == k)
    if not rows.height:
        raise ValueError(f"no map at k={k}; there are k={sorted(table['k'].unique())}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--counts", type=Path, required=True, help="the ingested download the runs were cleaned from, data/gbif/<key>: its marker gives the resolution")
    parser.add_argument("--run", default="native_only", help="the headline run, default native_only (DECISIONS, the circularity gate)")
    parser.add_argument("--base", default="off", help="the run whose boundary the headline's is set against, default off")
    parser.add_argument("--config", default="mt080", help="the Infomap configuration, default mt080 (DECISIONS, Infomap's scale)")
    parser.add_argument("--configs", type=Path, default=Path("explore/configs.toml"))
    parser.add_argument("--no-floor-ks", type=int, nargs="*", default=[19, 30], help="maps without the size floor, for M4 only")
    parser.add_argument("--region", default="Northern America", choices=sorted(REGION_VIEW), help="the TDWG level 1 region of the boundary figure")
    parser.add_argument("--tolerance", type=float, default=0.95, help="the share of its units inside a realm that makes a species endemic")
    parser.add_argument("--results", type=Path, default=Path("results/explore"))
    parser.add_argument("--reference", type=Path, default=Path("data/reference"))
    args = parser.parse_args()

    try:
        resolution = resolution_of(args.counts)
    except (IncompleteIngest, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error
    root = args.results / f"res{resolution}"
    out = root / "review"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(out / "run.log", mode="w")],
    )
    try:
        run(args, resolution, root, out)
    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(f"FAILED {error}") from error


def run(args: argparse.Namespace, resolution: int, root: Path, out: Path) -> None:
    started = time.perf_counter()
    timings: dict[str, float] = {}

    def lap(name: str, since: float) -> float:
        timings[name] = round(time.perf_counter() - since, 1)
        log.info(f"{name}: {timings[name]} s, peak {peak_gb():.1f} GB")
        return time.perf_counter()

    gpkg = manifest_path(args.reference, f"isea3h_{resolution}")
    land = manifest_path(args.reference, "natural_earth_land")
    level3 = manifest_path(args.reference, "wgsrpd_level3")
    log.info(f"---- review: resolution {resolution} (from {args.counts}), run {args.run} against {args.base}, configuration {args.config}")
    log.info(f"grid {gpkg}, land {land}, regions {level3}")
    run_dir = root / args.run
    head = read_selected(run_dir, args.config, args.configs)
    base = read_selected(root / args.base, args.config, args.configs)
    params = head.stamp["inputs"]["params"]
    k_sel = head.k
    ks = [k for k in (k_sel - 1, k_sel, k_sel + 1) if k in set(head.parts["k"].to_list())]
    log.info(f"{args.run}: selected k={k_sel} (fell short: {head.fell_short}), candidates k={sorted(head.parts['k'].unique())}; drawn k={ks}")
    log.info(f"{args.base}: selected k={base.k} (fell short: {base.fell_short})")
    units = read_units(run_dir)
    cells = pl.read_parquet(run_dir / UNITS / "cells.parquet")
    areas = part_areas(cells)
    t = lap("read the runs", time.perf_counter())

    wanted = set(head.parts["part"]) | set(base.parts["part"])
    geoms = part_geometries(gpkg, land, wanted)
    edges = edge_lengths(gpkg)
    adjacency = part_neighbours(wanted, edge_adjacency(edges))
    empty = sorted(geoms["part"][geoms.geometry.is_empty])
    log.info(f"{len(geoms)} part geometries, {len(edges)} edges; {len(empty)} parts have no area on the map: "
             + ", ".join(f"{p} ({areas[p]:.1f} km2 in cells.parquet)" for p in empty))
    t = lap("part geometries and adjacency", t)

    # ---- the selected map and its neighbours, one colour per realm across them
    labels = {k: dict(at(head.parts, k).select("part", "realm").iter_rows()) for k in ks}
    neighbours = {k: realm_neighbours(labels[k], adjacency) for k in ks}
    colours = {k_sel: distinct_colours(labels[k_sel], areas)}
    for k in ks:
        if k != k_sel:
            colours[k] = match_colours(labels[k_sel], colours[k_sel], labels[k], areas, neighbours[k])
    for k in ks:
        log.info(f"k={k} colours: " + ", ".join(f"{r} {c}" for r, c in sorted(colours[k].items())))
        rows = at(head.parts, k)
        subs = dict(at(head.candidates, k).filter(pl.col("sub_realm").is_not_null()).select("unit", "sub_realm").iter_rows())
        sub_realms = {p: subs[u] for p, u in rows.select("part", "unit").iter_rows() if u in subs}
        imputed = set(rows.filter(pl.col("imputed"))["part"])
        tag = "selected" if k == k_sel else "candidate"
        title = f"{args.run}, {args.config}, k = {k} ({tag}), resolution {resolution}: realms numbered, sub-realms outlined, imputed parts hatched"
        plot_parts(geoms[geoms["part"].isin(list(labels[k]))], labels[k], out / f"k{k}", title, sub_realms=sub_realms, imputed=imputed, adjacency=adjacency, colours=colours[k])
        t = lap(f"map k={k}", t)

    # ---- the matrix the select stage scored, the floor checked against the candidates, and the maps without the floor
    ct = network_matrix(run_dir, units, head.merged)
    partition = head.partition
    scored = assigned_only(partition, ct, units)
    del ct
    log.info(f"matrix: {len(scored.ct.cells)} units, {len(scored.ct.taxa)} taxa, {len(partition.unassigned)} unassigned; {len(set(partition.leaf_labels().values()))} modules")
    t = lap("matrix", t)

    resolved, no_floor = unfloored_maps(partition, scored.ct, scored.adjacency, args.no_floor_ks)
    t = lap(f"unfloored_maps k={args.no_floor_ks} from {len(set(resolved.values()))} modules", t)
    missing = sorted(set(args.no_floor_ks) - set(no_floor))
    if missing:
        log.info(f"no map without the floor at k={missing}: there are only {len(set(resolved.values()))} modules")
    floored, sub_of = apply_size_floor(resolved, scored.areas, scored.ct, params["floor_km2"])
    again = merge_sequence(floored, scored.ct, params["k_max"], params["k_min"])
    for k, found in again.items():
        direct = dict(at(head.candidates, k).select("unit", "realm").iter_rows())
        if found != direct:
            raise ValueError(f"the modules floored again do not give candidates.parquet at k={k}")
    direct_subs = dict(head.candidates.filter((pl.col("k") == k_sel) & pl.col("sub_realm").is_not_null()).select("unit", "sub_realm").iter_rows())
    if direct_subs != {u: m for u, m in resolved.items() if m in sub_of}:
        raise ValueError("the modules floored again do not give the sub-realms of candidates.parquet")
    log.info(f"check: the resolved modules, floored at {params['floor_km2']:,.0f} km2 and merged, give candidates.parquet at k={sorted(again)} and its {len(sub_of)} sub-realms")
    t = lap("floor check", t)

    no_floor_parts: pl.DataFrame | None = None
    if no_floor:
        frames = part_rows(no_floor, units, read_light_grid(run_dir), None)
        no_floor_parts = pl.concat(frames).drop("selected")
        no_floor_parts.write_parquet(out / "no_floor_parts.parquet")
        cells_sorted = sorted(resolved)
        pl.concat(
            [pl.DataFrame({"unit": cells_sorted, "k": k, "realm": [m[c] for c in cells_sorted], "module": [resolved[c] for c in cells_sorted]}) for k, m in sorted(no_floor.items())]
        ).write_parquet(out / "no_floor_units.parquet")
        centres = cells.select("cell", "lat", "lon")
        no_floor_realms = (
            no_floor_parts.with_columns(pl.col("part").replace_strict(areas, return_dtype=pl.Float64).alias("km2"), pl.col("part").str.replace(r":[^:]*$", "").alias("cell"))
            .join(centres, on="cell")
            .group_by("k", "realm")
            .agg(pl.col("km2").sum().alias("area_km2"), pl.len().alias("parts"), pl.col("part").str.ends_with(":land").sum().alias("land_parts"),
                 pl.col("imputed").sum().alias("imputed_parts"), pl.col("lat").mean().alias("mean_lat"),
                 # the mean direction, so that a realm across the antimeridian is not put at longitude 0
                 pl.arctan2(pl.col("lon").radians().sin().mean(), pl.col("lon").radians().cos().mean()).degrees().alias("mean_lon"))
            .sort(["k", "area_km2"], descending=[False, True])
        )
        no_floor_realms.write_csv(out / "no_floor_realms.csv")
        # colours carried on from the finest map drawn, so the large realms keep theirs; the rest only differ from their neighbours
        previous, previous_colours = labels[ks[-1]], colours[ks[-1]]
        for k in sorted(no_floor):
            rows = at(no_floor_parts, k)
            part_labels = dict(rows.select("part", "realm").iter_rows())
            colours[k] = match_colours(previous, previous_colours, part_labels, areas, realm_neighbours(part_labels, adjacency))
            previous, previous_colours = part_labels, colours[k]
            log.info(f"k={k} without the floor, colours: " + ", ".join(f"{r} {c}" for r, c in sorted(colours[k].items())))
            title = f"{args.run}, {args.config}, k = {k} WITHOUT THE SIZE FLOOR (for M4, never used to choose), resolution {resolution}: imputed parts hatched"
            plot_parts(geoms[geoms["part"].isin(list(part_labels))], part_labels, out / f"k{k}-no-floor", title,
                       imputed=set(rows.filter(pl.col("imputed"))["part"]), adjacency=adjacency, colours=colours[k])
            t = lap(f"map k={k} without the floor", t)

    # ---- the per-realm tables at the selected k
    unit_labels = dict(at(head.candidates, k_sel).select("unit", "realm").iter_rows())
    endemic = endemic_species(unit_labels, scored.ct, args.tolerance)
    t = lap("endemic species", t)
    realms = realm_table(at(head.parts, k_sel), at(head.candidates, k_sel), areas, scored.areas, endemic)
    realms.write_csv(out / f"realms-k{k_sel}.csv")
    subs = sub_realm_table(at(head.candidates, k_sel), scored.areas)
    subs.write_csv(out / f"sub_realms-k{k_sel}.csv")
    total = realms["area_km2"].sum()
    log.info(f"k={k_sel}: {realms.height} realms, {total:,.0f} km2 of parts ({realms['land_km2'].sum():,.0f} land, {realms['sea_km2'].sum():,.0f} sea); {subs.height} sub-realms")
    for r in realms.iter_rows(named=True):
        log.info(f"realm {r['realm']}: {r['area_km2'] / 1e6:.2f}M km2 ({r['area_share']:.1%}), {r['units']} units, {r['imputed_parts']} of {r['parts']} parts imputed, "
                 f"{r['endemic_species']:,} of {r['species']:,} species endemic, {r['sub_realms']} sub-realms ({r['sub_realm_km2'] / 1e6:.2f}M km2)")

    def sizes_of(part_rows: pl.DataFrame, unit_realm: dict[str, int]) -> dict:
        """The smallest realm by its parts (imputed included), and the largest realm's share of the map by parts and by units alone."""
        by_part = part_rows.with_columns(pl.col("part").replace_strict(areas, return_dtype=pl.Float64).alias("km2")).group_by("realm", maintain_order=True).agg(pl.col("km2").sum())["km2"]
        by_unit: dict[int, float] = {}
        for unit, realm in unit_realm.items():
            by_unit[realm] = by_unit.get(realm, 0.0) + scored.areas[unit]
        return {"smallest_realm_km2": by_part.min(), "largest_realm_share": by_part.max() / by_part.sum(),
                "largest_realm_unit_share": max(by_unit.values()) / sum(by_unit.values())}

    summary = []
    for k in ks:
        row = at(head.selection, k).row(0, named=True)
        summary.append({"k": k, "floor": "with the size floor", "realms": len(set(labels[k].values())), "stability": row["stability"], "excess": row["excess"],
                        "selected": k == k_sel, **sizes_of(at(head.parts, k), dict(at(head.candidates, k).select("unit", "realm").iter_rows()))})
    if no_floor_parts is not None:
        for k in sorted(no_floor):
            summary.append({"k": k, "floor": "without the size floor", "realms": len(set(no_floor[k].values())), "stability": None, "excess": None,
                            "selected": False, **sizes_of(at(no_floor_parts, k), no_floor[k])})
    summary_table = pl.DataFrame(summary, schema={"k": pl.Int64, "floor": pl.String, "realms": pl.Int64, "stability": pl.Float64, "excess": pl.Float64,
                                                  "selected": pl.Boolean, "smallest_realm_km2": pl.Float64, "largest_realm_share": pl.Float64,
                                                  "largest_realm_unit_share": pl.Float64})
    summary_table.write_csv(out / "summary.csv")

    # ---- each realm of the selected map, its sub-realms and imputed parts
    top = max(head.candidates["k"].to_list())
    above = dict(at(head.candidates, top).select("unit", "realm").iter_rows())
    for r in realms["realm"].to_list():
        row = realms.filter(pl.col("realm") == r).row(0, named=True)
        title = (f"Realm {r} of k = {k_sel} ({args.run}, {args.config}): {row['area_km2'] / 1e6:.1f}M km2 of parts, imputed included "
                 f"({row['unit_km2'] / 1e6:.1f}M km2 of units), {row['units']:,} units, {row['sub_realms']} sub-realms, {row['imputed_parts']:,} imputed parts")
        drawn = plot_realm(geoms, at(head.parts, k_sel), at(head.candidates, k_sel), r, colours[k_sel][r], out / f"realm-{r}-k{k_sel}", title, floored=above)
        log.info(f"realm {r} panel: {drawn}")
    t = lap("realm panels", t)

    # ---- the land boundary of the run against the base, in one region
    regions = gpd.read_file(level3).sort_values(TDWG_LEVEL3_CODE, kind="stable")
    region_of = centre_regions(read_centres(gpkg), regions)
    base_labels = dict(at(base.parts, k_sel).select("part", "realm").iter_rows())
    changes = attribute(edge_changes(boundary(base_labels, "land", edges), boundary(labels[k_sel], "land", edges)), region_of)
    here = changes.filter(pl.col("region") == args.region)
    slug = args.region.lower().replace(" ", "-")
    here.write_csv(out / f"boundary-{slug}-k{k_sel}.csv")
    for change, km, n in here.group_by("change").agg(pl.col("km").sum(), pl.len()).sort("change").iter_rows():
        log.info(f"{args.region}, k={k_sel}, {args.run} against {args.base}: {change} {km:,.0f} km ({n} edges)")
    lines = edge_lines(gpkg, list(zip(here["cell_a"], here["cell_b"], strict=True)))
    bounds, crs = REGION_VIEW[args.region]
    outside = sum(1 for line in lines if not (bounds[0] <= line.centroid.x <= bounds[2]))
    log.info(f"{outside} of {len(lines)} edges lie outside the view {bounds} (across the antimeridian)")
    title = (f"Land boundary in {args.region} (TDWG level 1), k = {k_sel}: {args.run} against {args.base} ({args.config}); "
             f"parts in {args.run}'s realm colours, TDWG level 3 outlines")
    inside = regions[[TDWG_LEVEL1.get(str(c), str(c)) == args.region for c in regions[TDWG_LEVEL1_CODE]]]
    area = shapely.make_valid(shapely.union_all(inside.geometry.values))
    focus = {p for p, g in zip(geoms["part"], geoms.geometry, strict=True) if p.endswith(":land") and g.intersects(area)}
    log.info(f"{len(focus)} land parts meet the TDWG level 3 regions of {args.region} ({len(inside)} of them)")
    plot_boundary_change(geoms, labels[k_sel], colours[k_sel], here, lines, regions, crs, bounds, out / f"boundary-{slug}-k{k_sel}", title, focus=focus)
    t = lap("boundary figure", t)

    # ---- the markdown
    shown = realms.with_columns(
        (pl.col("area_km2") / 1e6).alias("area Mkm2"), (pl.col("land_km2") / 1e6).alias("land Mkm2"), (pl.col("sea_km2") / 1e6).alias("sea Mkm2"),
        (pl.col("area_share") * 100).alias("share %"), (pl.col("unit_km2") / 1e6).alias("unit Mkm2"), (pl.col("unit_share") * 100).alias("unit share %"),
        (pl.col("imputed_km2") / 1e6).alias("imputed Mkm2"), (pl.col("sub_realm_km2") / 1e3).alias("sub-realm thousand km2"),
    ).select("realm", "area Mkm2", "land Mkm2", "sea Mkm2", "share %", "unit Mkm2", "unit share %", "imputed Mkm2", "units", "land_units", "sea_units",
             "parts", "imputed_parts", "species", "endemic_species", "sub_realms", "sub_realm_units", "sub-realm thousand km2")
    two = {c: precise(".2f") for c in ("area Mkm2", "land Mkm2", "sea Mkm2", "unit Mkm2", "imputed Mkm2")} | {
        "share %": ".1f", "unit share %": ".1f", "sub-realm thousand km2": precise(",.1f"),
    } | {c: "," for c in ("units", "land_units", "sea_units", "parts", "imputed_parts", "species", "endemic_species", "sub_realm_units")}
    unit_total, part_total, imputed_total = realms["unit_km2"].sum(), realms["area_km2"].sum(), realms["imputed_km2"].sum()
    lines_md = [
        f"## Realms of the selected map, {args.run} {args.config} k = {k_sel} (resolution {resolution})\n",
        f"Area is the sum of the realm's land and sea parts (imputed ones included). Endemic species: at least {args.tolerance:.0%} of the species' units "
        f"in the realm, species columns of the matrix the select stage scored. Sub-realms are modules the {params['floor_km2'] / 1e6:g}M km2 size floor merged into the realm.\n",
        markdown(shown, two),
        f"\nTwo bases of area: the parts, imputed ones included ({part_total / 1e6:.2f}M km2, `area` and `share`), and the units alone ({unit_total / 1e6:.2f}M km2, "
        f"`unit` and `unit share`; the largest realm {realms['unit_share'].max():.1%}), which is the basis of the select stage's scores and of the table in DECISIONS "
        f"(\"Infomap's scale on half A\"). "
        f"The difference is the imputed parts, {imputed_total / 1e6:.2f}M km2, as a part that is not imputed has its unit's area.\n",
        "\n### Sub-realms, the five largest of each realm\n",
    ]
    largest = subs.group_by("realm", maintain_order=True).head(5).select("realm", "sub_realm", "units", "area_km2")
    counts = subs.group_by("realm").agg(pl.len().alias("n"), pl.col("area_km2").median().alias("median"))
    lines_md.append(markdown(largest, {"area_km2": precise(",.1f")}))
    lines_md.append("\nAll sub-realms, by realm (count, median area in km2): "
                    + "; ".join(f"realm {r} {n} ({precise(',.1f')(m)})" for r, n, m in counts.sort("realm").iter_rows()) + f". Every one is in `sub_realms-k{k_sel}.csv`.\n")
    lines_md.append("\nOn the maps, the pale globe shows where no part is mapped: parts with no record left after cleaning, which are neither in a unit nor imputed.\n")
    lines_md.append("\n## The maps drawn\n")
    lines_md.append(markdown(summary_table.with_columns(
        (pl.col("smallest_realm_km2") / 1e3).alias("smallest realm, thousand km2 of parts"), (pl.col("largest_realm_share") * 100).alias("largest realm, % of parts"),
        (pl.col("largest_realm_unit_share") * 100).alias("largest realm, % of units"),
    ).select("k", "floor", "realms", "stability", "excess", "selected", "smallest realm, thousand km2 of parts", "largest realm, % of parts", "largest realm, % of units"),
        {"stability": ".3f", "excess": ".4f", "smallest realm, thousand km2 of parts": precise(",.1f"), "largest realm, % of parts": ".1f", "largest realm, % of units": ".1f"}))
    lines_md.append(f"\nThe maps without the size floor are `candidates.merge_sequence` on Infomap's {len(set(resolved.values()))} modules (separate pieces resolved), "
                    "with no floor. They were not bootstrapped or scored, are kept for M4's comparison and are never used to choose.\n")
    by_change = here.group_by("change").agg(pl.col("km").sum(), pl.len().alias("edges")).sort("change")
    lines_md.append(f"\n## Land boundary in {args.region}, {args.run} against {args.base}, k = {k_sel}\n")
    lines_md.append(markdown(by_change, {"km": ",.0f"}))
    (out / "tables.md").write_text("\n".join(lines_md))
    log.info(f"wrote {out / 'tables.md'}; timings {timings}; total {time.perf_counter() - started:.0f} s, peak {peak_gb():.1f} GB")


if __name__ == "__main__":
    main()
