"""Tables and figures of the review package (`scripts/review_maps.py`): realm colours kept across k, the per-realm table, endemic species,
the maps without the size floor, a realm's sub-realms and imputed parts, and a boundary of one run against another.

A map here is `parts.parquet` at one k: every part (`<cell>:land`, `<cell>:sea`) with its realm, its unit (none for an imputed part) and
whether it was imputed. The units and their sub-realms are `candidates.parquet` at that k.
"""

from collections import defaultdict
from collections.abc import Collection, Mapping
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import shapely  # noqa: E402
from matplotlib.colors import to_hex, to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.patheffects import withStroke  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from shapely import MultiLineString  # noqa: E402
from shapely.ops import polylabel  # noqa: E402

from realms.boundaries import _cells  # noqa: E402
from realms.candidates import merge_sequence, resolve_pieces  # noqa: E402
from realms.cluster import Partition  # noqa: E402
from realms.endemism import Labels  # noqa: E402
from realms.grid import SNAP, _shifted  # noqa: E402
from realms.hexmaps import DPI, EQUAL_EARTH, GLOBE, OUTLINE, WGS84, _densified, _dissolve, _polygons, _projected, _repaired  # noqa: E402
from realms.maps import HUES, SURFACE, TEXT  # noqa: E402
from realms.matrix import RANKS, CellTaxa  # noqa: E402

OTHER = "#b9b6ad"  # parts of other realms, for context; darker than the globe, which shows where no part is mapped
MEET_M = 1.0  # metres on the map within which the outlines of two groups of parts are taken to be the line where they meet
CHANGE_COLOURS = {"kept": "#3d3d3a", "added": "#d6332e", "removed": "#2a78d6"}
CHANGES = ("kept", "added", "removed")
REALM_TABLE_SCHEMA = {
    "realm": pl.Int64, "area_km2": pl.Float64, "land_km2": pl.Float64, "sea_km2": pl.Float64, "area_share": pl.Float64,
    "units": pl.Int64, "land_units": pl.Int64, "sea_units": pl.Int64, "parts": pl.Int64, "imputed_parts": pl.Int64,
    "species": pl.Int64, "endemic_species": pl.Int64, "sub_realms": pl.Int64, "sub_realm_units": pl.Int64, "sub_realm_km2": pl.Float64,
}


# ---------------------------------------------------------------------------------------------------------------------------- colours


def part_neighbours(parts: Collection[str], cell_adjacency: Mapping[str, Collection[str]]) -> dict[str, set[str]]:
    """Which of `parts` are next to which: the same kind in a neighbouring cell and the other kind in the part's own cell, as
    `grid.part_adjacency` has it, from a cell adjacency such as `boundaries.edge_adjacency`."""
    names = set(parts)
    near: dict[str, set[str]] = {p: set() for p in names}
    for part in names:
        cell, _, kind = part.rpartition(":")
        other = f"{cell}:{'sea' if kind == 'land' else 'land'}"
        found = [other] + [f"{c}:{kind}" for c in cell_adjacency.get(cell, ())]
        near[part].update(p for p in found if p in names)
    return near


def realm_neighbours(labels: Mapping[str, int], adjacency: Mapping[str, Collection[str]]) -> dict[int, set[int]]:
    """Realm to the realms that one of its parts is next to."""
    neighbours: dict[int, set[int]] = {r: set() for r in set(labels.values())}
    for part, near in adjacency.items():
        if part not in labels:
            continue
        for other in near:
            if other in labels and labels[other] != labels[part]:
                neighbours[labels[part]].add(labels[other])
    return neighbours


def overlap(a: Mapping[str, int], b: Mapping[str, int], areas: Mapping[str, float]) -> dict[tuple[int, int], float]:
    """Area shared by realm `ra` of map `a` and realm `rb` of map `b`, over the parts both maps hold, keyed `(ra, rb)`."""
    shared: dict[tuple[int, int], float] = defaultdict(float)
    for part in a.keys() & b.keys():
        shared[(a[part], b[part])] += areas[part]
    return dict(shared)


def distinct_colours(labels: Mapping[str, int], areas: Mapping[str, float], palette: list[str] = HUES) -> dict[int, str]:
    """A colour of its own for each realm, the palette's first going to the largest realm (a tie to the lower realm), for a map with no
    more realms than the palette has colours."""
    size: dict[int, float] = defaultdict(float)
    for part, realm in labels.items():
        size[realm] += areas[part]
    if len(size) > len(palette):
        raise ValueError(f"{len(size)} realms and only {len(palette)} colours")
    return {realm: palette[i] for i, realm in enumerate(sorted(size, key=lambda r: (-size[r], r)))}


def match_colours(
    reference: Mapping[str, int],
    colours: Mapping[int, str],
    other: Mapping[str, int],
    areas: Mapping[str, float],
    neighbours: Mapping[int, Collection[int]] | None = None,
    palette: list[str] = HUES,
) -> dict[int, str]:
    """A colour for each realm of `other` that follows `reference`'s colours: realms are paired one to one by the largest shared area first
    (a tie goes to the lower realms), so a realm takes the colour of its largest-overlap parent or child.

    A realm left unpaired (the extra realm of a finer map) takes the first colour of `palette` that no paired realm and none of its
    `neighbours` has; failing that, the first its neighbours lack; failing that, the first of the palette.
    """
    shared = overlap(reference, other, areas)
    pairs = sorted(((area, ra, rb) for (ra, rb), area in shared.items() if area > 0), key=lambda t: (-t[0], t[1], t[2]))
    result: dict[int, str] = {}
    taken = set()
    for _, ra, rb in pairs:
        if rb in result or ra in taken:
            continue
        result[rb] = colours[ra]
        taken.add(ra)
    neighbours = neighbours or {}
    for realm in sorted(set(other.values()) - set(result)):
        beside = {result[n] for n in neighbours.get(realm, ()) if n in result}
        used = set(result.values())
        result[realm] = next(
            (c for c in palette if c not in used and c not in beside), next((c for c in palette if c not in beside), palette[0])
        )
    return result


# ---------------------------------------------------------------------------------------------------------------------------- tables


def part_areas(cells: pl.DataFrame) -> dict[str, float]:
    """km2 of every land and sea part of the cells (`cell`, `land_km2`, `sea_km2`, as the units stage's `cells.parquet`)."""
    areas = {}
    for cell, land, sea in cells.select("cell", "land_km2", "sea_km2").iter_rows():
        areas[f"{cell}:land"] = land
        areas[f"{cell}:sea"] = sea
    return areas


def endemic_species(labels: Labels, ct: CellTaxa, tolerance: float = 0.95) -> dict[int, tuple[int, int]]:
    """Realm to (endemic species, species present), by `endemism.realm_share`'s rule on the species columns only.

    `labels` maps units (rows of `ct`) to realms. A species is endemic to a realm when at least `tolerance` of all the units that hold it,
    in the matrix, are in that realm.
    """
    species = ct.select(["species"]) if (ct.rank != RANKS.index("species")).any() else ct
    cells = sorted(labels)
    realms = sorted(set(labels.values()))
    index = {r: i for i, r in enumerate(realms)}
    table = species.rows(cells)
    membership = coo_matrix(
        (np.ones(len(cells), dtype=np.int32), (np.arange(len(cells)), [index[labels[c]] for c in cells])), shape=(len(cells), len(realms))
    ).tocsr()
    counts = (table.xt @ membership).tocoo()
    present = np.bincount(counts.col, minlength=len(realms))
    endemic = np.bincount(counts.col[counts.data / table.totals[counts.row] >= tolerance], minlength=len(realms))
    return {r: (int(endemic[i]), int(present[i])) for r, i in index.items()}


def sub_realm_table(candidates: pl.DataFrame, unit_areas: Mapping[str, float]) -> pl.DataFrame:
    """`realm`, `sub_realm`, `units`, `area_km2` of every sub-realm (a module the size floor merged into a realm) of one k of
    `candidates.parquet`, largest first within each realm."""
    rows = candidates.filter(pl.col("sub_realm").is_not_null()).with_columns(
        pl.col("unit").replace_strict(dict(unit_areas), return_dtype=pl.Float64).alias("area_km2")
    )
    return (
        rows.group_by("realm", "sub_realm")
        .agg(pl.len().cast(pl.Int64).alias("units"), pl.col("area_km2").sum())
        .sort(["realm", "area_km2", "sub_realm"], descending=[False, True, False])
    )


def realm_table(
    parts: pl.DataFrame,
    candidates: pl.DataFrame,
    areas: Mapping[str, float],
    unit_areas: Mapping[str, float],
    endemic: Mapping[int, tuple[int, int]],
) -> pl.DataFrame:
    """One row per realm of one k: its area from its parts (`area_km2`, `land_km2`, `sea_km2`, `area_share` of the map), its units (all,
    land and sea), its parts and imputed parts, its species and endemic species (`endemic_species`), and its sub-realms with their units and area.

    `parts` is `parts.parquet` and `candidates` `candidates.parquet`, both at that k. Largest realm first.
    """
    for name, table in (("parts", parts), ("candidates", candidates)):
        if table["k"].n_unique() != 1:
            raise ValueError(f"the {name} must be one k, not {sorted(table['k'].unique())}")
    if parts["k"][0] != candidates["k"][0]:
        raise ValueError(f"the parts are k={parts['k'][0]} and the candidates k={candidates['k'][0]}")
    mapped = dict(parts.filter(pl.col("unit").is_not_null()).select("unit", "realm").unique().iter_rows())
    direct = dict(candidates.select("unit", "realm").iter_rows())
    if mapped != {u: direct.get(u) for u in mapped} or set(mapped) != set(direct):
        raise ValueError("the realms of the units in the parts are not those of the candidates")
    sized = parts.with_columns(
        pl.col("part").replace_strict(dict(areas), return_dtype=pl.Float64).alias("km2"), pl.col("part").str.ends_with(":land").alias("land")
    )
    by_part = sized.group_by("realm").agg(
        pl.col("km2").sum().alias("area_km2"),
        pl.col("km2").filter(pl.col("land")).sum().alias("land_km2"),
        pl.col("km2").filter(~pl.col("land")).sum().alias("sea_km2"),
        pl.len().cast(pl.Int64).alias("parts"),
        pl.col("imputed").sum().cast(pl.Int64).alias("imputed_parts"),
    )
    by_unit = candidates.with_columns(pl.col("unit").str.ends_with(":land").alias("land")).group_by("realm").agg(
        pl.len().cast(pl.Int64).alias("units"),
        pl.col("land").sum().cast(pl.Int64).alias("land_units"),
        (~pl.col("land")).sum().cast(pl.Int64).alias("sea_units"),
    )
    subs = sub_realm_table(candidates, unit_areas).group_by("realm").agg(
        pl.len().cast(pl.Int64).alias("sub_realms"), pl.col("units").sum().cast(pl.Int64).alias("sub_realm_units"), pl.col("area_km2").sum().alias("sub_realm_km2")
    )
    species = pl.DataFrame(
        [(r, present, count) for r, (count, present) in endemic.items()], schema={"realm": pl.Int64, "species": pl.Int64, "endemic_species": pl.Int64}, orient="row"
    )
    table = (
        by_part.join(by_unit, on="realm", how="left")
        .join(species, on="realm", how="left")
        .join(subs, on="realm", how="left")
        .with_columns(
            (pl.col("area_km2") / pl.col("area_km2").sum()).alias("area_share"),
            pl.col("units", "land_units", "sea_units", "species", "endemic_species", "sub_realms", "sub_realm_units").fill_null(0),
            pl.col("sub_realm_km2").fill_null(0.0),
        )
    )
    return table.select(list(REALM_TABLE_SCHEMA)).cast(REALM_TABLE_SCHEMA).sort(["area_km2", "realm"], descending=[True, False])


def markdown(table: pl.DataFrame, formats: Mapping[str, str] | None = None) -> str:
    """A pipe table of `table`, each column formatted with `formats[column]` (a format spec) if given, an empty cell for a null."""
    formats = formats or {}

    def cell(column: str, value) -> str:
        if value is None:
            return ""
        return format(value, formats[column]) if column in formats else str(value)

    lines = ["| " + " | ".join(table.columns) + " |", "|" + "|".join("---" for _ in table.columns) + "|"]
    lines += ["| " + " | ".join(cell(c, v) for c, v in zip(table.columns, row, strict=True)) + " |" for row in table.iter_rows()]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------------------------------------------------------- without the floor


def unfloored_maps(partition: Partition, ct: CellTaxa, adjacency: Mapping[str, set[str]], ks: Collection[int]) -> tuple[Labels, dict[int, Labels]]:
    """Maps at `ks` realms without the size floor: Infomap's finest modules with their separate pieces resolved (`candidates.resolve_pieces`,
    as `candidate_maps` does first), then `candidates.merge_sequence` straight on those, with no floor in between.

    Returns the units before any merge and the map at each of `ks` that the number of modules allows. `ct` and `adjacency` are those the
    select stage scored (`explore.assigned_only`).
    """
    units = resolve_pieces(partition.leaf_labels(), dict(adjacency), ct)
    wanted = sorted(k for k in ks if k <= len(set(units.values())))
    if not wanted:
        return units, {}
    maps = merge_sequence(units, ct, k_max=max(wanted), k_min=min(wanted))
    return units, {k: maps[k] for k in wanted}


# ---------------------------------------------------------------------------------------------------------------------------- figures


def tint(colour: str, share: float) -> str:
    """`colour` moved `share` of the way to white, as hex."""
    rgb = np.array(to_rgb(colour))
    return to_hex(rgb + (1.0 - rgb) * share)


def _flat(geoms: gpd.GeoDataFrame, crs: str = EQUAL_EARTH) -> gpd.GeoDataFrame:
    shapes = [_projected(s) for s in geoms.geometry] if crs == EQUAL_EARTH else [_to_crs(s, crs) for s in geoms.geometry]
    vanished = [p for p, g, s in zip(geoms["part"], geoms.geometry, shapes, strict=True) if not g.is_empty and s.is_empty]
    if vanished:
        raise ValueError(f"{len(vanished)} parts vanish when projected, such as {vanished[0]!r}")
    return gpd.GeoDataFrame(geoms.drop(columns="geometry"), geometry=shapes, crs=crs)


def _to_crs(shape, crs: str):
    flat = gpd.GeoSeries([_densified(shape)], crs=WGS84).to_crs(crs).iloc[0]
    return flat if flat.is_valid else _repaired(flat)


def _label(ax, shape, text: str, size: float = 9, clip: bool = False) -> None:
    biggest = max(shapely.get_parts(shape), key=lambda p: p.area)
    spot = polylabel(biggest, tolerance=max(1.0, biggest.length / 2000))
    ax.text(spot.x, spot.y, text, ha="center", va="center", fontsize=size, fontweight="bold", color=TEXT, zorder=6, clip_on=clip,
            path_effects=[withStroke(linewidth=2.5, foreground=SURFACE)])


def view(flat: gpd.GeoDataFrame, share: float = 0.98) -> tuple[float, float, float, float]:
    """West, south, east, north of the window that holds `share` of the area of the projected parts: the area-weighted quantiles of their
    representative points, `(1 - share) / 2` in from either side on each axis, widened to the whole of every part whose point is inside."""
    points = flat.geometry.representative_point()
    weight = flat.geometry.area.to_numpy()
    cut = (1.0 - share) / 2.0

    def quantiles(values: np.ndarray) -> tuple[float, float]:
        order = np.argsort(values)
        cumulative = np.cumsum(weight[order]) / weight.sum()
        low = values[order][np.searchsorted(cumulative, cut, side="left")]
        high = values[order][min(len(values) - 1, np.searchsorted(cumulative, 1.0 - cut, side="left"))]
        return float(low), float(high)

    x0, x1 = quantiles(points.x.to_numpy())
    y0, y1 = quantiles(points.y.to_numpy())
    inside = points.x.between(x0, x1) & points.y.between(y0, y1)
    return tuple(float(v) for v in flat[inside].total_bounds)


def plot_realm(
    geoms: gpd.GeoDataFrame,
    parts: pl.DataFrame,
    candidates: pl.DataFrame,
    realm: int,
    colour: str,
    out: Path,
    title: str,
    floored: Mapping[str, int] | None = None,
    label_sub_realms: int = 6,
    view_share: float = 0.995,
) -> dict:
    """One realm of a map, zoomed to it, with every other realm in grey: its parts in `colour`, those of a sub-realm lighter and outlined,
    imputed parts hatched, and, with `floored` (unit to the unit it was above the floor in, such as the candidates' largest k), the line
    between the units above the floor that the realm merges. The largest `label_sub_realms` sub-realms carry their number. The view is `view`'s
    window on `view_share` of the realm's area, and the legend counts the parts left outside it. Written to `<out>.png`.

    `parts` and `candidates` are one k of `parts.parquet` and `candidates.parquet`. Returns what was drawn, for the log.
    """
    sub_of = dict(candidates.filter(pl.col("sub_realm").is_not_null()).select("unit", "sub_realm").iter_rows())
    info = {p: (r, u, i) for p, r, u, i in parts.select("part", "realm", "unit", "imputed").iter_rows()}
    frame = geoms[geoms["part"].isin(list(info))].reset_index(drop=True)
    if len(frame) != len(info):
        raise ValueError(f"{len(info) - len(frame)} parts of the map have no geometry")
    frame = frame[~frame.geometry.is_empty].reset_index(drop=True)
    frame["realm"] = [info[p][0] for p in frame["part"]]
    frame["unit"] = [info[p][1] for p in frame["part"]]
    frame["imputed"] = [bool(info[p][2]) for p in frame["part"]]
    frame["sub_realm"] = [sub_of.get(u) if u is not None else None for u in frame["unit"]]
    flat = _flat(frame)
    mine = flat[flat["realm"] == realm]
    # dissolved before projecting, as plot_parts does: cells that disagree by a hair at a shared vertex are snapped together in degrees,
    # and projected first they would leave seams along the grid's own lines
    lonlat = frame[frame["realm"] == realm]

    def dissolved(rows) -> object:
        return _projected(_dissolve(lonlat.geometry.loc[rows.index]))

    if mine.empty:
        raise ValueError(f"realm {realm} has no parts")

    fig, ax = plt.subplots(figsize=(12, 7.4), facecolor=SURFACE)
    globe = _projected(shapely.box(-180, -90, 180, 90))
    ax.fill(*_polygons(globe)[0].exterior.xy, facecolor=GLOBE, edgecolor=OUTLINE, linewidth=0.5, zorder=0)
    rest = flat[flat["realm"] != realm]
    rest.plot(ax=ax, color=OTHER, edgecolor=OTHER, linewidth=0.2, aspect=None, zorder=1)
    in_sub = mine["sub_realm"].notna()
    fills = [tint(colour, 0.55) if s else colour for s in in_sub]
    mine.plot(ax=ax, color=fills, edgecolor=fills, linewidth=0.2, aspect=None, zorder=2)
    imputed = mine[mine["imputed"]]
    if len(imputed):
        with plt.rc_context({"hatch.linewidth": 0.5}):
            imputed.plot(ax=ax, facecolor="none", edgecolor="#0b0b0bb0", hatch="//////", linewidth=0, aspect=None, zorder=3)
    outline = dissolved(mine)
    gpd.GeoSeries([outline], crs=EQUAL_EARTH).boundary.plot(ax=ax, color=TEXT, linewidth=0.8, aspect=None, zorder=5)
    subs = mine[in_sub]
    sub_shapes = {int(s): dissolved(g) for s, g in subs.groupby("sub_realm")}
    if sub_shapes:
        gpd.GeoSeries(list(sub_shapes.values()), crs=EQUAL_EARTH).boundary.plot(ax=ax, color=TEXT, linewidth=0.35, aspect=None, zorder=4)
    sizes = dict(subs.assign(a=subs.geometry.area).groupby("sub_realm")["a"].sum())
    for s in sorted(sizes, key=lambda s: (-sizes[s], s))[:label_sub_realms]:
        _label(ax, sub_shapes[int(s)], str(int(s)), size=7, clip=True)  # a zoomed view can leave a sub-realm outside it
    merged_units = []
    if floored:
        with_unit = mine[mine["unit"].notna()]
        groups = with_unit.groupby(with_unit["unit"].map(lambda u: floored[u]))
        merged_units = sorted(int(g) for g in groups.groups)
        if len(merged_units) > 1:
            # only where two of them meet: a group's own outline also runs along coasts, holes and imputed parts
            outlines = [dissolved(g).boundary for _, g in groups]
            lines = [
                shapely.intersection(a, shapely.buffer(b, MEET_M))
                for i, a in enumerate(outlines) for b in outlines[i + 1 :]
            ]
            gpd.GeoSeries(lines, crs=EQUAL_EARTH).plot(ax=ax, color="#ffffff", linewidth=1.4, linestyle=(0, (3, 2)), aspect=None, zorder=4.5)

    west, south, east, north = view(mine, share=view_share)
    centres = mine.geometry.representative_point()
    outside = int((~centres.x.between(west, east) | ~centres.y.between(south, north)).sum())
    gw, gs, ge, gn = globe.bounds
    pad_x, pad_y = 0.06 * (east - west) + 2e5, 0.06 * (north - south) + 2e5
    ax.set_xlim(max(gw, west - pad_x), min(ge, east + pad_x))
    ax.set_ylim(max(gs, south - pad_y), min(gn, north + pad_y))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left", color=TEXT, fontsize=12)
    keys = [
        Patch(facecolor=colour, edgecolor="none", label=f"realm {realm}, units above the floor"),
        Patch(facecolor=tint(colour, 0.55), edgecolor=TEXT, linewidth=0.35, label=f"sub-realms ({len(sub_shapes)}), outlined; the largest numbered"),
        Patch(facecolor="none", edgecolor="#0b0b0bb0", hatch="//////", label=f"imputed parts ({len(imputed)})"),
        Patch(facecolor=OTHER, edgecolor="none", label="other realms"),
        Patch(facecolor=GLOBE, edgecolor=OUTLINE, linewidth=0.3, label="not in the map (no records)"),
    ]
    if len(merged_units) > 1:
        keys.append(Line2D([], [], color="#7a7a7a", linewidth=1.4, linestyle=(0, (3, 2)), label=f"between the units above the floor ({', '.join(map(str, merged_units))})"))
    if outside:
        keys.append(Patch(facecolor="none", edgecolor="none", label=f"zoomed to {view_share:.1%} of its area: {outside:,} parts lie outside the view"))
    ax.legend(handles=keys, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8, labelcolor=TEXT)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.parent / f"{out.name}.png", dpi=DPI, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
    return {"realm": realm, "parts": len(mine), "sub_realms": len(sub_shapes), "imputed": len(imputed), "units_above_floor": merged_units, "outside_view": outside}


def edge_changes(base: pl.DataFrame, other: pl.DataFrame) -> pl.DataFrame:
    """`cell_a`, `cell_b`, `km`, `change` of two boundaries (`boundaries.boundary`): `kept` in both, `added` only in `other`, `removed` only in `base`."""
    key = ["cell_a", "cell_b"]
    kept = base.join(other.select(key), on=key, how="semi").with_columns(pl.lit("kept").alias("change"))
    added = other.join(base.select(key), on=key, how="anti").with_columns(pl.lit("added").alias("change"))
    removed = base.join(other.select(key), on=key, how="anti").with_columns(pl.lit("removed").alias("change"))
    return pl.concat([kept, added, removed]).sort(key)


def edge_lines(gpkg: Path, pairs: list[tuple[str, str]]) -> list:
    """The edge each pair of neighbouring cells shares, as a line in longitude and latitude moved to centre within -180 and 180.

    It is found as `boundaries.edge_lengths` finds it: the boundary of one cell inside a band of `SNAP` degrees around the other's, with
    the first cell moved a turn of longitude either way, keeping the longest piece. A pair that shares nothing raises.
    """
    cells, shapes, _ = _cells(gpkg)
    index = {c: i for i, c in enumerate(cells)}
    result = []
    for a, b in pairs:
        line_a, band_b = shapes[index[a]].boundary, shapely.buffer(shapes[index[b]].boundary, SNAP)
        found = max((shapely.intersection(_shifted(np.array([line_a]), turn)[0], band_b) for turn in (-360.0, 0.0, 360.0)), key=shapely.length)
        pieces = [p for p in shapely.get_parts(shapely.line_merge(found)) if p.geom_type == "LineString"] if found.length > 0 else []
        if not pieces:
            raise ValueError(f"cells {a} and {b} share no edge")
        line = MultiLineString(pieces) if len(pieces) > 1 else pieces[0]
        middle = shapely.centroid(line).x
        result.append(shapely.transform(line, lambda xy, t=360.0 * np.round(middle / 360.0): xy - [t, 0.0]))
    return result


def plot_boundary_change(
    geoms: gpd.GeoDataFrame,
    labels: Mapping[str, int],
    colours: Mapping[int, str],
    changes: pl.DataFrame,
    lines: list,
    regions: gpd.GeoDataFrame,
    crs: str,
    bounds: tuple[float, float, float, float],
    out: Path,
    title: str,
    focus: Collection[str] | None = None,
) -> None:
    """The edges of `changes` (`edge_changes`, one line each from `edge_lines`) drawn over the land parts of `geoms`, lightly in the colour of
    their realm in `labels`, with the outlines of `regions`, in the projection `crs`, cut to `bounds` (west, south, east, north, degrees).
    With `focus`, only those parts are coloured and numbered, and the other land is grey. The legend gives the km of each kind of change.
    Written to `<out>.png`."""
    box = shapely.box(*bounds)
    land = geoms[geoms["part"].str.endswith(":land") & geoms.geometry.intersects(box)]
    land = land[land["part"].isin(list(labels))]
    flat = _flat(land.reset_index(drop=True), crs)
    shown = flat["part"].isin(list(focus)) if focus is not None else flat["part"].notna()
    fig, ax = plt.subplots(figsize=(12, 9), facecolor=SURFACE)
    flat[~shown].plot(ax=ax, color=OTHER, edgecolor="none", linewidth=0, aspect=None, zorder=1)
    flat = flat[shown]
    flat.plot(ax=ax, color=[tint(colours[labels[p]], 0.6) for p in flat["part"]], edgecolor="none", linewidth=0, aspect=None, zorder=1)
    outlines = regions[regions.geometry.intersects(box)].copy()
    outlines["geometry"] = [shapely.intersection(shapely.make_valid(g), shapely.box(-180, -90, 180, 90)) for g in outlines.geometry]
    outlines.to_crs(crs).boundary.plot(ax=ax, color="#8f8d86", linewidth=0.35, aspect=None, zorder=2)
    keys = []
    for change in CHANGES:
        chosen = [line for line, c in zip(lines, changes["change"], strict=True) if c == change]
        km = changes.filter(pl.col("change") == change)["km"].sum()
        if chosen:
            gpd.GeoSeries(chosen, crs=WGS84).to_crs(crs).plot(ax=ax, color=CHANGE_COLOURS[change], linewidth=1.6 if change != "kept" else 1.2, aspect=None, zorder=3 + CHANGES.index(change))
        keys.append(Line2D([], [], color=CHANGE_COLOURS[change], linewidth=1.6, label=f"{change}: {km:,.0f} km, {len(chosen)} edges"))
    for realm in sorted(set(labels[p] for p in flat["part"])):
        shape = _dissolve(flat.geometry[flat["part"].map(labels) == realm])
        if not shape.is_empty:
            _label(ax, shape, str(realm), size=10)
    keys.append(Line2D([], [], color="#8f8d86", linewidth=0.6, label="TDWG level 3"))
    corners = gpd.GeoSeries([_densified(box, 1.0)], crs=WGS84).to_crs(crs).total_bounds
    ax.set_xlim(corners[0], corners[2])
    ax.set_ylim(corners[1], corners[3])
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left", color=TEXT, fontsize=12)
    ax.legend(handles=keys, loc="lower left", frameon=True, facecolor=SURFACE, edgecolor="none", fontsize=9, labelcolor=TEXT, framealpha=0.9)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.parent / f"{out.name}.png", dpi=DPI, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
