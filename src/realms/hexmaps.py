"""Equal-area maps of the land and sea parts of ISEA3H cells.

A part is drawn from the cell as `grid.py` draws it (`_planar`, around the cell's centre longitude) and from the same Natural Earth land
as `grid._land_km2`: a land part is the cell where it meets the land and a sea part is the rest of the cell. A cell that crosses the
antimeridian reaches past 180 degrees in that drawing, so for the map it is cut there and each piece moved back, and a cell over a pole
is already whole to the pole (`_planar` closes it there), so it draws as a cap. The map is Equal Earth (EPSG:8857).
"""

import logging
from collections import defaultdict
from collections.abc import Collection, Iterable, Mapping
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shapely  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.patheffects import withStroke  # noqa: E402
from shapely import MultiPolygon, Polygon  # noqa: E402
from shapely.ops import polylabel  # noqa: E402

from realms.grid import SNAP, ZERO_KM2, _km2, _land_overlaps, _neighbours, _planar, part_id  # noqa: E402
from realms.maps import SURFACE, TEXT, _separating_colours  # noqa: E402

log = logging.getLogger(__name__)
WGS84 = "EPSG:4326"
EQUAL_EARTH = "EPSG:8857"
# Degrees. Edges are straight in longitude and latitude, which are curves on the map, so they are cut into short ones before projecting.
SEGMENT = 0.25
GLOBE = "#eceae4"
OUTLINE = "#52514e"
DPI = 200
SLIVER_M2 = 1e6  # m2 on the map: a piece the repair of a projected part leaves that is smaller, and thin, is dropped
SLIVER_SHAPE = 0.05  # 4 pi area / perimeter^2 under which a small piece counts as thin


def _polygons(geometry) -> list[Polygon]:
    """The polygons of a geometry, however deeply collected, leaving out the lines and points that an intersection or a repair can leave."""
    found, todo = [], [geometry]
    while todo:
        shape = todo.pop()
        if shape.geom_type == "Polygon":
            if not shape.is_empty:
                found.append(shape)
        elif shape.geom_type in ("MultiPolygon", "GeometryCollection"):
            todo += list(shape.geoms)
    return found


def _part_of(polygons: list[Polygon]) -> MultiPolygon:
    """The polygons as one geometry, empty if there is no area worth the name."""
    shape = MultiPolygon(polygons)
    return shape if _km2(shape) >= ZERO_KM2 else MultiPolygon()


def _without_slivers(shape) -> MultiPolygon:
    """The polygons of a part, within -180 and 180, that have an area worth the name, as one geometry.

    Taking the land from a cell leaves slivers of no area along the cell's edges where the coast crosses them, which the map would draw as
    short ticks of the sea's colour reaching into the land. Measured after the cut at the antimeridian, where the geodesic area is that of
    the shape (beyond 180 degrees a sliver can measure as a hundred km2).
    """
    return MultiPolygon([p for p in _polygons(shape) if _km2(p) >= ZERO_KM2])


def _cut_at_antimeridian(shape) -> MultiPolygon:
    """The shape, drawn within 180 degrees of the centre longitude, as the pieces either side of -180 and 180, each moved back into range."""
    pieces = []
    west, _, east, _ = shape.bounds
    for turn in (-1, 0, 1):
        if east > -180.0 + 360.0 * turn and west < 180.0 + 360.0 * turn:
            # not clip_by_rect, which fails with "Invalid number of points in LinearRing" on some valid parts at the coast
            inside = shapely.intersection(shape, shapely.box(-180.0 + 360.0 * turn, -90.0, 180.0 + 360.0 * turn, 90.0))
            pieces += [shapely.transform(p, lambda xy, turn=turn: xy - [360.0 * turn, 0.0]) for p in _polygons(inside) if p.area > 0]
    return MultiPolygon(pieces)


def _split_part(part: str) -> tuple[str, str]:
    cell, colon, kind = part.rpartition(":")
    if not colon:
        raise ValueError(f"{part!r} is not a part id, which is <cell>:land or <cell>:sea")
    part_id(cell, kind)
    return cell, kind


def part_geometries(gpkg: Path, land_zip: Path, parts: Iterable[str]) -> gpd.GeoDataFrame:
    """`part` and `geometry` (EPSG:4326) of the given parts, sorted by id, each a MultiPolygon that no edge of runs round the globe.

    A land part is its cell where it meets Natural Earth land and a sea part is the cell without that land. A part with less than
    `ZERO_KM2` of area, such as the sea part of a cell wholly on land, is empty, and so is any polygon of a part with less than that. A part
    that has area as drawn around its cell but none once cut at the antimeridian (a sliver measured beyond 180 degrees) is empty too, and
    logged as a warning.
    """
    wanted = {}
    for part in sorted(set(parts)):
        wanted.setdefault(_split_part(part)[0], []).append(part)

    frame = gpd.read_file(gpkg)
    if frame["geogrid_id"].dtype.kind == "f":
        raise ValueError(f"geogrid_id was read as floating point, which loses digits: {gpkg}")
    frame = frame.assign(cell=[str(i) for i in frame["geogrid_id"]])
    frame = frame[frame["cell"].isin(wanted)]
    missing = sorted(set(wanted) - set(frame["cell"]))
    if missing:
        raise ValueError(f"cell {missing[0]!r} of part {wanted[missing[0]][0]!r} is not in the grid {gpkg}")

    cells = list(frame["cell"])
    shapes = np.array([_planar(cell, shape, lon) for cell, shape, lon in zip(cells, frame.geometry, frame["lon_c"], strict=True)], dtype=object)
    where, overlap = _land_overlaps(shapes, land_zip)
    on_land: dict[int, list[Polygon]] = defaultdict(list)
    for cell, piece in zip(where.tolist(), overlap, strict=True):
        on_land[cell] += _polygons(piece)

    rows = {}
    for index, (cell, shape) in enumerate(zip(cells, shapes, strict=True)):
        land = shapely.union_all(on_land[index]) if on_land[index] else Polygon()
        rows[part_id(cell, "land")] = _part_of(_polygons(land))
        rows[part_id(cell, "sea")] = _part_of(_polygons(shapely.difference(shape, land)))
    chosen = sorted(part for group in wanted.values() for part in group)
    drawn = [_without_slivers(_cut_at_antimeridian(rows[part])) for part in chosen]
    lost = [part for part, shape in zip(chosen, drawn, strict=True) if shape.is_empty and not rows[part].is_empty]
    if lost:
        log.warning(
            f"{len(lost)} parts have area as the grid measures them but none once cut at the antimeridian and rid of slivers, so they are "
            f"empty on the map: {', '.join(lost)}"
        )
    return gpd.GeoDataFrame({"part": chosen}, geometry=drawn, crs=WGS84)


def _realm_neighbours(realm_of: Mapping[str, int], adjacency: Mapping[str, Collection[str]]) -> dict[int, set[int]]:
    neighbours: dict[int, set[int]] = defaultdict(set)
    for part, near in adjacency.items():
        if part not in realm_of:
            continue
        for other in near:
            if other in realm_of and realm_of[other] != realm_of[part]:
                neighbours[realm_of[part]].add(realm_of[other])
                neighbours[realm_of[other]].add(realm_of[part])
    return neighbours


def _touching(geoms: gpd.GeoDataFrame) -> dict[str, set[str]]:
    """Parts that share boundary on the map, for when the caller has no `part_adjacency`. Land and sea of one cell share the coast."""
    names = geoms["part"].tolist()
    shapes = np.array(geoms.geometry.values, dtype=object)
    near: dict[str, set[str]] = {name: set() for name in names}
    nonempty = np.flatnonzero(~shapely.is_empty(shapes))
    for a, b in _neighbours(shapes[nonempty]):
        near[names[nonempty[a]]].add(names[nonempty[b]])
        near[names[nonempty[b]]].add(names[nonempty[a]])
    return near


def _dissolve(shapes) -> MultiPolygon | Polygon:
    # Cells disagree by up to 1e-8 degrees where they should share a vertex, which would leave seams
    return shapely.union_all(list(shapes), grid_size=SNAP)


def _densified(shape, step: float = SEGMENT) -> MultiPolygon:
    """The shape with no edge longer than `step` degrees in either direction. (`shapely.segmentize` returns an empty polygon for some parts, GEOS 3.14.)"""

    def ring(coords) -> np.ndarray:
        xy = np.asarray(coords)
        pieces = np.maximum(1, np.ceil(np.abs(np.diff(xy, axis=0)).max(axis=1) / step)).astype(int)
        start = np.repeat(np.arange(len(pieces)), pieces)
        along = (np.arange(len(start)) - np.repeat(np.cumsum(pieces) - pieces, pieces)) / pieces[start]
        return np.vstack([xy[start] + along[:, None] * (xy[start + 1] - xy[start]), xy[-1:]])

    return MultiPolygon(
        [Polygon(ring(p.exterior.coords), [ring(h.coords) for h in p.interiors]) for p in _polygons(shape)]
    )


def _repaired(flat) -> MultiPolygon:
    """A projected shape made valid, without the slivers the repair cuts off.

    Where a coast runs along a cell's edge, the two are cut into short edges at different points, and the projected chords of one curve
    cross each other. `make_valid` then leaves thin slivers between them, metres wide, which the map draws as lines of the part's colour
    into its neighbour. A piece under `SLIVER_M2` that is thin (`SLIVER_SHAPE`, 1 for a circle) is such a sliver and is dropped, unless it
    is the largest piece: a part that is itself a thin strip of sea along a coast keeps that strip.
    """
    pieces = _polygons(shapely.make_valid(flat))
    if not pieces:
        return MultiPolygon()
    largest = max(pieces, key=lambda p: p.area)
    return MultiPolygon([p for p in pieces if p is largest or p.area >= SLIVER_M2 or 4 * np.pi * p.area / p.length**2 >= SLIVER_SHAPE])


def _projected(shape):
    """The shape on the map. Projecting bends edges that lay along each other, such as a coast that follows the cell's edge, across each other."""
    flat = gpd.GeoSeries([_densified(shape)], crs=WGS84).to_crs(EQUAL_EARTH).iloc[0]
    if flat.is_valid:
        return flat
    return _repaired(flat)


def _flag(imputed: Collection[str] | Mapping[str, bool] | None) -> set[str]:
    if imputed is None:
        return set()
    return {part for part, yes in imputed.items() if yes} if isinstance(imputed, Mapping) else set(imputed)


def plot_parts(
    geoms: gpd.GeoDataFrame,
    labels: Mapping[str, int],
    out: Path,
    title: str,
    sub_realms: Mapping[str, int] | None = None,
    imputed: Collection[str] | Mapping[str, bool] | None = None,
    adjacency: Mapping[str, Collection[str]] | None = None,
    colours: Mapping[int, str] | None = None,
) -> None:
    """An Equal Earth map of `geoms` (from `part_geometries`) coloured by the realm in `labels`, written to `<out>.png`, with `<out>.gpkg`.

    - `labels` maps every part to its realm. Colours only tell neighbouring realms apart, so each realm's number is printed on it.
    - `sub_realms` maps parts to a sub-realm, which is outlined; parts it leaves out have none.
    - `imputed` is the parts, or a mapping of part to flag, whose realm was imputed rather than found; they are hatched.
    - `adjacency` is which parts are next to which, as `grid.part_adjacency` gives it; by default the parts that touch on the map.
    - `colours` gives each realm its colour, to keep one realm's colour across maps; by default they only tell neighbours apart.

    A part with area that comes out empty on the map raises `ValueError`, before anything is written.

    `<out>.gpkg` holds a layer `parts` (part, realm, sub_realm, imputed) and a layer `realms` (realm, the parts dissolved into one
    geometry, and `dissolved`, how many parts that was), both in EPSG:4326.
    """
    parts = geoms["part"].tolist()
    unlabelled = [p for p in parts if p not in labels]
    if unlabelled:
        raise ValueError(f"{len(unlabelled)} parts have no realm in the labels, such as {unlabelled[0]!r}")
    realm_of = {p: int(labels[p]) for p in parts}
    sub_of = {p: int(sub_realms[p]) for p in parts if sub_realms and p in sub_realms}
    flagged = _flag(imputed)

    table = gpd.GeoDataFrame(
        {
            "part": parts,
            "realm": [realm_of[p] for p in parts],
            "sub_realm": pd.array([sub_of.get(p) for p in parts], dtype="Int64"),
            "imputed": [p in flagged for p in parts],
        },
        geometry=geoms.geometry.values,
        crs=WGS84,
    )
    table = table[~table.geometry.is_empty].reset_index(drop=True)
    # projected before anything is written: a part with area that comes out empty would be a hole in the map that nothing points to
    projected = [_projected(s) for s in table.geometry]
    vanished = [part for part, shape in zip(table["part"], projected, strict=True) if shape.is_empty]
    if vanished:
        raise ValueError(f"{len(vanished)} parts vanish when projected to {EQUAL_EARTH}, such as {vanished[0]!r}: the map would have holes")
    realms = sorted(set(table["realm"]))
    dissolved = gpd.GeoDataFrame(
        {
            "realm": realms,
            "dissolved": [int((table["realm"] == r).sum()) for r in realms],
        },
        geometry=[_dissolve(table.geometry[table["realm"] == r]) for r in realms],
        crs=WGS84,
    )

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    gpkg = out.parent / f"{out.name}.gpkg"
    gpkg.unlink(missing_ok=True)
    table.to_file(gpkg, layer="parts", driver="GPKG")
    dissolved[["realm", "geometry", "dissolved"]].to_file(gpkg, layer="realms", driver="GPKG")

    if colours is not None:
        uncoloured = [r for r in realms if r not in colours]
        if uncoloured:
            raise ValueError(f"{len(uncoloured)} realms have no colour, such as {uncoloured[0]}")
    else:
        near = adjacency if adjacency is not None else _touching(table)
        colours = _separating_colours(realms, defaultdict(set, _realm_neighbours(realm_of, near)))

    fig, ax = plt.subplots(figsize=(14, 7.4), facecolor=SURFACE)
    globe = _projected(shapely.box(-180, -90, 180, 90))
    ax.fill(*_polygons(globe)[0].exterior.xy, facecolor=GLOBE, edgecolor=OUTLINE, linewidth=0.5, zorder=0)

    flat = gpd.GeoDataFrame(table.drop(columns="geometry"), geometry=projected, crs=EQUAL_EARTH)
    # the edge in the fill's colour hides the hairline gaps between neighbours that antialiasing leaves
    flat.plot(ax=ax, color=[colours[r] for r in flat["realm"]], edgecolor=[colours[r] for r in flat["realm"]], linewidth=0.3, aspect=None, zorder=1)
    if flat["imputed"].any():
        with plt.rc_context({"hatch.linewidth": 0.5}):
            flat[flat["imputed"]].plot(ax=ax, facecolor="none", edgecolor="#0b0b0bb0", hatch="//////", linewidth=0, aspect=None, zorder=2)

    outlines = gpd.GeoSeries([_projected(s) for s in dissolved.geometry], crs=EQUAL_EARTH)
    outlines.boundary.plot(ax=ax, color=SURFACE, linewidth=0.6, aspect=None, zorder=3)
    if sub_of:
        keyed = table[table["sub_realm"].notna()]
        # one plot call: geopandas draws the whole figure again after every call, which is minutes for hundreds of sub-realms
        subs = [_projected(_dissolve(group.geometry)) for _, group in keyed.groupby(["realm", "sub_realm"])]
        gpd.GeoSeries(subs, crs=EQUAL_EARTH).boundary.plot(ax=ax, color=TEXT, linewidth=0.9, aspect=None, zorder=4)

    for realm, shape in zip(realms, outlines, strict=True):
        biggest = max(shapely.get_parts(shape), key=lambda p: p.area)
        spot = polylabel(biggest, tolerance=5e4)
        ax.text(
            spot.x, spot.y, str(realm), ha="center", va="center", fontsize=10, fontweight="bold", color=TEXT, zorder=5,
            path_effects=[withStroke(linewidth=3, foreground=SURFACE)],
        )

    keys = [Patch(facecolor=GLOBE, edgecolor=OUTLINE, linewidth=0.3, label="not in the map")]
    if sub_of:
        keys.append(Line2D([], [], color=TEXT, linewidth=0.9, label="sub-realm outline"))
    if flat["imputed"].any():
        keys.append(Patch(facecolor="none", edgecolor="#0b0b0bb0", hatch="//////", label="imputed"))
    ax.legend(handles=keys, loc="lower left", frameon=False, fontsize=9, labelcolor=TEXT)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left", color=TEXT, fontsize=13)
    fig.savefig(out.parent / f"{out.name}.png", dpi=DPI, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
