"""The ISEA3H grid, with each cell split into the part on land and the part in the sea.

Cells come from GBIF's grid files, whose `geogrid_id` is what `GBIF_ISEA3HCode` returns for a point. Ids are long integers that a
float can't hold, so they are text from the moment they are read. A part is named `<cell>:land` or `<cell>:sea`: the coast is
a natural boundary for species, so a coastal cell's two halves are separate nodes, joined to each other and to the part of the same
kind in each neighbouring cell.

Areas are geodesic, because an equal-area projection gets the cells at the poles wrong.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import numpy as np
import polars as pl
import pyproj
import shapely
from shapely import Polygon, STRtree

GEOD = pyproj.Geod(ellps="WGS84")
KINDS = ("land", "sea")
# A part with no area still gets float noise: a cell wholly on land at resolution 6 has a sea_km2 of 3e-8. The smallest real parts on the
# grids are 0.006 km2 of land and 0.26 km2 of sea, so anything below 1000 m2 is none.
ZERO_KM2 = 1e-3
# Degrees, about 10 cm. Cells disagree by up to 1.4e-8 where they should share a vertex, and no two vertices lie between that and 1 apart.
SNAP = 1e-6
SHARED = 1e-4  # degrees of boundary two cells must share to be neighbours, not just touch at a corner
# An edge spanning this much longitude runs over a pole. Longitudes of polar vertices are only good to about 2e-6 degrees.
POLE_EDGE = 180 - 1e-4


@dataclass(frozen=True)
class Grid:
    """The cells of one resolution. `cells` has one row per cell, in the order of the file."""

    resolution: int
    cells: pl.DataFrame  # cell (the geogrid_id, as text), lat, lon (the centre), area_km2, land_km2, sea_km2
    adjacency: dict[str, set[str]]  # cells sharing a boundary segment, every cell a key


def _km2(geometry) -> float:
    """Geodesic area on WGS84. Lines and points left over from an intersection have none."""
    polygons = (p for p in shapely.get_parts(geometry) if p.geom_type == "Polygon")
    return sum(abs(GEOD.geometry_area_perimeter(p)[0]) for p in polygons) / 1e6


def _planar(cell: str, shape, lon_c: float) -> Polygon:
    """One simple polygon for the cell, drawn in the 360 degrees of longitude around its centre.

    The files cut a cell that crosses the antimeridian into pieces at -180 and 180, and draw some edges the long way round it. Near
    the poles an edge between two vertices of equal latitude and opposite longitude is drawn along the parallel, but runs over the
    pole. Measuring every longitude from the centre joins the pieces and makes each edge take the short way, and an edge that then
    spans 180 degrees is replaced by the detour along the meridians to the pole and back, so that the cells either side of a pole
    share a boundary and the polygon covers what the cell does.
    """
    pieces = []
    for piece in shapely.get_parts(shape):
        if len(piece.interiors):
            raise ValueError(f"cell {cell} has a hole")
        xy = np.asarray(piece.exterior.coords)
        # whole turns only, so that two pieces put the cut at exactly the same longitude
        pieces.append(Polygon(np.column_stack([xy[:, 0] + 360.0 * np.round((lon_c - xy[:, 0]) / 360.0), xy[:, 1]])))
    merged = pieces[0]
    if len(pieces) > 1:
        # the cut left vertices on edges, which would hide a pole edge that was cut in two
        merged = shapely.simplify(shapely.union_all(pieces), SNAP)
    if merged.geom_type != "Polygon":
        raise ValueError(f"cell {cell} does not make one polygon but a {merged.geom_type}")

    xy = np.asarray(merged.exterior.coords)
    ring = []
    for (x, y), (next_x, _) in zip(xy[:-1], xy[1:], strict=True):
        ring.append((x, y))
        if abs(next_x - x) > POLE_EDGE:
            pole = 90.0 if y > 0 else -90.0
            ring += [(x, pole), (next_x, pole)]
    return Polygon(ring)


def _shifted(geometries: np.ndarray, degrees: float) -> np.ndarray:
    return shapely.transform(geometries, lambda xy: xy + [degrees, 0.0]) if degrees else geometries


def _neighbours(shapes: np.ndarray) -> set[tuple[int, int]]:
    """Index pairs of cells that share at least `SHARED` degrees of boundary, up to `SNAP`, with 180 and -180 the same line.

    Cells are drawn around their own centres, so a neighbour may be a whole turn of longitude away: each is tried 360 degrees to
    either side as well.
    """
    boundary = shapely.boundary(shapes)
    band = shapely.buffer(boundary, SNAP)
    tree = STRtree(band)
    pairs = set()
    for degrees in (-360.0, 0.0, 360.0):
        moved = _shifted(boundary, degrees)
        i, j = tree.query(moved, predicate="intersects")
        shared = shapely.length(shapely.intersection(moved[i], band[j]))
        keep = (shared > SHARED) & (i != j)
        pairs.update((min(a, b), max(a, b)) for a, b in zip(i[keep].tolist(), j[keep].tolist(), strict=True))
    return pairs


def _land_km2(shapes: np.ndarray, land: Path) -> np.ndarray:
    """Geodesic area of each cell's overlap with the land.

    The land is drawn between -180 and 180, so it is also tried a turn of longitude to either side. The continents are single
    polygons of tens of thousands of vertices, so a cell is only tested against polygons that meet it, and then only against the
    part of each inside the cell's bounding box.
    """
    polygons = shapely.get_parts(shapely.union_all(gpd.read_file(land).geometry.values))
    copies = np.concatenate([_shifted(polygons, degrees) for degrees in (-360.0, 0.0, 360.0)])
    shapely.prepare(copies)
    i, j = STRtree(copies).query(shapes)
    meets = shapely.intersects(copies[j], shapes[i])
    i, j = i[meets], j[meets]
    inside = [shapely.clip_by_rect(land_piece, *box) for land_piece, box in zip(copies[j], shapely.bounds(shapes[i]), strict=True)]
    overlap = shapely.intersection(shapes[i], inside)
    result = np.zeros(len(shapes))
    for cell, piece in zip(i, overlap, strict=True):
        result[cell] += _km2(piece)
    return result


def load_grid(gpkg: Path, land_zip: Path) -> Grid:
    """Read a GBIF ISEA3H file and measure each cell against Natural Earth land, from `land_zip` or any file geopandas reads.

    `land_km2` is the geodesic area of the cell's overlap with the land, kept within the cell's own area, and `sea_km2` the rest.
    """
    frame = gpd.read_file(gpkg)
    ids = frame["geogrid_id"]
    if ids.dtype.kind == "f":
        raise ValueError(f"geogrid_id was read as floating point, which loses digits: {gpkg}")
    cells = [str(i) for i in ids]
    if len(set(cells)) != len(cells):
        raise ValueError(f"geogrid_id is not unique in {gpkg}")
    resolutions = frame["res"].unique()
    if len(resolutions) != 1:
        raise ValueError(f"{gpkg} holds several resolutions: {sorted(resolutions)}")

    drawn = zip(cells, frame.geometry, frame["lon_c"], strict=True)
    shapes = np.array([_planar(cell, shape, lon) for cell, shape, lon in drawn], dtype=object)
    area = np.array([_km2(shape) for shape in shapes])
    land = np.clip(_land_km2(shapes, land_zip), 0, area)

    adjacency: dict[str, set[str]] = {cell: set() for cell in cells}
    for a, b in _neighbours(shapes):
        adjacency[cells[a]].add(cells[b])
        adjacency[cells[b]].add(cells[a])
    table = pl.DataFrame(
        {
            "cell": cells,
            "lat": frame["lat_c"].to_numpy(),
            "lon": frame["lon_c"].to_numpy(),
            "area_km2": area,
            "land_km2": land,
            "sea_km2": area - land,
        }
    )
    return Grid(int(resolutions[0]), table, adjacency)


def part_id(cell: str, part: str) -> str:
    if part not in KINDS:
        raise ValueError(f"a part is one of {KINDS}, not {part!r}")
    return f"{cell}:{part}"


def cell_of(part: pl.Expr) -> pl.Expr:
    """The cell of a column of part ids (`<cell>:land` or `<cell>:sea`, see `part_id`), without checking them against a grid."""
    return part.str.replace(r":[^:]*$", "")


def _split(grid: Grid, present: Iterable[str]) -> list[tuple[str, str, str]]:
    """The distinct part ids as (part, cell, kind), sorted, each checked against the grid."""
    result = []
    for part in sorted(set(present)):
        cell, colon, kind = part.rpartition(":")
        if not colon:
            raise ValueError(f"{part!r} is not a part id, which is <cell>:land or <cell>:sea")
        part_id(cell, kind)
        if cell not in grid.adjacency:
            raise ValueError(f"cell {cell!r} of part {part!r} is not in the grid")
        result.append((part, cell, kind))
    return result


def parts(grid: Grid, present: Iterable[str]) -> pl.DataFrame:
    """The parts that hold records, from their ids: `part`, `cell`, `kind` ("land" or "sea") and the part's `area_km2`."""
    area = {cell: {"land": land, "sea": sea} for cell, land, sea in grid.cells.select("cell", "land_km2", "sea_km2").iter_rows()}
    rows = [(part, cell, kind, area[cell][kind]) for part, cell, kind in _split(grid, present)]
    return pl.DataFrame(rows, schema={"part": pl.String, "cell": pl.String, "kind": pl.String, "area_km2": pl.Float64}, orient="row")


def part_adjacency(grid: Grid, present: Iterable[str]) -> dict[str, set[str]]:
    """For each part that holds records, the parts that do too and are next to it.

    That is the same kind in a neighbouring cell, and the other kind in its own cell.
    """
    names = {part for part, _, _ in _split(grid, present)}
    near: dict[str, set[str]] = {part: set() for part in names}
    for part in names:
        cell, _, kind = part.rpartition(":")
        candidates = [part_id(cell, "sea" if kind == "land" else "land")] + [part_id(other, kind) for other in grid.adjacency[cell]]
        near[part].update(c for c in candidates if c in names)
    return near
