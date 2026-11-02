"""Realm boundaries on the ISEA3H grid, and a gate for whether they follow the edges of reference regions more than chance.

A realm boundary runs between two neighbouring cells whose parts of one kind (land or sea) carry different realms; its length is that of
the edge the cells share, in geodesic km. The land and sea parts of one cell are never a boundary: the coast is not a realm's edge.

Cleaning layer 3 drops plant records outside their WCVP native range, which is given as TDWG level 3 regions, so it could draw realm
boundaries along TDWG edges by construction. `gate` asks whether a WCVP mode puts more of its land boundary near TDWG edges than the map
without layer 3 does, beyond what randomly rotated copies of the regions give. The test is paired: each rotated copy is measured against
both maps, so that whatever a rotation does to the share near the regions (rotated edges fall on the sea, where no land boundary is) cancels.

Cells are drawn as `grid.py` draws them (`_planar`) and are neighbours where `grid._neighbours` says so, so the edges here are those of the
adjacency the rest of the analysis uses. An edge is the geodesic between two vertices of the grid file, as `grid.py` measures areas.
"""

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
import polars as pl
import shapely
from shapely import MultiPolygon, Polygon, STRtree

from realms.grid import GEOD, SNAP, _neighbours, _planar, _shifted
from realms.hexmaps import _cut_at_antimeridian, _densified, _polygons

ROOT = Path(__file__).resolve().parents[2]  # the repository, so that the defaults below do not depend on the working directory
LAND = ROOT / "data" / "reference" / "natural_earth" / "ne_10m_land.zip"
CACHE = ROOT / "results" / "boundaries"
KINDS = ("land", "sea")
# Degrees. Region edges are straight in longitude and latitude; they are cut this short before rotating, so that they keep their shape.
ROTATION_STEP = 0.25
# Degrees. The region boundaries are sampled at segments this long to measure how much of them lies over land (about 5 km or less).
SAMPLE_STEP = 0.05
# Degrees from a pole within which a vertex is at the pole, where its longitude means nothing.
AT_POLE = 1e-9
# Longitude of a vertex on the line where the files cut a region at the antimeridian, which is no boundary.
ANTIMERIDIAN = 180 - 1e-6
EDGE_SCHEMA = {"cell_a": pl.String, "cell_b": pl.String, "km": pl.Float64}
GATE_SCHEMA = {"mode": pl.String, "share": pl.Float64, "raise": pl.Float64, "null_q95": pl.Float64, "passes": pl.Boolean}


# ---------------------------------------------------------------------------------------------------------------------------- cells


@lru_cache(maxsize=4)
def _cells_of(path: str, size: int, mtime_ns: int) -> tuple[tuple[str, ...], np.ndarray, STRtree]:
    frame = gpd.read_file(path)
    if frame["geogrid_id"].dtype.kind == "f":
        raise ValueError(f"geogrid_id was read as floating point, which loses digits: {path}")
    cells = tuple(str(i) for i in frame["geogrid_id"])
    if len(set(cells)) != len(cells):
        raise ValueError(f"geogrid_id is not unique in {path}")
    shapes = np.array([_planar(c, s, lon) for c, s, lon in zip(cells, frame.geometry, frame["lon_c"], strict=True)], dtype=object)
    return cells, shapes, STRtree(shapes)


def _cells(gpkg: Path) -> tuple[tuple[str, ...], np.ndarray, STRtree]:
    """The cells of a grid file, drawn as `load_grid` draws them, with a tree of them. Kept while the file is unchanged."""
    path = Path(gpkg).resolve()
    stat = path.stat()
    return _cells_of(str(path), stat.st_size, stat.st_mtime_ns)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _line_km(lines: np.ndarray) -> np.ndarray:
    """Geodesic km of each line (or collection of lines and points), each segment a geodesic between its two vertices.

    A segment along a pole, between two vertices at latitude 90 of different longitude, is one point and has no length, so the detour that
    `_planar` draws over a pole measures as the geodesic over it.
    """
    parts, owner = shapely.get_parts(lines, return_index=True)
    coords, index = shapely.get_coordinates(parts, return_index=True)
    result = np.zeros(len(lines))
    if len(coords) < 2:
        return result
    same = index[1:] == index[:-1]
    start = np.flatnonzero(same)
    _, _, metres = GEOD.inv(coords[start, 0], coords[start, 1], coords[start + 1, 0], coords[start + 1, 1])
    np.add.at(result, owner[index[start]], metres / 1e3)
    return result


def _measure_edges(cells: tuple[str, ...], shapes: np.ndarray) -> dict[tuple[str, str], float]:
    """The shared edge of every pair that `grid._neighbours` finds, as `edge_lengths` gives it.

    The edge two cells share is the boundary of one inside a band of `SNAP` degrees around the boundary of the other. Cells are drawn
    around their own centres, so a neighbour may be a whole turn of longitude away, and a cell over a pole shares one meridian with its
    neighbour at one turn and the other meridian at another: every turn is tried and the pieces add up. Each piece is found from the cell
    of lower index only, so none is counted twice.
    """
    boundary = shapely.boundary(shapes)
    band = shapely.buffer(boundary, SNAP)
    tree = STRtree(band)
    km = np.zeros(0)
    pairs = []
    for degrees in (-360.0, 0.0, 360.0):
        moved = _shifted(boundary, degrees)
        i, j = tree.query(moved, predicate="intersects")
        keep = i < j
        i, j = i[keep], j[keep]
        pairs.append(np.column_stack([i, j]))
        km = np.concatenate([km, _line_km(shapely.intersection(moved[i], band[j]))])
    found: dict[tuple[int, int], float] = {}
    for (a, b), length in zip(np.concatenate(pairs).tolist(), km.tolist(), strict=True):
        found[(a, b)] = found.get((a, b), 0.0) + length

    result = {}
    for a, b in sorted(_neighbours(shapes)):
        length = found.get((a, b), 0.0)
        if length <= 0:
            raise ValueError(f"cells {cells[a]} and {cells[b]} are neighbours but no shared edge was measured")
        first, second = sorted((cells[a], cells[b]))
        result[(first, second)] = length
    return result


def edge_lengths(gpkg: Path, *, cache_dir: Path = CACHE) -> dict[tuple[str, str], float]:
    """Geodesic km (WGS84) of the edge each pair of neighbouring cells shares, keyed `(cell_a, cell_b)` with `cell_a < cell_b` as text.

    The pairs are those of `grid.load_grid`'s adjacency. An edge is the geodesic between the file's vertices, measured exactly, so no
    densifying is needed; across the antimeridian and over a pole it is the edge on the globe. The result is cached as
    `<cache_dir>/edges-<file stem>.parquet`, with the sha256 of the grid file next to it, and reused while that agrees.
    """
    gpkg = Path(gpkg)
    out = Path(cache_dir) / f"edges-{gpkg.stem}.parquet"
    stamp = out.with_suffix(".json")
    digest = _sha256(gpkg)
    if out.exists() and stamp.exists() and json.loads(stamp.read_text()).get("sha256") == digest:
        table = pl.read_parquet(out)
        return {(a, b): km for a, b, km in table.iter_rows()}
    cells, shapes, _ = _cells(gpkg)
    edges = _measure_edges(cells, shapes)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = [(a, b, km) for (a, b), km in edges.items()]
    pl.DataFrame(rows, schema=EDGE_SCHEMA, orient="row").write_parquet(out)
    stamp.write_text(json.dumps({"gpkg": str(gpkg), "sha256": digest, "pairs": len(rows)}, indent=1) + "\n")
    return edges


def edge_adjacency(edges: Iterable[tuple[str, str]]) -> dict[str, set[str]]:
    """Cell to neighbouring cells, from the keys of `edge_lengths`. A cell with no neighbour is not a key."""
    adjacency: dict[str, set[str]] = {}
    for a, b in edges:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


# ------------------------------------------------------------------------------------------------------------------------ boundaries


def boundary(labels: Mapping[str, int], kind: str, edges: Mapping[tuple[str, str], float]) -> pl.DataFrame:
    """`cell_a`, `cell_b`, `km`: the neighbouring cells whose parts of `kind` both carry a realm in `labels`, a different one.

    `labels` maps part ids (`<cell>:land`, `<cell>:sea`) to realms, such as one k of `parts.parquet`. Only parts of `kind` count, so the
    two parts of one cell are never a boundary. Sorted by `cell_a`, `cell_b`.
    """
    if kind not in KINDS:
        raise ValueError(f"a kind is one of {KINDS}, not {kind!r}")
    suffix = f":{kind}"
    realm = {part[: -len(suffix)]: label for part, label in labels.items() if part.endswith(suffix)}
    rows = [(a, b, km) for (a, b), km in edges.items() if a in realm and b in realm and realm[a] != realm[b]]
    return pl.DataFrame(sorted(rows), schema=EDGE_SCHEMA, orient="row")


def edge_share(boundary: pl.DataFrame, near: Iterable[str]) -> float:
    """The share of the boundary's km with either cell in `near`, NaN for a map with no boundary."""
    total = boundary["km"].sum()
    if not total > 0:
        return math.nan
    near = list(set(near))
    close = boundary.filter(pl.col("cell_a").is_in(near) | pl.col("cell_b").is_in(near))["km"].sum()
    return float(close / total)


def _pieces(regions: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    """The polygons of the regions, and the row of the region each belongs to."""
    shapes = shapely.make_valid(np.asarray(regions.geometry.values, dtype=object))
    polygons, owner = [], []
    for row, shape in enumerate(shapes):
        found = _polygons(shape)
        polygons += found
        owner += [row] * len(found)
    return np.array(polygons, dtype=object), np.array(owner, dtype=np.int64)


def near_cells(regions: gpd.GeoDataFrame, gpkg: Path, adjacency: Mapping[str, Iterable[str]]) -> set[str]:
    """The cells within one cell of an edge between regions: those whose polygon intersects two or more regions, and their neighbours.

    `regions` has one row per region (EPSG:4326); a region may be several polygons, as one cut at the antimeridian is. Cells are drawn
    around their own centres, so each region is also tried a turn of longitude to either side. A cell that only touches a region along an
    edge intersects it. `adjacency` maps a cell to its neighbours, as `edge_adjacency` or `Grid.adjacency` give them.
    """
    cells, _, tree = _cells(gpkg)
    polygons, owner = _pieces(regions)
    if not len(polygons):
        return set()
    copies = np.concatenate([_shifted(polygons, degrees) for degrees in (-360.0, 0.0, 360.0)])
    owners = np.tile(owner, 3)
    piece, cell = tree.query(copies, predicate="intersects")
    distinct = np.unique(np.column_stack([cell, owners[piece]]), axis=0)
    counts = np.bincount(distinct[:, 0], minlength=len(cells))
    straddle = {cells[i] for i in np.flatnonzero(counts >= 2)}
    near = set(straddle)
    for c in straddle:
        near.update(adjacency.get(c, ()))
    return near


# ------------------------------------------------------------------------------------------------------------------------- rotations


def _xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon, lat = np.radians(lon), np.radians(lat)
    return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def _lonlat(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])), np.degrees(np.arcsin(np.clip(xyz[:, 2], -1.0, 1.0)))


def _matrix(shift_deg: float, axis_lon_deg: float, tilt_deg: float) -> np.ndarray:
    """A turn of `shift_deg` about the polar axis, then a tilt of `tilt_deg` about the horizontal axis through the equator at `axis_lon_deg`."""
    s = math.radians(shift_deg)
    turn = np.array([[math.cos(s), -math.sin(s), 0.0], [math.sin(s), math.cos(s), 0.0], [0.0, 0.0, 1.0]])
    a, t = math.radians(axis_lon_deg), math.radians(tilt_deg)
    u = np.array([math.cos(a), math.sin(a), 0.0])
    cross = np.array([[0.0, -u[2], u[1]], [u[2], 0.0, -u[0]], [-u[1], u[0], 0.0]])
    tilt = math.cos(t) * np.eye(3) + math.sin(t) * cross + (1 - math.cos(t)) * np.outer(u, u)
    return tilt @ turn


def _rotate_points(matrix: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return _lonlat(_xyz(np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)) @ matrix.T)


def _inside_point(ring: np.ndarray) -> tuple[float, float]:
    """A point inside the polygon a ring of the region file bounds, as drawn in longitude and latitude (for Antarctica, down to the pole)."""
    shape = shapely.make_valid(Polygon(ring))
    pieces = _polygons(shape)
    if not pieces:
        raise ValueError("a ring of the regions bounds no area")
    point = max(pieces, key=lambda p: p.area).representative_point()
    return point.x, point.y


def _without_spikes(xyz: np.ndarray) -> np.ndarray:
    """The vertices of an open ring on the unit sphere without repeats and without spikes, where the ring runs out and back the same way.

    The region file cuts a region over a pole, such as Antarctica, at the antimeridian and draws it down to the pole: on the globe its
    ring runs down the antimeridian to the pole and back up the same line. Rotated, that is a line out of the region and back, which would
    not bound anything.
    """
    def same(p, q) -> bool:
        return bool(np.abs(p - q).max() < 1e-10)

    n = len(xyz)
    if n >= 3:
        step = np.abs(np.diff(np.vstack([xyz, xyz[:2]]), axis=0)).max(axis=1)
        back = np.abs(np.vstack([xyz[2:], xyz[:2]]) - xyz).max(axis=1)
        if step.min() >= 1e-10 and back.min() >= 1e-10:
            return xyz
    kept: list[np.ndarray] = []
    for p in xyz:
        if kept and same(p, kept[-1]):
            continue
        if len(kept) >= 2 and same(p, kept[-2]):
            kept.pop()
            continue
        kept.append(p)
    changed = True
    while changed and len(kept) >= 3:
        changed = False
        if same(kept[0], kept[-1]):
            kept.pop()
            changed = True
        elif same(kept[-1], kept[1]):  # out to kept[0] and back
            del kept[:2]
            changed = True
        elif same(kept[-2], kept[0]):  # out to kept[-1] and back
            kept.pop()
            del kept[0]
            changed = True
    return np.array(kept) if len(kept) >= 3 else np.zeros((0, 3))


def _rotated_ring(matrix: np.ndarray, ring: np.ndarray) -> list[Polygon]:
    """What a ring of the region file bounds, rotated, as polygons within -180 and 180 degrees of longitude.

    The ring is turned on the unit sphere and its longitudes unwrapped, so that no step between vertices exceeds 180 degrees. A ring that
    then winds once round the globe encloses a pole, and is closed along it; the pole is the one on the side of a point inside the ring.
    The result is cut at the antimeridian.
    """
    xyz = _without_spikes(_xyz(ring[:-1, 0], ring[:-1, 1]) @ matrix.T)
    if len(xyz) < 3:
        return []
    lon, lat = _lonlat(xyz)
    at_pole = np.flatnonzero(90 - np.abs(lat) < AT_POLE)
    if len(at_pole):
        # a vertex at a pole takes the longitude of the vertex before it, and is repeated with that of the one after it
        xs, ys, n = [], [], len(lon)
        for k in range(n):
            if k in at_pole:
                pole = math.copysign(90.0, lat[k])
                xs += [lon[(k - 1) % n], lon[(k + 1) % n]]
                ys += [pole, pole]
            else:
                xs.append(lon[k])
                ys.append(lat[k])
        lon, lat = np.array(xs), np.array(ys)
    lon, lat = np.append(lon, lon[0]), np.append(lat, lat[0])
    steps = np.diff(lon)
    steps = (steps + 180.0) % 360.0 - 180.0
    unwrapped = np.concatenate([[lon[0]], lon[0] + np.cumsum(steps)])
    winds = round((unwrapped[-1] - unwrapped[0]) / 360.0)
    unwrapped -= 360.0 * math.floor((unwrapped.min() + 180.0) / 360.0)
    coords = np.column_stack([unwrapped, lat])
    if winds == 0:
        shape = Polygon(coords)
    elif abs(winds) == 1:
        px, py = _inside_point(ring)
        px, py = _rotate_points(matrix, [px], [py])
        candidates = {}
        for pole in (90.0, -90.0):
            closed = np.vstack([coords, [[coords[-1, 0], pole], [coords[0, 0], pole]]])
            candidates[pole] = shapely.make_valid(Polygon(closed))
        tried = shapely.points(px[0] + 360.0 * np.arange(-2, 3), np.full(5, py[0]))
        north = shapely.intersects(candidates[90.0], tried).any()
        shape = candidates[90.0 if north else -90.0]
    else:
        raise ValueError(f"a rotated ring winds {winds} times round the globe")
    if not shape.is_valid:
        shape = shapely.make_valid(shape)
    return _polygons(_cut_at_antimeridian(MultiPolygon(_polygons(shape))))


@dataclass(frozen=True)
class Rotation:
    """A rotation of the sphere: a turn of `shift_deg` about the polar axis, then a tilt of `tilt_deg` about a horizontal axis.

    The tilt axis runs through the equator at `axis_lon_deg`. `matrix` is the 3x3 matrix, rows of a tuple, applied to unit vectors
    (x towards 0 E on the equator, z the north pole). `land_km` is the region boundary over land after rotating, where `rotations` measured it.
    """

    shift_deg: float
    axis_lon_deg: float
    tilt_deg: float
    matrix: tuple[tuple[float, float, float], ...]
    land_km: float = math.nan

    @classmethod
    def of(cls, shift_deg: float, axis_lon_deg: float, tilt_deg: float) -> "Rotation":
        m = _matrix(shift_deg, axis_lon_deg, tilt_deg)
        return cls(float(shift_deg), float(axis_lon_deg), float(tilt_deg), tuple(tuple(float(v) for v in row) for row in m))

    @property
    def array(self) -> np.ndarray:
        return np.array(self.matrix)

    def points(self, lon, lat) -> tuple[np.ndarray, np.ndarray]:
        """Rotated longitudes and latitudes, in degrees."""
        return _rotate_points(self.array, lon, lat)

    def geometry(self, shape) -> MultiPolygon:
        """A polygon or multipolygon of longitudes and latitudes, rotated, as valid polygons cut at the antimeridian.

        Edges are taken as straight in longitude and latitude, as the region file draws them, and cut to `ROTATION_STEP` degrees before
        rotating. A hole is rotated as a ring of its own and taken out.
        """
        matrix = self.array
        outer, holes = [], []
        for polygon in _polygons(shapely.make_valid(shape)):
            for dense in _densified(polygon, ROTATION_STEP).geoms:
                outer += _rotated_ring(matrix, np.asarray(dense.exterior.coords))
                for hole in dense.interiors:
                    holes += _rotated_ring(matrix, np.asarray(hole.coords))
        if not outer:
            return MultiPolygon()
        result = shapely.union_all(outer)
        if holes:
            result = shapely.difference(result, shapely.union_all(holes))
        if not result.is_valid:
            result = shapely.make_valid(result)
        return MultiPolygon(_polygons(result))

    def apply(self, regions: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """The regions, one row each as before, with their geometry rotated (`geometry`)."""
        rotated = [self.geometry(shape) for shape in regions.geometry]
        return regions.set_geometry(gpd.GeoSeries(rotated, index=regions.index, crs=regions.crs))


@lru_cache(maxsize=2)
def _land_from_file(path: str):
    land = shapely.union_all(shapely.make_valid(gpd.read_file(path).geometry.values))
    shapely.prepare(land)
    return land


def _land(land) -> object:
    if isinstance(land, (str, Path)):
        return _land_from_file(str(Path(land).resolve()))
    shapely.prepare(land)
    return land


def _boundary_samples(regions: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Midpoints (longitude, latitude) and geodesic km of short pieces of the union of the regions' boundaries.

    Edges shared by two regions are in the union once. The cut that the file makes at the antimeridian is no boundary and is left out.
    """
    lines = shapely.union_all(shapely.boundary(shapely.make_valid(np.asarray(regions.geometry.values, dtype=object))), grid_size=SNAP)
    lines = [line for line in shapely.get_parts(lines) if line.geom_type == "LineString" and len(line.coords) > 1]
    xs, ys, km = [], [], []
    for line in lines:
        xy = np.asarray(line.coords)
        pieces = np.maximum(1, np.ceil(np.abs(np.diff(xy, axis=0)).max(axis=1) / SAMPLE_STEP)).astype(int)
        start = np.repeat(np.arange(len(pieces)), pieces)
        along = (np.arange(len(start)) - np.repeat(np.cumsum(pieces) - pieces, pieces)) / pieces[start]
        dense = np.vstack([xy[start] + along[:, None] * (xy[start + 1] - xy[start]), xy[-1:]])
        a, b = dense[:-1], dense[1:]
        cut = (np.abs(a[:, 0]) > ANTIMERIDIAN) & (np.abs(b[:, 0]) > ANTIMERIDIAN)
        a, b = a[~cut], b[~cut]
        _, _, metres = GEOD.inv(a[:, 0], a[:, 1], b[:, 0], b[:, 1])
        middle = _xyz(a[:, 0], a[:, 1]) + _xyz(b[:, 0], b[:, 1])
        mid = _lonlat(middle / np.linalg.norm(middle, axis=1, keepdims=True))
        xs.append(mid[0])
        ys.append(mid[1])
        km.append(metres / 1e3)
    if not km:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(km)


def _over_land(samples: tuple[np.ndarray, np.ndarray, np.ndarray], land, rotation: Rotation | None) -> float:
    lon, lat, km = samples
    if rotation is not None:
        lon, lat = rotation.points(lon, lat)
    return float(km[shapely.contains_xy(land, lon, lat)].sum())


def boundary_km_over_land(regions: gpd.GeoDataFrame, rotation: Rotation | None = None, *, land=LAND) -> float:
    """Geodesic km of the union of the regions' boundaries that lies over land, after `rotation` if one is given.

    The boundaries are cut into pieces of at most `SAMPLE_STEP` degrees, and a piece counts as over land when its midpoint is. `land` is a
    file geopandas reads or a geometry, in longitude and latitude.
    """
    return _over_land(_boundary_samples(regions), _land(land), rotation)


def rotations(
    regions: gpd.GeoDataFrame,
    n: int,
    seed: int,
    max_tilt_deg: float = 30.0,
    land_tolerance: float | None = 0.5,
    *,
    land=LAND,
    max_attempts: int | None = None,
) -> list[Rotation]:
    """`n` random rotations of the sphere that keep the regions' boundary over land within `land_tolerance` of what it is unrotated.

    A draw is a longitude shift uniform in [0, 360), then a tilt about a horizontal axis through the equator at a uniform longitude, by an
    angle that puts the rotated pole uniformly on the cap within `max_tilt_deg` of where it was. It is kept when
    `|rotated - unrotated| <= land_tolerance * unrotated`, with the lengths of `boundary_km_over_land`. Draws stop when `n` are kept, or
    after `max_attempts` (by default 50 * n), when too few were kept and that is an error. The same seed gives the same rotations.

    The starting tolerance is 0.5: the regions are land, and rotated their boundary over land is 0.31 to 0.83 of what it is unrotated
    (median 0.49), so 0.5 keeps the more land-like half of the draws. A rotation with little land shrinks the null raises and makes the
    gate stricter, so loosening the tolerance does not make a mode pass more easily. `land_tolerance=None` keeps every draw, leaving the
    land bias to the pairing alone; the land is then not read, and `land_km` is NaN.
    """
    if n < 1:
        raise ValueError(f"n is at least 1, not {n}")
    rng = np.random.default_rng(seed)
    lowest = math.cos(math.radians(max_tilt_deg))

    def draw() -> Rotation:
        shift, axis, cos_tilt = rng.uniform(0.0, 360.0), rng.uniform(0.0, 360.0), rng.uniform(lowest, 1.0)
        return Rotation.of(shift, axis, math.degrees(math.acos(min(1.0, cos_tilt))))

    if land_tolerance is None:
        return [draw() for _ in range(n)]
    samples = _boundary_samples(regions)
    where = _land(land)
    unrotated = _over_land(samples, where, None)
    if not unrotated > 0:
        raise ValueError("no boundary of the regions lies over land")
    attempts = max_attempts if max_attempts is not None else 50 * n
    kept: list[Rotation] = []
    for _ in range(attempts):
        rotation = draw()
        km = _over_land(samples, where, rotation)
        if abs(km - unrotated) <= land_tolerance * unrotated:
            kept.append(replace(rotation, land_km=km))
            if len(kept) == n:
                return kept
    raise ValueError(
        f"only {len(kept)} of {attempts} rotations kept the boundary over land within {land_tolerance:.0%} of {unrotated:,.0f} km, "
        f"where {n} were wanted (max_tilt_deg={max_tilt_deg})"
    )


# ------------------------------------------------------------------------------------------------------------------------------ gate


def gate(
    maps: Mapping[str, Mapping[str, int]],
    base: str,
    regions: gpd.GeoDataFrame,
    rotations: Iterable[Rotation],
    *,
    edges: Mapping[tuple[str, str], float],
    gpkg: Path,
    adjacency: Mapping[str, Iterable[str]],
) -> pl.DataFrame:
    """`mode`, `share`, `raise`, `null_q95`, `passes`: whether each map follows the region edges more than chance, against `maps[base]`.

    `maps` is mode to labels (part id to realm, as `boundary` takes them); `base` names the map without layer 3. For a map m and a rotation
    r, `s_m(r)` is the share of m's land boundary near the regions rotated by r (`near_cells`, `edge_share`), and r = 0 leaves them as they
    are. `share` is `s_m(0)`, `raise` is `s_m(0) - s_base(0)`, `null_q95` is the 95th percentile (linear interpolation) of the null
    raises `s_m(r) - s_base(r)`, and a mode `passes` when its raise is at most that. Only land boundaries count, as the regions are of
    land. One row per mode other than `base`, in the order of `maps`.
    """
    if base not in maps:
        raise ValueError(f"the base mode {base!r} is not among the maps {sorted(maps)}")
    rotations = list(rotations)
    if not rotations:
        raise ValueError("the gate needs at least one rotation")
    lines = {mode: boundary(labels, "land", edges) for mode, labels in maps.items()}
    others = [mode for mode in maps if mode != base]

    def shares(near: set[str]) -> dict[str, float]:
        return {mode: edge_share(line, near) for mode, line in lines.items()}

    observed = shares(near_cells(regions, gpkg, adjacency))
    null: dict[str, list[float]] = {mode: [] for mode in others}
    for rotation in rotations:
        s = shares(near_cells(rotation.apply(regions), gpkg, adjacency))
        for mode in others:
            null[mode].append(s[mode] - s[base])

    rows = []
    for mode in others:
        lift = observed[mode] - observed[base]
        q95 = float(np.percentile(null[mode], 95))
        rows.append((mode, observed[mode], lift, q95, bool(lift <= q95)))
    return pl.DataFrame(rows, schema=GATE_SCHEMA, orient="row")
