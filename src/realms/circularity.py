"""Comparing the maps of cleaning runs part by part, and where their boundaries move, for the circularity gate (`scripts/circularity.py`).

Runs that clean differently have different network units, so maps are compared by part (`<cell>:land`, `<cell>:sea`) through `parts.parquet`,
never by unit. A boundary that moves between two runs is an edge in the boundary of one and not of the other (`boundaries.boundary`), and a
moved land edge is put in the TDWG level 1 region that holds the centre of its first cell, else of its second, else in `NONE`.
"""

from collections import Counter
from collections.abc import Iterable, Mapping

import geopandas as gpd
import numpy as np
import polars as pl
import shapely
from shapely import STRtree

from realms.boundaries import EDGE_SCHEMA
from realms.score import v_measure
from realms.natives import TDWG_LEVEL3_CODE
from realms.units import TDWG_LEVEL1, TDWG_LEVEL1_CODE

NONE = "none"  # a moved edge with neither cell's centre in a region
WCVP_ORDER = ("native_only", "introduced")  # the WCVP modes the headline may take, strongest first; `off` when none passes
BASE = "off"
REGION_SCHEMA = {"region": pl.String, "km": pl.Float64, "edges": pl.Int64}


def labels_at(parts: pl.DataFrame, k: int) -> dict[str, int]:
    """Part id to realm at one k of a `parts.parquet` table (imputed parts included)."""
    rows = parts.filter(pl.col("k") == k)
    if not rows.height:
        raise ValueError(f"the parts hold no map at k={k}; they have k={sorted(parts['k'].unique())}")
    return dict(rows.select("part", "realm").iter_rows())


def compare(a: Mapping[str, int], b: Mapping[str, int]) -> dict:
    """`v_measure` over the parts both maps hold, `shared` (how many), and `only_a` and `only_b`, the parts that one map holds and the other not."""
    shared = a.keys() & b.keys()
    return {"v_measure": v_measure(dict(a), dict(b)), "shared": len(shared), "only_a": len(a.keys() - shared), "only_b": len(b.keys() - shared)}


def common_k(candidates: Mapping[str, Iterable[int]], selected: Mapping[str, int]) -> int:
    """The k at which every run has a map that the most runs selected; a tie, or no selected k that all have, goes to the coarsest such k."""
    shared = set.intersection(*(set(ks) for ks in candidates.values()))
    if not shared:
        raise ValueError(f"no k is a candidate of every run: {({run: sorted(ks) for run, ks in candidates.items()})}")
    votes = Counter(k for k in selected.values() if k in shared)
    return min(shared, key=lambda k: (-votes[k], k))


def moved(a: pl.DataFrame, b: pl.DataFrame) -> pl.DataFrame:
    """The edges in exactly one of two boundaries (`boundaries.boundary` of one kind), with their km, sorted by `cell_a`, `cell_b`.

    An edge that is a boundary in both does not move, even where the realms on either side were renamed or changed. A boundary needs both
    of its parts mapped, so a part that only one run maps can put an edge in that run's boundary alone: such parts move edges too.
    """
    key = ["cell_a", "cell_b"]
    only_a = a.join(b.select(key), on=key, how="anti")
    only_b = b.join(a.select(key), on=key, how="anti")
    return pl.concat([only_a, only_b]).cast(EDGE_SCHEMA).sort(key)


def moves_by_region(edges: pl.DataFrame, region_of: Mapping[str, str | None]) -> pl.DataFrame:
    """`region`, `km`, `edges`: moved edges by the region of the centre of their first cell (`cell_a`), else their second, else `NONE`.

    `region_of` is cell to region, as `centre_regions` gives it; a cell missing from it is in no region. Sorted by km, largest first.
    """
    regions = [region_of.get(a) or region_of.get(b) or NONE for a, b in zip(edges["cell_a"], edges["cell_b"], strict=True)]
    return (
        edges.with_columns(pl.Series("region", regions, dtype=pl.String))
        .group_by("region")
        .agg(pl.col("km").sum(), pl.len().cast(pl.Int64).alias("edges"))
        .cast(REGION_SCHEMA)
        .select(list(REGION_SCHEMA))
        .sort(["km", "region"], descending=[True, False])
    )


def centre_regions(centres: pl.DataFrame, regions: gpd.GeoDataFrame) -> dict[str, str | None]:
    """Cell to the name of the TDWG level 1 region whose level 3 polygons hold its centre (`cell`, `lon`, `lat`), or None if none does.

    Unlike `units.part_regions`, a centre outside every region (on the sea, or on a coast WGSRPD draws elsewhere) takes no nearest region.
    A centre on an edge goes to the region of the lowest level 3 code; each centre is tried a turn of longitude to either side as well.
    """
    ordered = regions.sort_values(TDWG_LEVEL3_CODE, kind="stable")
    shapes = shapely.make_valid(np.asarray(ordered.geometry.values, dtype=object))
    names = [TDWG_LEVEL1.get(str(c), str(c)) for c in ordered[TDWG_LEVEL1_CODE]]
    lon, lat = centres["lon"].to_numpy(), centres["lat"].to_numpy()
    tree = STRtree(shapes)
    found = np.full(len(lon), len(shapes), dtype=np.int64)
    for degrees in (-360.0, 0.0, 360.0):
        point, shape = tree.query(shapely.points(lon + degrees, lat), predicate="intersects")
        np.minimum.at(found, point, shape)
    return {cell: names[i] if i < len(shapes) else None for cell, i in zip(centres["cell"], found.tolist(), strict=True)}


def headline(passes: Mapping[str, bool]) -> str:
    """The headline WCVP mode: the first of `WCVP_ORDER` whose gate passes (at the default land tolerance), else `BASE`, which is WCVP off.

    `passes` is mode to the gate's verdict; other modes in it (the layer variants) are information and are never chosen.
    """
    missing = [mode for mode in WCVP_ORDER if mode not in passes]
    if missing:
        raise ValueError(f"the gate has no verdict for {missing}")
    return next((mode for mode in WCVP_ORDER if passes[mode]), BASE)
