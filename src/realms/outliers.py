"""Layer 5 of the cleaning: the rare-outlier heuristics. A single wrong record is weak evidence on its own, so these judge it against the species'
own distribution and against how heavily its cell is sampled. They remove what no flag, list or coast could catch: a zoo animal in a city, a
vagrant, a misidentification, a wrong coordinate.

**Grain.** The rules judge cells, not parts. A species' records in a cell are summed over its parts and countries, and a rule that flags a
(species, cell) takes all the species' rows in that cell. Distances are great-circle arcs between cell centres, in kilometres, and the scale
of the grid is measured from the grid (`neighbour_spacing`), never taken from a resolution.

**Core.** A species' core cells are its top cells by records, up to and including the cell that takes the running total to `core_share` of its
records, and all its cells under `heuristic_min_records` records. Cells with equal records are ordered by cell id (as text), so among tied cells at
the line which are core is reproducible but arbitrary. A core cell is never removed by `detached_satellites` or `under_reported`.
`thin_long_range` is the exception: it takes a whole species out.

**The rules** (the thresholds are the `Settings` of `realms.clean`):

- `detached_satellites`: for a species with at least `heuristic_min_records` records, the occupied cells are joined into patches when their
  centres are at most `satellite_gap_km` plus the median neighbour spacing apart. A patch with under `satellite_record_share` of the species'
  records and under `satellite_cell_share` of its cells goes, unless it has at least `satellite_min_datasets` datasets and `satellite_min_years`
  years, and never if it holds a core cell. A patch's datasets and years are the largest over its cells, which is a lower bound (the distinct
  counts of a row do not add up).
- `under_reported`: the species' share of its group's records in a cell is its reporting rate. A cell with at least `group_min_records` records
  of the group, where the species has at most `under_max_records` and at least `heuristic_min_records` in all, goes if its rate is under
  `under_rate_factor` times the species' median rate over its core cells. A species with no group is not judged.
- `thin_long_range`: a species with under `thin_max_records` records whose cells are over `thin_km` apart goes whole.

Every rule judges the rows layers 1 to 4 left, and not one another's removals, so what a rule flags does not depend on the order of the rules.
A cell that several flag goes with the first of `RULE_ORDER`, and `Settings.rules_off` leaves a rule out.

`suspect_cells` only lists: cells where a far higher share of species were removed than in the cells next to them.
"""

from typing import TYPE_CHECKING

import numpy as np
import polars as pl
from scipy import sparse
from scipy.sparse import csgraph
from scipy.spatial import cKDTree

from realms.grid import Grid, cell_of

if TYPE_CHECKING:
    from realms.clean import Cleaned, Settings

EARTH_KM = 6371.0088  # the mean radius of the earth
EPS = 1e-9  # a share compared with a count is exact to this: 0.05 * 60 is 3.0000000000000004, and 3 cells of 60 are not under 5%
EXPANSIONS = 20_000_000  # neighbour pairs that the patch search holds at a time
RULE_ORDER = ("thin_long_range", "detached_satellite", "under_reported")
SUSPECT_RATIO = 3.0
SUSPECT_POINTS = 0.10
NEEDED = ("species_id", "part", "records", "n_datasets", "n_years")
FLAGGED_SCHEMA = {"species_id": pl.String, "cell": pl.String, "rule": pl.String}
SUSPECT_SCHEMA = {
    "cell": pl.String,
    "lat": pl.Float64,
    "lon": pl.Float64,
    "species": pl.Int64,
    "removed": pl.Int64,
    "share": pl.Float64,
    "neighbours": pl.Int64,
    "neighbour_median": pl.Float64,
}


# ---- geometry

def _vectors(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Unit vectors of points on the sphere, so that nearness does not depend on the antimeridian or the poles."""
    phi, lam = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])


def _arc_km(lat1, lon1, lat2, lon2) -> np.ndarray:
    """Great-circle distance in kilometres (haversine), between arrays of points."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * EARTH_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def neighbour_spacing(grid: Grid) -> float:
    """The median distance in kilometres between the centres of cells that share a boundary: the scale of the grid, 0 where no cell has a neighbour."""
    index = {cell: i for i, cell in enumerate(grid.cells["cell"])}
    pairs = [(index[cell], index[other]) for cell, near in grid.adjacency.items() for other in near if index[cell] < index[other]]
    if not pairs:
        return 0.0
    a, b = np.array(pairs).T
    lat, lon = grid.cells["lat"].to_numpy(), grid.cells["lon"].to_numpy()
    return float(np.median(_arc_km(lat[a], lon[a], lat[b], lon[b])))


def cell_graph(grid: Grid, km: float) -> sparse.csr_matrix:
    """The cells of the grid, in its order, joined where their centres are at most `km` apart: a symmetric matrix with no entry on the diagonal."""
    xyz = _vectors(grid.cells["lat"].to_numpy(), grid.cells["lon"].to_numpy())
    chord = 2 * np.sin(min(km, np.pi * EARTH_KM) / (2 * EARTH_KM))
    pairs = cKDTree(xyz).query_pairs(chord * (1 + 1e-9), output_type="ndarray")
    n = len(xyz)
    a, b = pairs[:, 0], pairs[:, 1]
    return sparse.csr_matrix((np.ones(2 * len(pairs), dtype=np.int8), (np.concatenate([a, b]), np.concatenate([b, a]))), shape=(n, n))


# ---- the cells of a species

def by_cell(presence: pl.DataFrame) -> pl.DataFrame:
    """`species_id`, `cell`, `records` (summed), `n_datasets` and `n_years` (the largest, a lower bound) per species and cell, from rows by part."""
    missing = [c for c in NEEDED if c not in presence.columns]
    if missing:
        raise ValueError(f"presence needs the columns {missing}, and has {presence.columns}")
    return presence.group_by("species_id", cell_of(pl.col("part")).alias("cell")).agg(
        pl.col("records").sum(), pl.col("n_datasets").max(), pl.col("n_years").max()
    )


def _marked(cells: pl.DataFrame, settings: "Settings") -> pl.DataFrame:
    """`cells` (from `by_cell`), sorted by species and then by records, with `species_records`, `species_cells` and whether the cell is `core`."""
    ordered = cells.sort(["species_id", "records", "cell"], descending=[False, True, False])
    total = pl.col("records").sum().over("species_id")
    before = pl.col("records").cum_sum().over("species_id") - pl.col("records")  # what the larger cells hold
    ordered = ordered.with_columns(total.alias("species_records"), pl.len().over("species_id").alias("species_cells"), before.alias("before"))
    core = (pl.col("species_records") < settings.heuristic_min_records) | (pl.col("before") + EPS < settings.core_share * pl.col("species_records"))
    return ordered.with_columns(core.alias("core")).drop("before")


def core_parts(presence: pl.DataFrame, settings: "Settings") -> pl.DataFrame:
    """The parts of each species' core cells: `species_id`, `part` and `cell`. `presence` is rows by species and part (`Cleaned.presence`)."""
    core = _marked(by_cell(presence), settings).filter(pl.col("core")).select("species_id", "cell")
    parts = presence.select("species_id", "part", cell_of(pl.col("part")).alias("cell")).unique()
    return parts.join(core, on=["species_id", "cell"], how="semi").sort("species_id", "part")


# ---- the rules

def _satellites(cells: pl.DataFrame, grid: Grid, settings: "Settings") -> pl.DataFrame:
    """The cells of the patches that `detached_satellites` takes out, from the cells `_marked` has marked."""
    schema = {"species_id": pl.String, "cell": pl.String, "records": pl.Int64, **dict.fromkeys(("patch_records", "patch_cells", "patch_datasets", "patch_years"), pl.Int64)}
    big = cells.filter(pl.col("species_records") >= settings.heuristic_min_records)
    if big.is_empty():
        return pl.DataFrame(schema=schema)
    index = grid.cells.select("cell").with_row_index("cell_idx").with_columns(pl.col("cell_idx").cast(pl.Int64))
    big = big.join(index, on="cell", how="left")
    if big["cell_idx"].null_count():
        raise ValueError(f"cells that are not in the grid: {big.filter(pl.col('cell_idx').is_null())['cell'].unique().head(5).to_list()}")
    big = big.with_columns(pl.col("species_id").rank("dense").cast(pl.Int64).alias("species_idx")).sort("species_idx", "cell_idx")
    graph = cell_graph(grid, settings.satellite_gap_km + neighbour_spacing(grid))
    labels = _patches(big["species_idx"].to_numpy(), big["cell_idx"].to_numpy(), graph)
    patches = (
        big.with_columns(pl.Series("patch", labels))
        .with_columns(
            pl.col("records").sum().over("patch").alias("patch_records"),
            pl.len().over("patch").cast(pl.Int64).alias("patch_cells"),
            pl.col("n_datasets").max().over("patch").alias("patch_datasets"),
            pl.col("n_years").max().over("patch").alias("patch_years"),
            pl.col("core").any().over("patch").alias("patch_core"),
        )
    )
    return patches.filter(
        ~pl.col("patch_core"),
        pl.col("patch_records") + EPS < settings.satellite_record_share * pl.col("species_records"),
        pl.col("patch_cells") + EPS < settings.satellite_cell_share * pl.col("species_cells"),
        ~((pl.col("patch_datasets") >= settings.satellite_min_datasets) & (pl.col("patch_years") >= settings.satellite_min_years)),
    ).select(*schema)


def _patches(species: np.ndarray, cells: np.ndarray, graph: sparse.csr_matrix) -> np.ndarray:
    """The patch of each (species, cell) node: the connected components of a species' cells on `graph`.

    The nodes are sorted by species and then cell, `species` and `cells` are integers, and `graph` is over the cells. Every node's neighbours on
    the graph are looked up among the nodes of its own species, a chunk of nodes at a time, so that the species are never joined to one another.
    """
    n, n_cells = len(species), graph.shape[0]
    keys = species * n_cells + cells
    indptr, indices = graph.indptr.astype(np.int64), graph.indices.astype(np.int64)
    degree = indptr[cells + 1] - indptr[cells]
    ends = np.cumsum(degree)
    first, second = [], []
    start = 0
    while start < n:
        base = int(ends[start - 1]) if start else 0
        stop = max(int(np.searchsorted(ends, base + EXPANSIONS, side="right")), start + 1)
        d = degree[start:stop]
        total = int(d.sum())
        if total:
            node = np.repeat(np.arange(start, stop, dtype=np.int64), d)
            offset = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(d) - d, d)
            near = indices[np.repeat(indptr[cells[start:stop]], d) + offset]
            want = species[node] * n_cells + near
            found = np.minimum(np.searchsorted(keys, want), n - 1)
            hit = (keys[found] == want) & (node < found)  # each pair once
            first.append(node[hit].astype(np.int32))
            second.append(found[hit].astype(np.int32))
        start = stop
    a = np.concatenate(first) if first else np.empty(0, dtype=np.int32)
    b = np.concatenate(second) if second else np.empty(0, dtype=np.int32)
    links = sparse.coo_matrix((np.ones(len(a), dtype=np.int8), (a, b)), shape=(n, n))
    return csgraph.connected_components(links, directed=False)[1]


def detached_satellites(presence: pl.DataFrame, grid: Grid, settings: "Settings") -> pl.DataFrame:
    """The cells of the patches that a species has apart from its main range, by the rule `detached_satellite`: `species_id`, `cell`, the cell's `records`,
    and its patch's `patch_records`, `patch_cells`, `patch_datasets` and `patch_years`. `presence` is rows by species and part, with the counts."""
    return _satellites(_marked(by_cell(presence), settings), grid, settings)


def _under_reported(cells: pl.DataFrame, groups: pl.DataFrame, settings: "Settings") -> pl.DataFrame:
    schema = {"species_id": pl.String, "cell": pl.String, "records": pl.Int64, "group_records": pl.Int64, "rate": pl.Float64, "core_rate": pl.Float64}
    group = groups.drop_nulls("tax_group").sort("species_id", "tax_group").unique("species_id", keep="first", maintain_order=True)
    grouped = cells.join(group.select("species_id", "tax_group"), on="species_id", how="inner")
    grouped = grouped.with_columns(pl.col("records").sum().over("tax_group", "cell").alias("group_records"))
    candidates = grouped.filter(
        ~pl.col("core"),
        pl.col("records") <= settings.under_max_records,
        pl.col("species_records") >= settings.heuristic_min_records,
        pl.col("group_records") >= settings.group_min_records,
    )
    if candidates.is_empty():
        return pl.DataFrame(schema=schema)
    rate = (pl.col("records") / pl.col("group_records")).alias("rate")
    core = (
        grouped.filter(pl.col("core"))
        .join(candidates.select("species_id").unique(), on="species_id", how="semi")
        .group_by("species_id")
        .agg(rate.median().alias("core_rate"))
    )
    return (
        candidates.join(core, on="species_id", how="inner")
        .with_columns(rate)
        .filter(pl.col("rate") < settings.under_rate_factor * pl.col("core_rate"))
        .select(*schema)
    )


def under_reported(presence: pl.DataFrame, groups: pl.DataFrame, settings: "Settings") -> pl.DataFrame:
    """The cells that a species reports far less in than its group does, by the rule `under_reported`: `species_id`, `cell`, the species' `records` there,
    the `group_records` there, its `rate` (records over group records) and its median `core_rate`. `groups` is `species_id` and `tax_group`."""
    return _under_reported(_marked(by_cell(presence), settings), groups, settings)


def _thin(cells: pl.DataFrame, grid: Grid, settings: "Settings") -> pl.DataFrame:
    """All the cells of the thin species whose cells are over `thin_km` apart, with `km`, the greatest distance between two of them."""
    schema = {"species_id": pl.String, "cell": pl.String, "records": pl.Int64, "km": pl.Float64}
    thin = cells.filter(pl.col("species_records") < settings.thin_max_records, pl.col("species_cells") >= 2)
    if thin.is_empty():
        return pl.DataFrame(schema=schema)
    thin = thin.join(grid.cells.select("cell", "lat", "lon"), on="cell", how="left")
    if thin["lat"].null_count():
        raise ValueError(f"cells that are not in the grid: {thin.filter(pl.col('lat').is_null())['cell'].unique().head(5).to_list()}")
    pairs = thin.join(thin, on="species_id", suffix="_b").filter(pl.col("cell") < pl.col("cell_b"))
    km = _arc_km(pairs["lat"].to_numpy(), pairs["lon"].to_numpy(), pairs["lat_b"].to_numpy(), pairs["lon_b"].to_numpy())
    far = pairs.select("species_id").with_columns(pl.Series("km", km)).group_by("species_id").agg(pl.col("km").max()).filter(pl.col("km") > settings.thin_km)
    return thin.join(far, on="species_id", how="inner").select(*schema)


def thin_long_range(presence: pl.DataFrame, grid: Grid, settings: "Settings") -> pl.DataFrame:
    """The species that `thin_long_range` takes out whole: `species_id`, `cell` (every cell of the species), `records` and `km`, the greatest distance
    between its cells. The core guardrail does not apply: a species with under `thin_max_records` records is all core."""
    return _thin(_marked(by_cell(presence), settings), grid, settings)


def flag_outliers(presence: pl.DataFrame, groups: pl.DataFrame, grid: Grid, settings: "Settings") -> pl.DataFrame:
    """The (species, cell) pairs that layer 5 takes out, with the `rule` of each: `species_id`, `cell` and `rule`.

    Each rule is judged on `presence` as it is, and a pair that several flag goes with the first of `RULE_ORDER`. Rules in `settings.rules_off`
    are not run.
    """
    cells = _marked(by_cell(presence), settings)
    runs = {
        "thin_long_range": lambda: _thin(cells, grid, settings),
        "detached_satellite": lambda: _satellites(cells, grid, settings),
        "under_reported": lambda: _under_reported(cells, groups, settings),
    }
    found = [runs[rule]().select("species_id", "cell", pl.lit(rule).alias("rule")) for rule in RULE_ORDER if rule not in settings.rules_off]
    if not found:
        return pl.DataFrame(schema=FLAGGED_SCHEMA)
    return pl.concat(found).unique(["species_id", "cell"], keep="first", maintain_order=True)


# ---- the report

def suspect_cells(cleaned: "Cleaned", grid: Grid) -> pl.DataFrame:
    """The cells where layer 5 took a far higher share of the species than it did next door. Lists, and removes nothing.

    A cell's share is the species that layer 5 took out of it over those it had before layer 5, which is what `cleaned.presence` holds there plus
    what layer 5 removed. The cell is suspect when its share is at least 3 times the median share of its neighbours (the cells that share a boundary
    with it, as `grid.adjacency` has them) and at least 10 percentage points above it. A neighbour with no species has no share and is left out of
    the median, and a cell with no neighbour that has a share is not judged. Sorted by share, highest first, with the cell's centre.

    `cell`, `lat`, `lon`, `species` (before layer 5), `removed`, `share`, `neighbours` (those with a share) and `neighbour_median`.
    """
    removed = cleaned.removed.filter(pl.col("layer") == 5).select("species_id", cell_of(pl.col("part")).alias("cell")).unique()
    left = cleaned.presence.select("species_id", cell_of(pl.col("part")).alias("cell")).unique()

    def count(frame: pl.DataFrame, name: str) -> pl.DataFrame:
        return frame.group_by("cell").agg(pl.len().cast(pl.Int64).alias(name))

    cells = (
        count(left, "left")
        .join(count(removed, "removed"), on="cell", how="full", coalesce=True)
        .with_columns(pl.col("left", "removed").fill_null(0))
        .with_columns((pl.col("left") + pl.col("removed")).alias("species"))
        .with_columns((pl.col("removed") / pl.col("species")).alias("share"))
    )
    adjacent = pl.DataFrame(
        [(cell, other) for cell, near in grid.adjacency.items() for other in near], schema={"cell": pl.String, "neighbour": pl.String}, orient="row"
    )
    around = (
        adjacent.join(cells.select(pl.col("cell").alias("neighbour"), pl.col("share").alias("neighbour_share")), on="neighbour", how="inner")
        .group_by("cell")
        .agg(pl.len().cast(pl.Int64).alias("neighbours"), pl.col("neighbour_share").median().alias("neighbour_median"))
    )
    found = (
        cells.join(around, on="cell", how="inner")
        .filter(
            pl.col("removed") > 0,
            pl.col("share") + EPS >= SUSPECT_RATIO * pl.col("neighbour_median"),
            pl.col("share") - pl.col("neighbour_median") + EPS >= SUSPECT_POINTS,
        )
        .join(grid.cells.select("cell", "lat", "lon"), on="cell", how="left")
    )
    return found.select(*SUSPECT_SCHEMA).sort(["share", "species", "cell"], descending=[True, True, False])
