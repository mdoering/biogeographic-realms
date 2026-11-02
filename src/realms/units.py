"""The units of the network: which cell parts enter it, and as what.

A part with too little sampling says little about which realm it belongs to, so it is not a node of its own:

1. A part whose Good-Turing coverage or species count is too low joins the coarse unit of its coarse cell and kind, `r5:<cell>:<kind>`
   (the resolution is the setting's). The coarse cell is the one that holds the centre of the part's own cell. All sparse parts of one
   coarse cell and kind merge, and the merge counts each species once. Parts that are not sparse stay as they are.
2. A coarse unit that is still sparse is left out of the network, as `sparse`.
3. A unit with too few major groups sampled is left out of the network, as `composition`.

A part that is left out is only in `Units.excluded`. `impute` then gives it the realm of its nearest network parts of the same kind,
and says it did, so that no map shows an imputed part as if it were clustered.

The records here are the clean records, the `records` of `Cleaned.presence`.
"""

import logging
import re
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd
import numpy as np
import polars as pl
import shapely
from scipy.spatial import cKDTree

from realms.grid import KINDS, ZERO_KM2, Grid, locate, part_adjacency, parts
from realms.natives import TDWG_LEVEL3_CODE, _cell_shapes, _gpkg

log = logging.getLogger(__name__)

# COL's `taxGroup`, grouped into the major groups the composition gate counts. A group that COL gives no major group (below) is not counted.
_GROUPS = {
    "plants": ("angiosperms", "gymnosperms", "pteridophytes", "bryophytes", "plants"),
    "insects": ("coleoptera", "lepidoptera", "diptera", "hymenoptera", "hemiptera", "orthoptera", "trichoptera", "otherinsects", "insects"),
    "other arthropods": ("arachnids", "crustacean", "otherarthropods", "arthropods"),
    "molluscs": ("gastropods", "bivalves", "othermolluscs", "molluscs"),
    "fungi": ("ascomycetes", "basidiomycetes", "otherfungi", "fungi"),
    "other invertebrates": ("otheranimals", "animals"),
}
MAJOR_GROUPS = {group: major for major, groups in _GROUPS.items() for group in groups}
NOT_COUNTED = ("algae", "protists", "bacteria", "archaea", "viruses", "eukaryotes", "pseudofungi", "prokaryotes")
CHORDATES = "chordates"  # COL's one `taxGroup` for the chordates, which `CHORDATE_CLASSES` splits
# The classes of the 92,927 chordates of COL release 315557 (20 of them, and some chordates with no class), as the extended export gives them. COL's
# `taxGroup` has one group for the chordates, but the gate needs land vertebrates, fish and the invertebrate chordates (tunicates and lancelets)
# apart. A chordate is a species of phylum Chordata whatever its group: 149 are grouped as `animals`, which would make fish and reptiles into
# other invertebrates. A chordate of a class in none of these lists, or of no class, is not counted.
CHORDATE_CLASSES = {
    "land vertebrates": ("Mammalia", "Aves", "Reptilia", "Amphibia"),
    "fish": (
        "Teleostei", "Elasmobranchii", "Chondrostei", "Holocephali", "Myxini", "Holostei", "Petromyzonti", "Coelacanthi", "Cladistii", "Dipneusti",
        "Cephalaspidomorphi",
    ),
    "other invertebrates": ("Ascidiacea", "Thaliacea", "Appendicularia", "Leptocardii"),
}
CHORDATE_MAJOR = {cls: major for major, classes in CHORDATE_CLASSES.items() for cls in classes}
# The WGSRPD level 1 regions, by the `LEVEL1_COD` of the polygon file
TDWG_LEVEL1 = {
    "1": "Europe", "2": "Africa", "3": "Asia-Temperate", "4": "Asia-Tropical", "5": "Australasia", "6": "Pacific",
    "7": "Northern America", "8": "Southern America", "9": "Antarctic",
}
TDWG_LEVEL1_CODE = "LEVEL1_COD"
SEA = "sea"
NETWORK = "network"
UNPLACED = "unplaced"  # a land part that no region of the polygon file could be found for
COARSE = re.compile(r"^r\d+:")

STATS_SCHEMA = {
    "part": pl.String, "kind": pl.String, "records": pl.Int64, "species": pl.Int64, "singletons": pl.Int64, "coverage": pl.Float64, "groups_sampled": pl.Int64,
}
EXCLUDED_SCHEMA = {"part": pl.String, "reason": pl.String, "tried_unit": pl.String, "records": pl.Int64, "species": pl.Int64}
PRESENCE_SCHEMA = {"species_id": pl.String, "unit": pl.String, "records": pl.Int64}


@dataclass(frozen=True)
class UnitSettings:
    """Starting values, to be tuned on half A only.

    `min_coverage` and `min_species` are the design's, as is the coarse resolution. `min_group_species` (a group is sampled where it has this
    many species) and `min_groups` (a unit needs this many sampled groups) make the design's composition gate concrete. `impute_k` is how many
    network parts vote for a part that was left out.
    """

    min_coverage: float = 0.5
    min_species: int = 50
    coarse_resolution: int = 5
    min_group_species: int = 5
    min_groups: int = 2
    impute_k: int = 6


@dataclass(frozen=True)
class Units:
    """What enters the network, and what does not.

    `members`: each part that is in a unit, to its unit, which is the part's own id or a coarse id.
    `presence`: `species_id`, `unit`, `records`, summed over each unit's parts, for the units of the network only.
    `areas`: km2 of each unit, the sum of its parts' land or sea area (the whole cell of a coarse unit is not counted, only the parts that joined it).
    `adjacency`: for each unit, the units that hold a part next to one of its parts.
    `excluded`: one row per part left out: `part`, `reason` (`sparse` or `composition`), `tried_unit` (the unit it was tested as: the coarse unit of a sparse part,
    the part itself otherwise) and its `records` and `species`, so that what the gate took out can be counted. Such a part is in no other field.
    `stats`: every unit that was tested, `unit`, `kind`, the `part_stats` columns, `parts` (how many parts it holds) and `status` (`network`, or the reason it was left out).
    """

    members: dict[str, str]
    presence: pl.DataFrame
    areas: dict[str, float]
    adjacency: dict[str, set[str]]
    excluded: pl.DataFrame
    stats: pl.DataFrame
    settings: UnitSettings = field(default_factory=UnitSettings)


def coarse_id(resolution: int, cell: str, kind: str) -> str:
    """The id of the unit that the sparse parts of one coarse cell and kind merge into: `r<resolution>:<cell>:<kind>`.

    Read as a part id (`grid.cell_of`), its cell is `r<resolution>:<cell>`, which no cell of the grid is: cells are numbers, and
    `network_units` refuses a grid where one starts with `r<digits>:`. A coarse unit and a fine part therefore never share an id.
    """
    if kind not in KINDS:
        raise ValueError(f"a unit is one of {KINDS}, not {kind!r}")
    return f"r{resolution}:{cell}:{kind}"


def is_coarse(unit: str) -> bool:
    """Whether a unit id is a coarse unit's (see `coarse_id`) and not a part's."""
    return COARSE.match(unit) is not None


def _kind_of(unit: pl.Expr) -> pl.Expr:
    return unit.str.extract(r":([^:]*)$")


def _rows(presence: pl.DataFrame) -> pl.DataFrame:
    """`species_id`, `part` and the clean `records`, one row for each species and part (rows that repeat a pair add up).

    The clean count is `records`, or else `n - n_human_text`, as `clean` leaves them.
    """
    columns = presence.columns
    if "records" in columns:
        records = pl.col("records")
    elif "n" in columns and "n_human_text" in columns:
        records = pl.col("n") - pl.col("n_human_text")
    else:
        raise ValueError(f"the presence needs `records`, or `n` and `n_human_text`, and has {columns}")
    return (
        presence.select("species_id", "part", records.cast(pl.Int64).alias("records"))
        .group_by("part", "species_id")
        .agg(pl.col("records").sum())
        .sort("part", "species_id")
    )


def _stats(rows: pl.DataFrame, major: pl.DataFrame, settings: UnitSettings) -> pl.DataFrame:
    """`part_stats` on rows as `_rows` gives them. A part whose rows hold no clean record is still a part, with nothing in it."""
    live = rows.filter(pl.col("records") > 0)
    base = rows.group_by("part").agg(pl.col("records").sum().alias("records"))
    counted = live.group_by("part").agg(pl.len().alias("species"), (pl.col("records") == 1).sum().alias("singletons"))
    per_group = (
        live.join(major.drop_nulls("major_group").unique("species_id"), on="species_id")
        .group_by("part", "major_group")
        .agg(pl.len().alias("n"))
        .filter(pl.col("n") >= settings.min_group_species)
        .group_by("part")
        .agg(pl.len().alias("groups_sampled"))
    )
    return (
        base.join(counted, on="part", how="left")
        .join(per_group, on="part", how="left")
        .with_columns(pl.col("species", "singletons", "groups_sampled").fill_null(0).cast(pl.Int64))
        .with_columns(
            _kind_of(pl.col("part")).alias("kind"),
            pl.when(pl.col("records") > 0).then(1 - pl.col("singletons") / pl.col("records")).otherwise(0.0).alias("coverage"),
        )
        .select(*STATS_SCHEMA)
        .sort("part")
    )


def part_stats(presence: pl.DataFrame, major: pl.DataFrame, settings: UnitSettings = UnitSettings()) -> pl.DataFrame:
    """`part`, `kind`, `records`, `species`, `singletons`, `coverage` and `groups_sampled` for each part of a presence table.

    `presence` has `species_id` and `part`, and the clean records as `records` (or `n` and `n_human_text`). `major` is `species_id`,
    `major_group`, as `major_groups` gives it. Species, singletons and coverage come from the rows with clean records: a singleton is a
    species with exactly one in the part, and coverage is Good-Turing's, `1 - singletons / records`, or 0 for a part with none.
    `groups_sampled` is how many major groups have at least `settings.min_group_species` species there. The same function serves the units
    after coarsening, whose ids are in the `part` column.
    """
    return _stats(_rows(presence), major, settings)


def major_groups(taxa: pl.DataFrame, groups: pl.DataFrame) -> pl.DataFrame:
    """`species_id`, `major_group` for the species that count towards the composition gate, from COL's `taxGroup` and, for chordates, `class`.

    `taxa` is `References.taxa` (it needs `species_id`, `phylum` and `class`) and `groups` is `References.groups` (`species_id`, `tax_group`).
    A species of no group, or of a group that is not counted (algae, protists, bacteria...), is left out, as is a chordate whose class is in
    none of `CHORDATE_CLASSES`. A `taxGroup` that is in neither list is an error rather than a species quietly not counted.
    """
    joined = taxa.select("species_id", "phylum", "class").join(groups.select("species_id", "tax_group"), on="species_id", how="left")
    known = {*MAJOR_GROUPS, *NOT_COUNTED, CHORDATES}
    unknown = sorted(set(joined["tax_group"].drop_nulls().unique().to_list()) - known)
    if unknown:
        raise ValueError(f"COL taxGroups that are not in the major groups: {unknown}, add them to units.MAJOR_GROUPS or NOT_COUNTED")
    chordate = (pl.col("tax_group") == CHORDATES) | (pl.col("phylum") == "Chordata")
    return (
        joined.with_columns(
            pl.when(chordate)
            .then(pl.col("class").replace_strict(CHORDATE_MAJOR, default=None, return_dtype=pl.String))
            .otherwise(pl.col("tax_group").replace_strict(MAJOR_GROUPS, default=None, return_dtype=pl.String))
            .alias("major_group")
        )
        .drop_nulls("major_group")
        .select("species_id", "major_group")
    )


def _coarse_cells(grid: Grid, coarse: Grid, coarse_gpkg: Path | None, cells: list[str]) -> dict[str, str]:
    """For each of `cells` of the fine grid, the coarse cell that holds its centre."""
    if not cells:
        return {}
    shapes = _cell_shapes(coarse, _gpkg(coarse, coarse_gpkg))
    ids = coarse.cells["cell"].to_numpy()
    order = np.argsort(ids, kind="stable")  # so that a centre on an edge between cells goes to the lowest id
    centres = grid.cells.filter(pl.col("cell").is_in(cells)).sort("cell")
    found = locate(centres["lon"].to_numpy(), centres["lat"].to_numpy(), shapes[order])
    return dict(zip(centres["cell"].to_list(), ids[order][found].tolist(), strict=True))


def network_units(
    presence: pl.DataFrame, grid: Grid, coarse: Grid, major: pl.DataFrame, settings: UnitSettings = UnitSettings(), *, coarse_gpkg: Path | None = None
) -> Units:
    """Coarsen the sparse parts, apply the composition gate, and say what is left out.

    `presence` is `Cleaned.presence` (`species_id`, `part`, and the counts), `grid` the grid it is on and `coarse` the grid of
    `settings.coarse_resolution`, whose file is `coarse_gpkg` (`data/reference/grids/ISEA3H-<resolution>.gpkg` if not given). `major` is
    `species_id`, `major_group` (see `major_groups`). A part with no clean record is a part, and sparse.
    """
    if coarse.resolution != settings.coarse_resolution:
        raise ValueError(f"the coarse grid is at resolution {coarse.resolution}, and the settings ask for resolution {settings.coarse_resolution}")
    rows = _rows(presence)
    table = parts(grid, rows["part"].unique())  # checks that every part is in the grid
    clash = [p for p in table["part"] if is_coarse(p)]
    if clash:
        raise ValueError(f"part ids that read as a coarse unit's, which would collide with one: {clash[:3]}")

    fine = _stats(rows, major, settings)
    fails = (pl.col("coverage") < settings.min_coverage) | (pl.col("species") < settings.min_species)
    sparse = set(fine.filter(fails)["part"])
    by_cell = _coarse_cells(grid, coarse, coarse_gpkg, table.filter(pl.col("part").is_in(sparse))["cell"].unique().to_list())
    unit_of = {
        part: coarse_id(coarse.resolution, by_cell[cell], kind) if part in sparse else part for part, cell, kind in table.select("part", "cell", "kind").iter_rows()
    }
    parts_of = pl.DataFrame({"part": list(unit_of), "unit": list(unit_of.values())}, schema={"part": pl.String, "unit": pl.String})

    unit_rows = rows.join(parts_of, on="part").group_by("unit", "species_id").agg(pl.col("records").sum()).rename({"unit": "part"})
    tested = (
        _stats(unit_rows, major, settings)
        .join(parts_of.group_by("unit").agg(pl.len().cast(pl.Int64).alias("parts")), left_on="part", right_on="unit")
        .with_columns(
            pl.when(fails).then(pl.lit("sparse")).when(pl.col("groups_sampled") < settings.min_groups).then(pl.lit("composition")).otherwise(pl.lit(NETWORK)).alias("status")
        )
        .rename({"part": "unit"})
        .select("unit", "kind", "records", "species", "singletons", "coverage", "groups_sampled", "parts", "status")
        .sort("unit")
    )
    judged = parts_of.join(tested.select("unit", "status"), on="unit")
    kept = judged.filter(pl.col("status") == NETWORK)
    excluded = (
        judged.filter(pl.col("status") != NETWORK)
        .join(fine.select("part", "records", "species"), on="part")
        .select("part", pl.col("status").alias("reason"), pl.col("unit").alias("tried_unit"), "records", "species")
        .sort("part")
        .cast(EXCLUDED_SCHEMA)
    )

    members = dict(sorted(kept.select("part", "unit").iter_rows()))
    area = dict(table.select("part", "area_km2").iter_rows())
    areas: dict[str, float] = {}
    for part, unit in members.items():
        areas[unit] = areas.get(unit, 0.0) + area[part]
    near: dict[str, set[str]] = {unit: set() for unit in areas}
    for part, neighbours in part_adjacency(grid, members).items():
        near[members[part]].update(members[n] for n in neighbours if members[n] != members[part])
    log.info(
        "units: %d parts, %d of them sparse, %d units in the network (%d coarse), %d parts left out",
        len(unit_of), len(sparse), len(areas), sum(is_coarse(u) for u in areas), excluded.height,
    )
    return Units(
        members=members,
        presence=unit_rows.filter(pl.col("part").is_in(set(areas)), pl.col("records") > 0).rename({"part": "unit"}).select(*PRESENCE_SCHEMA).sort("unit", "species_id"),
        areas=dict(sorted(areas.items())),
        adjacency=dict(sorted(near.items())),
        excluded=excluded,
        stats=tested,
        settings=settings,
    )


def _vectors(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Points on the unit sphere. Their straight-line distances are in the order of their great-circle distances."""
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def impute(labels: Mapping[str, int], units: Units, grid: Grid) -> tuple[dict[str, int], set[str]]:
    """The realm of every part, and which parts only got theirs by imputation.

    `labels` gives the realm of each unit of the network (and of no other). A part of a unit has its unit's realm. An excluded part gets the
    majority realm of its `units.settings.impute_k` nearest parts of the network of the same kind, nearest by great-circle distance between
    the cells' centres, and every part of a coarse unit is one such part with its own cell's centre. A tie of votes goes to the lowest realm
    id, and a tie in distance at the edge of the k to the lower part id. The set is the parts that were imputed.
    """
    network = set(units.areas)
    if set(labels) != network:
        missing, extra = sorted(network - set(labels)), sorted(set(labels) - network)
        raise ValueError(f"labels must be those of the units of the network: missing {missing[:3]}, not units {extra[:3]}")
    realms = {part: labels[unit] for part, unit in units.members.items()}
    excluded = sorted(units.excluded["part"].to_list())
    if not excluded:
        return realms, set()

    centre = {cell: (lat, lon) for cell, lat, lon in grid.cells.select("cell", "lat", "lon").iter_rows()}

    def at(names: list[str]) -> np.ndarray:
        found = [centre[n.rpartition(":")[0]] for n in names]
        return _vectors(np.array([lat for lat, _ in found]), np.array([lon for _, lon in found]))

    for kind in KINDS:
        wanted = [p for p in excluded if p.rpartition(":")[2] == kind]
        if not wanted:
            continue
        donors = sorted(p for p in units.members if p.rpartition(":")[2] == kind)
        if not donors:
            raise ValueError(f"no part of kind {kind!r} is in the network to impute {wanted[:3]} from")
        votes = np.array([realms[p] for p in donors])
        where = at(donors)
        tree = cKDTree(where)
        k = min(units.settings.impute_k, len(donors))
        there = at(wanted)
        reach, nearest = tree.query(there, k=k)
        reach, nearest = reach.reshape(len(wanted), k), nearest.reshape(len(wanted), k)
        for n, part in enumerate(wanted):
            chosen = nearest[n]
            # parts as far away as the k-th are as near as it, and the ones that fit are taken by their id
            tied = tree.query_ball_point(there[n], reach[n, -1] * (1 + 1e-9) + 1e-12)
            if len(tied) > k:
                gap = np.round(np.linalg.norm(where[tied] - there[n], axis=1), 12)
                ranked = sorted(zip(gap.tolist(), tied, strict=True), key=lambda g: (g[0], donors[g[1]]))
                chosen = [d for _, d in ranked[:k]]
            count = Counter(votes[chosen].tolist())
            top = max(count.values())
            realms[part] = min(realm for realm, c in count.items() if c == top)
    return realms, set(excluded)


def part_regions(grid: Grid, level3: Path) -> pl.DataFrame:
    """`part`, `region`: the WGSRPD level 1 region of the centre of the cell of each land part of the grid (a part with land to speak of).

    A centre that lies in no region, as a coastal cell's can where WGSRPD draws its own coastline, takes the nearest region. `level3` is the
    level 3 polygon file, whose `LEVEL1_COD` gives the level 1 region.
    """
    regions = gpd.read_file(level3).sort_values(TDWG_LEVEL3_CODE, kind="stable")
    shapes = shapely.make_valid(regions.geometry.to_numpy())
    land = grid.cells.filter(pl.col("land_km2") > ZERO_KM2)
    found = locate(land["lon"].to_numpy(), land["lat"].to_numpy(), shapes)
    codes = [str(c) for c in regions[TDWG_LEVEL1_CODE].to_numpy()[found]]
    return pl.DataFrame(
        {"part": [f"{c}:land" for c in land["cell"]], "region": [TDWG_LEVEL1.get(c, c) for c in codes]}, schema={"part": pl.String, "region": pl.String}
    )


def removed_by_region(units: Units, regions: pl.DataFrame) -> pl.DataFrame:
    """`reason`, `region`, `parts` and `records`: what was left out of the network, by reason and where.

    A land part's region is the one `part_regions` gives (`regions` is its result), `unplaced` if there is none, and a sea part's is `sea`.
    `records` are the clean records of the parts.
    """
    return (
        units.excluded.join(regions, on="part", how="left")
        .with_columns(
            pl.when(pl.col("part").str.ends_with(":sea")).then(pl.lit(SEA)).otherwise(pl.col("region").fill_null(UNPLACED)).alias("region")
        )
        .group_by("reason", "region")
        .agg(pl.len().cast(pl.Int64).alias("parts"), pl.col("records").sum().cast(pl.Int64).alias("records"))
        .sort("reason", "region")
    )
