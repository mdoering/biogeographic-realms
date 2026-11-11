"""Native ranges and alien lists, joined to COL ids.

WCVP, GRIIS and WRiMS know species by name, and their own ids mean nothing in COL (a GRIIS `col:ID` is a GBIF URL). So every join
here goes through `name_index`: the names of COL's accepted species, and of the synonyms that point to them. Reads existing exports
only, and never asks ChecklistBank for a new one.

A name that fits several species of one standing is ambiguous and is never guessed. Nothing joins to it, and the counts say how many.
"""

import hashlib
import json
import logging
import re
import shutil
import tempfile
import time
import unicodedata
import zipfile
from collections.abc import Iterable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd
import numpy as np
import polars as pl
import shapely
from shapely import STRtree

from realms.gbif import gbif_match
from realms.grid import ZERO_KM2, Grid, _planar, _shifted, part_id
from realms.lists import resolve_no_wild_range
from realms.reference import ACCEPTED, GRIIS_ALIAS, load_tax_groups, load_taxa

log = logging.getLogger(__name__)

INDEX_COLUMNS = ["col:ID", "col:parentID", "col:status", "col:rank", "col:scientificName", "col:kingdom"]
# The kingdoms whose disagreement counts as a clash, as in the occurrence query: verbatim names like `Metazoa`, and real disagreements
# like Chromista against Plantae, aren't wrong matches, whereas a plant matched to an animal homonym is.
CLASH_KINGDOMS = ("Animalia", "Plantae", "Fungi")
WCVP_KINGDOM = "Plantae"  # WCVP's own `col:kingdom` is empty, and it lists vascular plants only
MATCHED, AMBIGUOUS, UNMATCHED = "matched", "ambiguous", "unmatched"
SYNONYM = "synonym"
STATUS_ORDER = {"accepted": 0, "provisionally accepted": 1, SYNONYM: 2}  # which of two usages with one name COL would have us read
TDWG3 = r"^[A-Z]{3}$"  # level 3 codes, as against the level 1 and 2 codes (digits) that WCVP sometimes gives
WCVP_MEANS = (None, "introduced", "uncertain")
# What a GRIIS list says an entry is, as `griis_names` gives it. `unstated` is an entry whose distribution has no `establishmentMeans` (or has no
# row at all): a GRIIS list is a list of aliens, so that is introduced. `name_only` is an entry of a list whose export has no distribution table,
# which says nothing. In the order that settles an entry with several: native wins, as in WCVP, and uncertain only if nothing else is said.
GRIIS_MEANS = ("native", "introduced", "unstated", "uncertain")
GRIIS_ALIEN = ("introduced", "unstated")
# ChecklistBank values that are a kind of one of GRIIS_MEANS in Darwin Core's establishmentMeans vocabulary, read as that one. `introduced assisted
# colonisation` is how the reimport of 2026-11-11 gives the US registers' `introduced: assisted colonization`, species released on purpose (DECISIONS)
GRIIS_MEANS_KIND = {"introduced assisted colonisation": "introduced"}
NAME_ONLY = "name_only"
GRIIS_LIST_SCHEMA = {
    "alias": pl.String, "country": pl.String, "subdivision": pl.String, "table": pl.Boolean,
    "entries": pl.UInt32, "introduced": pl.UInt32, "unstated": pl.UInt32, "native": pl.UInt32, "uncertain": pl.UInt32,
}
TDWG_LEVEL3_CODE = "LEVEL3_COD"
FALLBACK_FILE = "gbif_fallback.parquet"  # what `scripts/check_natives.py` got from GBIF for the names that failed
# the names Natural Earth gives an admin-1 polygon, in the languages it holds
NAME_COLUMNS = ("name", "name_alt", "name_local", "woe_name", "woe_label", "gn_name", "gns_name")
# `GRIIS-US-United` lists the conterminous United States, which has no admin-1 polygon of its own: the country less its two outliers
WHOLE_COUNTRY_EXCEPT = {("US", "United"): ("Alaska", "Hawaii")}
ADMIN1_SCHEMA = {"cell": pl.String, "country": pl.String, "subdivision": pl.String, "via": pl.String}
GRIIS_PAIR_SCHEMA = {"species_id": pl.String, "country": pl.String, "subdivision": pl.String}


@contextmanager
def _member(path: Path, member: str) -> Iterator[Path]:
    """`member` of a COLDP zip (or of an unzipped export) as a file Polars can scan, which is extracted to a temporary directory."""
    path = Path(path)
    if path.is_dir():
        file = path / member
        if not file.is_file():
            raise ValueError(f"no {member} in {path}")
        yield file
        return
    with zipfile.ZipFile(path) as archive:
        name = next((n for n in archive.namelist() if n.rsplit("/", 1)[-1] == member), None)
        if name is None:
            raise ValueError(f"no {member} in {path}")
        with tempfile.TemporaryDirectory(prefix="realms-") as tmp:
            target = Path(tmp) / member
            with archive.open(name) as source, target.open("wb") as out:
                shutil.copyfileobj(source, out, 1 << 20)
            yield target


def has_member(path: Path, member: str) -> bool:
    """Whether the COLDP zip (or unzipped export) at `path` holds `member`."""
    path = Path(path)
    if path.is_dir():
        return (path / member).is_file()
    with zipfile.ZipFile(path) as archive:
        return any(n.rsplit("/", 1)[-1] == member for n in archive.namelist())


def read_tsv(path: Path, member: str, columns: Sequence[str], optional: Sequence[str] = ()) -> pl.DataFrame:
    """Only the named columns of a COLDP file, all as text with empty values as null. The file is scanned, never held whole."""
    with _member(path, member) as file:
        lazy = pl.scan_csv(file, separator="\t", quote_char=None, infer_schema=False, null_values=[""])
        header = lazy.collect_schema().names()
        missing = [c for c in columns if c not in header]
        if missing:
            raise ValueError(f"missing {missing} in {member} of {path}, whose columns are {header}")
        return lazy.select([*columns, *(c for c in optional if c in header)]).collect(engine="streaming")


def name_index(coldp_zip: Path) -> pl.DataFrame:
    """Names that identify a COL species: `name`, `kingdom`, `species_id` and `status` (that of the usage the name belongs to).

    That is every accepted and provisionally accepted species, and the synonyms of species rank that point to one of them. A
    synonym of a subspecies or variety rolls up to that infraspecific taxon's species, as GBIF counts its records for the species.
    Ambiguous synonyms and misapplied names are left out, as a pro parte name could belong to several species. A synonym takes the
    kingdom of its accepted species, as the export gives synonyms none. `name` is `col:scientificName`, with no authorship.
    """
    usages = read_tsv(coldp_zip, "NameUsage.tsv", INDEX_COLUMNS)
    species = usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED)).select(
        pl.col("col:scientificName").alias("name"), pl.col("col:kingdom").alias("kingdom"), pl.col("col:ID").alias("species_id"), pl.col("col:status").alias("status")
    )
    accepted = usages.filter(pl.col("col:status").is_in(ACCEPTED)).select(
        pl.col("col:ID").alias("parent_id"), pl.col("col:rank").alias("parent_rank"), pl.col("col:parentID").alias("grandparent_id")
    )
    synonyms = (
        usages.filter(pl.col("col:rank") == "species", pl.col("col:status") == "synonym")
        .select(pl.col("col:scientificName").alias("name"), pl.col("col:parentID").alias("parent_id"))
        .join(accepted, on="parent_id", how="inner")
        .with_columns(
            pl.when(pl.col("parent_rank") == "species").then(pl.col("parent_id")).otherwise(pl.col("grandparent_id")).alias("species_id")
        )
        .select("name", "species_id")
        .join(species.select(pl.col("species_id"), pl.col("kingdom")), on="species_id", how="inner")  # a parent that isn't a species falls out here
    )
    return (
        pl.concat([species, synonyms.select("name", "kingdom", "species_id").with_columns(pl.lit(SYNONYM).alias("status"))])
        .drop_nulls("name")
        .with_columns(pl.col("status").replace_strict(STATUS_ORDER).alias("_order"))
        .sort("_order", maintain_order=True)  # an accepted name keeps its row over the same species' synonym of that name
        .unique(["name", "kingdom", "species_id"], keep="first", maintain_order=True)
        .select("name", "kingdom", "species_id", "status")
    )


def match_names(names: pl.DataFrame, index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> pl.DataFrame:
    """Match `names` (a `name` column, and `kingdom` if known) to species by exact name, adding `species_id`, `status`, `n_species` and `via`.

    `status` is `matched`, `ambiguous` (several species of one standing fit, so `species_id` stays null and `n_species` says how
    many) or `unmatched`. With a kingdom given, a name matches within it, or failing that outside a clash: a name only fails on its
    kingdom where the two differ and both are animals, plants or fungi. Where one name fits several species, an accepted species
    wins over a provisionally accepted one and that over a synonym's, as COL and GBIF rank them. A tie within a rank stays
    ambiguous.

    `fallback` is `gbif_fallback`'s answer for names that stayed unmatched: those it lists become `matched`, with `via` `gbif`
    where `via` is `name` for the rest. Other columns of `names`, and their order, are kept.
    """
    if "name" not in names.columns:
        raise ValueError(f"names needs a `name` column, has {names.columns}")
    rows = names.drop("species_id", "status", "n_species", "via", strict=False).with_row_index("_row")
    kingdom = pl.col("kingdom") if "kingdom" in rows.columns else pl.lit(None, dtype=pl.String)
    wanted = rows.select("_row", "name", kingdom.alias("_kingdom")).join(
        index.select("name", pl.col("kingdom").alias("_index_kingdom"), "species_id", pl.col("status").replace_strict(STATUS_ORDER).alias("_order")),
        on="name",
        how="inner",
    )
    clash = pl.col("_kingdom").is_in(CLASH_KINGDOMS) & pl.col("_index_kingdom").is_in(CLASH_KINGDOMS) & (pl.col("_kingdom") != pl.col("_index_kingdom"))
    same = (pl.col("_kingdom") == pl.col("_index_kingdom")).fill_null(False)
    fitting = (
        wanted.with_columns(same.alias("_same"), (~clash.fill_null(False)).alias("_fits"))
        .with_columns(pl.col("_same").any().over("_row").alias("_any_same"))
        .filter(pl.col("_kingdom").is_null() | pl.col("_same") | (~pl.col("_any_same") & pl.col("_fits")))
    )
    candidates = (
        fitting.filter(pl.col("_order") == pl.col("_order").min().over("_row"))
        .group_by("_row")
        .agg(pl.col("species_id").n_unique().alias("n_species"), pl.col("species_id").first())
    )
    matched = rows.join(candidates, on="_row", how="left").with_columns(pl.col("n_species").fill_null(0).cast(pl.UInt32))
    result = matched.with_columns(
        pl.when(pl.col("n_species") == 1).then(pl.col("species_id")).alias("species_id"),
        pl.when(pl.col("n_species") == 1).then(pl.lit(MATCHED)).when(pl.col("n_species") > 1).then(pl.lit(AMBIGUOUS)).otherwise(pl.lit(UNMATCHED)).alias("status"),
        pl.when(pl.col("n_species") == 1).then(pl.lit("name")).alias("via"),
    )
    if fallback is not None and fallback.height:
        keys = [c for c in ("name", "kingdom") if c in result.columns]
        found = fallback.select(*keys, pl.col("species_id").alias("_fallback")).unique(keys)
        result = result.join(found, on=keys, how="left", nulls_equal=True).with_columns(
            pl.when((pl.col("status") == UNMATCHED) & pl.col("_fallback").is_not_null()).then(pl.col("_fallback")).otherwise(pl.col("species_id")).alias("species_id"),
            pl.when((pl.col("status") == UNMATCHED) & pl.col("_fallback").is_not_null()).then(pl.lit(MATCHED)).otherwise(pl.col("status")).alias("status"),
            pl.when((pl.col("status") == UNMATCHED) & pl.col("_fallback").is_not_null()).then(pl.lit("gbif")).otherwise(pl.col("via")).alias("via"),
        ).drop("_fallback")
    return result.sort("_row").drop("_row")


def match_stats(matched: pl.DataFrame) -> dict:
    """How many names `match_names` matched, left ambiguous and could not place, counting each distinct name once."""
    keys = [c for c in ("name", "kingdom") if c in matched.columns]
    counts = dict(matched.unique(keys).group_by("status").len().iter_rows())
    names = sum(counts.values())
    return {
        "names": names,
        MATCHED: counts.get(MATCHED, 0),
        AMBIGUOUS: counts.get(AMBIGUOUS, 0),
        UNMATCHED: counts.get(UNMATCHED, 0),
        "matched_share": counts.get(MATCHED, 0) / names if names else 0.0,
    }


def wcvp_names(wcvp_zip: Path) -> pl.DataFrame:
    """WCVP's accepted species, as `wcvp_id`, `name` and `kingdom` (always Plantae, as WCVP is a list of plants)."""
    usages = read_tsv(wcvp_zip, "NameUsage.tsv", ["col:ID", "col:status", "col:rank", "col:scientificName"])
    return usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED)).select(
        pl.col("col:ID").alias("wcvp_id"), pl.col("col:scientificName").alias("name"), pl.lit(WCVP_KINGDOM).alias("kingdom")
    )


def wcvp_distributions(wcvp_zip: Path) -> pl.DataFrame:
    """WCVP's TDWG distributions as `wcvp_id`, `tdwg3` and `native`. Empty means native, `introduced` is not, and `uncertain` is left out.

    Only level 3 codes are kept: the few rows with a level 1 or 2 code (digits) can't say which region of it, and are dropped.
    """
    rows = read_tsv(wcvp_zip, "Distribution.tsv", ["col:taxonID", "col:areaID", "col:gazetteer", "col:establishmentMeans"])
    unknown = sorted(set(rows["col:establishmentMeans"].unique().to_list()) - set(WCVP_MEANS))
    if unknown:
        raise ValueError(f"WCVP establishmentMeans has values nobody has judged, {unknown}: native, introduced or uncertain?")
    gazetteers = rows["col:gazetteer"].unique().to_list()
    if gazetteers != ["tdwg"]:
        raise ValueError(f"WCVP distributions are expected in the tdwg gazetteer alone, not {gazetteers}")
    return (
        rows.filter(pl.col("col:areaID").str.contains(TDWG3), pl.col("col:establishmentMeans").is_null() | (pl.col("col:establishmentMeans") == "introduced"))
        .select(pl.col("col:taxonID").alias("wcvp_id"), pl.col("col:areaID").alias("tdwg3"), pl.col("col:establishmentMeans").is_null().alias("native"))
        .unique(maintain_order=True)
    )


def wcvp_regions(wcvp_zip: Path, index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> pl.DataFrame:
    """`species_id`, `tdwg3` and `native` for each COL species WCVP places in a TDWG level 3 region.

    A species that WCVP holds under several names, or that has a region both native and introduced, is native in the region. A WCVP
    name that is ambiguous or unmatched has no row. Regions of infraspecific taxa aren't added: of 371,856 species, 3 have a
    distribution of a subspecies or variety and none of their own.
    """
    matched = match_names(wcvp_names(wcvp_zip), index, fallback).filter(pl.col("status") == MATCHED)
    return (
        wcvp_distributions(wcvp_zip)
        .join(matched.select("wcvp_id", "species_id"), on="wcvp_id", how="inner")
        .group_by("species_id", "tdwg3")
        .agg(pl.col("native").any())
        .sort("species_id", "tdwg3")
    )


def _kingdoms(usages: pl.DataFrame) -> pl.Series:
    """The name of the kingdom above each usage, found by walking `col:parentID`, null where the list has no kingdom row."""
    parent = dict(zip(usages["col:ID"].to_list(), usages["col:parentID"].to_list(), strict=True))
    name = {i: n for i, r, n in zip(usages["col:ID"], usages["col:rank"], usages["col:scientificName"], strict=True) if r == "kingdom"}
    found = []
    for usage in usages["col:ID"].to_list():
        seen = set()
        while usage is not None and usage not in name and usage not in seen:
            seen.add(usage)
            usage = parent.get(usage)
        found.append(name.get(usage))
    return pl.Series("kingdom", found, dtype=pl.String)


def _griis_means(path: Path) -> pl.DataFrame | None:
    """`col:ID` and `means` of each taxon the distribution table of a GRIIS export gives an establishment means, or None for an export with no table.

    An entry's means is that of its distribution rows, and where they differ native wins and then introduced, then an empty value (`unstated`), and
    `uncertain` only where nothing else is said: a species native somewhere in a list is no alien of it. A value other than `native`,
    `introduced` and `uncertain` is refused, as nobody has judged it, except a kind of one of them (`GRIIS_MEANS_KIND`), which is read as that one.
    """
    if not has_member(path, "Distribution.tsv"):
        return None
    rows = read_tsv(path, "Distribution.tsv", ["col:taxonID", "col:establishmentMeans"])
    rows = rows.with_columns(pl.col("col:establishmentMeans").replace(GRIIS_MEANS_KIND))
    unknown = sorted(set(rows["col:establishmentMeans"].drop_nulls().unique().to_list()) - set(GRIIS_MEANS))
    if unknown:
        raise ValueError(f"{path} has establishmentMeans values nobody has judged, {unknown}: native, introduced or uncertain?")
    order = {m: i for i, m in enumerate(GRIIS_MEANS)}
    return (
        rows.select(pl.col("col:taxonID").alias("col:ID"), pl.col("col:establishmentMeans").fill_null("unstated").alias("means"))
        .with_columns(pl.col("means").replace_strict(order, return_dtype=pl.Int64).alias("_order"))
        .sort("_order")
        .unique("col:ID", keep="first", maintain_order=True)
        .select("col:ID", "means")
    )


def griis_names(griis_dirs: Sequence[Path], species_only: bool = True) -> pl.DataFrame:
    """The entries of GRIIS checklists, as `alias`, `country`, `subdivision`, `name`, `kingdom`, `rank` and `means`.

    Country and subdivision come from the file name (`GRIIS-US-Hawaii.zip`). A plain export has no kingdom column, so the kingdom is
    that of the kingdom row above the entry. Only accepted entries are alien taxa: a synonym is another name of one. With
    `species_only` False the other ranks stay, which is how to count what the default drops.

    `means` is what the list says the entry is, from its distribution table (see `_griis_means`): `native`, `introduced`, `uncertain` or `unstated`,
    which is an empty value or no row in a list that has a table. An export with no table, the plain exports, says nothing and its entries
    are `name_only`.
    """
    frames = []
    for path in griis_dirs:
        path = Path(path)
        alias = path.name.removesuffix(".zip")
        match = GRIIS_ALIAS.match(alias)
        if not match:
            raise ValueError(f"{path} is not a GRIIS checklist, which is named GRIIS-<country>[-<subdivision>]")
        usages = read_tsv(path, "NameUsage.tsv", ["col:ID", "col:parentID", "col:status", "col:rank", "col:scientificName"], optional=["col:kingdom"])
        kingdom = usages["col:kingdom"] if "col:kingdom" in usages.columns else pl.Series([None] * usages.height, dtype=pl.String)
        # most of an extended export's own kingdom is empty too, so those entries fall back on the kingdom row above them
        usages = usages.with_columns(pl.coalesce(kingdom, _kingdoms(usages)).alias("kingdom"))
        means = _griis_means(path)
        if means is None:
            usages = usages.with_columns(pl.lit(NAME_ONLY).alias("means"))
        else:
            usages = usages.join(means, on="col:ID", how="left").with_columns(pl.col("means").fill_null("unstated"))
        entries = usages.filter(pl.col("col:status").is_in(ACCEPTED))
        if species_only:
            entries = entries.filter(pl.col("col:rank") == "species")
        frames.append(
            entries.select(
                pl.lit(alias).alias("alias"),
                pl.lit(match.group(1)).alias("country"),
                pl.lit(match.group(2), dtype=pl.String).alias("subdivision"),
                pl.col("col:scientificName").alias("name"),
                "kingdom",
                pl.col("col:rank").alias("rank"),
                "means",
            )
        )
    if not frames:
        return pl.DataFrame(
            schema={"alias": pl.String, "country": pl.String, "subdivision": pl.String, "name": pl.String, "kingdom": pl.String, "rank": pl.String, "means": pl.String}
        )
    return pl.concat(frames).drop_nulls("name")


def griis_lists(names: pl.DataFrame) -> pl.DataFrame:
    """What each GRIIS list holds, from `griis_names`: `alias`, `country`, `subdivision`, `table`, `entries`, and the entries that are `introduced`,
    `unstated` (an empty value, counted as introduced), `native` and `uncertain`, which `griis_pairs` leaves out.

    `table` is whether the list's export has a distribution table. A list without one says nothing of its entries, so its counts of the kinds are 0
    and `entries` is all there is to say. Entries are those of the list (accepted species, as named), whether COL knows them or not.
    """
    kinds = ["introduced", "unstated", "native", "uncertain"]
    return (
        names.group_by("alias", "country", "subdivision")
        .agg(
            (pl.col("means") != NAME_ONLY).all().alias("table"),
            pl.len().cast(pl.UInt32).alias("entries"),
            *(pl.col("means").eq(kind).sum().cast(pl.UInt32).alias(kind) for kind in kinds),
        )
        .select(*GRIIS_LIST_SCHEMA)
        .sort("alias")
    )


def griis_entries(names: pl.DataFrame, index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> pl.DataFrame:
    """`species_id`, `country`, `subdivision` and `means` for each COL species the entries of `names` (`griis_names`) fit, one row per species and list.

    Two entries of one list that fit one species (a name and its synonym) are one row, whose means is the one that wins (see `_griis_means`).
    """
    order = {m: i for i, m in enumerate([*GRIIS_MEANS, NAME_ONLY])}
    matched = match_names(names, index, fallback).filter(pl.col("status") == MATCHED)
    return (
        matched.with_columns(pl.col("means").replace_strict(order, return_dtype=pl.Int64).alias("_order"))
        .sort("_order")
        .unique(["species_id", "country", "subdivision"], keep="first", maintain_order=True)
        .select("species_id", "country", "subdivision", "means")
    )


def _pairs(entries: pl.DataFrame, means: Sequence[str]) -> pl.DataFrame:
    keep = entries.filter(pl.col("means").is_in(list(means))).select("species_id", "country", "subdivision")
    return keep.unique(maintain_order=True).sort("species_id", "country", "subdivision", nulls_last=True)


def griis_split(names: pl.DataFrame, index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The pairs of `griis_pairs` and of `griis_name_only_pairs`, from the entries `names` (`griis_names`): the one place the two are made, which
    `load_references` calls, so that what the tests see is what the references hold."""
    entries = griis_entries(names, index, fallback)
    return _pairs(entries, GRIIS_ALIEN), _pairs(entries, (NAME_ONLY,))


def griis_pairs(griis_dirs: Sequence[Path], index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> pl.DataFrame:
    """`species_id`, `country` and `subdivision` for each COL species an alien list names as an alien. Ranks above species (and below) are dropped.

    An entry is an alien where its distribution says `introduced`, or says nothing (`unstated`: a GRIIS list is a list of aliens). `native` and `uncertain`
    entries are left out, and so are the entries of a list whose export has no distribution table (`griis_name_only_pairs`).
    """
    return griis_split(griis_names(griis_dirs), index, fallback)[0]


def griis_name_only_pairs(griis_dirs: Sequence[Path], index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> pl.DataFrame:
    """The pairs of `griis_pairs`' shape for the lists whose export has no distribution table, which say nothing of what their entries are."""
    return griis_split(griis_names(griis_dirs), index, fallback)[1]


def _cell_shapes(grid: Grid, gpkg: Path) -> np.ndarray:
    """The polygon of every cell of `grid`, in the order of `grid.cells`, drawn as `load_grid` draws them."""
    frame = gpd.read_file(gpkg)
    ids = [str(i) for i in frame["geogrid_id"]]
    drawn = {cell: _planar(cell, shape, lon) for cell, shape, lon in zip(ids, frame.geometry, frame["lon_c"], strict=True)}
    cells = grid.cells["cell"].to_list()
    if set(cells) != set(drawn):
        raise ValueError(f"{gpkg} holds {len(drawn)} cells, which are not the {len(cells)} of the grid at resolution {grid.resolution}")
    return np.array([drawn[c] for c in cells], dtype=object)


def _overlaps(cells: np.ndarray, regions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Index pairs (cell, region) whose polygons share some area, with the regions also tried a turn of longitude to either side.

    Cells are drawn around their own centre, so a region at the antimeridian (Chukotka, Fiji, the Aleutians) may sit 360 degrees
    from the cell that holds it. Touching along an edge is not overlapping.
    """
    copies = np.concatenate([_shifted(regions, degrees) for degrees in (-360.0, 0.0, 360.0)])
    shapely.prepare(copies)
    i, j = STRtree(copies).query(cells, predicate="intersects")
    shared = shapely.area(shapely.intersection(cells[i], copies[j])) > 0
    return i[shared], j[shared] % len(regions)


def tdwg3_parts(grid: Grid, level3: Path, gpkg: Path | None = None) -> pl.DataFrame:
    """`part` and `tdwg3`: each land part, with every TDWG level 3 region that its cell's polygon intersects.

    The cell's polygon is the whole cell, not the land of it, as WGSRPD draws its own coastline. A cell on a border has several
    regions, and a land part no region reaches (a speck of coast the scheme leaves out) has none. `gpkg` is the grid's file, which
    holds the polygons `Grid` doesn't keep, and defaults to `data/reference/grids/ISEA3H-<resolution>.gpkg`.
    """
    regions = gpd.read_file(level3)
    codes = regions[TDWG_LEVEL3_CODE].to_numpy()
    cells = grid.cells.with_row_index("_i").filter(pl.col("land_km2") > ZERO_KM2)
    polygons = _cell_shapes(grid, _gpkg(grid, gpkg))[cells["_i"].to_numpy()]
    i, j = _overlaps(polygons, shapely.make_valid(regions.geometry.to_numpy()))
    found = pl.DataFrame({"cell": cells["cell"].to_numpy()[i], "tdwg3": codes[j]}, schema={"cell": pl.String, "tdwg3": pl.String})
    return (
        cells.select("cell")
        .join(found.unique(), on="cell", how="left")
        .select(pl.col("cell").map_elements(lambda c: part_id(c, "land"), return_dtype=pl.String).alias("part"), "tdwg3")
        .sort("part", "tdwg3", nulls_last=True)
    )


def _gpkg(grid: Grid, gpkg: Path | None) -> Path:
    return Path(gpkg) if gpkg else Path("data/reference/grids") / f"ISEA3H-{grid.resolution}.gpkg"


def _fold(text: object) -> str:
    """A place name as bare lower-case letters and digits, so that `Ha'apai`, `HAAPAI` and `Haapai` are one name."""
    plain = unicodedata.normalize("NFKD", str(text))
    return re.sub(r"[^a-z0-9]", "", "".join(c for c in plain if not unicodedata.combining(c)).lower())


def admin1_places(admin1: Path, subdivisions: Iterable[tuple[str, str]]) -> gpd.GeoDataFrame:
    """The Natural Earth admin-1 polygons that sub-national GRIIS lists name, as `country`, `subdivision` and `geometry`.

    A GRIIS subdivision (`Hawaii`, `StHelena`, `Rapa-Nui-2`) is only a name. It is compared with every name Natural Earth gives a
    polygon, folded to letters and digits. Where several polygons carry the name, those of the list's own country are kept, as
    `Madeira` is a place in Portugal and nowhere else. Natural Earth sometimes puts a place in another country than GRIIS does (Saba
    is `NL` there and `BQ` for GRIIS), so a name that one polygon alone carries is taken whatever its country. A subdivision no
    polygon carries gets no row. `GRIIS-US-United` is the conterminous United States, so it takes every US polygon but Alaska's and
    Hawaii's.
    """
    places = gpd.read_file(admin1)
    columns = [c for c in NAME_COLUMNS if c in places.columns]
    names = (
        pl.DataFrame({c: [v if isinstance(v, str) else None for v in places[c]] for c in columns}, schema={c: pl.String for c in columns})
        .with_row_index("polygon")
        .unpivot(columns, index="polygon", value_name="alias")
        .drop_nulls("alias")
        .with_columns(pl.col("alias").str.split("|"))
        .explode("alias")
        .select("polygon", pl.col("alias").map_elements(_fold, return_dtype=pl.String).alias("name"))
        .join(pl.DataFrame({"polygon": range(places.shape[0]), "iso": [v if isinstance(v, str) else None for v in places["iso_a2"]]}, schema={"polygon": pl.UInt32, "iso": pl.String}), on="polygon")
        .unique()
    )
    found = []
    for country, subdivision in sorted(set(subdivisions)):
        if (country, subdivision) in WHOLE_COUNTRY_EXCEPT:
            apart = [_fold(n) for n in WHOLE_COUNTRY_EXCEPT[(country, subdivision)]]
            left_out = names.filter(pl.col("iso") == country, pl.col("name").is_in(apart))["polygon"]
            hits = names.filter(pl.col("iso") == country, ~pl.col("polygon").is_in(left_out.implode()))
        else:
            hits = names.filter(pl.col("name") == _fold(subdivision))
            if hits["polygon"].n_unique() > 1:
                hits = hits.filter(pl.col("iso") == country)
        found += [(country, subdivision, polygon) for polygon in hits["polygon"].unique().to_list()]
    chosen = pl.DataFrame(found, schema={"country": pl.String, "subdivision": pl.String, "polygon": pl.Int64}, orient="row")
    return gpd.GeoDataFrame(
        {"country": chosen["country"].to_list(), "subdivision": chosen["subdivision"].to_list()},
        geometry=places.geometry.to_numpy()[chosen["polygon"].to_numpy()],
        crs=places.crs,
    )


def admin1_cells(grid: Grid, admin1: Path, subdivisions: Iterable[tuple[str, str]], gpkg: Path | None = None) -> pl.DataFrame:
    """`cell`, `country`, `subdivision` and `via`: the cells that belong to a sub-national GRIIS list's admin-1 polygons.

    A cell belongs to a subdivision when its polygon overlaps the subdivision's with some area, so coastal and border cells count, and
    one whose centre is at sea. The cells are big (7,292 at resolution 6) and the places small, and by centre alone a coastal cell, a
    border cell or an island fell out (the starling's US presences by the sea were left alone). `via` is `centre` where the cell's centre
    is in the polygon and `overlap` where the overlap alone puts it there. A cell can belong to several subdivisions, and a country's border
    cell to a subdivision and to a neighbour's country, which `clean` tells apart by the country code of the rows.
    A subdivision Natural Earth has no polygon for has no rows, and the caller finds those by comparing with `subdivisions`.
    `gpkg` is the grid's file, as for `tdwg3_parts`.
    """
    places = admin1_places(admin1, subdivisions)
    if places.empty:
        return pl.DataFrame(schema=ADMIN1_SCHEMA)
    centres = gpd.GeoDataFrame({"cell": grid.cells["cell"].to_list()}, geometry=gpd.points_from_xy(grid.cells["lon"], grid.cells["lat"]), crs=places.crs)
    inside = gpd.sjoin(centres, places, predicate="within")
    by_centre = pl.DataFrame(
        {"cell": inside["cell"].to_list(), "country": inside["country"].to_list(), "subdivision": inside["subdivision"].to_list(), "via": ["centre"] * len(inside)},
        schema=ADMIN1_SCHEMA,
    ).unique()
    polygons = _cell_shapes(grid, _gpkg(grid, gpkg))
    i, j = _overlaps(polygons, shapely.make_valid(places.geometry.to_numpy()))
    cells = grid.cells["cell"].to_numpy()
    by_overlap = pl.DataFrame(
        {"cell": cells[i], "country": places["country"].to_numpy()[j], "subdivision": places["subdivision"].to_numpy()[j], "via": ["overlap"] * len(i)},
        schema=ADMIN1_SCHEMA,
    ).unique()
    key = ["cell", "country", "subdivision"]
    return pl.concat([by_centre, by_overlap]).unique(key, keep="first", maintain_order=True).sort("country", "subdivision", "cell")


def gbif_fallback(
    names: pl.DataFrame, index: pl.DataFrame, per_second: float = 5.0, match=gbif_match, sleep=time.sleep, clock=time.monotonic
) -> pl.DataFrame:
    """Ask GBIF's match API (v2, on the COL checklist) for the names `match_names` could not place: `name`, `kingdom`, `species_id`, `match_type`.

    A name is taken only for an `EXACT` or `VARIANT` match whose accepted usage is a species of this release (not a synonym, a
    subspecies or a higher rank) and whose kingdom doesn't clash with the one given. Requests are sequential and at most
    `per_second` a second. This is the one place that asks GBIF for anything, so `names` should be what failed and no more.
    """
    own = index.filter(pl.col("status") != SYNONYM).select("species_id", pl.col("kingdom").alias("_found")).unique("species_id")
    known = dict(zip(own["species_id"], own["_found"], strict=True))
    interval = 1.0 / per_second
    last = None
    found = []
    todo = names.select([c for c in ("name", "kingdom") if c in names.columns]).unique(maintain_order=True)
    for row in todo.iter_rows(named=True):
        if last is not None and (wait := interval - (clock() - last)) > 0:
            sleep(wait)
        last = clock()
        result = match(row["name"])
        species, rank = (result.accepted_id, result.accepted_rank) if result.accepted_id else (result.usage_id, result.rank)
        if result.match_type not in ("EXACT", "VARIANT") or rank != "SPECIES" or species not in known:
            continue
        given, got = row.get("kingdom"), known[species]
        if given in CLASH_KINGDOMS and got in CLASH_KINGDOMS and given != got:
            continue
        found.append({"name": row["name"], "kingdom": given, "species_id": species, "match_type": result.match_type})
    return pl.DataFrame(found, schema={"name": pl.String, "kingdom": pl.String, "species_id": pl.String, "match_type": pl.String})


def write_fallback(found: pl.DataFrame, path: Path, release: int | str) -> None:
    """Save what `gbif_fallback` placed, with the COL release it was placed in: a species id means nothing in another release."""
    found.write_parquet(path, metadata={"release": str(release)})


def read_fallback(path: Path, release: int | str | None, index: pl.DataFrame) -> pl.DataFrame:
    """The saved fallback, if it was made for `release` (the COL release of `index`), without any species the index doesn't hold.

    A file with no release, or another one, raises: its species ids could point at different species. A species id that this index
    lacks is dropped, with a warning that counts them.
    """
    saved = pl.read_parquet_metadata(path).get("release")
    if saved is None or release is None or str(saved) != str(release):
        raise ValueError(f"{path} was made for COL release {saved}, not {release}: run scripts/check_natives.py again, or delete it")
    fallback = pl.read_parquet(path)
    known = index["species_id"].unique().implode()
    kept = fallback.filter(pl.col("species_id").is_in(known))
    if kept.height < fallback.height:
        log.warning("%d of the %d names GBIF placed in %s are on species that release %s lacks, and are dropped", fallback.height - kept.height, fallback.height, path, release)
    return kept


# GBIF's countries, by the adjective and the name Marine Regions uses for them in "<Italian> part of the <Tyrrhenian Sea>",
# "<Italian> Exclusive Economic Zone" and "Italy - <Western Mediterranean>". `Dominican` and `Congolese` could be two countries each
# and are left out, with the overlapping claims.
MARINE_COUNTRIES = {
    "AG": ("Antiguan and Barbudan", "Antigua and Barbuda"), "AL": ("Albanian", "Albania"), "AO": ("Angolan", "Angola"), "AR": ("Argentinian", "Argentina"),
    "AU": ("Australian", "Australia"), "BD": ("Bangladeshi", "Bangladesh"), "BE": ("Belgian", "Belgium"), "BG": ("Bulgarian", "Bulgaria"),
    "BA": ("Bosnian and Herzegovinian", "Bosnia and Herzegovina"), "BJ": ("Beninese", "Benin"), "BR": ("Brazilian", "Brazil"), "BS": ("Bahamian", "Bahamas"),
    "BZ": ("Belizean", "Belize"), "CA": ("Canadian", "Canada"), "CL": ("Chilean", "Chile"), "CM": ("Cameroonian", "Cameroon"), "CN": ("Chinese", "China"),
    "CO": ("Colombian", "Colombia"), "CR": ("Costa Rican", "Costa Rica"), "CU": ("Cuban", "Cuba"), "CV": ("Cape Verdean", "Cape Verde"),
    "CY": ("Cypriote", "Cyprus"), "DE": ("German", "Germany"), "DK": ("Danish", "Denmark"), "DZ": ("Algerian", "Algeria"), "EC": ("Ecuadorian", "Ecuador"),
    "EE": ("Estonian", "Estonia"), "EG": ("Egyptian", "Egypt"), "ES": ("Spanish", "Spain"), "FI": ("Finnish", "Finland"), "FJ": ("Fijian", "Fiji"),
    "FM": ("Micronesian", "Micronesia"), "FR": ("French", "France"), "GB": ("British", "United Kingdom"), "GE": ("Georgian", "Georgia"),
    "GH": ("Ghanaian", "Ghana"), "GM": ("Gambian", "Gambia"), "GR": ("Greek", "Greece"), "HN": ("Honduran", "Honduras"), "HR": ("Croatian", "Croatia"),
    "HT": ("Haitian", "Haiti"), "ID": ("Indonesian", "Indonesia"), "IE": ("Irish", "Ireland"), "IL": ("Israeli", "Israel"), "IN": ("Indian", "India"),
    "IQ": ("Iraqi", "Iraq"), "IR": ("Iranian", "Iran"), "IS": ("Icelandic", "Iceland"), "IT": ("Italian", "Italy"), "JM": ("Jamaican", "Jamaica"),
    "JP": ("Japanese", "Japan"), "KE": ("Kenyan", "Kenya"), "KI": ("Kiribati",), "KN": ("Kittitian and Nevisian", "Saint Kitts and Nevis"),
    "KR": ("South Korean", "South Korea"), "KZ": ("Kazakh", "Kazakhstan"), "LB": ("Lebanese", "Lebanon"), "LK": ("Sri Lankan", "Sri Lanka"),
    "LT": ("Lithuanian", "Lithuania"), "LV": ("Latvian", "Latvia"), "LY": ("Libyan", "Libya"), "MA": ("Moroccan", "Morocco"), "ME": ("Montenegrin", "Montenegro"),
    "MH": ("Marshallese", "Marshall Islands"), "MT": ("Maltese", "Malta"), "MX": ("Mexican", "Mexico"), "MY": ("Malaysian", "Malaysia"),
    "MZ": ("Mozambican", "Mozambique"), "NA": ("Namibian", "Namibia"), "NG": ("Nigerian", "Nigeria"), "NI": ("Nicaraguan", "Nicaragua"),
    "NL": ("Dutch", "Netherlands"), "NO": ("Norwegian", "Norway"), "NZ": ("New Zealand",), "OM": ("Omani", "Oman"), "PA": ("Panamanian", "Panama"),
    "PE": ("Peruvian", "Peru"), "PH": ("Philippine", "Philippines"), "PK": ("Pakistani", "Pakistan"), "PL": ("Polish", "Poland"), "PT": ("Portuguese", "Portugal"),
    "PW": ("Palauan", "Palau"), "RO": ("Romanian", "Romania"), "RU": ("Russian", "Russia"), "SA": ("Saudi Arabian", "Saudi Arabia"), "SD": ("Sudanese", "Sudan"),
    "SE": ("Swedish", "Sweden"), "SG": ("Singaporean", "Singapore"), "SI": ("Slovenian", "Slovenia"), "SL": ("Sierra Leonean", "Sierra Leone"),
    "SN": ("Senegalese", "Senegal"), "SO": ("Somali", "Somalia"), "SY": ("Syrian", "Syria"), "TH": ("Thai", "Thailand"), "TN": ("Tunisian", "Tunisia"),
    "TR": ("Turkish", "Turkey", "Türkiye"), "TT": ("Trinidadian and Tobagonian", "Trinidad and Tobago"), "TZ": ("Tanzanian", "Tanzania"),
    "UA": ("Ukrainian", "Ukraine"), "US": ("United States",), "UY": ("Uruguayan", "Uruguay"), "VE": ("Venezuelan", "Venezuela"), "VN": ("Vietnamese", "Vietnam"),
    "WS": ("Samoan", "Samoa"), "YE": ("Yemeni", "Yemen"), "ZA": ("South African", "South Africa"),
}
_MARINE_PHRASE = {phrase: code for code, phrases in MARINE_COUNTRIES.items() for phrase in phrases}
_MARINE_LEAD = re.compile(r"^(.+?)(?: part of the .+| Exclusive Economic Zone.*| - .+)?$")


def marine_country(area: str | None) -> str | None:
    """The country a Marine Regions area name names (`Greek part of the Aegean Sea`, `Tunisian Exclusive Economic Zone`, `Spain - ...`), or None."""
    match = _MARINE_LEAD.match(area or "")
    return _MARINE_PHRASE.get(match.group(1)) if match else None


def wrims_names(wrims_zip: Path) -> pl.DataFrame:
    """WRiMS's accepted species, as `wrims_id`, `name` and `kingdom`."""
    usages = read_tsv(wrims_zip, "NameUsage.tsv", ["col:ID", "col:status", "col:rank", "col:scientificName"], optional=["col:kingdom"])
    kingdom = pl.col("col:kingdom") if "col:kingdom" in usages.columns else pl.lit(None, dtype=pl.String)
    return usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED)).select(
        pl.col("col:ID").alias("wrims_id"), pl.col("col:scientificName").alias("name"), kingdom.alias("kingdom")
    )


def wrims_distributions(wrims_zip: Path) -> pl.DataFrame:
    """The areas WRiMS lists a taxon as introduced in: `wrims_id`, `area_id` (a Marine Regions MRGID) and `area`.

    `native`, `uncertain` and rows with no establishment means at all are left out.
    """
    rows = read_tsv(wrims_zip, "Distribution.tsv", ["col:taxonID", "col:areaID", "col:area", "col:establishmentMeans"])
    return (
        rows.filter(pl.col("col:establishmentMeans") == "introduced", pl.col("col:areaID").is_not_null())
        .select(pl.col("col:taxonID").alias("wrims_id"), pl.col("col:areaID").alias("area_id"), pl.col("col:area").alias("area"))
        .unique(maintain_order=True)
    )


def wrims_pairs(wrims_zip: Path, index: pl.DataFrame, fallback: pl.DataFrame | None = None) -> pl.DataFrame:
    """`species_id`, `area_id`, `area` and `country` for each COL species WRiMS lists as introduced in a marine area.

    The areas are Marine Regions places with no geometry in the export, and mostly the part of a sea inside a country's waters, so
    `country` is what the name says, and null for an open sea, a bay or a harbour. Mapping them to sea parts needs their polygons.
    """
    matched = match_names(wrims_names(wrims_zip), index, fallback).filter(pl.col("status") == MATCHED)
    return (
        wrims_distributions(wrims_zip)
        .join(matched.select("wrims_id", "species_id"), on="wrims_id", how="inner")
        .select("species_id", "area_id", "area", pl.col("area").map_elements(marine_country, return_dtype=pl.String).alias("country"))
        .unique(maintain_order=True)
        .sort("species_id", "area_id")
    )


@dataclass(frozen=True)
class References:
    """Everything the cleaning of occurrences needs to know about where taxa are native, alien or not wild at all."""

    taxa: pl.DataFrame  # `load_taxa`: accepted species with their classification, `extinct` and `environment`
    groups: pl.DataFrame  # `species_id`, `tax_group`
    wcvp: pl.DataFrame  # `species_id`, `tdwg3`, `native`: native and introduced regions both, as the cleaning can use either
    griis: pl.DataFrame  # `species_id`, `country`, `subdivision`: the species a list with distributions names as introduced, for a country or part of one
    wrims: pl.DataFrame  # `species_id`, `area_id`, `area`, `country`: the marine areas WRiMS lists a species as introduced in
    tdwg3_parts: pl.DataFrame  # `part`, `tdwg3`: each land part, with the WGSRPD regions its cell overlaps
    admin1_cells: pl.DataFrame  # `cell`, `country`, `subdivision`, `via`: the cells that overlap a sub-national GRIIS list's subdivision
    griis_no_cells: pl.DataFrame  # the `griis` rows of a subdivision with no admin-1 cells, which can't be placed on the grid
    wrims_no_country: pl.DataFrame  # the `wrims` rows whose area names no country: open seas, bays, harbours
    no_wild_range: frozenset[str]  # COL ids of humans, livestock and crops that have no wild range
    missing: tuple[str, ...] = ()  # manifest items with no export, which the references go without
    # the pairs of the lists whose export has no distribution table (a plain export), which say nothing of what their entries are. `clean` applies them
    # only under `Settings.griis_without_status="apply"`, and `griis_name_only_no_cells` is their `griis_no_cells`
    griis_name_only: pl.DataFrame = field(default_factory=lambda: pl.DataFrame(schema=GRIIS_PAIR_SCHEMA))
    griis_name_only_no_cells: pl.DataFrame = field(default_factory=lambda: pl.DataFrame(schema=GRIIS_PAIR_SCHEMA))
    griis_lists: pl.DataFrame = field(default_factory=lambda: pl.DataFrame(schema=GRIIS_LIST_SCHEMA))  # `griis_lists`: what each list holds, native entries and all
    # `species_id`, `worms_environment`: WoRMS's environment of the COL species it matched (`scripts/worms_environment.py`), empty if the manifest has none
    worms_environment: pl.DataFrame = field(default_factory=lambda: pl.DataFrame(schema=WORMS_ENVIRONMENT_SCHEMA))


WORMS_ENVIRONMENT_SCHEMA = {"species_id": pl.String, "worms_environment": pl.String}


def load_worms_environment(manifest: dict) -> pl.DataFrame:
    """The manifest's `worms_environment` table (`species_id`, `worms_environment`), checked against its sha256, or an empty one if it has none."""
    path = _manifest_path(manifest, "worms_environment")
    if path is None:
        return pl.DataFrame(schema=WORMS_ENVIRONMENT_SCHEMA)
    expected = manifest["worms_environment"].get("sha256")
    if expected and hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"{path} does not have the sha256 the manifest gives it: run scripts/worms_environment.py again")
    return pl.read_parquet(path).select(list(WORMS_ENVIRONMENT_SCHEMA))


def griis_without_cells(griis: pl.DataFrame, cells: pl.DataFrame) -> pl.DataFrame:
    """The `griis` pairs of a sub-national list whose subdivision has no row in `admin1_cells`."""
    placed = cells.select("country", "subdivision").unique()
    return griis.filter(pl.col("subdivision").is_not_null()).join(placed, on=["country", "subdivision"], how="anti")


def unknown_tdwg3(wcvp: pl.DataFrame, level3: Path) -> pl.DataFrame:
    """The level 3 codes of `wcvp` that no region of the polygon file has: `tdwg3`, `rows` and `species`, the most rows first."""
    codes = set(gpd.read_file(level3)[TDWG_LEVEL3_CODE])
    return (
        wcvp.filter(~pl.col("tdwg3").is_in(sorted(codes)))
        .group_by("tdwg3")
        .agg(pl.len().alias("rows"), pl.col("species_id").n_unique().alias("species"))
        .sort("rows", "tdwg3", descending=[True, False])
    )


def _manifest_path(manifest: dict, item: str) -> Path | None:
    entry = manifest.get(item)
    return Path(entry["path"]) if entry and entry.get("path") else None


def load_admin1_cells(grid: Grid, reference_dir: Path = Path("data/reference")) -> pl.DataFrame:
    """The `admin1_cells` of `load_references` alone: the cells of every sub-national GRIIS list that the manifest in `reference_dir` has an export for.

    It needs the manifest's admin-1 polygons and the grid's file, and nothing of COL, WCVP or WRiMS, so that `scripts/griis_codes.py` need not load the
    references (16 GB) to decide which country code each sub-national list applies to. `load_references` calls it, so the two agree.
    """
    reference_dir = Path(reference_dir)
    manifest = json.loads((reference_dir / "manifest.json").read_text())
    needed = {item: _manifest_path(manifest, item) for item in ("natural_earth_admin1", f"isea3h_{grid.resolution}")}
    for item, path in needed.items():
        if path is None:
            raise ValueError(f"{reference_dir / 'manifest.json'} has no {item}: run scripts/fetch_reference.py")
    griis = [k for k in sorted(manifest) if k.startswith("griis/") and _manifest_path(manifest, k) is not None]
    subdivisions = sorted({(m.group(1), m.group(2)) for k in griis if (m := GRIIS_ALIAS.match(k.removeprefix("griis/"))) and m.group(2)})
    return admin1_cells(grid, needed["natural_earth_admin1"], subdivisions, needed[f"isea3h_{grid.resolution}"])


def load_references(grid: Grid, reference_dir: Path = Path("data/reference")) -> References:
    """Read the reference data `fetch_reference.py` put in `reference_dir`, as listed in its manifest, and join it to COL ids.

    A WRiMS or GRIIS export the manifest lists as missing is skipped and named in `missing`. COL, WCVP, the TDWG polygons, the
    admin-1 polygons and the grid are needed, and a manifest without them raises. GRIIS sub-national lists are placed on `grid`.
    """
    reference_dir = Path(reference_dir)
    manifest = json.loads((reference_dir / "manifest.json").read_text())
    needed = {}
    for item in ("col", "col_groups", "wcvp", "natural_earth_admin1", "wgsrpd_level3", f"isea3h_{grid.resolution}"):
        needed[item] = _manifest_path(manifest, item)
        if needed[item] is None:
            raise ValueError(f"{reference_dir / 'manifest.json'} has no {item}: run scripts/fetch_reference.py")
    griis = {k: _manifest_path(manifest, k) for k in sorted(manifest) if k.startswith("griis/")}
    wrims = _manifest_path(manifest, "wrims")
    missing = tuple(sorted([k for k, v in griis.items() if v is None] + ([] if wrims else ["wrims"])))
    for item in missing:
        log.warning("%s has no export, and the references go without it", item)
    griis = {k: v for k, v in griis.items() if v is not None}

    index = name_index(needed["col"])
    fallback = None
    if (cache := reference_dir / FALLBACK_FILE).exists():
        fallback = read_fallback(cache, (manifest.get("col") or {}).get("release"), index)
    gpkg = needed[f"isea3h_{grid.resolution}"]
    wcvp = wcvp_regions(needed["wcvp"], index, fallback)
    unknown = unknown_tdwg3(wcvp, needed["wgsrpd_level3"])
    if unknown.height:
        shown = ", ".join(f"{r['tdwg3']} ({r['rows']} rows)" for r in unknown.head(5).iter_rows(named=True))
        log.warning(
            "%d WCVP level 3 codes are not in %s of %s, so no land part has them (%d rows): %s",
            unknown.height, TDWG_LEVEL3_CODE, needed["wgsrpd_level3"], unknown["rows"].sum(), shown,
        )
    names = griis_names(list(griis.values()))
    griis_joined, griis_name_only = griis_split(names, index, fallback)
    lists = griis_lists(names)
    del names
    wrims_ = wrims_pairs(wrims, index, fallback) if wrims else pl.DataFrame(schema={"species_id": pl.String, "area_id": pl.String, "area": pl.String, "country": pl.String})
    cells = load_admin1_cells(grid, reference_dir)
    griis_no_cells = griis_without_cells(griis_joined, cells)
    griis_name_only_no_cells = griis_without_cells(griis_name_only, cells)
    wrims_no_country = wrims_.filter(pl.col("country").is_null())
    if griis_no_cells.height:
        log.warning(
            "%d GRIIS pairs of %d sub-national lists have no admin-1 cells on the grid", griis_no_cells.height, griis_no_cells.select("country", "subdivision").unique().height
        )
    if griis_name_only_no_cells.height:
        log.warning(
            "%d GRIIS pairs of %d sub-national lists with no distributions have no admin-1 cells on the grid",
            griis_name_only_no_cells.height, griis_name_only_no_cells.select("country", "subdivision").unique().height,
        )
    if wrims_no_country.height:
        log.warning("%d of %d WRiMS pairs are in areas that name no country", wrims_no_country.height, wrims_.height)
    no_wild_range, _ = resolve_no_wild_range(index)
    del index
    return References(
        taxa=load_taxa(needed["col"]),
        groups=load_tax_groups(needed["col_groups"]),
        wcvp=wcvp,
        griis=griis_joined,
        wrims=wrims_,
        tdwg3_parts=tdwg3_parts(grid, needed["wgsrpd_level3"], gpkg),
        admin1_cells=cells,
        griis_no_cells=griis_no_cells,
        wrims_no_country=wrims_no_country,
        no_wild_range=no_wild_range,
        missing=missing,
        griis_name_only=griis_name_only,
        griis_name_only_no_cells=griis_name_only_no_cells,
        griis_lists=lists,
        worms_environment=load_worms_environment(manifest),
    )
