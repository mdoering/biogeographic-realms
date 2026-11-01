"""The cleaning of GBIF occurrence counts, before anything is clustered: layers 1 to 4 of the design.

`clean` takes the counts per species, country and part of a cell (`load_country_counts`) and returns what is left per species and part,
and a log of every row it took out or moved, with the layer and the rule that did. Nothing here is silent: a removal that cannot be
traced to a layer and a rule does not happen, and every rule can be switched off, either with its layer or on its own, so that a later
run can measure what each source of knowledge does to the map.

Order and attribution:

- **Zero-area parts** (layer 0, rule `zero_area_part`) fold first and always. A part with less than `ZERO_KM2` of area cannot carry a
  realm (`excess_table` divides by its area), so its records move into the other part of the same cell, if that has area. This is a
  repair of the data, not a cleaning layer, and no setting turns it off.
- **Layer 1** counts the records of a row: `n - n_human_text`, or `n_raw` with the layer off. A row with none is taken out.
- **Layer 2** takes out COL's extinct species and `NO_WILD_RANGE`.
- **Layer 3** takes out what the alien lists and WCVP say is not wild there. The rules run in this order, and a row goes with the first
  rule that catches it: the country's GRIIS list, the sub-national lists, WRiMS, then WCVP last. A row that is both on a list and outside
  WCVP's native regions is therefore the list's, and what `wcvp_*` reports is what no list had already removed. To measure one source
  alone, switch the others off (`Settings.rules_off`, `wcvp_mode`).
- **Layer 4** folds a clearly land or clearly sea species' records on the wrong side to the right side. COL's environment is the surer
  word of the two ways to be clear, so where it disagrees with the records `Settings.side_conflict` lets it decide, or leaves the species be.
- **Layer 5**, the rare-outlier heuristics, is accepted in `layers` and does nothing yet.

The records of a row are `records` throughout: the clean records the first layer leaves, or `n_raw` without it. Other counts (datasets,
years, licence and source counts) are kept as the download has them, summed as `load_presence` sums them over countries, folds and
parts, with the distinct counts taking the maximum, which is a lower bound.
"""

import logging
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import polars as pl

from realms.grid import KINDS, ZERO_KM2, Grid, part_id, parts
from realms.natives import SYNONYM, WCVP_KINGDOM, References
from realms.occurrences import COUNTS, DISTINCT, SUMMED

log = logging.getLogger(__name__)

LAYERS = (1, 2, 3, 4, 5)
WCVP_MODES = ("native_only", "introduced", "off")
SIDE_CONFLICTS = ("environment", "leave")
# the rules that can be switched off one at a time. WCVP has `wcvp_mode`, layer 1 and layer 4 have their layers, and the fold of
# zero-area parts has no switch.
SWITCHABLE = ("extinct", "no_wild_range", "griis", "griis_subnational", "wrims")
RULES = {
    "zero_area_part": 0,
    "failed_flags": 1,
    "human_text": 1,
    "extinct": 2,
    "no_wild_range": 2,
    "griis": 3,
    "griis_subnational": 3,
    "wrims": 3,
    "wcvp_outside_native": 3,
    "wcvp_introduced": 3,
    "land_sea_fold": 4,
    "land_sea_neighbour": 4,
    "land_sea_dropped": 4,
}
KEY = ["species_id", "countrycode", "part"]
ROW_COLUMNS = [*KEY, *COUNTS, "records"]
REMOVED_SCHEMA = {
    "layer": pl.Int8,
    "rule": pl.String,
    "species_id": pl.String,
    "part": pl.String,
    "countrycode": pl.String,
    "records": pl.Int64,
    "to_part": pl.String,
    "presence_lost": pl.Boolean,
}
UNPLACED_SCHEMA = {"rule": pl.String, "species_id": pl.String, "country": pl.String, "place": pl.String, "records": pl.Int64}
SUMMARY_KEYS = ("tax_group", "countrycode")
UNKNOWN_GROUP = "unknown"
# COL's environment is a sorted, comma-separated list ("brackish,freshwater,marine,terrestrial")
FRESH_OR_BRACKISH = r"(^|,)(freshwater|brackish)(,|$)"


@dataclass(frozen=True)
class Settings:
    """The thresholds of layers 4 and 5, at the design's values, and the choices between sources. All are starting values until the freeze.

    `land_share` is the share of a species' records on one side that makes it a land or a sea species, whichever it is. `wcvp_mode` is how
    layer 3 uses WCVP (see `clean`): `native_only` is the design's rule, and which mode the analysis takes is decided in M1.
    `side_conflict` is what layer 4 does where COL's environment and the records disagree about a species' side, such as a species
    that is exactly terrestrial with nearly all its records in a sea part: `environment` lets COL decide, `leave` leaves the species
    alone. `rules_off` names rules of `SWITCHABLE` to leave out one at a time.
    """

    land_share: float = 0.95
    satellite_gap_km: int = 500
    satellite_record_share: float = 0.01
    satellite_cell_share: float = 0.05
    satellite_min_datasets: int = 2
    satellite_min_years: int = 3
    group_min_records: int = 1000
    under_max_records: int = 5
    under_rate_factor: float = 0.01
    thin_max_records: int = 5
    thin_km: int = 1000
    core_share: float = 0.9
    heuristic_min_records: int = 20
    wcvp_mode: Literal["native_only", "introduced", "off"] = "native_only"
    side_conflict: Literal["environment", "leave"] = "environment"
    rules_off: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        if self.wcvp_mode not in WCVP_MODES:
            raise ValueError(f"wcvp_mode is one of {WCVP_MODES}, not {self.wcvp_mode!r}")
        if self.side_conflict not in SIDE_CONFLICTS:
            raise ValueError(f"side_conflict is one of {SIDE_CONFLICTS}, not {self.side_conflict!r}")
        if not 0.5 < self.land_share <= 1:
            raise ValueError(f"land_share is a share above one half, so that a species has one side at most, not {self.land_share}")
        object.__setattr__(self, "rules_off", frozenset(self.rules_off))
        unknown = sorted(self.rules_off - set(SWITCHABLE))
        if unknown:
            raise ValueError(f"rules_off can name {SWITCHABLE}, not {unknown}: layers go off with `layers`, WCVP with wcvp_mode, and zero_area_part is not a switch")


@dataclass(frozen=True)
class Cleaned:
    """What `clean` leaves, and what it did.

    `presence`: `species_id`, `part`, the counts, and `records`, summed over countries.
    `removed`: one row per row of the country counts that a layer took out or moved: `layer` (0 for the zero-area fold), `rule`,
    `species_id`, `part` (where the row was), `countrycode`, `records`, `to_part` (where a fold put it, null for a removal) and
    `presence_lost`, whether the species had nothing left at that part afterwards.
    `unplaced`: the pairs Task 10 could not place on the grid, for the species that have records: `rule` (the one that would
    have applied them), `species_id`, `country`, `place` (the subdivision or the marine area) and the species' `records`.
    """

    presence: pl.DataFrame
    removed: pl.DataFrame
    unplaced: pl.DataFrame = field(default_factory=lambda: pl.DataFrame(schema=UNPLACED_SCHEMA))


def _other_kind(kind: str) -> str:
    return "sea" if kind == "land" else "land"


def _part_table(grid: Grid) -> pl.DataFrame:
    """Every part of the grid: `part`, `cell`, `kind`, `area_km2`, the id and kind of the other part of its cell, and whether each has area."""
    table = parts(grid, [part_id(cell, kind) for cell in grid.cells["cell"] for kind in KINDS])
    other = [part_id(cell, _other_kind(kind)) for cell, kind in zip(table["cell"], table["kind"], strict=True)]
    table = table.with_columns(pl.Series("other", other, dtype=pl.String), pl.col("kind").replace_strict({"land": "sea", "sea": "land"}).alias("other_kind"))
    return (
        table.join(table.select(pl.col("part").alias("other"), pl.col("area_km2").alias("other_area")), on="other")
        .with_columns((pl.col("area_km2") >= ZERO_KM2).alias("has_area"), (pl.col("other_area") >= ZERO_KM2).alias("other_has_area"))
        .drop("other_area")
    )


def _sum_rows(rows: pl.DataFrame, by: Sequence[str] = KEY) -> pl.DataFrame:
    """Rows summed over `by`, the counts as `load_presence` has them: records add up, distinct counts take the larger."""
    return rows.group_by(list(by)).agg(*(pl.col(c).sum() for c in (*SUMMED, "records")), *(pl.col(c).max() for c in DISTINCT)).select(*by, *COUNTS, "records")


def _cell_key(cell: str) -> tuple[int, int, str]:
    """Cell ids are long integers kept as text, and the lowest is the lowest as a number: `9` comes before `100`."""
    try:
        return (0, int(cell), "")
    except ValueError:
        return (1, 0, cell)


class _Run:
    """The rows as the layers have left them, and the log of what they took."""

    def __init__(self, rows: pl.DataFrame) -> None:
        self.rows = rows
        self.removed: list[pl.DataFrame] = []

    def _log(self, hit: pl.DataFrame, layer: int, records: pl.Expr | None = None) -> None:
        """Log `hit`, rows with `rule` and `to_part`, after `self.rows` no longer has them (where they were taken out)."""
        if hit.is_empty():
            return
        pairs = hit.select("species_id", "part").unique()
        alive = pairs.join(self.rows.select("species_id", "part"), on=["species_id", "part"], how="semi")
        lost = pairs.join(alive, on=["species_id", "part"], how="anti").with_columns(pl.lit(True).alias("presence_lost"))
        entry = (
            hit.select(
                pl.lit(layer, dtype=pl.Int8).alias("layer"),
                "rule",
                "species_id",
                "part",
                "countrycode",
                (pl.col("records") if records is None else records).alias("records"),
                "to_part",
            )
            .join(lost, on=["species_id", "part"], how="left")
            .with_columns(pl.col("presence_lost").fill_null(False))
        )
        self.removed.append(entry.select(*REMOVED_SCHEMA))

    def drop(self, hit: pl.DataFrame, layer: int, rule: str | pl.Expr, records: pl.Expr | None = None, log: pl.Expr | None = None) -> None:
        """Take `hit`, a subset of the rows, out, and log it under `rule` (a name, or an expression for a column of `hit`)."""
        if hit.is_empty():
            return
        self.rows = self.rows.join(hit.select(KEY), on=KEY, how="anti", nulls_equal=True)
        rule_column = pl.lit(rule) if isinstance(rule, str) else rule
        logged = hit.with_columns(rule_column.alias("rule"), pl.lit(None, dtype=pl.String).alias("to_part"))
        self._log(logged if log is None else logged.filter(log), layer, records)

    def move(self, moves: pl.DataFrame, on: list[str], layer: int) -> None:
        """Move the rows that `moves` names (by `on`) to its `to_part`, or take them out where that is null, and log each under its `rule`.

        A moved row is added to the row of the same species and country at the part it goes to, if there is one. Rows with no records are
        moved, but not logged: a row that counts for nothing is not a presence.
        """
        hit = self.rows.join(moves, on=on, how="inner", nulls_equal=True)
        if hit.is_empty():
            return
        rest = self.rows.join(moves.select(on).unique(), on=on, how="anti", nulls_equal=True)
        moved = hit.filter(pl.col("to_part").is_not_null()).with_columns(pl.col("to_part").alias("part")).select(rest.columns)
        self.rows = _sum_rows(pl.concat([rest, moved]))
        self._log(hit.filter(pl.col("records") > 0), layer)


def _layer_1(run: _Run, layers: frozenset[int]) -> None:
    empty = run.rows.filter(pl.col("records") <= 0)
    if 1 not in layers:
        run.rows = run.rows.filter(pl.col("records") > 0)
        return
    # the records of a row that nothing counts for are all lost with it, so it is logged with `n_raw`. A row of no record at all is no row.
    rule = pl.when(pl.col("n") == 0).then(pl.lit("failed_flags")).otherwise(pl.lit("human_text"))
    run.drop(empty, 1, rule, records=pl.col("n_raw"), log=pl.col("n_raw") > 0)


def _layer_2(run: _Run, refs: References, settings: Settings) -> None:
    extinct = refs.taxa.filter(pl.col("extinct").fill_null(False)).select("species_id")
    no_wild_range = pl.DataFrame({"species_id": sorted(refs.no_wild_range)}, schema={"species_id": pl.String})
    for rule, species in (("extinct", extinct), ("no_wild_range", no_wild_range)):
        if rule not in settings.rules_off:
            run.drop(run.rows.join(species, on="species_id", how="semi"), 2, rule)


def _unplaced(rows: pl.DataFrame, refs: References, settings: Settings) -> pl.DataFrame:
    """The pairs that Task 10 could not place and layer 3 therefore cannot apply, of the species that have records, with their records."""
    found = []
    if "griis_subnational" not in settings.rules_off:
        found.append(
            refs.griis_no_cells.select(pl.lit("griis_subnational").alias("rule"), "species_id", "country", pl.col("subdivision").alias("place"))
        )
    if "wrims" not in settings.rules_off:
        found.append(
            refs.wrims_no_country.select(pl.lit("wrims").alias("rule"), "species_id", pl.lit(None, dtype=pl.String).alias("country"), pl.col("area").alias("place"))
        )
    if not found:
        return pl.DataFrame(schema=UNPLACED_SCHEMA)
    held = rows.group_by("species_id").agg(pl.col("records").sum())
    return (
        pl.concat(found)
        .unique()
        .join(held, on="species_id", how="inner")
        .select(*UNPLACED_SCHEMA)
        .sort("rule", "species_id", "country", "place", nulls_last=True)
    )


def _wcvp(rows: pl.DataFrame, refs: References, mode: str) -> pl.DataFrame:
    """The rows WCVP takes out: of plants, on land parts that have a TDWG region, by the mode.

    A part touches a region when `tdwg3_parts` lists it for the region. A sea part has none, and neither has a land part that the
    scheme leaves out, so WCVP says nothing of them and they stay. A plant WCVP does not list is not touched either.
    """
    plants = refs.taxa.filter(pl.col("kingdom") == WCVP_KINGDOM).select("species_id")
    region = refs.tdwg3_parts.drop_nulls("tdwg3").select("part", "tdwg3")
    pairs = (
        rows.select("species_id", "part")
        .unique()
        .join(plants, on="species_id", how="semi")
        .join(region.select("part").unique(), on="part", how="semi")
    )
    native_species = refs.wcvp.filter(pl.col("native")).select("species_id").unique().with_columns(pl.lit(True).alias("has_native"))
    touch = (
        pairs.join(region, on="part")
        .join(refs.wcvp.select("species_id", "tdwg3", "native"), on=["species_id", "tdwg3"], how="left")
        .group_by("species_id", "part")
        .agg(pl.col("native").any().alias("native"), (~pl.col("native")).any().alias("introduced"))
        .join(native_species, on="species_id", how="left")
        .with_columns(pl.col("has_native").fill_null(False))
    )
    if mode == "native_only":
        outside = touch.filter(pl.col("has_native") & ~pl.col("native"))
    else:
        outside = touch.filter(pl.col("introduced") & ~pl.col("native"))
    return rows.join(outside.select("species_id", "part"), on=["species_id", "part"], how="semi")


def _layer_3(run: _Run, refs: References, meta: pl.DataFrame, settings: Settings) -> pl.DataFrame:
    """Returns the pairs it could not place."""
    off = settings.rules_off
    unplaced = _unplaced(run.rows, refs, settings)
    places = meta.select("part", "cell", "kind")

    def located() -> pl.DataFrame:
        return run.rows.join(places, on="part", how="left")

    if "griis" not in off:
        # a country's list is for its rows by country code, except in the cells of its sub-national lists, which have lists of their own
        country = refs.griis.filter(pl.col("subdivision").is_null()).select("species_id", pl.col("country").alias("countrycode")).unique()
        covered = refs.admin1_cells.select("cell", pl.col("country").alias("countrycode")).unique()
        hit = located().join(country, on=["species_id", "countrycode"], how="semi").join(covered, on=["cell", "countrycode"], how="anti")
        run.drop(hit, 3, "griis")
    if "griis_subnational" not in off:
        # by the cell alone: an overseas place carries the sovereign's country code in its list and GBIF's own in its records
        listed = (
            refs.griis.filter(pl.col("subdivision").is_not_null())
            .join(refs.admin1_cells.select("cell", "country", "subdivision"), on=["country", "subdivision"])
            .select("species_id", "cell")
            .unique()
        )
        run.drop(located().join(listed, on=["species_id", "cell"], how="semi"), 3, "griis_subnational")
    if "wrims" not in off:
        introduced = refs.wrims.filter(pl.col("country").is_not_null()).select("species_id", pl.col("country").alias("countrycode")).unique()
        run.drop(located().filter(pl.col("kind") == "sea").join(introduced, on=["species_id", "countrycode"], how="semi"), 3, "wrims")
    if settings.wcvp_mode != "off":
        rule = "wcvp_outside_native" if settings.wcvp_mode == "native_only" else "wcvp_introduced"
        run.drop(_wcvp(run.rows, refs, settings.wcvp_mode), 3, rule)
    return unplaced


def _sides(rows: pl.DataFrame, taxa: pl.DataFrame, meta: pl.DataFrame, settings: Settings) -> pl.DataFrame:
    """`species_id` and `side` ("land", "sea" or null) of each species, from COL's environment and the side of its records.

    A land species is exactly terrestrial in COL or has at least `land_share` of its records on land parts, and a sea species is the
    same the other way round. Freshwater and brackish species and all of Aves are exempt. Where the two disagree, as for a species that
    is exactly marine with nearly all its records in land parts, `settings.side_conflict` says: COL's environment decides, or the
    species is left alone.
    """
    share = (
        rows.join(meta.select("part", "kind"), on="part", how="left")
        .group_by("species_id")
        .agg(
            pl.col("records").filter(pl.col("kind") == "land").sum().alias("land"),
            pl.col("records").filter(pl.col("kind") == "sea").sum().alias("sea"),
        )
        .with_columns((pl.col("land") + pl.col("sea")).alias("total"))
    )
    environment = pl.col("environment")
    exempt = environment.str.contains(FRESH_OR_BRACKISH).fill_null(False) | (pl.col("class") == "Aves").fill_null(False)
    land_env, sea_env = (environment == "terrestrial").fill_null(False), (environment == "marine").fill_null(False)
    land_records, sea_records = pl.col("land") / pl.col("total") >= settings.land_share, pl.col("sea") / pl.col("total") >= settings.land_share
    none = pl.lit(None, dtype=pl.String)
    found = (
        share.join(taxa.select("species_id", "environment", "class").unique("species_id"), on="species_id", how="left")
        .with_columns(exempt.alias("exempt"), land_env.alias("land_env"), sea_env.alias("sea_env"), land_records.alias("land_records"), sea_records.alias("sea_records"))
        .with_columns((pl.col("land_env") | pl.col("land_records")).alias("land_side"), (pl.col("sea_env") | pl.col("sea_records")).alias("sea_side"))
        .with_columns((~pl.col("exempt") & pl.col("land_side") & pl.col("sea_side")).alias("conflict"))
        .with_columns(
            pl.when(pl.col("exempt") | (pl.col("conflict") & (settings.side_conflict == "leave")))
            .then(none)
            .when(pl.col("land_env"))
            .then(pl.lit("land"))
            .when(pl.col("sea_env"))
            .then(pl.lit("sea"))
            .when(pl.col("land_records"))
            .then(pl.lit("land"))
            .when(pl.col("sea_records"))
            .then(pl.lit("sea"))
            .alias("side")
        )
    )
    both = found.filter(pl.col("conflict"))
    log.info(
        "layer 4: %d species, of which %d land, %d sea, %d exempt, and %d where COL's environment and the records disagree (%d exactly terrestrial "
        "with at least %s of the records at sea, %d exactly marine with that on land), which side_conflict=%s %s",
        found.height,
        found.filter(pl.col("side") == "land").height,
        found.filter(pl.col("side") == "sea").height,
        found.filter(pl.col("exempt")).height,
        both.height,
        both.filter(pl.col("land_env") & pl.col("sea_records")).height,
        settings.land_share,
        both.filter(pl.col("sea_env") & pl.col("land_records")).height,
        settings.side_conflict,
        "decides by the environment" if settings.side_conflict == "environment" else "leaves alone",
    )
    return found.select("species_id", "side")


def _layer_4(run: _Run, refs: References, grid: Grid, meta: pl.DataFrame, settings: Settings) -> None:
    sides = _sides(run.rows, refs.taxa, meta, settings).drop_nulls("side")
    places = meta.select("part", "cell", "kind", "other", "other_kind", "other_has_area")
    wrong = (
        run.rows.select("species_id", "part")
        .unique()
        .join(places, on="part", how="left")
        .join(sides, on="species_id")
        .filter(pl.col("side") != pl.col("kind"))
    )
    if wrong.is_empty():
        return
    same = wrong.filter(pl.col("other_has_area")).select("species_id", "part", pl.col("other").alias("to_part"), pl.lit("land_sea_fold").alias("rule"))
    far = wrong.filter(~pl.col("other_has_area")).select("species_id", "part", "cell", pl.col("other_kind").alias("right_kind"))

    # no such part in the cell: the right-side part of the neighbouring cell that holds most of the species' records, the lowest id on a tie
    adjacent = [(cell, other) for cell, near in grid.adjacency.items() for other in near]
    neighbours = pl.DataFrame(adjacent, schema={"cell": pl.String, "neighbour": pl.String}, orient="row")
    order = sorted(grid.adjacency, key=_cell_key)
    rank = pl.DataFrame({"neighbour": order, "rank": range(len(order))}, schema={"neighbour": pl.String, "rank": pl.Int64})
    right = meta.filter(pl.col("has_area")).select(pl.col("cell").alias("neighbour"), pl.col("kind").alias("right_kind"), pl.col("part").alias("to_part"))
    held = run.rows.group_by("species_id", "part").agg(pl.col("records").sum().alias("held")).rename({"part": "to_part"})
    chosen = (
        far.join(neighbours, on="cell")
        .join(right, on=["neighbour", "right_kind"])
        .join(held, on=["species_id", "to_part"])
        .join(rank, on="neighbour")
        .sort(["species_id", "part", "held", "rank"], descending=[False, False, True, False])
        .unique(["species_id", "part"], keep="first", maintain_order=True)
        .select("species_id", "part", "to_part", pl.lit("land_sea_neighbour").alias("rule"))
    )
    dropped = (
        far.select("species_id", "part")
        .join(chosen.select("species_id", "part"), on=["species_id", "part"], how="anti")
        .select("species_id", "part", pl.lit(None, dtype=pl.String).alias("to_part"), pl.lit("land_sea_dropped").alias("rule"))
    )
    run.move(pl.concat([same, chosen, dropped]), ["species_id", "part"], 4)


def clean(
    counts: pl.DataFrame,
    refs: References,
    grid: Grid,
    layers: frozenset[int] = frozenset(LAYERS),
    settings: Settings = Settings(),
) -> Cleaned:
    """Clean the country counts of `load_country_counts` through the `layers` that are on, with the rules of `settings`.

    See the module's docstring for the order, and how a row is attributed to a rule. WCVP is used by `settings.wcvp_mode`, for plants
    only (`kingdom` Plantae in `refs.taxa`), on land parts with a TDWG region:

    - `native_only`: a plant with a native region loses the parts touching none of them (`wcvp_outside_native`).
    - `introduced`: a plant loses the parts touching an introduced region and no native one (`wcvp_introduced`).
    - `off`: WCVP removes nothing, and GRIIS and WRiMS still apply.

    A country's GRIIS list applies to rows by `countrycode`, except in the cells of its sub-national lists, which apply by cell whatever
    the code (`griis_subnational`). WRiMS applies to sea parts, by `countrycode` (`wrims`). What `refs` could not place is not applied,
    and is listed in `Cleaned.unplaced`.

    Layer 4 folds a land species' row on a sea part, and a sea species' on a land part, into the other part of its cell if that has
    area (`land_sea_fold`), or else into the right-side part of the neighbouring cell that holds most of the species' records, the
    lowest cell id first on a tie (`land_sea_neighbour`), or else out (`land_sea_dropped`). A neighbour counts only if the species has
    records on that part. Where COL's environment and the records disagree about a species' side, `settings.side_conflict` decides.
    """
    layers = frozenset(layers)
    unknown = sorted(layers - set(LAYERS))
    if unknown:
        raise ValueError(f"layers are of {LAYERS}, not {unknown}")
    missing = [c for c in (*KEY, *COUNTS) if c not in counts.columns]
    if missing:
        raise ValueError(f"counts needs the columns {missing}, and has {counts.columns}")
    parts(grid, counts["part"].unique())  # refuses a part that is not on the grid, before anything is joined to the grid
    meta = _part_table(grid)

    records = (pl.col("n") - pl.col("n_human_text")) if 1 in layers else pl.col("n_raw")
    run = _Run(_sum_rows(counts.select(*KEY, *COUNTS).with_columns(records.alias("records"))))

    zero = meta.filter(pl.col("area_km2") < ZERO_KM2, pl.col("other_has_area")).select("part", pl.col("other").alias("to_part"), pl.lit("zero_area_part").alias("rule"))
    run.move(zero, ["part"], 0)
    _layer_1(run, layers)
    if 2 in layers:
        _layer_2(run, refs, settings)
    unplaced = pl.DataFrame(schema=UNPLACED_SCHEMA)
    if 3 in layers:
        unplaced = _layer_3(run, refs, meta, settings)
    if 4 in layers:
        _layer_4(run, refs, grid, meta, settings)

    presence = _sum_rows(run.rows, ["species_id", "part"]).sort("species_id", "part")
    removed = pl.concat(run.removed) if run.removed else pl.DataFrame(schema=REMOVED_SCHEMA)
    removed = removed.sort("layer", "rule", "species_id", "part", "countrycode", nulls_last=True)
    return Cleaned(presence, removed, unplaced)


def removal_summary(cleaned: Cleaned, refs: References, by: Sequence[str] = SUMMARY_KEYS) -> pl.DataFrame:
    """What was removed per `layer`, `rule` and the columns of `by` (`tax_group`, from `refs.groups`, and `countrycode`).

    `rows` are the rows of the country counts, `records` their records and `presences` the (species, part) pairs that were left with
    nothing by it. Rows of folds are counted as removed from the part they left. A country's presences do not add up over countries,
    as one part's presence can be in several, so use `by=("tax_group",)` or `by=()` for totals.

    `applied` is False for what the lists could not place (`unplaced_griis`, `unplaced_wrims`), which are not removed: `rows` are the
    pairs of species that have records, `records` the records of those species in all, as an upper bound on what the gap leaves, and
    `countrycode` the country of the list, null for WRiMS. Species with no group are `unknown`.
    """
    bad = [c for c in by if c not in SUMMARY_KEYS]
    if bad:
        raise ValueError(f"by is made of {SUMMARY_KEYS}, not {bad}")
    keys = list(by)
    groups = refs.groups.unique("species_id")
    group = pl.col("tax_group").fill_null(UNKNOWN_GROUP)
    done = (
        cleaned.removed.join(groups, on="species_id", how="left")
        .with_columns(group)
        .group_by("layer", "rule", *keys)
        .agg(
            pl.len().cast(pl.Int64).alias("rows"),
            pl.struct("species_id", "part").filter(pl.col("presence_lost")).n_unique().cast(pl.Int64).alias("presences"),
            pl.col("records").sum().cast(pl.Int64),
        )
        .with_columns(pl.col("layer").cast(pl.Int64), pl.lit(True).alias("applied"))
    )
    gap = (
        cleaned.unplaced.join(groups, on="species_id", how="left")
        .with_columns(
            group,
            pl.col("country").alias("countrycode"),
            pl.lit(3, dtype=pl.Int64).alias("layer"),
            pl.when(pl.col("rule") == "wrims").then(pl.lit("unplaced_wrims")).otherwise(pl.lit("unplaced_griis")).alias("rule"),
        )
        .group_by("layer", "rule", *keys)
        .agg(pl.len().cast(pl.Int64).alias("rows"), pl.lit(None, dtype=pl.Int64).alias("presences"), pl.col("records").sum().cast(pl.Int64))
        .with_columns(pl.lit(False).alias("applied"))
    )
    columns = ["layer", "rule", "applied", *keys, "rows", "presences", "records"]
    return pl.concat([done.select(columns), gap.select(columns)]).sort("layer", "rule", *keys, nulls_last=True)


def provisional_twins(index: pl.DataFrame) -> pl.DataFrame:
    """The provisionally accepted species whose name also belongs to an accepted species of the same kingdom, from `natives.name_index`.

    `species_id`, `name`, `kingdom` and `accepted_ids`, the number of accepted species with that name. Task 10 joins such a name to the
    accepted twin (and leaves it ambiguous where there are several), so a source's record on the name never reaches the provisional id.
    """
    species = index.filter(pl.col("status") != SYNONYM)
    accepted = species.filter(pl.col("status") == "accepted").group_by("name", "kingdom").agg(pl.len().alias("accepted_ids"))
    return (
        species.filter(pl.col("status") == "provisionally accepted")
        .select("species_id", "name", "kingdom")
        .join(accepted, on=["name", "kingdom"], how="inner", nulls_equal=True)
        .select("species_id", "name", "kingdom", "accepted_ids")
        .sort("name", "species_id")
    )


def provisional_twin_records(table: pl.DataFrame, twins: pl.DataFrame) -> dict:
    """How many species, rows and records (`n`, `n_raw` and `records` where `table` has them) of `table` sit on the ids of `twins`.

    `ambiguous_species` are those whose name has several accepted twins, which Task 10 does not join at all. Reports only.
    """
    hit = table.join(twins.select("species_id", "accepted_ids"), on="species_id", how="inner")
    return {
        "species": hit["species_id"].n_unique(),
        "rows": hit.height,
        **{c: int(hit[c].sum() or 0) for c in ("n", "n_raw", "records") if c in hit.columns},
        "ambiguous_species": hit.filter(pl.col("accepted_ids") > 1)["species_id"].n_unique(),
    }
