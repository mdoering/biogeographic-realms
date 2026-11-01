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
  alone, switch the others off (`Settings.rules_off`, `wcvp_mode`). A GRIIS list applies by what it says its entries are: introduced ones (an
  empty value counts), never native or uncertain ones, and not at all a list with no distribution table, unless `griis_without_status` is `apply`.
- **Layer 4** folds a clearly land or clearly sea species' records on the wrong side to the right side. COL's environment is the surer
  word of the two ways to be clear, so where it disagrees with the records `Settings.side_conflict` lets it decide, or leaves the species be.
- **Layer 5**, the rare-outlier heuristics of `realms.outliers`, judges cells and not parts: what one of its rules flags in a species' cell goes with
  every row of the species there, under the rule's name. The rules are `detached_satellite`, `under_reported` and `thin_long_range`; each judges
  the rows layers 1 to 4 left, and a cell that several flag goes with the first of `realms.outliers.RULE_ORDER`.

The records of a row are `records` throughout: the clean records `n - n_human_text` that layer 1 leaves, or `n_raw` without it. The one
exception is the `records` that layer 1 logs for the rows it empties, which is their `n_raw` (see `Cleaned`). Other counts (datasets,
years, licence and source counts) are kept as the download has them, summed as `load_presence` sums them over countries, folds and
parts, with the distinct counts taking the maximum, which is a lower bound.

The summaries of a `Cleaned`, and the report on provisionally accepted twins, are in `realms.clean_report`.
"""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Literal

import polars as pl

from realms.grid import KINDS, ZERO_KM2, Grid, cell_of, part_id, parts
from realms.natives import WCVP_KINGDOM, References
from realms.occurrences import COUNTS, DISTINCT, SUMMED
from realms.outliers import flag_outliers

LAYERS = (1, 2, 3, 4, 5)
WCVP_MODES = ("native_only", "introduced", "off")
GRIIS_WITHOUT_STATUS = ("skip", "apply")
SIDE_CONFLICTS = ("environment", "leave")
# the rules that can be switched off one at a time. WCVP has `wcvp_mode`, layer 1 and layer 4 have their layers, and the fold of
# zero-area parts has no switch.
SWITCHABLE = ("extinct", "no_wild_range", "griis", "griis_subnational", "wrims", "detached_satellite", "under_reported", "thin_long_range")
# Authoritative: every logged row's layer is looked up here from its rule, and a rule that is not in it raises.
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
    "detached_satellite": 5,
    "under_reported": 5,
    "thin_long_range": 5,
}
KEY = ["species_id", "countrycode", "part"]
OTHER = dict(zip(KINDS, reversed(KINDS), strict=True))  # the other part of a cell: land for sea and sea for land
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
# COL's environment is a sorted, comma-separated list ("brackish,freshwater,marine,terrestrial")
FRESH_OR_BRACKISH = r"(^|,)(freshwater|brackish)(,|$)"


@dataclass(frozen=True)
class Settings:
    """The thresholds of layers 4 and 5, at the design's values, and the choices between sources. All are starting values until the freeze.

    The layer 5 thresholds are those of `realms.outliers`: `satellite_*` for `detached_satellite`, `group_min_records`, `under_*` for
    `under_reported`, `thin_*` for `thin_long_range`, and `core_share` and `heuristic_min_records` for the core cells every rule spares (but the
    thin rule). `satellite_gap_km` is added to the grid's median neighbour spacing.

    `land_share` is the share of a species' records on one side that makes it a land or a sea species, whichever it is. `wcvp_mode` is how
    layer 3 uses WCVP (see `clean`): `native_only` is the design's rule, and which mode the analysis takes is decided in M1.
    `griis_without_status` is what layer 3 does with a GRIIS list whose export has no distribution table, so that nothing says which of its
    entries are introduced: `skip` (the rule, and the default) does not apply it, and `apply` applies it name-only, as every list was before
    GRIIS was read by its own establishment means. `apply` is an interim setting until the 250 lists with a plain export have an extended one.
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
    griis_without_status: Literal["skip", "apply"] = "skip"
    side_conflict: Literal["environment", "leave"] = "environment"
    rules_off: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        if self.wcvp_mode not in WCVP_MODES:
            raise ValueError(f"wcvp_mode is one of {WCVP_MODES}, not {self.wcvp_mode!r}")
        if self.griis_without_status not in GRIIS_WITHOUT_STATUS:
            raise ValueError(f"griis_without_status is one of {GRIIS_WITHOUT_STATUS}, not {self.griis_without_status!r}")
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
    `removed`: one row per row of the country counts that a layer took out or moved: `layer` (0 for the zero-area fold, from `RULES`), `rule`,
    `species_id`, `part` (where the row was), `countrycode`, `records`, `to_part` (where a fold put it, null for a removal) and
    `presence_lost`, whether the species had nothing left at that part afterwards.
    **The unit of `removed.records` differs by layer:** it is `n_raw` for layer 1 (everything the emptied row held, flagged records
    included) and the row's `records` for every other layer, which is the clean records `n - n_human_text`, or `n_raw` where layer 1 is
    off. So do not sum it across layer 1 and the rest. Count is conserved from layer 2 on: `presence.records` plus the `records` of the
    removals of layers 2 to 5 (those with a null `to_part`) are the input's `n - n_human_text`.
    `unplaced`: the pairs Task 10 could not place on the grid, for the species that have records: `rule` (the one that would
    have applied them), `species_id`, `country`, `place` (the subdivision or the marine area) and the species' `records`.
    `sides`: how layer 4 sorted the species (`species`, `land`, `sea`, `exempt`, `conflict`, and the two kinds of conflict), empty if
    layer 4 was off.
    """

    presence: pl.DataFrame
    removed: pl.DataFrame
    unplaced: pl.DataFrame = field(default_factory=lambda: pl.DataFrame(schema=UNPLACED_SCHEMA))
    sides: dict[str, int] = field(default_factory=dict)


def _part_table(grid: Grid) -> pl.DataFrame:
    """Every part of the grid: `part`, `cell`, `kind`, `area_km2`, the id and kind of the other part of its cell, and whether each has area."""
    table = parts(grid, [part_id(cell, kind) for cell in grid.cells["cell"] for kind in KINDS])
    other = [part_id(cell, OTHER[kind]) for cell, kind in zip(table["cell"], table["kind"], strict=True)]
    table = table.with_columns(pl.Series("other", other, dtype=pl.String), pl.col("kind").replace_strict(OTHER).alias("other_kind"))
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


def _check_rules(rules: Iterable[str]) -> None:
    unknown = sorted(set(rules) - set(RULES))
    if unknown:
        raise ValueError(f"{unknown} are not rules of RULES, which says which layer a rule is of: {sorted(RULES)}")


class _Run:
    """The rows as the layers have left them, the log of what they took, and what the layers have to report besides."""

    def __init__(self, rows: pl.DataFrame) -> None:
        self.rows = rows
        self.removed: list[pl.DataFrame] = []
        self.unplaced = pl.DataFrame(schema=UNPLACED_SCHEMA)  # set by layer 3
        self.sides: dict[str, int] = {}  # set by layer 4

    def _log(self, hit: pl.DataFrame, records: pl.Expr | None = None) -> None:
        """Log `hit`, rows with `rule` and `to_part`, after `self.rows` no longer has them (where they were taken out).

        The layer comes from `RULES`. `records` is an expression for the logged records where that is not the row's `records`.
        """
        if hit.is_empty():
            return
        _check_rules(hit["rule"].unique())
        pairs = hit.select("species_id", "part").unique()
        alive = pairs.join(self.rows.select("species_id", "part"), on=["species_id", "part"], how="semi")
        lost = pairs.join(alive, on=["species_id", "part"], how="anti").with_columns(pl.lit(True).alias("presence_lost"))
        entry = (
            hit.select(
                pl.col("rule").replace_strict(RULES, return_dtype=pl.Int8).alias("layer"),
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

    def drop(self, hit: pl.DataFrame, rule: str | pl.Expr, records: pl.Expr | None = None, log: pl.Expr | None = None) -> None:
        """Take `hit`, a subset of the rows, out, and log it under `rule` (a name, or an expression for a column of `hit`)."""
        if hit.is_empty():
            return
        rule_column = pl.lit(rule) if isinstance(rule, str) else rule
        logged = hit.with_columns(rule_column.alias("rule"), pl.lit(None, dtype=pl.String).alias("to_part"))
        _check_rules(logged["rule"].unique())  # before anything is changed
        self.rows = self.rows.join(hit.select(KEY), on=KEY, how="anti", nulls_equal=True)
        self._log(logged if log is None else logged.filter(log), records)

    def move(self, moves: pl.DataFrame, on: list[str]) -> None:
        """Move the rows that `moves` names (by `on`) to its `to_part`, or take them out where that is null, and log each under its `rule`.

        A moved row is added to the row of the same species and country at the part it goes to, if there is one. Rows with no records are
        moved, but not logged: a row that counts for nothing is not a presence.
        """
        _check_rules(moves["rule"].unique())
        hit = self.rows.join(moves, on=on, how="inner", nulls_equal=True)
        if hit.is_empty():
            return
        rest = self.rows.join(moves.select(on).unique(), on=on, how="anti", nulls_equal=True)
        moved = hit.filter(pl.col("to_part").is_not_null()).with_columns(pl.col("to_part").alias("part")).select(rest.columns)
        self.rows = _sum_rows(pl.concat([rest, moved]))
        self._log(hit.filter(pl.col("records") > 0))


def _layer_1(run: _Run) -> None:
    """Take out the rows with no clean record. The records of such a row are all lost with it, so it is logged with `n_raw`.

    A row with no record at all is no row, and goes unlogged.
    """
    empty = run.rows.filter(pl.col("records") <= 0)
    rule = pl.when(pl.col("n") == 0).then(pl.lit("failed_flags")).otherwise(pl.lit("human_text"))
    run.drop(empty, rule, records=pl.col("n_raw"), log=pl.col("n_raw") > 0)


def _layer_2(run: _Run, refs: References, settings: Settings) -> None:
    extinct = refs.taxa.filter(pl.col("extinct").fill_null(False)).select("species_id")
    no_wild_range = pl.DataFrame({"species_id": sorted(refs.no_wild_range)}, schema={"species_id": pl.String})
    for rule, species in (("extinct", extinct), ("no_wild_range", no_wild_range)):
        if rule not in settings.rules_off:
            run.drop(run.rows.join(species, on="species_id", how="semi"), rule)


def _griis(refs: References, settings: Settings) -> tuple[pl.DataFrame, pl.DataFrame]:
    """The GRIIS pairs layer 3 applies, and those of them that could not be placed on the grid: the lists with distributions, and under
    `griis_without_status="apply"` the name-only pairs of the lists without."""
    if settings.griis_without_status == "apply":
        return pl.concat([refs.griis, refs.griis_name_only]).unique(), pl.concat([refs.griis_no_cells, refs.griis_name_only_no_cells]).unique()
    return refs.griis, refs.griis_no_cells


def _unplaced(rows: pl.DataFrame, refs: References, settings: Settings) -> pl.DataFrame:
    """The pairs that Task 10 could not place and layer 3 therefore cannot apply, of the species that have records, with their records."""
    found = []
    if "griis_subnational" not in settings.rules_off:
        found.append(
            _griis(refs, settings)[1].select(pl.lit("griis_subnational").alias("rule"), "species_id", "country", pl.col("subdivision").alias("place"))
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


def _layer_3(run: _Run, refs: References, meta: pl.DataFrame, settings: Settings, seen: pl.DataFrame) -> None:
    """Layer 3, and `run.unplaced`: the pairs it could not place. `seen` are the `countrycode` and `cell` of every row of the counts, which says
    where GBIF uses a country code."""
    off = settings.rules_off
    run.unplaced = _unplaced(run.rows, refs, settings)
    places = meta.select("part", "cell", "kind")
    griis, _ = _griis(refs, settings)

    def located() -> pl.DataFrame:
        return run.rows.join(places, on="part", how="left")

    if "griis" not in off:
        # a country's list is for its rows by country code, except in the cells of its sub-national lists, which have lists of their own
        country = griis.filter(pl.col("subdivision").is_null()).select("species_id", pl.col("country").alias("countrycode")).unique()
        covered = refs.admin1_cells.select("cell", pl.col("country").alias("countrycode")).unique()
        hit = located().join(country, on=["species_id", "countrycode"], how="semi").join(covered, on=["cell", "countrycode"], how="anti")
        run.drop(hit, "griis")
    if "griis_subnational" not in off:
        # a sub-national list is for the rows in the cells that overlap its subdivision, and the cells are bigger than the places, so a border cell holds
        # a neighbour country's rows too. Where GBIF uses the list's country code in those cells only the rows with that code count, so a CA row in a cell
        # that overlaps a US state is spared. Where it does not, the cells alone decide: the Canary Islands are `IC` in the list and `ES` in GBIF, and an
        # overseas place carries the sovereign's code in its list (`FR`) and its own in the records (`RE`, `GP`).
        cells = refs.admin1_cells.select("cell", "country", "subdivision")
        uses_code = (
            cells.join(seen.rename({"countrycode": "country"}), on=["cell", "country"])
            .select("country", "subdivision")
            .unique()
            .with_columns(pl.lit(True).alias("by_code"))
        )
        listed = (
            griis.filter(pl.col("subdivision").is_not_null())
            .join(cells, on=["country", "subdivision"])
            .join(uses_code, on=["country", "subdivision"], how="left")
            .with_columns(pl.col("by_code").fill_null(False))
        )
        by_code = listed.filter(pl.col("by_code")).select("species_id", "cell", pl.col("country").alias("countrycode")).unique()
        by_cell = listed.filter(~pl.col("by_code")).select("species_id", "cell").unique()
        rows = located()
        hit = pl.concat(
            [rows.join(by_code, on=["species_id", "cell", "countrycode"], how="semi"), rows.join(by_cell, on=["species_id", "cell"], how="semi")]
        ).unique(KEY, maintain_order=True)
        run.drop(hit, "griis_subnational")
    if "wrims" not in off:
        introduced = refs.wrims.filter(pl.col("country").is_not_null()).select("species_id", pl.col("country").alias("countrycode")).unique()
        run.drop(located().filter(pl.col("kind") == "sea").join(introduced, on=["species_id", "countrycode"], how="semi"), "wrims")
    if settings.wcvp_mode != "off":
        rule = "wcvp_outside_native" if settings.wcvp_mode == "native_only" else "wcvp_introduced"
        run.drop(_wcvp(run.rows, refs, settings.wcvp_mode), rule)


def _sides(rows: pl.DataFrame, taxa: pl.DataFrame, meta: pl.DataFrame, settings: Settings) -> tuple[pl.DataFrame, dict[str, int]]:
    """`species_id` and `side` ("land", "sea" or null) of each species, from COL's environment and the side of its records, and the counts.

    A land species is exactly terrestrial in COL or has at least `land_share` of its records on land parts, and a sea species is the
    same the other way round. Freshwater and brackish species and all of Aves are exempt. Where the two disagree, as for a species that
    is exactly marine with nearly all its records in land parts, `settings.side_conflict` says: COL's environment decides, or the
    species is left alone. The counts are `Cleaned.sides`.
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
    counts = {
        "species": found.height,
        "land": found.filter(pl.col("side") == "land").height,
        "sea": found.filter(pl.col("side") == "sea").height,
        "exempt": found.filter(pl.col("exempt")).height,
        "conflict": both.height,
        "conflict_terrestrial_at_sea": both.filter(pl.col("land_env") & pl.col("sea_records")).height,
        "conflict_marine_on_land": both.filter(pl.col("sea_env") & pl.col("land_records")).height,
    }
    return found.select("species_id", "side"), counts


def _layer_4(run: _Run, refs: References, grid: Grid, meta: pl.DataFrame, settings: Settings) -> None:
    sides, run.sides = _sides(run.rows, refs.taxa, meta, settings)
    sides = sides.drop_nulls("side")
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

    # no such part in the cell: the right-side part of the neighbouring cell that holds most of the species' records, the lowest id on a tie.
    # `held` is what the species has there before this pass moves anything, the same-cell folds of `same` included.
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
    run.move(pl.concat([same, chosen, dropped]), ["species_id", "part"])


def outlier_rows(rows: pl.DataFrame, refs: References, grid: Grid, settings: Settings) -> pl.DataFrame:
    """The rows of `rows` that layer 5 takes out, with the `rule` of each: every row of a flagged species in a flagged cell.

    `rows` is what layers 1 to 4 leave, by species and part or also by country (`Cleaned.presence`, or the rows `clean` holds). Layer 5 is the
    last layer and judges species and cells, so the rows it takes from `Cleaned.presence` are the presences that cleaning with it would
    lose. That lets the calibration run layers 1 to 4 once and layer 5 for every setting.
    """
    flagged = flag_outliers(rows, refs.groups, grid, settings)
    return rows.with_columns(cell_of(pl.col("part")).alias("cell")).join(flagged, on=["species_id", "cell"], how="inner").drop("cell")


def _layer_5(run: _Run, refs: References, grid: Grid, settings: Settings) -> None:
    """Take out what `realms.outliers` flags: every row of a flagged species in a flagged cell, under the rule that flagged it."""
    run.drop(outlier_rows(run.rows, refs, grid, settings), pl.col("rule"))


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

    A country's GRIIS list applies to rows by `countrycode`, except in the cells of its sub-national lists (`griis`). A sub-national list applies to the
    rows in the cells that overlap its subdivision (`griis_subnational`), and where GBIF uses the list's country code in those cells, to the rows with that
    code only, so that a neighbour country's row in a border cell stays. Where the code does not occur among the rows of the cells, as for `IC` (the
    Canary Islands, `ES` in GBIF) or `FR` for Reunion (`RE`), the cells alone decide. Which GRIIS entries apply is `Settings.griis_without_status`
    and what the lists say (`References.griis`). WRiMS applies to sea parts, by `countrycode` (`wrims`). What `refs` could not place is not applied,
    and is listed in `Cleaned.unplaced`.

    Layer 4 folds a land species' row on a sea part, and a sea species' on a land part, into the other part of its cell if that has
    area (`land_sea_fold`), or else into the right-side part of the neighbouring cell that holds most of the species' records, the
    lowest cell id first on a tie (`land_sea_neighbour`), or else out (`land_sea_dropped`). A neighbour counts only if the species has
    records on that part. Where COL's environment and the records disagree about a species' side, `settings.side_conflict` decides.

    Layer 5 takes what `realms.outliers.flag_outliers` flags, with the three rules of that module (see there), each of which `settings.rules_off`
    can leave out. It needs `refs.groups` for the reporting rate, and the grid for distances.
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
    run.move(zero, ["part"])
    if 1 in layers:
        _layer_1(run)
    # a row with no record at all is no row. With layer 1 on that is gone already, and with it off it is a row of `n_raw` 0.
    run.rows = run.rows.filter(pl.col("records") > 0)
    if 2 in layers:
        _layer_2(run, refs, settings)
    if 3 in layers:
        # which country codes GBIF uses in which cells is read from the counts as they came, before layer 1 takes rows out
        seen = counts.select("countrycode", cell_of(pl.col("part")).alias("cell")).drop_nulls("countrycode").unique()
        _layer_3(run, refs, meta, settings, seen)
    if 4 in layers:
        _layer_4(run, refs, grid, meta, settings)
    if 5 in layers:
        _layer_5(run, refs, grid, settings)

    presence = _sum_rows(run.rows, ["species_id", "part"]).sort("species_id", "part")
    removed = pl.concat(run.removed) if run.removed else pl.DataFrame(schema=REMOVED_SCHEMA)
    removed = removed.sort("layer", "rule", "species_id", "part", "countrycode", nulls_last=True)
    return Cleaned(presence, removed, run.unplaced, run.sides)
