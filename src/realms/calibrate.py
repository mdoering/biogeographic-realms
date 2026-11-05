"""The calibration of layer 5 on WCVP plants, and the known cases of the design's Verification section.

**What the calibration measures.** WCVP says where a plant is native, and a plant record outside those regions is either a mistake or an
introduction, which layer 3 would take out. With layer 3 off, layer 5's heuristics have to catch such presences by their own evidence, so
the plants tell how well they do and what they cost. For a table of presences (species and part), `native_status` marks each plant presence
`native`, where its part touches one of the species' native TDWG level 3 regions, or `outside`, where it touches none. Only plants with a
native region get a status, and only land parts with a TDWG region: a sea part, or a land part the scheme leaves out, says nothing.
`rates` are the two numbers that matter: the share of outside presences that the heuristics remove, and the share of native presences they
lose. Both count presences, which is what the network is built from, and not records.

**Three readings of the 5% limit.** The thresholds only change layer 5, so layers 1, 2 and 4 run once (`baseline`) and layer 5 runs for each setting
(`clean.outlier_rows`, the same code that `clean` runs for layer 5). What counts as lost can be read three ways, and `readings` gives all three:

- `layer5`: presences that layer 5 takes, out of those that layers 1, 2 and 4 leave. What the four thresholds do, and nothing else.
- `run_all`: the whole of layers 1, 2, 4 and 5, from the raw presences (`n_raw` above 0), counting every raw presence that is not in the result.
- `run_drops`: the same, but **only drops count**: a presence whose rows a fold moved to another part (layer 0's zero-area fold, or layer 4, a row of the
  removal log with a `to_part`) is moved, not lost, even if it is gone from where it was. The shares are over all the raw presences of the kind.

The plan words the limit over the whole run, and Markus's ruling of 2026-11-01 (reading c) is `run_drops`: the 5% is read over the whole run of layers 1, 2, 4 and 5,
counting removals only, and a moved presence is not lost. `choose` reads it by default, and `calibrate` has its columns as `run_*`. The `layer5` columns keep the first reading.

**The species** are the plants WCVP gives a native region and every other species of their taxonomic groups. Layer 5 judges a species against its
group's records in a cell (`under_reported`), and a group can hold species of other kingdoms (in the data, `eukaryotes` does, and so those species are in the set), so the whole group is kept and the
group totals are what they are in the real run, which also makes a score on half B like for like with half A's. Everything else in layers 1 to 5 is per
species. `kingdom_counts` says how many of a set of species are not Plantae, from the taxonomy alone.

**The setting** is chosen by `choose`, by a knee rule on the convex hull (Markus, 2026-11-01). The settings that lose at most 5% of the native presences are kept, `pareto` builds their
frontier on native presences lost and outside presences removed, and `hull` takes its upper convex hull (a dent of the frontier and collinear points go). `knee` walks the hull from the lowest
native loss: a step to the next vertex is accepted while it costs at most `MAX_MARGINAL` (1, break-even) native presences lost for each extra outside presence removed, and the walk stops at the
first step that costs more. The setting chosen is the last vertex accepted. Where no setting is within 5%, the one that loses least is taken, and the caller is told. The first version of the knee
walked the frontier itself, from one setting to the next, and stopped at a dent of the 144-setting grid; it is in `DECISIONS.md`. `choose(rule="most_removed")` is the rule before the knee, the most
outside presences removed within 5%, with ties to the lower native loss and then to the grid order, kept to report what it would have chosen.

**The animal cross-check against GRIIS** (`griis_cross_check`) sets layer 5 against what GRIIS says is introduced. The heuristics have no range to go by,
so the animals that GRIIS lists as introduced (by its own establishment means, from the lists with distributions) tell how much of that layer 5 finds on its own,
set against what it takes of the other animals. Reports only.

**Known cases** (`CASES`) are checked on the fully cleaned half A with `known_cases`. An area is the cells whose land touches a TDWG level 3 region (or
the regions of a level 1 region), so both parts of a coastal cell count. The cases of a family in half B are not run (the half is looked up in the split
file, which holds no occurrences), and are returned with nothing but the family and the half.
"""

import itertools
import json
import logging
from collections.abc import Collection, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

import polars as pl

from realms.clean import RULES, SWITCHABLE, Settings, clean, outlier_rows
from realms.clean_report import GRIIS_RULES
from realms.grid import Grid, cell_of
from realms.natives import WCVP_KINGDOM, References, match_names

log = logging.getLogger(__name__)

# The design's layer 5 thresholds (docs/design.md), the second value of every threshold of the grid (setting 64). They were the defaults of `Settings` until
# Markus chose setting 132 of the grid (2026-11-01) and it became them. The calibration does not lean on those defaults, so that they can move without the calibration
# moving: the grid sets all four thresholds itself, and what takes a setting by default here is this one.
# Layer 4 as the calibration ran it, before 2026-11-05: 95% of the records, and no class taken as terrestrial (DECISIONS, "Terrestrial unless COL says marine").
DESIGN = Settings(satellite_record_share=0.01, satellite_cell_share=0.05, under_rate_factor=0.01, group_min_records=1000, land_share=0.95, land_classes=frozenset())
NATIVE_LOSS_LIMIT = 0.05
# the knee rule's price: at most this many native presences lost for each extra outside presence removed, at each step along the hull: break-even, one for one
MAX_MARGINAL = 1.0
EPS = 1e-9  # a share compared with a count: 0.05 * 60 is 3.0000000000000004
LEFTOVER = 10  # the regions listed for the presences that break a case
# the thresholds to try, in the order of the grid (the first varies slowest): 3 x 4 x 4 x 3 = 144 settings, the design's value the second of each
GRID = {
    "satellite_record_share": (0.005, 0.01, 0.02),
    "satellite_cell_share": (0.02, 0.05, 0.10, 0.20),
    "under_rate_factor": (0.005, 0.01, 0.015, 0.02),
    "group_min_records": (300, 1000, 3000),
}
RULES_OF_CHOICE = ("hull_knee", "most_removed")
BEFORE_5 = frozenset({1, 2, 4})  # layer 3 is off in the calibration, and layer 5 is what is measured
# the rules that layers 1, 2 and 4 read from `Settings.rules_off`, and that so have to agree across the grid: from `RULES`, which says each rule's layer
EARLY_RULES = frozenset(rule for rule in SWITCHABLE if RULES[rule] in BEFORE_5)
STATUS_SCHEMA = {"species_id": pl.String, "part": pl.String, "status": pl.String}
COUNT_KEYS = ("outside", "outside_removed", "native", "native_lost", "outside_rate", "native_loss")
SCORE = (*COUNT_KEYS, *(f"run_{k}" for k in COUNT_KEYS))
# the columns `choose` reads for each reading: the outside presences removed, the native presences lost and the native presences
READINGS = {"layer5": ("outside_removed", "native_lost", "native"), "run": ("run_outside_removed", "run_native_lost", "run_native")}
TABLE_SCHEMA = {
    "setting": pl.Int64,
    **{k: pl.Float64 if k.endswith(("share", "factor")) else pl.Int64 for k in GRID},
    **{k: pl.Float64 if k.endswith(("rate", "loss")) else pl.Int64 for k in SCORE},
}
REGION_SCHEMA = {"tdwg3": pl.String, "name": pl.String, "level2": pl.Int64, "level1": pl.Int64}


# ---- native and outside presences

def native_status(presence: pl.DataFrame, wcvp: pl.DataFrame, tdwg3_parts: pl.DataFrame, taxa: pl.DataFrame) -> pl.DataFrame:
    """`species_id`, `part` and `status` (`native` or `outside`) of the plant presences in `presence`, a table with `species_id` and `part`.

    A plant is a species of kingdom Plantae in `taxa` (`species_id`, `kingdom`), as layer 3's WCVP rule has it, with at least one native row in `wcvp`
    (`species_id`, `tdwg3`, `native`). Its presence on a land part is `native` when the part touches a native region of the species, as `tdwg3_parts`
    (`part`, `tdwg3`) has it, and `outside` when it touches none: a region where WCVP lists the species as introduced, or does not list it, is not native.
    A part touching a native and an introduced region is native, as layer 3 softens it. Presences with no status are left out: other species (a species
    WCVP lists that COL puts in no kingdom or another one included), sea parts and land parts with no region. Each (species, part) is counted once.
    """
    plants = taxa.filter(pl.col("kingdom") == WCVP_KINGDOM).select("species_id").unique()
    native = wcvp.filter(pl.col("native")).join(plants, on="species_id", how="semi").select("species_id", "tdwg3").unique()
    region = tdwg3_parts.drop_nulls("tdwg3").select("part", "tdwg3")
    pairs = (
        presence.select("species_id", "part")
        .unique()
        .join(native.select("species_id").unique(), on="species_id", how="semi")
        .join(region.select("part").unique(), on="part", how="semi")
    )
    touching = (
        pairs.join(region, on="part")
        .join(native.with_columns(pl.lit(True).alias("is_native")), on=["species_id", "tdwg3"], how="left")
        .group_by("species_id", "part")
        .agg(pl.col("is_native").is_not_null().any().alias("native"))
    )
    status = pl.when(pl.col("native")).then(pl.lit("native")).otherwise(pl.lit("outside")).alias("status")
    return touching.select("species_id", "part", status).cast(STATUS_SCHEMA).sort("species_id", "part")


def tally(status: pl.DataFrame, removed: pl.DataFrame) -> dict[str, int]:
    """How many `outside` and `native` presences of `status` there are, and how many of each `removed` (a table with `species_id` and `part`) holds.

    A presence is a (species, part) and counts once however many rows of `removed` name it. Presences of `removed` with no status are not counted.
    """
    gone = removed.select("species_id", "part").unique().with_columns(pl.lit(True).alias("gone"))
    by = (
        status.join(gone, on=["species_id", "part"], how="left")
        .group_by("status")
        .agg(pl.len().alias("n"), pl.col("gone").is_not_null().sum().alias("lost"))
    )
    found = {row["status"]: (int(row["n"]), int(row["lost"])) for row in by.iter_rows(named=True)}
    outside, outside_removed = found.get("outside", (0, 0))
    native, native_lost = found.get("native", (0, 0))
    return {"outside": outside, "outside_removed": outside_removed, "native": native, "native_lost": native_lost}


def _share(part: int, whole: int) -> float:
    return part / whole if whole else float("nan")


def rates(status: pl.DataFrame, removed: pl.DataFrame) -> tuple[float, float]:
    """The share of the outside presences of `status` that `removed` holds, and the share of the native ones: what the heuristics got and what they cost.

    Not a number (`nan`) where there are no presences of the kind.
    """
    t = tally(status, removed)
    return _share(t["outside_removed"], t["outside"]), _share(t["native_lost"], t["native"])


# ---- the settings

def settings_grid(base: Settings = DESIGN) -> list[Settings]:
    """The 144 settings of `GRID`, with `base` for everything else, in the order of the grid: the first threshold varies slowest. The four thresholds of the
    grid are set for every setting, so `base`'s own do not matter."""
    return [replace(base, **dict(zip(GRID, values, strict=True))) for values in itertools.product(*GRID.values())]


def _reading(reading: str) -> tuple[str, str, str]:
    if reading not in READINGS:
        raise ValueError(f"a reading is one of {sorted(READINGS)}, not {reading!r}")
    return READINGS[reading]


def within_limit(limit: float = NATIVE_LOSS_LIMIT, reading: Literal["layer5", "run"] = "run") -> pl.Expr:
    """Whether a row of `calibrate`'s table loses at most `limit` of the native presences, by `reading` (`layer5`, or the whole `run` counting drops only)."""
    _, lost, native = _reading(reading)
    return pl.col(lost) <= limit * pl.col(native) + EPS


def pareto(table: pl.DataFrame, reading: Literal["layer5", "run"] = "run") -> pl.DataFrame:
    """The Pareto frontier of the settings of `table` on (native presences lost, outside presences removed), by `reading`: the settings that no other beats on both.

    The settings are sorted by native lost, lowest first, then by outside removed, highest first, then by `setting`, and a row is kept if it removes more outside
    presences than every row before it. So a setting is left out if another loses no more native presences and removes at least as many outside ones, and of
    two that tie on both the lower setting stays. The rows come in that order, native lost rising and outside removed rising with it, and are `setting`,
    `outside_removed`, `native_lost` and `native` (by the columns of `reading`).
    """
    removed, lost, native = _reading(reading)
    if table.is_empty():
        raise ValueError("there are no settings to choose from")
    ordered = table.sort([lost, removed, "setting"], descending=[False, True, False])
    best_before = pl.col(removed).cum_max().shift(1).fill_null(-1)
    kept = ordered.filter(pl.col(removed) > best_before)
    return kept.select("setting", pl.col(removed).alias("outside_removed"), pl.col(lost).alias("native_lost"), pl.col(native).alias("native"))


def with_marginals(points: pl.DataFrame) -> pl.DataFrame:
    """`points` (`pareto`'s rows, or a subset of them in the same order) with the step into each row from the one before it: `extra_outside` and `extra_native`
    (what the row adds to the one before) and `marginal`, `extra_native` over `extra_outside`: the native presences it costs for each extra outside presence removed.
    The first row has no step, and `None` in those three.
    """
    stepped = points.with_columns(pl.col("outside_removed").diff().alias("extra_outside"), pl.col("native_lost").diff().alias("extra_native"))
    return stepped.with_columns((pl.col("extra_native") / pl.col("extra_outside")).alias("marginal"))


def frontier(table: pl.DataFrame, reading: Literal["layer5", "run"] = "run") -> pl.DataFrame:
    """`pareto` with the step into each row and its marginal cost (`with_marginals`)."""
    return with_marginals(pareto(table, reading))


def hull(points: pl.DataFrame) -> pl.DataFrame:
    """The upper convex hull of `points` (`pareto`'s rows) in (native presences lost, outside presences removed): the rows to walk, as columns `pareto` has.

    The points are on a frontier, so they rise in both, and the hull is the chain from the first to the last on which each step buys fewer outside presences per
    native one than the step before. A point on or below the segment that joins its neighbours on the hull is left out (a dent of the frontier: a step that is dear,
    after which the cost falls again), and so is a point on that segment (collinear), as it adds no step of its own. The first and the last are always on it.
    The cross products are exact for counts.
    """
    rows = points.to_dicts()
    chain: list[dict] = []
    for point in rows:
        while len(chain) >= 2:
            o, a = chain[-2], chain[-1]
            cross = (a["native_lost"] - o["native_lost"]) * (point["outside_removed"] - o["outside_removed"]) - (a["outside_removed"] - o["outside_removed"]) * (point["native_lost"] - o["native_lost"])
            if cross < 0:  # a right turn: a is above the segment from o to the point, and stays
                break
            chain.pop()
        chain.append(point)
    return pl.DataFrame(chain, schema=points.schema)


@dataclass(frozen=True)
class Knee:
    """What `knee` found, with what it was asked (`limit`, `reading`, `max_marginal`) and the walk along the hull that led to the chosen setting.

    `setting` is the chosen one and `met` whether any setting was within the limit. `steps` is the whole Pareto frontier (`frontier`), with `within_limit` (the setting
    loses at most the limit of the native presences) and `on_hull` (it is a vertex of the hull of the settings within the limit) flags, so that a setting off the hull stays in
    sight. `hull` is that hull in walking order with each step into a vertex (`with_marginals`) and two flags: `affordable` (the step into it costs at most `max_marginal`;
    the first vertex, with no step, is) and `accepted` (the walk got that far, that is the vertex and every vertex before it are affordable: it stops at the first that is not).
    Where no setting is within the limit the hull is empty, and `setting` is the first of the frontier.
    """

    setting: int
    met: bool
    steps: pl.DataFrame
    hull: pl.DataFrame
    limit: float
    reading: str
    max_marginal: float

    def _at(self) -> int | None:
        ids = self.hull["setting"].to_list()
        return ids.index(self.setting) if self.setting in ids else None

    def step_into(self) -> dict | None:
        """The hull vertex of the chosen setting: its counts, and the step into it from the vertex before (`None` in the step's columns if it is the first). `None` if no setting is within the limit."""
        at = self._at()
        return None if at is None else self.hull.row(at, named=True)

    def step_after(self) -> dict | None:
        """The hull vertex after the chosen setting, whose step the walk refused, or `None` if the chosen setting is the last of the hull (or there is no hull)."""
        at = self._at()
        return self.hull.row(at + 1, named=True) if at is not None and at + 1 < self.hull.height else None


def knee(
    table: pl.DataFrame, limit: float = NATIVE_LOSS_LIMIT, reading: Literal["layer5", "run"] = "run", max_marginal: float = MAX_MARGINAL
) -> Knee:
    """Walk the upper convex hull of the Pareto frontier of the settings within `limit` from the lowest native loss, and take the last vertex the walk accepts.

    Only the settings of `table` that lose at most `limit` of the native presences (`within_limit`) are taken, their `pareto` frontier is built (same sort, same ties), and `hull`
    of it. A step to the next vertex is accepted while it costs at most `max_marginal` native presences lost for each extra outside presence removed (counts compared, not the
    ratio), and the walk stops at the first step that costs more. On a hull the steps only get dearer, so there is no cheaper step after the first dear one, and a dent of the
    frontier, which is off the hull, does not stop it. At `max_marginal` 1 a step is accepted while it removes at least as many outside presences as it loses native ones.
    If no setting is within the limit, the one that loses the fewest native presences is taken (the first of the frontier: of those that tie, the one that removes most outside
    presences), and `met` is false.
    """
    points = pareto(table, reading)
    inside = table.filter(within_limit(limit, reading))
    inside_ids = set(inside["setting"].to_list())
    allowed = hull(pareto(inside, reading)) if inside.height else points.head(0)
    walk = with_marginals(allowed)
    affordable = [m is None or n <= max_marginal * o + EPS for m, n, o in zip(walk["marginal"].to_list(), walk["extra_native"].to_list(), walk["extra_outside"].to_list(), strict=True)]
    accepted, going = [], True
    for ok in affordable:
        going = going and ok
        accepted.append(going)
    walk = walk.with_columns(pl.Series("affordable", affordable, pl.Boolean), pl.Series("accepted", accepted, pl.Boolean))
    on_hull = set(walk["setting"].to_list())
    steps = frontier(table, reading).with_columns(
        pl.col("setting").is_in(list(inside_ids)).alias("within_limit"), pl.col("setting").is_in(list(on_hull)).alias("on_hull")
    )
    taken = walk.filter(pl.col("accepted"))
    chosen, met = (int(taken["setting"][-1]), True) if taken.height else (int(points["setting"][0]), False)
    return Knee(chosen, met, steps, walk, limit, reading, max_marginal)


def choose(
    table: pl.DataFrame,
    limit: float = NATIVE_LOSS_LIMIT,
    reading: Literal["layer5", "run"] = "run",
    rule: Literal["hull_knee", "most_removed"] = "hull_knee",
    max_marginal: float = MAX_MARGINAL,
) -> tuple[int, bool]:
    """The `setting` of `table` (from `calibrate`) to take, and whether the `limit` on the native presences lost was met.

    `reading` is which columns count: `run`, the whole of layers 1, 2, 4 and 5 counting drops only (the default), or `layer5`. `rule` is how:

    - `hull_knee` (the default): `knee`, which takes the settings within the `limit`, walks the upper convex hull of their Pareto frontier from the lowest native
      loss while each step costs at most `max_marginal` native presences for each extra outside presence removed (break-even at the default), and takes the last vertex accepted.
    - `most_removed`: the rule of before, kept for comparison: the most outside presences removed among the settings within the `limit`. A tie goes to the fewer
      native presences lost, and then to the lower setting, which is the order of the grid. Where no setting is within the limit, the one that loses
      the fewest native presences is taken, the lower setting again winning a tie.

    The native presences are the same for every setting of one table, so counts are compared and not shares.
    """
    if table.is_empty():
        raise ValueError("there are no settings to choose from")
    if rule not in RULES_OF_CHOICE:
        raise ValueError(f"a rule is one of {list(RULES_OF_CHOICE)}, not {rule!r}")
    removed, lost, _ = _reading(reading)
    if rule == "hull_knee":
        found = knee(table, limit, reading, max_marginal)
        return found.setting, found.met
    within = table.filter(within_limit(limit, reading))
    if within.height:
        best = within.sort([removed, lost, "setting"], descending=[True, False, False]).row(0, named=True)
        return int(best["setting"]), True
    return int(table.sort([lost, "setting"]).row(0, named=True)["setting"]), False


# ---- the species and the calibration

def calibration_species(refs: References) -> pl.Series:
    """The species the calibration cleans: plants with a native region in WCVP, and every other species of their taxonomic groups."""
    plants = refs.wcvp.filter(pl.col("native")).select("species_id").unique()
    groups = refs.groups.drop_nulls("tax_group").select("species_id", "tax_group").unique()
    plant_groups = groups.join(plants, on="species_id", how="semi").select("tax_group").unique()
    mates = groups.join(plant_groups, on="tax_group", how="semi").select("species_id")
    return pl.concat([plants, mates]).unique().sort("species_id")["species_id"]


def calibration_counts(counts: pl.DataFrame, refs: References) -> pl.DataFrame:
    """The rows of `counts` (the country counts) of `calibration_species`."""
    return counts.join(calibration_species(refs).to_frame(), on="species_id", how="semi")


def _early(settings: Settings) -> tuple:
    """What layers 1, 2 and 4 take from the settings, which `calibrate` runs once: `land_share`, `land_classes`, `marine_exceptions`, `side_conflict` and the rules of `EARLY_RULES` that are off."""
    return settings.land_share, settings.land_classes, settings.marine_exceptions, settings.side_conflict, settings.rules_off & EARLY_RULES


@dataclass(frozen=True)
class Baseline:
    """Layers 1, 2 and 4 run on the calibration species with `settings`, and what the readings need of them.

    `presence`: what they leave, rows by species and part. `status`: the `native_status` of those presences, for the `layer5` reading.
    `raw_status`: the status of the raw presences (`n_raw` above 0), the base of the whole-run readings. `absent`: the raw presences (`species_id`, `part`)
    that `presence` does not have, moved ones included. `moved`: the (`species_id`, `part`) whose rows a fold moved to another part (a removal-log row with
    a `to_part`, from layer 0 or 4).
    """

    settings: Settings
    presence: pl.DataFrame
    status: pl.DataFrame
    raw_status: pl.DataFrame
    absent: pl.DataFrame
    moved: pl.DataFrame


def baseline(counts: pl.DataFrame, refs: References, grid: Grid, settings: Settings = DESIGN) -> Baseline:
    """Run layers 1, 2 and 4 on the `calibration_counts` of `counts` with `settings`, and mark the plant presences, raw and left."""
    plants = calibration_counts(counts, refs)
    cleaned = clean(plants, refs, grid, BEFORE_5, settings)
    raw = plants.filter(pl.col("n_raw") > 0).select("species_id", "part").unique()
    absent = raw.join(cleaned.presence.select("species_id", "part"), on=["species_id", "part"], how="anti")
    moved = cleaned.removed.filter(pl.col("to_part").is_not_null()).select("species_id", "part").unique()
    return Baseline(
        settings,
        cleaned.presence,
        native_status(cleaned.presence, refs.wcvp, refs.tdwg3_parts, refs.taxa),
        native_status(raw, refs.wcvp, refs.tdwg3_parts, refs.taxa),
        absent,
        moved,
    )


def _with_rates(t: dict[str, int]) -> dict:
    return {**t, "outside_rate": _share(t["outside_removed"], t["outside"]), "native_loss": _share(t["native_lost"], t["native"])}


def _run_lost(start: Baseline, removed: pl.DataFrame, drops_only: bool) -> pl.DataFrame:
    """The raw presences that the whole run loses: those layers 1, 2 and 4 left no presence of, and those layer 5 takes (`removed`, rows by species and part).

    With `drops_only`, a presence whose rows a fold moved is not lost. The result is a table with `species_id` and `part`, to be counted with `tally`.
    """
    lost = pl.concat([start.absent, removed.select("species_id", "part")]).unique()
    return lost.join(start.moved, on=["species_id", "part"], how="anti") if drops_only else lost


def readings(start: Baseline, refs: References, grid: Grid, setting: Settings) -> dict[str, dict]:
    """The three readings of what `setting` costs, each the counts of `tally` and the two rates: `layer5`, `run_all` and `run_drops` (see the module)."""
    removed = outlier_rows(start.presence, refs, grid, setting)
    return {
        "layer5": _with_rates(tally(start.status, removed)),
        "run_all": _with_rates(tally(start.raw_status, _run_lost(start, removed, False))),
        "run_drops": _with_rates(tally(start.raw_status, _run_lost(start, removed, True))),
    }


def calibrate(counts: pl.DataFrame, settings_grid: Sequence[Settings], refs: References, grid: Grid, start: Baseline | None = None) -> pl.DataFrame:
    """One row per setting of `settings_grid`: its thresholds, and what it does to the plants' native and outside presences.

    Runs on the `calibration_species` of `counts` (country counts, `load_country_counts`), with layer 3 off: layers 1, 2 and 4 once, and layer 5
    at each setting. `setting` is the index in `settings_grid`. `outside` and `native` are the presences of each kind that layer 5 judges, the same in
    every row, `outside_removed` and `native_lost` those it takes, and `outside_rate` and `native_loss` the shares (the `layer5` reading). The `run_*`
    columns are the same for the whole run counting drops only (`run_drops`), over the raw presences. The settings must agree on what layers 1, 2 and 4
    use (`land_share`, `side_conflict` and the rules of `EARLY_RULES`), as those run once. `start` is that run, if the caller has it (`baseline`):
    it must have been made with settings that agree with the grid's on those.
    """
    if not len(settings_grid):
        raise ValueError("there are no settings to calibrate")
    first = settings_grid[0]
    for setting in settings_grid[1:]:
        if _early(setting) != _early(first):
            raise ValueError("the settings differ in land_share, side_conflict or the layer 2 rules of rules_off, which layers 1, 2 and 4 use and which run once")
    if start is None:
        start = baseline(counts, refs, grid, first)
    elif _early(start.settings) != _early(first):
        raise ValueError("the baseline was made with other land_share, side_conflict or layer 2 rules of rules_off than the settings have, so it is not the layers 1, 2 and 4 of this grid")
    first_tally = {"layer5": tally(start.status, start.presence.head(0)), "run": tally(start.raw_status, start.presence.head(0))}
    for name, present in first_tally.items():
        if not present["outside"] or not present["native"]:
            raise ValueError(f"the calibration needs outside and native presences, and the {name} reading has {present['outside']} outside and {present['native']} native")
    table = []
    for i, setting in enumerate(settings_grid):
        removed = outlier_rows(start.presence, refs, grid, setting)
        layer5 = _with_rates(tally(start.status, removed))
        run = _with_rates(tally(start.raw_status, _run_lost(start, removed, True)))
        table.append({"setting": i, **{k: getattr(setting, k) for k in GRID}, **layer5, **{f"run_{k}": v for k, v in run.items()}})
        if (i + 1) % 9 == 0:
            log.info("setting %d of %d", i + 1, len(settings_grid))
    return pl.DataFrame(table, schema=TABLE_SCHEMA)


def score(counts: pl.DataFrame, setting: Settings, refs: References, grid: Grid) -> dict:
    """What one `setting` does to the calibration species of `counts` (the plants and the other species of their groups): the counts and rates of
    `calibrate`, for a setting chosen on other species.

    Only numbers come back, no species and no cell: the keys of `SCORE`.
    """
    row = calibrate(counts, [setting], refs, grid).row(0, named=True)
    return {k: row[k] for k in SCORE}


def kingdom_counts(species: Collection[str], taxa: pl.DataFrame) -> dict[str, int]:
    """How many of `species` there are, and how many of them are Plantae and not, by `taxa` (`species_id`, `kingdom`). A species `taxa` lacks is not Plantae.

    From the taxonomy alone: it reads no occurrences, so it is safe on the species of half B.
    """
    mine = pl.DataFrame({"species_id": list(species)}, schema={"species_id": pl.String}).unique()
    plantae = mine.join(taxa.filter(pl.col("kingdom") == "Plantae").select("species_id").unique(), on="species_id", how="semi").height
    return {"species": mine.height, "plantae": plantae, "not_plantae": mine.height - plantae}


# ---- regions

def tdwg_regions(level3: Path) -> pl.DataFrame:
    """`tdwg3`, `name`, `level2` and `level1` of the TDWG level 3 regions, from the WGSRPD GeoJSON (its properties only)."""
    features = json.loads(Path(level3).read_text())["features"]
    rows = [(p["LEVEL3_COD"], p["LEVEL3_NAM"], p["LEVEL2_COD"], p["LEVEL1_COD"]) for p in (f["properties"] for f in features)]
    return pl.DataFrame(rows, schema=REGION_SCHEMA, orient="row").sort("tdwg3")


def area_cells(tdwg3_parts: pl.DataFrame, regions: pl.DataFrame, area: str) -> set[str]:
    """The cells whose land touches the area: `level1:<n>` for the regions of a TDWG level 1 region, or `tdwg3:<code>` for one level 3 region.

    Both parts of such a cell are in the area, and so is a cell on the border of it, as the cells are large.
    """
    kind, _, value = area.partition(":")
    if kind == "level1":
        codes = regions.filter(pl.col("level1") == int(value))["tdwg3"]
    elif kind == "tdwg3":
        codes = regions.filter(pl.col("tdwg3") == value)["tdwg3"]
    else:
        raise ValueError(f"an area is level1:<n> or tdwg3:<code>, not {area!r}")
    if codes.is_empty():
        raise ValueError(f"no region for the area {area!r}")
    found = tdwg3_parts.filter(pl.col("tdwg3").is_in(codes.implode())).select(cell_of(pl.col("part")).alias("cell")).unique()
    return set(found["cell"].to_list())


def region_loss(status: pl.DataFrame, removed: pl.DataFrame, wcvp: pl.DataFrame, tdwg3_parts: pl.DataFrame, regions: pl.DataFrame) -> pl.DataFrame:
    """The native presences lost per TDWG level 3 region: `tdwg3`, `name`, `level1`, `parts` (land parts touching the region), `native`, `lost` and `loss`.

    A native presence is counted in each region its part touches and where the species is native, so a part on a border counts in both and the rows
    do not add up to the total. Sorted by `loss`, highest first, and then by `native`.
    """
    native_in = wcvp.filter(pl.col("native")).select("species_id", "tdwg3")
    region = tdwg3_parts.drop_nulls("tdwg3").select("part", "tdwg3")
    gone = removed.select("species_id", "part").unique().with_columns(pl.lit(True).alias("gone"))
    per = (
        status.filter(pl.col("status") == "native")
        .join(region, on="part")
        .join(native_in, on=["species_id", "tdwg3"], how="inner")
        .join(gone, on=["species_id", "part"], how="left")
        .group_by("tdwg3")
        .agg(pl.len().cast(pl.Int64).alias("native"), pl.col("gone").is_not_null().sum().cast(pl.Int64).alias("lost"))
    )
    size = region.group_by("tdwg3").agg(pl.col("part").n_unique().cast(pl.Int64).alias("parts"))
    return (
        per.join(size, on="tdwg3", how="left")
        .join(regions.select("tdwg3", "name", "level1"), on="tdwg3", how="left")
        .with_columns((pl.col("lost") / pl.col("native")).alias("loss"))
        .select("tdwg3", "name", "level1", "parts", "native", "lost", "loss")
        .sort(["loss", "native", "tdwg3"], descending=[True, True, False])
    )


def top_regions(losses: pl.DataFrame, min_native: int = 100, n: int = 20) -> pl.DataFrame:
    """The `n` regions of `region_loss` with the highest loss among those with at least `min_native` native presences."""
    return losses.filter(pl.col("native") >= min_native).head(n)


# ---- the cross-check against GRIIS

def _unit_counts(frame: pl.DataFrame) -> dict:
    """How many (species, cell) pairs, (species, part) presences and records `frame` holds, and how many layer 5 took, with the shares."""
    pairs = frame.group_by("species_id", "cell").agg(pl.col("removed").any(), pl.col("records").sum())
    out = {
        "pairs": pairs.height,
        "pairs_removed": int(pairs["removed"].sum() or 0),
        "presences": frame.height,
        "presences_removed": int(frame["removed"].sum() or 0),
        "records": int(frame["records"].sum() or 0),
        "records_removed": int(frame.filter(pl.col("removed"))["records"].sum() or 0),
        "species": frame["species_id"].n_unique(),
    }
    return {**out, **{f"{k}_share": _share(out[f"{k}_removed"], out[k]) for k in ("pairs", "presences", "records")}}


def griis_cross_check(counts: pl.DataFrame, refs: References, grid: Grid, settings: Settings = DESIGN, list_codes: pl.DataFrame | None = None) -> dict:
    """What layer 5 takes of the animals that GRIIS lists as introduced, against what it takes of all the other animals, at the thresholds of `settings`.

    Layers 1, 2 and 4 run on everything `counts` holds, so that layer 5 sees the group totals of the real run, with layer 3 off, and layer 5 is
    `clean.outlier_rows`, the code that `clean` runs. An animal presence (a species on a part) is **listed** where the species' rows in its cell are taken by a
    GRIIS rule of layer 3 (`griis` or `griis_subnational`), that is where GRIIS lists the species as introduced in that country or sub-national place, and
    **other** where not. **The unit that is listed is the (species, cell) pair**, as layer 3's sub-national rule and layer 5 both judge cells: another country's rows in a
    listed species' border cell (a CA row in a cell that overlaps a US state) count as listed with it, though no GRIIS rule took them. The lists are those with a distribution table, whatever `settings.griis_without_status` says, and an entry applies as in `clean`:
    introduced, or with no value, but not native or uncertain. Layer 5 judges cells, so the units are the (species, cell) pairs it judges, the presences (species, part)
    and the records; the shares are those layer 5 takes. A listed presence is one the real run's layer 3 takes first, so the share is how much of
    that the heuristics would find unaided. Reports only: nothing is changed, and nothing of half B is read (the caller passes half A's counts). `list_codes` is the
    frozen table that says which country code each sub-national list applies to (`clean`), or `None` for the counts' own.
    """
    listing = replace(settings, wcvp_mode="off", rules_off=frozenset({"wrims"}), griis_without_status="skip")
    taken = clean(counts, refs, grid, frozenset({1, 2, 3}), listing, list_codes).removed.filter(pl.col("rule").is_in(GRIIS_RULES))
    listed = taken.select("species_id", cell_of(pl.col("part")).alias("cell")).unique().with_columns(pl.lit(True).alias("listed"))
    early = clean(counts, refs, grid, BEFORE_5, settings)
    gone = outlier_rows(early.presence, refs, grid, settings).select("species_id", "part").unique().with_columns(pl.lit(True).alias("removed"))
    animals = refs.taxa.filter(pl.col("kingdom") == "Animalia").select("species_id")
    judged = (
        early.presence.join(animals, on="species_id", how="semi")
        .with_columns(cell_of(pl.col("part")).alias("cell"))
        .join(listed, on=["species_id", "cell"], how="left")
        .join(gone, on=["species_id", "part"], how="left")
        .with_columns(pl.col("listed").fill_null(False), pl.col("removed").fill_null(False))
    )
    return {
        "listed": _unit_counts(judged.filter(pl.col("listed"))),
        "other": _unit_counts(judged.filter(~pl.col("listed"))),
        "animal_presences": judged.height,
    }


# ---- the known cases

@dataclass(frozen=True)
class Case:
    """A known case of the design. `kind` says what is checked, in the area `area` (see `area_cells`) for the first, third and fourth:

    - `absent`: no presence in the area, and some left outside it.
    - `native_only`: no presence outside the species' WCVP native regions, and some inside.
    - `kept`: at most `NATIVE_LOSS_LIMIT` of its presences in the area lost, of at least one.
    - `only_in`: no presence outside the area, and some inside.

    A case is a species by `name` (matched in COL, so a synonym finds its accepted species), or all the species of a `family`.
    """

    key: str
    claim: str
    kind: Literal["absent", "native_only", "kept", "only_in"]
    kingdom: str
    name: str | None = None
    family: str | None = None
    area: str | None = None


CASES = (
    Case("tiger", "no tiger presence in Europe", "absent", "Animalia", name="Panthera tigris", area="level1:1"),
    Case("potato", "potato removed outside its native regions", "native_only", "Plantae", name="Solanum tuberosum"),
    Case("starling", "starling in North America removed by GRIIS", "absent", "Animalia", name="Sturnus vulgaris", area="level1:7"),
    Case("hawaiian_bird", "a Hawaiian endemic bird keeps its range", "kept", "Animalia", name="Drepanis coccinea", area="tdwg3:HAW"),
    Case("lemurs", "Lemuridae stays endemic to Madagascar", "only_in", "Animalia", family="Lemuridae", area="tdwg3:MDG"),
)
CRITERIA = {
    "absent": "no presence in the area, and some left outside it",
    "native_only": "no presence outside the native regions, and some inside",
    "kept": f"at least {1 - NATIVE_LOSS_LIMIT:.0%} of the presences in the area kept, of at least one",
    "only_in": "no presence outside the area, and some inside",
}


def resolve_species(case: Case, index: pl.DataFrame, taxa: pl.DataFrame) -> list[str]:
    """The COL ids of the species of a case: the one `case.name` names in `index` (`name_index`), or those of `case.family` in `taxa`, sorted."""
    if case.name is not None:
        found = match_names(pl.DataFrame({"name": [case.name], "kingdom": [case.kingdom]}), index).row(0, named=True)
        if found["status"] != "matched":
            raise ValueError(f"{case.name} is not found among COL's species ({found['status']})")
        return [found["species_id"]]
    species = taxa.filter((pl.col("family") == case.family) & (pl.col("kingdom") == case.kingdom))["species_id"]
    if species.is_empty():
        raise ValueError(f"there are no species in the family {case.family}")
    return sorted(species.to_list())


def _lost_by_rule(removed: pl.DataFrame, species: Collection[str], zone: pl.Expr | None) -> dict[str, int]:
    """The presences of `species` that cleaning lost or moved, per rule, in `zone` (an expression on `cell`) or anywhere."""
    lost = removed.filter(pl.col("presence_lost"), pl.col("species_id").is_in(list(species)))
    if zone is not None:
        lost = lost.with_columns(cell_of(pl.col("part")).alias("cell")).filter(zone)
    by = lost.select("rule", "species_id", "part").unique().group_by("rule").agg(pl.len().cast(pl.Int64).alias("presences")).sort("rule")
    return dict(by.iter_rows())


def _where(presences: pl.DataFrame, tdwg3_parts: pl.DataFrame, regions: pl.DataFrame, n: int = LEFTOVER) -> list[dict]:
    """The TDWG level 3 regions that `presences` (`species_id`, `part`) touch, with the presences of each, the most first and at most `n`.

    A presence on a part that touches several regions counts in each, and one on a part that touches none (a sea part, mostly) is counted as `(no region)`.
    """
    found = (
        presences.select("species_id", "part")
        .unique()
        .join(tdwg3_parts.drop_nulls("tdwg3").select("part", "tdwg3"), on="part", how="left")
        .join(regions.select("tdwg3", "name"), on="tdwg3", how="left")
        .with_columns(pl.col("name").fill_null("(no region)"))
        .group_by("name")
        .agg(pl.struct("species_id", "part").n_unique().cast(pl.Int64).alias("presences"))
        .sort(["presences", "name"], descending=[True, False])
    )
    return found.head(n).rename({"name": "region"}).to_dicts()


def evaluate_case(
    case: Case,
    species: Collection[str],
    before: pl.DataFrame,
    after: pl.DataFrame,
    removed: pl.DataFrame,
    cells: Collection[str],
    wcvp: pl.DataFrame | None = None,
    tdwg3_parts: pl.DataFrame | None = None,
    regions: pl.DataFrame | None = None,
    taxa: pl.DataFrame | None = None,
) -> dict:
    """Check a case against the presences of its `species` before cleaning (`before`) and after it (`after`), tables with `species_id` and `part`.

    `removed` is `Cleaned.removed`, for what took the presences in the area (outside it for `only_in`, anywhere for `native_only`). `cells` are the area's
    cells (`area_cells`), and `wcvp`, `tdwg3_parts` and `taxa` are needed for `native_only` (`native_status`). The result has `case`, `claim`, `kind`, `criterion`, `species`, the
    presences before and after in all and in the area or outside it or native and outside, `removed_by_rule` and `passed`. With `tdwg3_parts` and `regions`
    it also has `leftover_regions`, where the presences that break an `absent`, `only_in` or `native_only` case are (see `_where`), empty if none do.
    """
    mine = pl.DataFrame({"species_id": list(species)}, schema={"species_id": pl.String})

    def own(frame: pl.DataFrame) -> pl.DataFrame:
        return frame.join(mine, on="species_id", how="semi").select("species_id", "part").unique()

    b, a = own(before), own(after)
    area = list(cells)

    def inside(frame: pl.DataFrame) -> int:
        return frame.filter(cell_of(pl.col("part")).is_in(area)).height

    out: dict = {
        "case": case.key,
        "claim": case.claim,
        "kind": case.kind,
        "criterion": CRITERIA[case.kind],
        "species": len(set(species)),
        "presences_before": b.height,
        "presences_after": a.height,
    }
    in_area = pl.col("cell").is_in(area)
    in_area_part = cell_of(pl.col("part")).is_in(area)
    if case.kind == "native_only":
        if wcvp is None or tdwg3_parts is None or taxa is None:
            raise ValueError("native_only needs wcvp, tdwg3_parts and taxa")
        s_before, s_after = (native_status(frame, wcvp, tdwg3_parts, taxa)["status"].value_counts() for frame in (b, a))
        n_b, n_a = (dict(v.iter_rows()) for v in (s_before, s_after))
        out |= {
            "native_before": n_b.get("native", 0),
            "outside_before": n_b.get("outside", 0),
            "native_after": n_a.get("native", 0),
            "outside_after": n_a.get("outside", 0),
        }
        out["passed"] = out["outside_after"] == 0 and out["native_after"] > 0
        zone = None
        offenders = native_status(a, wcvp, tdwg3_parts, taxa).filter(pl.col("status") == "outside")
    else:
        out |= {"area_before": inside(b), "area_after": inside(a)}
        if case.kind == "absent":
            out["passed"] = out["area_after"] == 0 and a.height - out["area_after"] > 0
            zone = in_area
            offenders = a.filter(in_area_part)
        elif case.kind == "kept":
            kept_enough = out["area_after"] >= (1 - NATIVE_LOSS_LIMIT) * out["area_before"] - EPS
            out["passed"] = out["area_before"] > 0 and kept_enough
            zone = in_area
            offenders = a.head(0)
        else:
            out |= {"outside_before": b.height - out["area_before"], "outside_after": a.height - out["area_after"]}
            out["passed"] = out["outside_after"] == 0 and out["area_after"] > 0
            zone = ~in_area
            offenders = a.filter(~in_area_part)
    out["passed"] = bool(out["passed"])
    out["removed_by_rule"] = _lost_by_rule(removed, species, zone)
    if tdwg3_parts is not None and regions is not None:
        out["leftover_regions"] = _where(offenders, tdwg3_parts, regions)
    return out


def known_cases(cases: Collection[Case], index: pl.DataFrame, split: pl.DataFrame, counts: pl.DataFrame, cleaned, refs: References, regions: pl.DataFrame) -> list[dict]:
    """Run the `cases` whose species are in half A, and list the others as waiting for step 9.

    `split` is `species_id, lineage, half` and is where the half comes from, never the occurrences. `counts` are half A's country counts and `cleaned`
    their `Cleaned` result with every layer on. A case whose family is in half B is returned with its `case`, `claim`, `family`, `half` and `status` and
    nothing else. A case whose species are in both halves is an error, as the split keeps a family whole.
    """
    results = []
    for case in cases:
        species = resolve_species(case, index, refs.taxa)
        held = split.filter(pl.col("species_id").is_in(species))
        halves = sorted(held["half"].unique().to_list())
        if held.height != len(species) or len(halves) != 1:
            raise ValueError(f"{case.key}: the species of the case are in the halves {halves} of the split, with {held.height} of {len(species)} in it, and a case belongs to one half")
        family = case.family or refs.taxa.filter(pl.col("species_id") == species[0])["family"][0]
        head = {"case": case.key, "claim": case.claim, "family": family, "half": halves[0]}
        if halves[0] != "A":
            results.append({**head, "status": "waiting for step 9"})
            continue
        cells = area_cells(refs.tdwg3_parts, regions, case.area) if case.area else set()
        before = counts.filter(pl.col("species_id").is_in(species), pl.col("n_raw") > 0).select("species_id", "part")
        done = evaluate_case(case, species, before, cleaned.presence, cleaned.removed, cells, refs.wcvp, refs.tdwg3_parts, regions, refs.taxa)
        results.append({**head, "status": "run", **{k: v for k, v in done.items() if k not in head}})
    return results
