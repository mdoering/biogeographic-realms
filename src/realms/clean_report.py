"""What `realms.clean` did, as tables: the removals per layer, rule, group and country, the provisionally accepted species that share a name, and
what GRIIS took out: how much (`griis_totals`), where it changed (`griis_change`) and which of its removals look like a native species' feral or released forms
(`griis_feral_candidates`).

Reports only. Nothing here changes what `clean` returns.
"""

from collections.abc import Sequence

import polars as pl

from realms.clean import RULES, Cleaned
from realms.grid import Grid, cell_of
from realms.natives import SYNONYM, References

SUMMARY_KEYS = ("tax_group", "countrycode")
UNKNOWN_GROUP = "unknown"
# the rules of `Cleaned.unplaced` and the names their gap goes by in the summary. A pair that could not be placed is not removed.
UNPLACED_RULES = {"griis_subnational": "unplaced_griis", "wrims": "unplaced_wrims"}
GRIIS_RULES = ("griis", "griis_subnational")


def removal_summary(cleaned: Cleaned, refs: References, by: Sequence[str] = SUMMARY_KEYS) -> pl.DataFrame:
    """What was removed per `layer`, `rule` and the columns of `by` (`tax_group`, from `refs.groups`, and `countrycode`).

    `rows` are the rows of the country counts, `records` their records and `presences` the (species, part) pairs that were left with
    nothing by it. Rows of folds are counted as removed from the part they left. A country's presences do not add up over countries,
    as one part's presence can be in several, so use `by=("tax_group",)` or `by=()` for totals.

    **`records` has two units, as in `Cleaned.removed`:** `n_raw` for layer 1 and the clean records `n - n_human_text` for every other
    layer (`n_raw` where layer 1 is off), so do not sum it across layer 1 and the rest. Layer 0 and the folds of layer 4 move records
    rather than lose them.

    `unplaced` is True for what the lists could not place (rules `unplaced_griis` and `unplaced_wrims`), which is not removed: `rows` are the
    pairs of species that have records, `records` the records of those species in all, as an upper bound on what the gap leaves, and
    `countrycode` the country of the list, null for WRiMS. A consumer that wants removals only filters on `~unplaced`. Species with no
    group are `unknown`.
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
        .with_columns(pl.col("layer").cast(pl.Int64), pl.lit(False).alias("unplaced"))
    )
    gap = (
        cleaned.unplaced.join(groups, on="species_id", how="left")
        .with_columns(
            group,
            pl.col("country").alias("countrycode"),
            pl.col("rule").replace_strict(RULES, return_dtype=pl.Int64).alias("layer"),
            pl.col("rule").replace_strict(UNPLACED_RULES).alias("rule"),
        )
        .group_by("layer", "rule", *keys)
        .agg(pl.len().cast(pl.Int64).alias("rows"), pl.lit(None, dtype=pl.Int64).alias("presences"), pl.col("records").sum().cast(pl.Int64))
        .with_columns(pl.lit(True).alias("unplaced"))
    )
    columns = ["layer", "rule", "unplaced", *keys, "rows", "presences", "records"]
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


# ---- GRIIS

def griis_totals(removed: pl.DataFrame) -> dict:
    """What the GRIIS rules took, from `Cleaned.removed`: `rows`, `presences` (the (species, part) pairs they left with nothing), `records` and `species`,
    and the same for each rule. The two rules never take the same row, so rows and records add up, and presences are counted once across them."""
    taken = removed.filter(pl.col("rule").is_in(GRIIS_RULES))

    def count(frame: pl.DataFrame) -> dict:
        return {
            "rows": frame.height,
            "presences": frame.filter(pl.col("presence_lost")).select("species_id", "part").unique().height,
            "records": int(frame["records"].sum() or 0),
            "species": frame["species_id"].n_unique(),
        }

    return {**count(taken), **{rule: count(taken.filter(pl.col("rule") == rule)) for rule in GRIIS_RULES}}


def griis_change(before: pl.DataFrame, after: pl.DataFrame, n: int = 20) -> pl.DataFrame:
    """The `n` species whose GRIIS removals changed most in records between two runs (`Cleaned.removed` of each): `species_id`, `before`, `after`
    and `change` (`after - before`), the largest change in either direction first."""
    def per_species(removed: pl.DataFrame, name: str) -> pl.DataFrame:
        return removed.filter(pl.col("rule").is_in(GRIIS_RULES)).group_by("species_id").agg(pl.col("records").sum().alias(name))

    both = per_species(before, "before").join(per_species(after, "after"), on="species_id", how="full", coalesce=True).fill_null(0)
    return (
        both.with_columns((pl.col("after") - pl.col("before")).alias("change"))
        .sort(pl.col("change").abs(), "species_id", descending=[True, False])
        .head(n)
        .select(pl.col("species_id"), *(pl.col(c).cast(pl.Int64) for c in ("before", "after", "change")))
    )


def neighbour_countries(rows: pl.DataFrame, grid: Grid) -> pl.DataFrame:
    """`country` and `neighbour`: the pairs of countries whose rows (`countrycode` and `part`) are in the same cell or in adjacent cells of `grid`.

    Defined from the data and not from borders: two countries are neighbours where GBIF has records of both in cells that touch, which takes in
    sea neighbours (the countries of one strait) and overseas places next to another country's waters. A country is not its own neighbour.
    """
    where = rows.select("countrycode", cell_of(pl.col("part")).alias("cell")).drop_nulls("countrycode").unique()
    near = [(cell, other) for cell, others in grid.adjacency.items() for other in others] + [(cell, cell) for cell in grid.adjacency]
    touching = pl.DataFrame(near, schema={"cell": pl.String, "near": pl.String}, orient="row")
    return (
        where.join(touching, on="cell")
        .join(where.rename({"countrycode": "neighbour", "cell": "near"}), on="near")
        .filter(pl.col("countrycode") != pl.col("neighbour"))
        .select(pl.col("countrycode").alias("country"), "neighbour")
        .unique()
        .sort("country", "neighbour")
    )


def griis_feral_candidates(counts: pl.DataFrame, removed: pl.DataFrame, grid: Grid, n: int = 20, min_neighbours: int = 3, min_records: int = 1) -> pl.DataFrame:
    """The `n` largest GRIIS removals, by records, of a species in a country where the species also has clean, unlisted records in at least
    `min_neighbours` neighbouring countries (`neighbour_countries`, from all the rows of `counts`), of at least `min_records` records each (1 is any).
    Report only: GRIIS calls these introduced, so no rule changes.

    These are the likely feral or released forms of native species (the Mallard in Denmark and Sweden, the Rock Dove in India), which GRIIS lists
    wherever some population is introduced. `counts` are the country counts and `removed` is `Cleaned.removed`. A record is clean where `n - n_human_text`
    is above 0, and unlisted where no GRIIS rule took its row. The result has `species_id`, `country`, `records` and `rows` (what the GRIIS rules took
    of the species there), `neighbours` (how many neighbouring countries hold its unlisted records), `neighbour_records` and `neighbour_countries`
    (`CC:records`, the largest first). It is sorted by `records`, largest first.
    """
    clean_rows = counts.select("species_id", "countrycode", "part", (pl.col("n") - pl.col("n_human_text")).alias("records")).filter(pl.col("records") > 0).drop_nulls("countrycode")
    taken = removed.filter(pl.col("rule").is_in(GRIIS_RULES)).drop_nulls("countrycode")
    listed = taken.group_by("species_id", "countrycode").agg(pl.col("records").sum().alias("records"), pl.len().alias("rows"))
    unlisted = (
        clean_rows.join(taken.select("species_id", "countrycode", "part"), on=["species_id", "countrycode", "part"], how="anti")
        .group_by("species_id", pl.col("countrycode").alias("neighbour"))
        .agg(pl.col("records").sum().alias("neighbour_records"))
        .filter(pl.col("neighbour_records") >= min_records)
    )
    found = (
        listed.join(neighbour_countries(counts, grid).rename({"country": "countrycode"}), on="countrycode")
        .join(unlisted, on=["species_id", "neighbour"])
        .sort("neighbour_records", "neighbour", descending=[True, False])
        .group_by("species_id", "countrycode", "records", "rows", maintain_order=True)
        .agg(
            pl.col("neighbour").n_unique().alias("neighbours"),
            pl.col("neighbour_records").sum().alias("neighbour_records"),
            pl.format("{}:{}", pl.col("neighbour"), pl.col("neighbour_records")).alias("neighbour_countries"),
        )
        .filter(pl.col("neighbours") >= min_neighbours)
    )
    return (
        found.rename({"countrycode": "country"})
        .sort("records", "species_id", "country", descending=[True, False, False])
        .head(n)
        .select("species_id", "country", "records", "rows", "neighbours", "neighbour_records", "neighbour_countries")
    )
