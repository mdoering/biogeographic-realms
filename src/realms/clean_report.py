"""What `realms.clean` did, as tables: the removals per layer, rule, group and country, and the provisionally accepted species that share a name.

Reports only. Nothing here changes what `clean` returns.
"""

from collections.abc import Sequence

import polars as pl

from realms.clean import RULES, Cleaned
from realms.natives import SYNONYM, References

SUMMARY_KEYS = ("tax_group", "countrycode")
UNKNOWN_GROUP = "unknown"
# the rules of `Cleaned.unplaced` and the names their gap goes by in the summary. A pair that could not be placed is not removed.
UNPLACED_RULES = {"griis_subnational": "unplaced_griis", "wrims": "unplaced_wrims"}


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
