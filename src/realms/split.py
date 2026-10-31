"""Split species into two evolutionarily independent halves.

Whole families go to one half or the other, so relatives of an endemic family never sit on both sides of the confirmation
test. Half A is for exploring, half B stays untouched until the pipeline is frozen.
"""

import numpy as np
import polars as pl

SPLIT_SEED: int = 20261031
OTHER_GROUP = "other"


def _lineage(columns: list[str]) -> pl.Expr:
    """A family, or the next rank up if it has none.

    By name on purpose, unlike the network: COL holds many families twice under two ids, and those must not straddle the halves, while
    a true homonym going along does no harm. Names are qualified by kingdom where known.
    """
    kingdom = pl.col("kingdom").fill_null("?") + pl.lit(":") if "kingdom" in columns else pl.lit("")
    return pl.coalesce(
        pl.lit("family:") + kingdom + pl.col("family"),
        pl.lit("order:") + kingdom + pl.col("order"),
        pl.lit("class:") + kingdom + pl.col("class"),
        pl.lit("species:") + pl.col("species_id"),
    ).alias("lineage")


def _lineage_groups(df: pl.DataFrame) -> pl.DataFrame:
    """Each lineage's size and its most common group, so a family whose species carry several groups stays whole."""
    common = (
        df.filter(pl.col("group").is_not_null())
        .group_by("lineage", "group")
        .len()
        .sort(["lineage", "len", "group"], descending=[False, True, False])
        .unique("lineage", keep="first", maintain_order=True)
        .select("lineage", pl.col("group").alias("main_group"))
    )
    sizes = df.group_by("lineage").len()
    return sizes.join(common, on="lineage", how="left").with_columns(pl.col("main_group").fill_null(OTHER_GROUP)).sort("main_group", "lineage")


def assign_halves(taxa: pl.DataFrame, seed: int = SPLIT_SEED) -> pl.DataFrame:
    """Add `lineage` and `half` ("A" or "B") to a table with columns species_id, group, family, order, class.

    A lineage belongs to its most common group, and species with no group follow their lineage. Within each group, lineages are
    shuffled and then placed largest first into whichever half holds fewer species, so the halves stay balanced by species count
    and a lineage's half depends only on the seed and the table, never on the row order.
    """
    df = taxa.with_columns(_lineage(taxa.columns))
    rng = np.random.default_rng(seed)

    halves: dict[str, str] = {}
    for _, rows in _lineage_groups(df).group_by(["main_group"], maintain_order=True):
        lineages = rows["lineage"].to_numpy(allow_copy=True)
        counts = rows["len"].to_numpy()
        order = rng.permutation(len(lineages))
        largest_first = order[np.argsort(-counts[order], kind="stable")]
        totals = {"A": 0, "B": 0}
        for i in largest_first:
            if totals["A"] == totals["B"]:
                half = "A" if rng.random() < 0.5 else "B"
            else:
                half = "A" if totals["A"] < totals["B"] else "B"
            halves[lineages[i]] = half
            totals[half] += int(counts[i])

    assignment = pl.DataFrame({"lineage": list(halves), "half": list(halves.values())}, schema={"lineage": pl.String, "half": pl.String})
    return df.join(assignment, on="lineage", how="left", maintain_order="left")
