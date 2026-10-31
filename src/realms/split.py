"""Split species into two evolutionarily independent halves.

Whole families go to one half or the other, so relatives of an endemic family never sit on both sides of the confirmation
test. Half A is for exploring, half B stays untouched until the pipeline is frozen.
"""

import numpy as np
import polars as pl

SPLIT_SEED: int = 20261031


def _lineage() -> pl.Expr:
    return pl.coalesce(
        pl.lit("family:") + pl.col("family"),
        pl.lit("order:") + pl.col("order"),
        pl.lit("class:") + pl.col("class"),
        pl.lit("species:") + pl.col("species_id"),
    ).alias("lineage")


def assign_halves(taxa: pl.DataFrame, seed: int = SPLIT_SEED) -> pl.DataFrame:
    """Add `lineage` and `half` ("A" or "B") to a table with columns species_id, group, family, order, class.

    Within each group, lineages are shuffled and then placed largest first into whichever half holds fewer species, so the
    halves stay balanced by species count and a lineage's half depends only on the seed, never on the input order.
    """
    df = taxa.with_columns(_lineage())
    sizes = df.group_by("group", "lineage").len().sort("group", "lineage")
    rng = np.random.default_rng(seed)

    halves: dict[tuple[str, str], str] = {}
    for (group,), rows in sizes.group_by(["group"], maintain_order=True):
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
            halves[(group, lineages[i])] = half
            totals[half] += int(counts[i])

    assignment = pl.DataFrame(
        [(g, lineage, half) for (g, lineage), half in halves.items()],
        schema={"group": df.schema["group"], "lineage": pl.String, "half": pl.String},
        orient="row",
    )
    return df.join(assignment, on=["group", "lineage"], how="left", maintain_order="left")
