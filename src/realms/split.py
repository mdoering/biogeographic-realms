"""Split species into two evolutionarily independent halves.

Whole families go to one half or the other, so relatives of an endemic family never sit on both sides of the confirmation
test. Half A is for exploring, half B stays untouched until the pipeline is frozen.
"""

import hashlib
from pathlib import Path

import numpy as np
import polars as pl

SPLIT_SEED: int = 20261031
# The canonical hash of the split for COL release 315557, as made by scripts/write_split.py.
SPLIT_SHA256: str = "88f237b47d4964108c161e08acb00c9cf03601418dd2e0f16b49a43704845fa9"
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


def split_sha256(split: pl.DataFrame) -> str:
    """The sha256 of the split's content: the rows `species_id TAB lineage TAB half LF`, sorted by species_id, as UTF-8.

    The content is hashed, not the Parquet file, whose bytes carry the writer's version and would change with a Polars upgrade.
    """
    rows = split.select("species_id", "lineage", "half")
    if any(rows.null_count().row(0)):
        raise ValueError("the split has null values in species_id, lineage or half")
    lines = pl.concat_str("species_id", pl.lit("\t"), "lineage", pl.lit("\t"), "half", pl.lit("\n"))
    content = rows.sort("species_id").select(lines.str.join("")).item()
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def load_split(path: Path, sha256: str = SPLIT_SHA256) -> pl.DataFrame:
    """Read the frozen split (species_id, lineage, half), refusing a file whose content hash differs from `sha256`."""
    split = pl.read_parquet(path)
    found = split_sha256(split)
    if found != sha256:
        raise ValueError(f"{path} has content hash {found}, not the frozen {sha256}")
    return split
