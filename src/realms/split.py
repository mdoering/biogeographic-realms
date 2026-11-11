"""Split species into two evolutionarily independent halves.

Whole families go to one half or the other, so relatives of an endemic family never sit on both sides of the confirmation
test. Half A is for exploring, half B stays untouched until the pipeline is frozen.
"""

import hashlib
from pathlib import Path

import numpy as np
import polars as pl

SPLIT_SEED: int = 20261031
# The canonical hash of the split for COL release 315557, as made by scripts/write_split.py; scripts/carry_split.py carries it to later releases.
SPLIT_SHA256_315557: str = "88f237b47d4964108c161e08acb00c9cf03601418dd2e0f16b49a43704845fa9"
# The split the pipeline loads by default.
SPLIT_SHA256: str = SPLIT_SHA256_315557
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
        halves.update(_place(rows, rng, {"A": 0, "B": 0}))
    return _with_halves(df, halves)


def carry_over(taxa: pl.DataFrame, previous: pl.DataFrame, seed: int = SPLIT_SEED) -> pl.DataFrame:
    """Add `lineage` and `half` to the taxa of a new release, keeping each lineage of `previous` (species_id, lineage, half) in its half.

    A lineage is a family name (or the next rank up), so a species follows its family of the new release, and a species that changed
    family may change half. Lineages that `previous` lacks are placed as `assign_halves` places them, within their group, largest first
    into the half with fewer species, counting the carried-over lineages of that group from the start. The halves explored and held out
    under the earlier release so stay apart (DECISIONS, 2026-11-11, the move to COL 316691).
    """
    known = previous.select("lineage", "half").unique()
    doubled = sorted(known.group_by("lineage").len().filter(pl.col("len") > 1)["lineage"])
    if doubled:
        raise ValueError(f"{len(doubled)} lineages are in both halves of the previous split: {', '.join(doubled[:10])}")
    df = taxa.with_columns(_lineage(taxa.columns))
    groups = _lineage_groups(df).join(known, on="lineage", how="left", maintain_order="left")
    halves = dict(groups.filter(pl.col("half").is_not_null()).select("lineage", "half").iter_rows())

    rng = np.random.default_rng(seed)
    for _, rows in groups.group_by(["main_group"], maintain_order=True):
        new = rows.filter(pl.col("half").is_null())
        if new.height:
            carried = rows.filter(pl.col("half").is_not_null())
            totals = {h: int(carried.filter(pl.col("half") == h)["len"].sum()) for h in ("A", "B")}
            halves.update(_place(new, rng, totals))
    return _with_halves(df, halves)


def crossed(before: pl.DataFrame, after: pl.DataFrame) -> pl.DataFrame:
    """The species in both splits whose half changed: species_id, lineage_before, half_before, lineage, half, sorted by species_id."""
    old = before.select("species_id", pl.col("lineage").alias("lineage_before"), pl.col("half").alias("half_before"))
    both = old.join(after.select("species_id", "lineage", "half"), on="species_id", how="inner")
    return both.filter(pl.col("half") != pl.col("half_before")).sort("species_id")


def id_changes(before: pl.DataFrame, after: pl.DataFrame) -> pl.DataFrame:
    """Species whose id changed between releases: `species_id_before, species_id` for each name (with authorship) found once in each
    release under two different ids, sorted. A name held twice in either release is left out, as it can't be paired."""
    once = lambda df: df.filter(pl.col("name").is_unique()).select("species_id", "name")  # noqa: E731
    pairs = once(before).rename({"species_id": "species_id_before"}).join(once(after), on="name", how="inner")
    return pairs.filter(pl.col("species_id_before") != pl.col("species_id")).select("species_id_before", "species_id").sort("species_id_before")


def _place(rows: pl.DataFrame, rng: np.random.Generator, totals: dict[str, int]) -> dict[str, str]:
    """Each lineage of `rows` (lineage, len) in a half: shuffled, then largest first into whichever half holds fewer species so far."""
    lineages = rows["lineage"].to_numpy(allow_copy=True)
    counts = rows["len"].to_numpy()
    order = rng.permutation(len(lineages))
    largest_first = order[np.argsort(-counts[order], kind="stable")]
    halves = {}
    for i in largest_first:
        if totals["A"] == totals["B"]:
            half = "A" if rng.random() < 0.5 else "B"
        else:
            half = "A" if totals["A"] < totals["B"] else "B"
        halves[lineages[i]] = half
        totals[half] += int(counts[i])
    return halves


def _with_halves(df: pl.DataFrame, halves: dict[str, str]) -> pl.DataFrame:
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
