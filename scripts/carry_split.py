"""Carry the half A/B split of COL release 315557 over to the release in the manifest (Markus, 2026-11-11; DECISIONS, the move to COL 316691).

A fresh `write_split.py` on a new release would shuffle and balance a different list of families, so many families explored in half A would land in
half B. Here each lineage (a family name, or the next rank up) keeps its half from 315557, and only lineages new in the release are placed, as
`split.carry_over` says. This reads taxonomy only, never an occurrence.

Writes data/split/halves-<release>.parquet (species_id, lineage, half) and data/split/crossed-315557-<release>.parquet (the species of both releases
whose half changed, with their lineage before and after). Species ids churn between releases (16,548 accepted species ids of 315192 were gone in
315557), so it also pairs the ids of names found once in each release (`split.id_changes`) and counts the renamed species that changed half too. It
prints the split's canonical hash, which `SPLIT_SHA256` in src/realms/split.py must carry
once the pipeline moves to the release.

Usage: uv run python -I scripts/carry_split.py
"""

import io
import json
import zipfile
from pathlib import Path

import polars as pl

from realms.reference import ACCEPTED, load_tax_groups, load_taxa
from realms.split import SPLIT_SEED, SPLIT_SHA256_315557, carry_over, crossed, id_changes, load_split, split_sha256

BEFORE = 315557
ROOT = Path("data/reference")
SPLITS = Path("data/split")
# the COL export of 315557, kept with the rest of its reference data when the manifest moved on
BEFORE_COL = Path(f"data/reference-{BEFORE}/exports/extended/col.zip")


def species_names(coldp_zip: Path) -> pl.DataFrame:
    """The accepted species of a COLDP export as `species_id, name`, the name with its authorship."""
    with zipfile.ZipFile(coldp_zip) as archive:
        member = next(n for n in archive.namelist() if n.rsplit("/", 1)[-1] == "NameUsage.tsv")
        content = archive.read(member)
    columns = ["col:ID", "col:scientificName", "col:authorship", "col:rank", "col:status"]
    usages = pl.read_csv(io.BytesIO(content), separator="\t", quote_char=None, columns=columns, schema_overrides={c: pl.String for c in columns}, null_values=[""])
    return usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED)).select(
        pl.col("col:ID").alias("species_id"), pl.concat_str("col:scientificName", "col:authorship", separator=" ", ignore_nulls=True).alias("name")
    )


def main() -> None:
    manifest = json.loads((ROOT / "manifest.json").read_text())
    release = manifest["col"]["dataset"]
    if manifest["col_groups"]["dataset"] != release:
        raise SystemExit(f"manifest col is release {release} but col_groups is {manifest['col_groups']['dataset']}")
    if release == BEFORE:
        raise SystemExit(f"the manifest is still on release {BEFORE}: nothing to carry over")

    before = load_split(SPLITS / f"halves-{BEFORE}.parquet", sha256=SPLIT_SHA256_315557)
    taxa = load_taxa(Path(manifest["col"]["path"]))
    groups = load_tax_groups(Path(manifest["col_groups"]["path"])).rename({"tax_group": "group"})
    taxa = taxa.join(groups, on="species_id", how="left", maintain_order="left")
    if taxa["species_id"].n_unique() != taxa.height:
        raise SystemExit("species_id is not unique after joining the groups")

    split = carry_over(taxa, before.select("species_id", "lineage", "half"), SPLIT_SEED)
    out = SPLITS / f"halves-{release}.parquet"
    split.select("species_id", "lineage", "half").sort("species_id").write_parquet(out)
    moved = crossed(before, split)
    moved_out = SPLITS / f"crossed-{BEFORE}-{release}.parquet"
    moved.write_parquet(moved_out)

    known = set(before["lineage"])
    lineages = split.group_by("lineage").agg(pl.len().alias("species"), pl.first("half")).with_columns(pl.col("lineage").is_in(known).alias("carried"))
    shared = split.join(before.select("species_id"), on="species_id", how="semi").height
    print(f"release {BEFORE} -> {release}, seed {SPLIT_SEED}: {split.height} species in {lineages.height} lineages -> {out}")
    print(f"species: {shared} in both releases, {split.height - shared} new, {before.height - shared} of {BEFORE} gone")
    for carried, label in ((True, f"lineages of {BEFORE}"), (False, "new lineages")):
        part = lineages.filter(pl.col("carried") == carried)
        print(f"{label}: {part.height}, {part['species'].sum()} species")
    print(f"species whose half changed: {moved.height} ({moved.height / max(shared, 1):.3%} of those in both) -> {moved_out}")
    renamed = id_changes(species_names(BEFORE_COL), species_names(Path(manifest["col"]["path"])))
    halves_before = before.select(pl.col("species_id").alias("species_id_before"), pl.col("half").alias("half_before"))
    renamed = renamed.join(halves_before, on="species_id_before").join(split.select("species_id", "half"), on="species_id")
    print(f"species under a new id, paired by name: {renamed.height}, of which {renamed.filter(pl.col('half') != pl.col('half_before')).height} changed half")

    digest = split_sha256(split)
    print(f"canonical sha256 {digest}")

    with pl.Config(tbl_rows=200, tbl_hide_dataframe_shape=True, tbl_hide_column_data_types=True):
        print(split.group_by("half").agg(pl.len().alias("species"), pl.col("lineage").n_unique().alias("lineages")).sort("half"))
        per_group = (
            split.with_columns(pl.col("group").fill_null("(none)"))
            .group_by("group", "half")
            .len()
            .pivot(on="half", index="group", values="len")
            .pipe(lambda d: d.with_columns(pl.lit(0).alias(h) for h in ("A", "B") if h not in d.columns))
            .fill_null(0)
            .with_columns((pl.col("A") / (pl.col("A") + pl.col("B"))).round(3).alias("share_A"))
            .sort("group")
        )
        print(per_group)
        print(moved.group_by("half_before", "half").len().sort("half_before", "half"))
        print(moved.group_by("lineage_before", "lineage").len().sort("len", descending=True).head(20))


if __name__ == "__main__":
    main()
