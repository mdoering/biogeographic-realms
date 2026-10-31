"""Freeze the family split for the COL release, before anything is clustered on GBIF data.

Writes data/split/halves-<release>.parquet (species_id, lineage, half) and prints its canonical content hash, which
`SPLIT_SHA256` in src/realms/split.py must carry.
"""

import json
from pathlib import Path

import polars as pl

from realms.reference import load_tax_groups, load_taxa
from realms.split import SPLIT_SEED, SPLIT_SHA256, assign_halves, split_sha256

RELEASE = 315557
ROOT = Path("data/reference")
OUT = Path("data/split") / f"halves-{RELEASE}.parquet"


def main() -> None:
    manifest = json.loads((ROOT / "manifest.json").read_text())
    for key in ("col", "col_groups"):
        if manifest[key]["dataset"] != RELEASE:
            raise SystemExit(f"manifest {key} is dataset {manifest[key]['dataset']}, not release {RELEASE}")

    taxa = load_taxa(Path(manifest["col"]["path"]))
    groups = load_tax_groups(Path(manifest["col_groups"]["path"])).rename({"tax_group": "group"})
    taxa = taxa.join(groups, on="species_id", how="left", maintain_order="left")
    if taxa["species_id"].n_unique() != taxa.height:
        raise SystemExit("species_id is not unique after joining the groups")

    split = assign_halves(taxa, SPLIT_SEED)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    split.select("species_id", "lineage", "half").sort("species_id").write_parquet(OUT)

    digest = split_sha256(split)
    print(f"release {RELEASE}, seed {SPLIT_SEED}: {split.height} species, {split['lineage'].n_unique()} lineages -> {OUT}")
    print(f"canonical sha256 {digest}")
    print("SPLIT_SHA256 in split.py " + ("matches" if digest == SPLIT_SHA256 else "differs: set it to the hash above"))

    per_half = split.group_by("half").agg(pl.len().alias("species"), pl.col("lineage").n_unique().alias("lineages")).sort("half")
    per_group = (
        split.with_columns(pl.col("group").fill_null("(none)"))
        .group_by("group", "half")
        .agg(pl.len().alias("species"), pl.col("lineage").n_unique().alias("lineages"))
        .sort("group", "half")
    )
    with pl.Config(tbl_rows=200, tbl_hide_dataframe_shape=True, tbl_hide_column_data_types=True):
        print(per_half)
        print(per_group)


if __name__ == "__main__":
    main()
