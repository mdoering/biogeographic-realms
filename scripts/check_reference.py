"""Measure what the reference data can support: rank coverage in COL, and whether WCVP's empty value means native."""

import io
import json
import zipfile
from pathlib import Path

import numpy as np
import polars as pl

from realms.reference import load_tax_groups, load_taxa, rank_coverage, repeated_names

ROOT = Path("data/reference")
OUT = Path("results/coverage.json")
ANCHORS = ("Solanum tuberosum", "Zea mays", "Eucalyptus globulus", "Acer pseudoplatanus", "Lantana camara", "Solanum lycopersicum")


def read_tsv(archive: zipfile.ZipFile, name: str) -> pl.DataFrame:
    with archive.open(name) as f:
        return pl.read_csv(io.BytesIO(f.read()), separator="\t", quote_char=None, infer_schema_length=0, null_values=[""], truncate_ragged_lines=True)


def main() -> None:
    manifest = json.loads((ROOT / "manifest.json").read_text())
    taxa = load_taxa(Path(manifest["col"]["path"]))
    coverage = rank_coverage(taxa)
    coverage["release"] = manifest["col"]["release"]
    coverage["repeated_names"] = {rank: repeated_names(taxa, rank) for rank in ("genus", "family")}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(coverage, indent=2))

    n = coverage["species"]
    print(f"release {coverage['release']}: {n} accepted species")
    for key in ("with_genus", "with_family", "with_environment", "extinct"):
        print(f"  {key:18} {coverage[key]:>9}  {coverage[key] / n:.1%}")
    for kingdom, c in coverage["per_kingdom"].items():
        print(f"  {kingdom:18} {c['species']:>9}  genus {c['with_genus'] / max(c['species'], 1):.1%}  family {c['with_family'] / max(c['species'], 1):.1%}")
    print(f"names held by several ids within a kingdom: {coverage['repeated_names']}")

    groups = load_tax_groups(Path(manifest["col_groups"]["path"])).join(taxa.select("species_id"), on="species_id", how="right")
    print(f"\nspecies with a taxGroup: {groups['tax_group'].drop_nulls().len() / n:.1%}")
    with pl.Config(tbl_rows=40, tbl_hide_dataframe_shape=True, tbl_hide_column_data_types=True):
        print(groups.group_by("tax_group").len().sort("len", descending=True))

    with zipfile.ZipFile(manifest["wcvp"]["path"]) as archive:
        names = {n.rsplit("/", 1)[-1]: n for n in archive.namelist()}
        usages = read_tsv(archive, names["NameUsage.tsv"]).select("col:ID", "col:scientificName", "col:rank", "col:status")
        distribution = read_tsv(archive, names["Distribution.tsv"]).select("col:taxonID", "col:areaID", "col:gazetteer", "col:establishmentMeans")
    print("\nWCVP distribution values")
    print(distribution.group_by("col:gazetteer", "col:establishmentMeans").len().sort("len", descending=True))

    species = usages.filter(pl.col("col:rank") == "species", pl.col("col:status") == "accepted")
    rows = distribution.join(species, left_on="col:taxonID", right_on="col:ID").filter(pl.col("col:gazetteer") == "tdwg")

    def describe(name: str, group: pl.DataFrame) -> str:
        native = group.filter(pl.col("col:establishmentMeans").is_null())["col:areaID"].sort().to_list()
        introduced = group.filter(pl.col("col:establishmentMeans") == "introduced")["col:areaID"].sort().to_list()
        return f"{name[:38]:38} unflagged {len(native):>3} {','.join(native[:9]):42} introduced {len(introduced):>3} {','.join(introduced[:7])}"

    print("\nKnown cases")
    for anchor in ANCHORS:
        group = rows.filter(pl.col("col:scientificName") == anchor)
        print(describe(anchor, group))

    both = rows.group_by("col:scientificName").agg(
        pl.col("col:establishmentMeans").is_null().sum().alias("native"), (pl.col("col:establishmentMeans") == "introduced").sum().alias("introduced")
    ).filter(pl.col("native") > 0, pl.col("introduced") > 0)
    print(f"\n{both.height} of {rows['col:scientificName'].n_unique()} species have both unflagged and introduced regions. A random 50:")
    picks = sorted(both["col:scientificName"].to_list())
    for i in np.random.default_rng(3).choice(len(picks), size=50, replace=False):
        print(describe(picks[i], rows.filter(pl.col("col:scientificName") == picks[i])))


if __name__ == "__main__":
    main()
