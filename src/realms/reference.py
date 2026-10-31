"""Reference data from ChecklistBank. Reads existing exports only, and never asks for a new one."""

import io
import re
import urllib.parse
import zipfile
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from realms.gbif import CLB_API, get_json

SPECIES_COLUMNS = {
    "col:ID": "species_id",
    "col:genus": "genus",
    "col:family": "family",
    "col:order": "order",
    "col:class": "class",
    "col:phylum": "phylum",
    "col:kingdom": "kingdom",
    "col:environment": "environment",
    "col:extinct": "extinct",
}
USED_COLUMNS = ["col:status", "col:rank", *SPECIES_COLUMNS]
ACCEPTED = ("accepted", "provisionally accepted")
GRIIS_ALIAS = re.compile(r"^GRIIS-([A-Z]{2})(?:-(.+))?$")


@dataclass(frozen=True)
class ExportInfo:
    key: str
    url: str
    size: int
    created: str


@dataclass(frozen=True)
class GriisDataset:
    key: int
    alias: str
    country: str
    subdivision: str | None


def latest_export(dataset_key: int, *, classification: bool = False, extended: bool = False, api: str = CLB_API) -> ExportInfo | None:
    """The newest finished COLDP export with exactly these flags, or None. This lists exports and creates nothing."""
    _, body = get_json(f"{api}/export?datasetKey={dataset_key}&format=COLDP&limit=50")
    matching = [
        e
        for e in (body or {}).get("result", [])
        if e.get("status") == "finished"
        and (e.get("request") or {}).get("synonyms")
        and bool((e.get("request") or {}).get("classification")) == classification
        and bool((e.get("request") or {}).get("extended")) == extended
    ]
    if not matching:
        return None
    newest = max(matching, key=lambda e: e["created"])
    return ExportInfo(newest["key"], newest["download"], newest["size"], newest["created"])


def griis_datasets(api: str = CLB_API) -> list[GriisDataset]:
    """Country and sub-national GRIIS checklists, without the separate lists for protected areas."""
    query = urllib.parse.urlencode({"q": "GRIIS", "limit": 1000})
    _, body = get_json(f"{api}/dataset?{query}")
    found = []
    for dataset in (body or {}).get("result", []):
        match = GRIIS_ALIAS.match(dataset.get("alias") or "")
        if match and "Protected Areas" not in (dataset.get("title") or ""):
            found.append(GriisDataset(dataset["key"], dataset["alias"], match.group(1), match.group(2)))
    return found


def load_taxa(coldp_zip: Path) -> pl.DataFrame:
    """Accepted species from a COLDP export with a flat classification."""
    with zipfile.ZipFile(coldp_zip) as archive:
        name = next((n for n in archive.namelist() if n.rsplit("/", 1)[-1] == "NameUsage.tsv"), None)
        if name is None:
            raise ValueError(f"no NameUsage.tsv in {coldp_zip}")
        with archive.open(name) as f:
            content = f.read()
    header = content[: content.index(b"\n")].decode().rstrip("\r").split("\t")
    missing = [c for c in USED_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"missing {missing} in {coldp_zip}, whose columns are {header}")
    usages = pl.read_csv(
        io.BytesIO(content),
        separator="\t",
        quote_char=None,
        columns=USED_COLUMNS,
        schema_overrides={c: pl.String for c in USED_COLUMNS},
        null_values=[""],
    )
    return (
        usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED))
        .select(pl.col(c).alias(n) for c, n in SPECIES_COLUMNS.items())
        .with_columns((pl.col("extinct") == "true").fill_null(False).alias("extinct"))
    )


def load_tax_groups(classification_zip: Path) -> pl.DataFrame:
    """Accepted species with their `clb:taxGroup`, from a classification-only export (a single `dataset-<key>.tsv`)."""
    with zipfile.ZipFile(classification_zip) as archive:
        name = next((n for n in archive.namelist() if n.endswith(".tsv")), None)
        if name is None:
            raise ValueError(f"no tsv file in {classification_zip}")
        with archive.open(name) as f:
            content = f.read()
    columns = ["col:ID", "col:status", "col:rank", "clb:taxGroup"]
    usages = pl.read_csv(
        io.BytesIO(content),
        separator="\t",
        quote_char=None,
        columns=columns,
        schema_overrides={c: pl.String for c in columns},
        null_values=[""],
    )
    return usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED)).select(
        pl.col("col:ID").alias("species_id"), pl.col("clb:taxGroup").alias("tax_group")
    )


def _counts(taxa: pl.DataFrame) -> dict:
    return {
        "species": taxa.height,
        "with_genus": taxa["genus"].drop_nulls().len(),
        "with_family": taxa["family"].drop_nulls().len(),
    }


def rank_coverage(taxa: pl.DataFrame) -> dict:
    """How many species have the ranks the network needs, overall and per kingdom."""
    per_kingdom = {str(k): _counts(rows) for (k,), rows in taxa.group_by(["kingdom"])}
    return {
        **_counts(taxa),
        "with_environment": taxa["environment"].drop_nulls().len(),
        "extinct": int(taxa["extinct"].sum()),
        "per_kingdom": dict(sorted(per_kingdom.items())),
    }
