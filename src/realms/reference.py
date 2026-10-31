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
USED_COLUMNS = ["col:parentID", "col:status", "col:rank", *SPECIES_COLUMNS]
ACCEPTED = ("accepted", "provisionally accepted")
ID_RANKS = ("genus", "family")
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


def _ancestor_ids(usages: pl.DataFrame, species: list[str]) -> pl.DataFrame:
    """The id of each species' genus and family, walking up `parentID`, which is the only place the export keeps them."""
    parent = dict(zip(usages["col:ID"].to_list(), usages["col:parentID"].to_list()))
    rank = dict(zip(usages["col:ID"].to_list(), usages["col:rank"].to_list()))
    status = dict(zip(usages["col:ID"].to_list(), usages["col:status"].to_list()))
    found: dict[str, list[str | None]] = {r: [] for r in ID_RANKS}
    for species_id in species:
        ancestors = {}
        usage = parent.get(species_id)
        while usage is not None:
            if rank.get(usage) in ID_RANKS:
                if status[usage] not in ACCEPTED:
                    raise ValueError(f"accepted species {species_id} sits under {status[usage]} {rank[usage]} {usage}")
                ancestors.setdefault(rank[usage], usage)
            usage = parent.get(usage)
        for r in ID_RANKS:
            found[r].append(ancestors.get(r))
    return pl.DataFrame(
        {"col:ID": species, **{f"{r}_id": found[r] for r in ID_RANKS}},
        schema={"col:ID": pl.String, **{f"{r}_id": pl.String for r in ID_RANKS}},
    )


def load_taxa(coldp_zip: Path) -> pl.DataFrame:
    """Accepted species from a COLDP export, with the names of their flat classification and the ids of their genus and family."""
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
    species = usages.filter(pl.col("col:rank") == "species", pl.col("col:status").is_in(ACCEPTED))
    return (
        species.join(_ancestor_ids(usages, species["col:ID"].to_list()), on="col:ID", maintain_order="left")
        .select(*(pl.col(c).alias(n) for c, n in SPECIES_COLUMNS.items()), *(f"{r}_id" for r in ID_RANKS))
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
        "with_genus": taxa["genus_id"].drop_nulls().len(),
        "with_family": taxa["family_id"].drop_nulls().len(),
    }


def repeated_names(taxa: pl.DataFrame, rank: str) -> dict:
    """Names that several ids of a rank hold species under within one kingdom, the species below them, and those not under the name's largest id."""
    per_id = taxa.filter(pl.col(f"{rank}_id").is_not_null()).group_by("kingdom", rank, f"{rank}_id").len()
    repeated = per_id.filter(pl.len().over("kingdom", rank) > 1).group_by("kingdom", rank).agg(pl.col("len").sum().alias("species"), pl.col("len").max().alias("largest"))
    return {
        "names": repeated.height,
        "species": int(repeated["species"].sum()),
        "outside_largest": int((repeated["species"] - repeated["largest"]).sum()),
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
