"""WoRMS's environment for COL's species, to fill what COL's own environment leaves empty (DECISIONS, "WoRMS's environment", 2026-11-05).

`uv run python -I scripts/worms_environment.py`

Reads:
- the matching of WoRMS (ChecklistBank dataset 2011) against the COL release 315557, a ChecklistBank matching job (Markus, 2026-11-05, job
  59fcfcc1-565c-470f-963e-57be8814a05d, relaxed mode, synonyms on): `data/reference/exports/worms-match/match.zip`;
- WoRMS's extended ColDP export, which holds `col:environment`: `data/reference/exports/worms/worms-2011-coldp-extended.zip`;
- the accepted species of the release (`reference.load_taxa` of the manifest's `col`).

A WoRMS species usage takes its own environment if accepted (or provisionally), and its accepted name's if a synonym. It counts where it matched a COL
name of rank species exactly or as a variant, the COL name accepted (or provisionally), or a synonym whose accepted name is then the species; ambiguous
synonyms and misapplied names do not count. Where several WoRMS usages reach one COL species, their environments are joined (the union of the values,
sorted and comma-separated as COL writes them). Writes `data/reference/worms_environment.parquet` (`species_id`, `worms_environment`, `worms_ids`, the
number of WoRMS usages behind it) and enters it in the manifest as `worms_environment`, with the job, the exports and the file's size and sha256.
"""

import hashlib
import io
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from realms.reference import load_taxa

REFERENCE = Path("data/reference")
MATCH = REFERENCE / "exports" / "worms-match" / "match.zip"
WORMS = REFERENCE / "exports" / "worms" / "worms-2011-coldp-extended.zip"
OUT = REFERENCE / "worms_environment.parquet"
JOB = "59fcfcc1-565c-470f-963e-57be8814a05d"
ACCEPTED = ("accepted", "provisionally accepted")


def read_tsv(archive: Path, member: str, columns: list[str]) -> pl.DataFrame:
    with zipfile.ZipFile(archive) as z, z.open(member) as f:
        return pl.read_csv(io.BytesIO(f.read()), separator="\t", quote_char=None, infer_schema=False, columns=columns)


def main() -> None:
    manifest_path = REFERENCE / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    worms = read_tsv(WORMS, "NameUsage.tsv", ["col:ID", "col:parentID", "col:status", "col:environment"])
    worms = worms.rename({c: c.split(":")[1] for c in worms.columns})
    own = worms.filter(pl.col("status").is_in(ACCEPTED)).select(pl.col("ID").alias("accepted"), pl.col("environment").alias("env"))
    worms = worms.join(own, left_on=pl.when(pl.col("status").is_in(ACCEPTED)).then(pl.col("ID")).otherwise(pl.col("parentID")), right_on="accepted", how="left")
    worms = worms.select(pl.col("ID").alias("original_ID"), "env")
    print(f"WoRMS: {worms.height} usages, {worms['env'].is_not_null().sum()} with an environment (their own or their accepted name's)")

    match = read_tsv(MATCH, "match-dataset-2011.tsv", ["original_ID", "original_rank", "matchType", "ID", "rank", "status", "acceptedID"])
    usable = match.filter(
        pl.col("matchType").is_in(["exact", "variant"]) & (pl.col("rank") == "species") & pl.col("status").is_in(["accepted", "provisionally accepted", "synonym"])
    ).with_columns(pl.when(pl.col("status") == "synonym").then(pl.col("acceptedID")).otherwise(pl.col("ID")).alias("species_id"))
    print(f"match: {match.height} rows, {usable.height} usable (exact or variant, a COL species accepted or a plain synonym)")

    species = set(load_taxa(Path(manifest["col"]["path"]))["species_id"])
    joined = usable.join(worms, on="original_ID", how="inner").filter(pl.col("env").is_not_null() & pl.col("species_id").is_in(list(species)))
    table = (
        joined.with_columns(pl.col("env").str.split(","))
        .explode("env")
        .group_by("species_id")
        .agg(pl.col("env").unique().sort().str.join(",").alias("worms_environment"), pl.col("original_ID").n_unique().alias("worms_ids"))
        .sort("species_id")
    )
    table.write_parquet(OUT)
    with pl.Config(tbl_rows=12):
        print(f"{table.height} COL species with a WoRMS environment\n{table.group_by('worms_environment').len().sort('len', descending=True)}")
    print(f"{table.filter(pl.col('worms_ids') > 1).height} of them from more than one WoRMS usage")

    data = OUT.read_bytes()
    manifest["worms_environment"] = {
        "path": str(OUT), "size": len(data), "sha256": hashlib.sha256(data).hexdigest(), "species": table.height,
        "worms_dataset": 2011, "worms_export": (WORMS.parent / "source-url-extended.txt").read_text().strip(),
        "match_job": JOB, "match_url": (MATCH.parent / "source-url.txt").read_text().strip(), "release": manifest["col"]["release"],
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {OUT} and the manifest's worms_environment")


if __name__ == "__main__":
    main()
