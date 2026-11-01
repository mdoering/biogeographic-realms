"""GBIF occurrence counts per species, country and part of a cell: streamed from the SQL download into Parquet, checked, and read per half.

The download is one tab-separated member, named `.csv`, with one row per (species, country, cell, part) and the counts of its clean
records (`queries/occurrence_res<r>.sql`). Its rows hold no free text, so it can be cut at any newline. It is never written out as text:
`ingest` feeds it to the Parquet writer a chunk of lines at a time.

Half B of the species stays sealed until the pipeline is frozen (`docs/design.md`). `load_presence` and `load_country_counts` refuse it
unless the caller says `unsealed=True`, and each such call site is listed in `DECISIONS.md`. `qa` reports global numbers only.

An ingest is complete when its directory holds `_COMPLETE` (`mark_complete`), which `scripts/ingest.py` writes once the ingest and its QA have passed,
with the rows, the resolution and the zip's sha256. `ingest` removes it before it writes anything, so a run that is killed half way leaves parts and no
marker, and `load_presence` and `load_country_counts` refuse a directory without one (`IncompleteIngest`), or whose parts hold other rows than it says.
"""

import hashlib
import io
import itertools
import json
import os
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import polars as pl

from realms.grid import KINDS, ZERO_KM2, Grid

COUNTS = (
    "n_raw",
    "n",
    "n_datasets",
    "n_years",
    "n_nc",
    "n_open",
    "n_datasets_open",
    "n_years_open",
    "n_ebird",
    "n_inat",
    "n_synonym",
    "n_human_text",
)
# distinct counts, which don't add up over countries
DISTINCT = ("n_datasets", "n_years", "n_datasets_open", "n_years_open")
SUMMED = tuple(c for c in COUNTS if c not in DISTINCT)
# the columns of the download, in the order of the query, with the types they are stored as. Ids are text, as floats lose digits.
COLUMNS = {"specieskey": pl.String, "countrycode": pl.String, "cell": pl.String, "part": pl.String, **{c: pl.Int64 for c in COUNTS}}
HALVES = ("A", "B", "all")
BUFFER = 1 << 20
EXAMPLES = 10
COMPLETE = "_COMPLETE"  # the marker of a complete ingest, in the directory of its parts
MARKER_KEYS = ("rows", "resolution", "zip", "zip_sha256")


class SealedHalf(Exception):
    """Half B was asked for without `unsealed=True`."""


class IncompleteIngest(Exception):
    """A directory of parts has no `_COMPLETE`, or its parts are not what the marker says: an ingest that did not finish, or did not pass its QA."""


def _clear(out_dir: Path) -> None:
    """Remove the marker first and then the parts, so that no marker ever stands beside parts it was not written for."""
    (out_dir / COMPLETE).unlink(missing_ok=True)
    for path in out_dir.glob("part-*.parquet"):
        path.unlink()


def _parse(content: bytes, source: str) -> pl.DataFrame:
    """The rows of a chunk that starts with the header. An empty field is a null."""
    try:
        df = pl.read_csv(io.BytesIO(content), separator="\t", quote_char=None, schema=COLUMNS, null_values=[""])
    except pl.exceptions.PolarsError as error:
        raise ValueError(f"{source} has a row that is not {list(COLUMNS)}: {error}") from error
    unknown = df.filter(~pl.col("part").is_in(KINDS).fill_null(False))["part"].unique().to_list()
    if unknown:
        raise ValueError(f"part is one of {KINDS}, not {unknown}, in {source}")
    return df


def ingest(zip_path: Path, out_dir: Path, chunk_lines: int = 5_000_000) -> int:
    """Stream the zip's single member into `out_dir/part-00000.parquet` and onwards, `chunk_lines` lines to a part. Returns the rows.

    The marker and the parts of an earlier run are removed first, and the parts written by a run that fails are removed too. A run that is killed
    leaves its parts, but no marker: the loaders refuse them until `mark_complete` is written after the QA.
    """
    if chunk_lines < 1:
        raise ValueError(f"chunk_lines is at least 1, not {chunk_lines}")
    out_dir.mkdir(parents=True, exist_ok=True)
    _clear(out_dir)
    rows = 0
    try:
        with zipfile.ZipFile(zip_path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise ValueError(f"{zip_path} should hold one member, not {len(members)}: {members}")
            # a buffered reader reads lines in C, where the zip member's own readline is Python
            with io.BufferedReader(archive.open(members[0]), BUFFER) as lines:
                header = lines.readline()
                found = header.decode().rstrip("\r\n").split("\t")
                if found != list(COLUMNS):
                    raise ValueError(f"the columns of {zip_path} are {found}, not {list(COLUMNS)}")
                for i in itertools.count():
                    chunk = b"".join(itertools.islice(lines, chunk_lines))
                    if not chunk:
                        break
                    df = _parse(header + chunk, f"{zip_path.name}, chunk {i}")
                    df.write_parquet(out_dir / f"part-{i:05d}.parquet")
                    rows += df.height
        if not rows:
            raise ValueError(f"{zip_path} has no rows")
    except BaseException:
        _clear(out_dir)
        raise
    return rows


def _parts(out_dir: Path) -> list[Path]:
    parts = sorted(out_dir.glob("part-*.parquet"))
    if not parts:
        raise FileNotFoundError(f"no part-*.parquet in {out_dir}")
    return parts


def scan_counts(out_dir: Path) -> pl.LazyFrame:
    """The parts of `ingest`, lazily, with `specieskey` renamed `species_id`, whether the ingest is complete or not (`qa` reads them before the marker).

    It holds every species, half B included, and is not sealed. It has two readers: `qa`, whose numbers are global, and `_select`, which checks the
    seal and the marker first and then reads it per species, for the half it was asked for.
    """
    return pl.scan_parquet(_parts(out_dir)).rename({"specieskey": "species_id"})


def part_rows(out_dir: Path) -> int:
    """How many rows the parts of `out_dir` hold, from their Parquet metadata: a global number, which reads no row."""
    return pl.scan_parquet(_parts(out_dir)).select(pl.len()).collect().item()


def sha256(path: Path) -> str:
    """The sha256 of a file, read a megabyte at a time."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(BUFFER):
            digest.update(chunk)
    return digest.hexdigest()


def mark_complete(out_dir: Path, zip_path: Path, resolution: int, rows: int, how: str = "ingest") -> dict:
    """Write `out_dir/_COMPLETE`: the parts are a whole ingest of `zip_path` at `resolution`, and passed their QA. Returns what was written.

    `rows` is what the caller counted (the ingest's rows, or `qa.json`'s), and it must be what the parts hold. The marker holds `rows`, `resolution`, `zip`
    (the file's name), `zip_sha256`, `written` (the time) and `how` (`ingest`, or `mark-complete` for a marker added to an earlier ingest). It is written to
    a temporary file and renamed, so a marker is never half written.
    """
    found = part_rows(out_dir)
    if found != rows:
        raise IncompleteIngest(f"the parts of {out_dir} hold {found:,} rows, not the {rows:,} of the ingest: no {COMPLETE} is written")
    marker = {
        "rows": rows,
        "resolution": resolution,
        "zip": zip_path.name,
        "zip_sha256": sha256(zip_path),
        "written": datetime.now(UTC).isoformat(timespec="seconds"),
        "how": how,
    }
    temporary = out_dir / f"{COMPLETE}.tmp"
    temporary.write_text(json.dumps(marker, indent=2) + "\n")
    os.replace(temporary, out_dir / COMPLETE)
    return marker


def read_complete(out_dir: Path) -> dict:
    """The marker of a complete ingest in `out_dir` (`mark_complete`), after checking that its parts hold the rows it says.

    Raises `IncompleteIngest` if there is no marker, if it lacks a key, or if the parts hold other rows: an ingest that was killed, failed its QA or
    was changed since.
    """
    path = out_dir / COMPLETE
    if not path.exists():
        raise IncompleteIngest(
            f"{out_dir} has no {COMPLETE}: its ingest did not finish or did not pass its QA. Run scripts/ingest.py on the download again, "
            "or, for an ingest made before the marker existed, scripts/ingest.py --mark-complete"
        )
    marker = json.loads(path.read_text())
    missing = [k for k in MARKER_KEYS if k not in marker]
    if missing:
        raise IncompleteIngest(f"{path} lacks {missing}")
    try:
        found = part_rows(out_dir)
    except FileNotFoundError as error:
        raise IncompleteIngest(f"{path} is there, and no part is: {error}") from error
    if found != marker["rows"]:
        raise IncompleteIngest(f"the parts of {out_dir} hold {found:,} rows, and {COMPLETE} says {marker['rows']:,}: the parts changed after the ingest")
    return marker


def _tally(flag: pl.Expr, name: str) -> list[pl.Expr]:
    return [
        flag.sum().alias(f"{name}.rows"),
        pl.col("n_raw").filter(flag).sum().alias(f"{name}.n_raw"),
        pl.col("n").filter(flag).sum().alias(f"{name}.n"),
    ]


def qa(counts: pl.LazyFrame, taxa: pl.DataFrame, grid: Grid) -> dict:
    """Global numbers on the counts, and nothing about any taxon, so that half B is not looked at before the freeze.

    - `rows`, `n_raw`, `n`, `n_synonym` and `synonym_share`, which is `sum(n_synonym) / sum(n)` (None if no clean record)
    - `unaccepted`: rows and records of species that are not in `taxa`, and how many species that is
    - `unknown_cells`: rows and records of cell codes the grid doesn't have, how many codes, and a few as examples
    - `zero_area`: rows and records of land parts in cells with no land, and sea parts in cells with no sea (lakes and the like, or a
      coast the land polygons miss), in how many cells. An area under `ZERO_KM2` is none.

    Records are given as `n_raw` (every record) and `n` (the clean ones).
    """
    areas = grid.cells.select("cell", "land_km2", "sea_km2").lazy()
    area = pl.when(pl.col("part") == "land").then(pl.col("land_km2")).otherwise(pl.col("sea_km2"))
    accepted = pl.col("species_id").is_in(taxa["species_id"].implode()).fill_null(False)
    unknown = pl.col("land_km2").is_null()  # the join found no cell, and every cell of the grid has a land area
    zero = {kind: ((pl.col("part") == kind) & (area < ZERO_KM2)).fill_null(False) for kind in KINDS}

    row = (
        counts.select("species_id", "cell", "part", "n_raw", "n", "n_synonym")
        .join(areas, on="cell", how="left")
        .select(
            pl.len().alias("rows"),
            pl.col("n_raw").sum(),
            pl.col("n").sum(),
            pl.col("n_synonym").sum(),
            *_tally(~accepted, "unaccepted"),
            pl.col("species_id").filter(~accepted).n_unique().alias("unaccepted.species"),
            *_tally(unknown, "unknown"),
            pl.col("cell").filter(unknown).n_unique().alias("unknown.cells"),
            pl.col("cell").filter(unknown).drop_nulls().unique().sort().head(EXAMPLES).implode().alias("unknown.examples"),
            *(e for kind in KINDS for e in _tally(zero[kind], kind)),
            *(pl.col("cell").filter(zero[kind]).n_unique().alias(f"{kind}.cells") for kind in KINDS),
        )
        .collect()
        .row(0, named=True)
    )

    def group(prefix: str, distinct: str) -> dict:
        """`distinct` is what is counted besides rows and records: the species or the cells."""
        return {
            distinct: row[f"{prefix}.{distinct}"],
            "rows": row[f"{prefix}.rows"],
            "n_raw": row[f"{prefix}.n_raw"],
            "n": row[f"{prefix}.n"],
        }

    n, n_synonym = row["n"], row["n_synonym"]
    return {
        "rows": row["rows"],
        "n_raw": row["n_raw"],
        "n": n,
        "n_synonym": n_synonym,
        "synonym_share": n_synonym / n if n else None,
        "unaccepted": group("unaccepted", "species"),
        "unknown_cells": {**group("unknown", "cells"), "examples": row["unknown.examples"]},
        "zero_area": {kind: group(kind, "cells") for kind in KINDS},
    }


def resolution_of(out_dir: Path, stated: int | None = None) -> int:
    """The resolution of the complete ingest in `out_dir`, from its marker (`read_complete`), so that a script knows it before it loads anything large.

    `stated` is a resolution the caller was given, and it must be the marker's: a mismatch raises `ValueError`.
    """
    found = read_complete(out_dir)["resolution"]
    if stated is not None and stated != found:
        raise ValueError(f"{out_dir} is an ingest at resolution {found}, not {stated}")
    return found


def _select(out_dir: Path, split: pl.DataFrame, half: str, unsealed: bool, keep: tuple[str, ...]) -> pl.DataFrame:
    """The counts of the species of `half`, summed over what `keep` doesn't name. A species that is not in `split` is left out.

    The seal is checked before anything is read, and then the marker of a complete ingest (`read_complete`), before any row is.
    """
    if half not in HALVES:
        raise ValueError(f"half is one of {HALVES}, not {half!r}")
    if half != "A" and not unsealed:
        raise SealedHalf(
            f"half {half!r} includes half B, which is sealed until the pipeline is frozen. Pass unsealed=True only at a call site "
            "that DECISIONS.md lists."
        )
    if split["species_id"].n_unique() != split.height:
        raise ValueError("the split has a species_id more than once")
    read_complete(out_dir)
    species = split if half == "all" else split.filter(pl.col("half") == half)
    sums = [pl.col(c).sum() for c in SUMMED]
    # the distinct counts of the countries in one part overlap, so the sum would overcount and the maximum is a lower bound
    largest = [pl.col(c).max() for c in DISTINCT]
    part = (pl.col("cell") + pl.lit(":") + pl.col("part")).alias("part")
    return (
        scan_counts(out_dir)
        .join(species.select("species_id").lazy(), on="species_id", how="inner")
        .group_by(*keep, "cell", "part")
        .agg(*sums, *largest)
        .with_columns(part)
        .select(*keep, "part", *COUNTS)
        .sort(*keep, "part")
        .collect()
    )


def load_presence(out_dir: Path, split: pl.DataFrame, half: Literal["A", "B", "all"], *, unsealed: bool = False) -> pl.DataFrame:
    """`species_id`, `part` (a `part_id`) and the counts of the species of `half`, summed over countries.

    `n_datasets`, `n_years` and their `_open` forms take the maximum over the countries instead. A dataset or a year can count in
    two countries of a border cell, so the true number lies between the maximum and the sum, and the maximum is a lower bound.
    `split` is `species_id, lineage, half`. Species it doesn't have are left out. "B" and "all" raise `SealedHalf` unless `unsealed`, and a directory
    without the marker of a complete ingest raises `IncompleteIngest`.
    """
    return _select(out_dir, split, half, unsealed, ("species_id",))


def load_country_counts(out_dir: Path, split: pl.DataFrame, half: Literal["A", "B", "all"], *, unsealed: bool = False) -> pl.DataFrame:
    """The same as `load_presence`, with a row per country as well (`species_id`, `countrycode`, `part`, counts), for GRIIS."""
    return _select(out_dir, split, half, unsealed, ("species_id", "countrycode"))
