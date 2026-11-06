"""GBIF occurrence counts per species, country and part of a cell: streamed from the SQL download into Parquet, checked, and read per half.

The download is one tab-separated member, named `.csv`, with one row per (species, country, cell, part) and the counts of its clean
records (`queries/occurrence_res<r>.sql`). Its rows hold no free text, so it can be cut at any newline. It is never written out as text:
`ingest` feeds it to the Parquet writer a chunk of lines at a time.

Half B of the species stays sealed until the pipeline is frozen (`docs/design.md`). `load_presence` and `load_country_counts` refuse it
unless the caller says `unsealed=True`, and each such call site is listed in `DECISIONS.md`. `qa` reports global numbers only.

GBIF's code for a point and the grid file's id of the same cell can differ by a few hundred metres in a handful of cells (seven codes in 35 rows of the
resolution 7 download). `ingest_remapped` maps such a code to the cell of the grid that it is a near miss of, while it ingests (`Centres`), so that the parts
never hold a code that is not in the grid, unless it is no near miss of any cell: those stay, and the QA fails on them.

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
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import numpy as np
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
# the visits of a download made with `gbif.build_sql`'s `visit_resolution` (DECISIONS, "Species-days", 2026-11-06): distinct dates and fine cells, which
# add up over the countries of a part (a fine cell lies in one country but at a border), so they are summed
VISITS = ("n_visits", "n_visits_open")
VISIT_COLUMNS = {**COLUMNS, **{c: pl.Int64 for c in VISITS}}


def _schema(found: list[str], source: object) -> dict[str, pl.DataType]:
    """The columns of a download with this header: `COLUMNS`, or `VISIT_COLUMNS` for one with visits. Anything else raises."""
    for schema in (COLUMNS, VISIT_COLUMNS):
        if found == list(schema):
            return schema
    raise ValueError(f"the columns of {source} are {found}, not {list(COLUMNS)} (or with {list(VISITS)} after them)")
HALVES = ("A", "B", "all")
BUFFER = 1 << 20
EXAMPLES = 10
COMPLETE = "_COMPLETE"  # the marker of a complete ingest, in the directory of its parts
MARKER_KEYS = ("rows", "resolution", "zip", "zip_sha256")
KEY = ("specieskey", "countrycode", "cell", "part")  # what makes a row of the download: it has one row of each
# A cell code that is not in the grid is a near miss of a cell, and is remapped to it, when that cell is the only one with the same 2-digit code within
# REMAP_KM of the code's own centre and the next is more than REMAP_RATIO times as far. More than MAX_REMAPS such codes in a download is another grid.
REMAP_KM = 5.0
REMAP_RATIO = 3.0
MAX_REMAPS = 100
EARTH_KM = 6371.0088  # the mean radius of the earth, which gets a distance within 0.5% of the geodesic one: all that a limit of 5 km needs
HEMISPHERE_TOLERANCE = 1e-6  # degrees, the precision of a cell id: a centre this close to the equator or the meridian is on neither side


class SealedHalf(Exception):
    """Half B was asked for without `unsealed=True`."""


class IncompleteIngest(Exception):
    """A directory of parts has no `_COMPLETE`, or its parts are not what the marker says: an ingest that did not finish, or did not pass its QA."""


@dataclass(frozen=True)
class Remap:
    """One cell code of the download that is not in the grid and was mapped to the cell `target`, which is `km` from where the code puts it.

    `rows` and the records (`n_raw` all of them, `n` the clean ones) are those of the download's rows with the code.
    """

    code: str
    target: str
    rows: int
    n_raw: int
    n: int
    km: float


@dataclass(frozen=True)
class Ingested:
    """What `ingest_remapped` made: the `rows` of the parts, the `downloaded` rows of the download, and every cell code it remapped, sorted by code."""

    rows: int
    downloaded: int
    remaps: tuple[Remap, ...]


def remaps_report(done: Ingested) -> dict:
    """The remaps of an ingest as the QA reports them: their totals, the rows that were merged into others, and every one of them."""
    remaps = done.remaps
    return {
        "codes": len(remaps),
        "rows": sum(r.rows for r in remaps),
        "n_raw": sum(r.n_raw for r in remaps),
        "n": sum(r.n for r in remaps),
        "merged_rows": done.downloaded - done.rows,
        "list": [asdict(r) for r in remaps],
    }


def decode(code: str) -> tuple[str, float, float] | None:
    """The 2-digit code of a cell code and the latitude and longitude it gives, as absolute values, or None if it is no place.

    A `geogrid_id` without its sign and zero-padded to 19 digits is the code (2 digits, which hold the face and the hemisphere), then |lat| in millionths of a
    degree in 8 digits and |lon| in 9. The sign says nothing of the hemisphere: `Centres` takes that from the code, as the grid shows it.
    """
    digits = code[1:] if code.startswith("-") else code
    if not (digits.isascii() and digits.isdigit()) or len(digits) > 19:
        return None
    padded = digits.zfill(19)
    lat, lon = int(padded[2:10]) / 1e6, int(padded[10:]) / 1e6
    if lat > 90 or lon > 180:
        return None
    return padded[:2], lat, lon


def _km(lat: float, lon: float, lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """Great-circle distances in km from a point to arrays of points, in degrees (haversine)."""
    p, ps = np.radians(lat), np.radians(lats)
    a = np.sin((ps - p) / 2) ** 2 + np.cos(p) * np.cos(ps) * np.sin(np.radians(lons - lon) / 2) ** 2
    return 2 * EARTH_KM * np.arcsin(np.sqrt(np.minimum(a, 1.0)))


def _side(values: np.ndarray) -> int | None:
    """1 if no value of a group is below zero, -1 if none is above, None if there are both (the sign of the group is not known)."""
    if (values >= -HEMISPHERE_TOLERANCE).all():
        return 1
    if (values <= HEMISPHERE_TOLERANCE).all():
        return -1
    return None


@dataclass(frozen=True)
class _Group:
    """The cells of a grid with one 2-digit code, with the hemisphere that the code stands for."""

    ids: list[str]
    lat: np.ndarray
    lon: np.ndarray
    north: int | None
    east: int | None


class Centres:
    """The cells of a grid by the 2-digit code of their ids, to find the cell that a code that is not in the grid is a near miss of.

    A code is a near miss of a cell if the cell has the same 2-digit code, and is the only one within `REMAP_KM` of the code's own centre and the next is more than
    `REMAP_RATIO` times as far, and the code decodes to a place (`decode`). The code of a cell stands for a hemisphere, which is read from the centres of the
    cells that have it (a code whose cells lie on both sides of the equator or the prime meridian gives no place), so that the centre of a code is a place and
    not a pair of absolute values. `cells` is `cell` (text), `lat` and `lon` (the centres, in degrees).
    """

    def __init__(self, cells: pl.DataFrame):
        self.cells = cells["cell"]
        found: dict[str, tuple[list[str], list[float], list[float]]] = {}
        for cell, lat, lon in cells.select("cell", "lat", "lon").iter_rows():
            decoded = decode(cell)
            if decoded is not None:
                ids, lats, lons = found.setdefault(decoded[0], ([], [], []))
                ids.append(cell)
                lats.append(lat)
                lons.append(lon)
        self._groups = {}
        for code, (ids, lats, lons) in found.items():
            lat, lon = np.array(lats), np.array(lons)
            self._groups[code] = _Group(ids, lat, lon, _side(lat), _side(lon))
        self._nearest: dict[str, tuple[str, float] | None] = {}

    def in_grid(self) -> pl.Expr:
        """Whether the cell of a row is a cell of the grid, and not null: an expression on the column `cell`."""
        return pl.col("cell").is_in(self.cells.implode()).fill_null(False)

    def _locate(self, code: str) -> tuple[_Group, float, float] | None:
        decoded = decode(code)
        group = self._groups.get(decoded[0]) if decoded else None
        if group is None or group.north is None or group.east is None:
            return None
        return group, group.north * decoded[1], group.east * decoded[2]

    def place(self, code: str) -> tuple[float, float] | None:
        """The latitude and longitude of a cell code, or None if it is no place or the grid doesn't show its hemisphere."""
        located = self._locate(code)
        return None if located is None else located[1:]

    def nearest(self, code: str) -> tuple[str, float] | None:
        """The cell that `code` is a near miss of, and its distance in km (to a metre), or None if it is a near miss of no cell."""
        if code not in self._nearest:
            self._nearest[code] = self._find(code)
        return self._nearest[code]

    def _find(self, code: str) -> tuple[str, float] | None:
        located = self._locate(code)
        if located is None:
            return None
        group, lat, lon = located
        km = _km(lat, lon, group.lat, group.lon)
        order = np.argsort(km, kind="stable")
        first = float(km[order[0]])
        if first > REMAP_KM or (len(order) > 1 and not km[order[1]] > REMAP_RATIO * first):
            return None
        return group.ids[order[0]], round(first, 3)

    @property
    def remapped(self) -> int:
        """How many codes have been found to be near misses so far."""
        return sum(found is not None for found in self._nearest.values())


def _clear(out_dir: Path) -> None:
    """Remove the marker first and then the parts (and the half written copy of one that `_fold` was rewriting), so that no marker ever stands beside
    parts it was not written for."""
    (out_dir / COMPLETE).unlink(missing_ok=True)
    for path in out_dir.glob("part-*.parquet*"):
        path.unlink()


def _parse(content: bytes, source: str, schema: dict[str, pl.DataType] = COLUMNS) -> pl.DataFrame:
    """The rows of a chunk that starts with the header, in `schema`. An empty field is a null."""
    try:
        df = pl.read_csv(io.BytesIO(content), separator="\t", quote_char=None, schema=schema, null_values=[""])
    except pl.exceptions.PolarsError as error:
        raise ValueError(f"{source} has a row that is not {list(schema)}: {error}") from error
    unknown = df.filter(~pl.col("part").is_in(KINDS).fill_null(False))["part"].unique().to_list()
    if unknown:
        raise ValueError(f"part is one of {KINDS}, not {unknown}, in {source}")
    return df


def _set_aside(df: pl.DataFrame, centres: Centres, asides: list[pl.DataFrame]) -> pl.DataFrame:
    """`df` without the rows whose cell code is not in the grid but is a near miss of a cell (`Centres`), which are added to `asides` as they are.

    A row with no cell, and one whose code is no near miss, stays: the QA counts them and fails. Stops if the download has more than `MAX_REMAPS` near misses.
    """
    codes = df.filter(centres.in_grid().not_() & pl.col("cell").is_not_null())["cell"].unique().to_list()
    near = [code for code in codes if centres.nearest(code) is not None]
    if centres.remapped > MAX_REMAPS:
        raise ValueError(
            f"{centres.remapped} cell codes are not in the grid but near a cell of it, more than the {MAX_REMAPS} that a rounding difference of a few cells "
            "explains: this is another grid, and nothing is remapped"
        )
    if not near:
        return df
    moved = pl.col("cell").is_in(near).fill_null(False)
    asides.append(df.filter(moved))
    return df.filter(moved.not_())


def _fold(out_dir: Path, asides: list[pl.DataFrame], centres: Centres, index: int) -> tuple[tuple[Remap, ...], int]:
    """Give the rows that were set aside the cell they are a near miss of, and write them, merged, as the last part, `part-<index>.parquet`.

    Returns the remaps and the change in the rows of the parts. A row of the download is one (species, country, cell, part), so a remapped row can have the
    same key as one of the parts, or as another remapped row. All of them are merged as `load_presence` merges the rows of a country: the counts add up and
    the distinct counts take the larger. The rows of the earlier parts that it merges are taken out of them, so that no key is in two rows.
    """
    moved = pl.concat(asides)
    remaps = []
    for code, rows, n_raw, n in moved.group_by("cell").agg(pl.len(), pl.col("n_raw").sum(), pl.col("n").sum()).sort("cell").iter_rows():
        target, km = centres.nearest(code)
        remaps.append(Remap(code, target, rows, n_raw, n, km))
    moved = moved.with_columns(pl.col("cell").replace_strict({r.code: r.target for r in remaps}, return_dtype=pl.String))
    keys = moved.select(KEY).unique()
    targets = keys["cell"].unique().implode()
    pool, taken = [moved], 0
    for path in sorted(out_dir.glob("part-*.parquet")):
        hit = pl.scan_parquet(path).filter(pl.col("cell").is_in(targets)).join(keys.lazy(), on=KEY, how="semi", nulls_equal=True).collect()
        if hit.is_empty():
            continue
        # streamed to a file of its own and renamed, as the part holds millions of rows and a part is never half written
        kept = path.with_name(path.name + ".tmp")
        pl.scan_parquet(path).join(keys.lazy(), on=KEY, how="anti", nulls_equal=True).sink_parquet(kept)
        os.replace(kept, path)
        pool.append(hit)
        taken += hit.height
    visits = [c for c in VISITS if c in moved.columns]
    sums = [pl.col(c).sum() for c in (*SUMMED, *visits)]
    largest = [pl.col(c).max() for c in DISTINCT]
    merged = pl.concat(pool).group_by(*KEY).agg(*sums, *largest).select(*COLUMNS, *visits).sort(*KEY)
    merged.write_parquet(out_dir / f"part-{index:05d}.parquet")
    return tuple(remaps), merged.height - taken


def _ingest(zip_path: Path, out_dir: Path, chunk_lines: int, centres: Centres | None) -> Ingested:
    if chunk_lines < 1:
        raise ValueError(f"chunk_lines is at least 1, not {chunk_lines}")
    out_dir.mkdir(parents=True, exist_ok=True)
    _clear(out_dir)
    rows = downloaded = 0
    asides: list[pl.DataFrame] = []
    remaps: tuple[Remap, ...] = ()
    try:
        with zipfile.ZipFile(zip_path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/")]
            if len(members) != 1:
                raise ValueError(f"{zip_path} should hold one member, not {len(members)}: {members}")
            # a buffered reader reads lines in C, where the zip member's own readline is Python
            with io.BufferedReader(archive.open(members[0]), BUFFER) as lines:
                header = lines.readline()
                schema = _schema(header.decode().rstrip("\r\n").split("\t"), zip_path)
                for i in itertools.count():
                    chunk = b"".join(itertools.islice(lines, chunk_lines))
                    if not chunk:
                        break
                    df = _parse(header + chunk, f"{zip_path.name}, chunk {i}", schema)
                    downloaded += df.height
                    if centres is not None:
                        df = _set_aside(df, centres, asides)
                    df.write_parquet(out_dir / f"part-{i:05d}.parquet")
                    rows += df.height
        if asides:
            remaps, change = _fold(out_dir, asides, centres, i)
            rows += change
        if not rows:
            raise ValueError(f"{zip_path} has no rows")
    except BaseException:
        _clear(out_dir)
        raise
    return Ingested(rows, downloaded, remaps)


def ingest(zip_path: Path, out_dir: Path, chunk_lines: int = 5_000_000) -> int:
    """Stream the zip's single member into `out_dir/part-00000.parquet` and onwards, `chunk_lines` lines to a part. Returns the rows.

    The marker and the parts of an earlier run are removed first, and the parts written by a run that fails are removed too. A run that is killed
    leaves its parts, but no marker: the loaders refuse them until `mark_complete` is written after the QA. Cell codes are left as they are: see
    `ingest_remapped`.
    """
    return _ingest(zip_path, out_dir, chunk_lines, None).rows


def ingest_remapped(zip_path: Path, out_dir: Path, cells: pl.DataFrame, chunk_lines: int = 5_000_000) -> Ingested:
    """`ingest`, with the cell codes that are not in the grid but are near misses of a cell of it (`Centres`) mapped to that cell.

    `cells` is `cell`, `lat` and `lon` of the grid (`realms.grid.read_centres`). The rows of a remapped code are held back, merged with any row of the
    cell they go to (`_fold`), and written as the last part, so that the parts never hold such a code and no key is in two rows. Every remap is in the result.
    Any other code that the grid doesn't have is left for the QA to count: nothing is dropped.
    """
    return _ingest(zip_path, out_dir, chunk_lines, Centres(cells))


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
    - `unknown_cells`: rows and records of cell codes the grid doesn't have, how many codes, and a few as examples. The ingest has remapped the near misses of
      a cell (`ingest_remapped`), and the script reports those apart, under `remaps`: what is left here is a code that is no near miss.
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
    scan = scan_counts(out_dir)
    visits = [c for c in VISITS if c in scan.collect_schema().names()]  # only in a download with visits, and then after the other counts
    sums = [pl.col(c).sum() for c in (*SUMMED, *visits)]
    # the distinct counts of the countries in one part overlap, so the sum would overcount and the maximum is a lower bound
    largest = [pl.col(c).max() for c in DISTINCT]
    part = (pl.col("cell") + pl.lit(":") + pl.col("part")).alias("part")
    return (
        scan
        .join(species.select("species_id").lazy(), on="species_id", how="inner")
        .group_by(*keep, "cell", "part")
        .agg(*sums, *largest)
        .with_columns(part)
        .select(*keep, "part", *COUNTS, *visits)
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
