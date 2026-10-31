"""Check the columns of a GBIF SQL download against what the analysis assumes of them. This submits and fetches nothing.

The download is `queries/check_columns_dk.sql`, run on Denmark alone. It settles three things before a large download is trusted:

1. `taxonomicstatus` holds COL's status values, not the backbone's. The query's clean filter reads them. The API's own top-level status is
   the backbone's, so the comparison is with the facet of `checklistKey` COL, saved as `<download>.api-facet.json` when the download was
   submitted. The download is a fixed snapshot and a later API call drifts with new records, so the script never asks the API.
2. Every `specieskey` is an accepted species of the COL release in `data/reference`.
3. Every `cell6` and `cell7` is a `geogrid_id` of GBIF's ISEA3H file at that resolution.

Checks 1 to 3 are assertions, and the script exits with a message when one fails. Check 4 (records on a part of zero area) and check 5
(growth of the (species, cell) pairs from resolution 6 to 7) are reported only.

The download is untrusted input, so run it as `uv run python -I scripts/check_columns.py data/gbif/<key>.zip`.
"""

import argparse
import io
import json
import zipfile
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from realms.grid import ZERO_KM2, Grid, load_grid
from realms.reference import load_taxa

COLUMNS = {"specieskey": pl.String, "taxonomicstatus": pl.String, "cell6": pl.String, "cell7": pl.String, "part": pl.String, "n": pl.Int64}
PARTS = ("sea", "land")
RESOLUTIONS = (6, 7)
STATUS_TOLERANCE = 0.01  # of the API's count
EXAMPLES = 10


class CheckFailed(AssertionError):
    """A check that must hold doesn't."""


def read_download(path: Path) -> pl.DataFrame:
    """The rows of the download's one member. GBIF names it `.csv` but separates with tabs, and writes a null as nothing."""
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if not name.endswith("/")]
        if len(members) != 1:
            raise ValueError(f"{path} should hold one member, not {len(members)}: {members}")
        content = archive.read(members[0])
    header = content[: content.index(b"\n")].decode().rstrip("\r").split("\t") if b"\n" in content else []
    if header != list(COLUMNS):
        raise ValueError(f"the columns of {path} are {header}, not {list(COLUMNS)}")
    df = pl.read_csv(io.BytesIO(content), separator="\t", quote_char=None, infer_schema=False, null_values=[""])
    n = df["n"].cast(pl.Int64, strict=False)
    if n.null_count():
        raise ValueError(f"n is not a whole number in {n.null_count()} rows of {path}, for example {df.filter(n.is_null())['n'].head(3).to_list()}")
    df = df.with_columns(n)
    unknown = df.filter(~pl.col("part").is_in(PARTS))["part"].unique().to_list()
    if unknown:
        raise ValueError(f"part is {PARTS}, not {unknown}, in {path}")
    return df


def read_facet(path: Path) -> tuple[int, dict[str, int]]:
    """The total and the count per status of a saved occurrence search response, faceted on TAXONOMIC_STATUS."""
    try:
        response = json.loads(Path(path).read_text())
        facets = [f for f in response["facets"] if f["field"] == "TAXONOMIC_STATUS"]
        if len(facets) != 1:
            raise KeyError("TAXONOMIC_STATUS")
        return response["count"], {c["name"]: c["count"] for c in facets[0]["counts"]}
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{path} is not an occurrence search response with one TAXONOMIC_STATUS facet ({error!r})") from error


def check_status(df: pl.DataFrame, api: dict[str, int]) -> dict[str | None, tuple[int, int | None]]:
    """Check 1. The statuses are the API's COL statuses plus null, each count within 1% of the API's. Returns (sql, api) per status."""
    sql = dict(df.group_by("taxonomicstatus").agg(pl.col("n").sum()).iter_rows())
    rows = {status: (sql.get(status, 0), api.get(status)) for status in sorted(set(sql) | set(api), key=lambda s: (s is None, str(s)))}
    problems = []
    extra = sorted(s for s in sql if s is not None and s not in api)
    if extra:
        problems.append(f"SQL has statuses that COL's facet in the API doesn't: {', '.join(extra)}")
    absent = sorted(s for s in api if s not in sql)
    if absent:
        problems.append(f"SQL has none of the API's statuses: {', '.join(absent)}")
    if None not in sql:
        problems.append("SQL has no null status, though records COL doesn't cover should have one")
    for status, count in api.items():
        if status in sql and abs(sql[status] - count) > STATUS_TOLERANCE * count:
            problems.append(f"{status} is {sql[status]:,} in SQL against {count:,} in the API, {(sql[status] - count) / count:+.2%}")
    if problems:
        raise CheckFailed("check 1, the status values: " + "; ".join(problems))
    return rows


@dataclass(frozen=True)
class Species:
    species: int  # distinct specieskeys
    rows: int
    records_without_key: int


def check_species(df: pl.DataFrame, accepted: pl.Series) -> Species:
    """Check 2. No specieskey lies outside the accepted species."""
    keyed = df.filter(pl.col("specieskey").is_not_null())
    outside = keyed.filter(~pl.col("specieskey").is_in(accepted.implode()))
    if outside.height:
        ids = outside["specieskey"].unique().sort()
        raise CheckFailed(
            f"check 2, the species: {ids.len():,} species ({outside['n'].sum():,} records) are not accepted species of the release, "
            f"for example {', '.join(ids.head(EXAMPLES).to_list())}"
        )
    return Species(keyed["specieskey"].n_unique(), keyed.height, df.filter(pl.col("specieskey").is_null())["n"].sum())


def check_cells(df: pl.DataFrame, column: str, grid: Grid) -> int:
    """Check 3. Every cell is a geogrid_id of the grid, compared as text. Returns the number of distinct cells."""
    cells = df[column]
    if cells.null_count():
        raise CheckFailed(f"check 3, the cells: {cells.null_count():,} rows have a null {column}, though every record has coordinates")
    codes = cells.unique().sort()
    missing = codes.filter(~codes.is_in(grid.cells["cell"].implode()))
    if missing.len():
        raise CheckFailed(
            f"check 3, the cells: {missing.len():,} of {codes.len():,} codes in {column} are not geogrid_ids of the resolution "
            f"{grid.resolution} grid, "
            f"for example {', '.join(missing.head(EXAMPLES).to_list())}, where the grid has {', '.join(grid.cells['cell'].head(3).to_list())}"
        )
    return codes.len()


@dataclass(frozen=True)
class ZeroArea:
    records: int
    cells: int
    of: int  # the records of this part, in cells of any area


def zero_area_parts(df: pl.DataFrame, column: str, grid: Grid) -> dict[str, ZeroArea]:
    """Check 4. The records of each part (sea, land) that lie in a cell where the grid has none of that part."""
    areas = grid.cells.select(pl.col("cell").alias(column), pl.col("sea_km2").alias("sea"), pl.col("land_km2").alias("land"))
    joined = df.select(column, "part", "n").join(areas, on=column, how="left")
    result = {}
    for part in PARTS:
        of_part = joined.filter(pl.col("part") == part)
        none = of_part.filter(pl.col(part) < ZERO_KM2)
        result[part] = ZeroArea(none["n"].sum(), none[column].n_unique(), of_part["n"].sum())
    return result


def pair_counts(df: pl.DataFrame) -> tuple[int, int]:
    """Check 5. The distinct (species, cell) pairs at resolution 6 and at 7."""
    keyed = df.filter(pl.col("specieskey").is_not_null())
    return tuple(keyed.select("specieskey", f"cell{r}").unique().height for r in RESOLUTIONS)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("download", type=Path, help="the SQL download, data/gbif/<key>.zip")
    parser.add_argument("--facet", type=Path, help="the API's COL status facet saved with the download, default <download>.api-facet.json beside it")
    parser.add_argument("--reference", type=Path, default=Path("data/reference"), help="the reference data, default data/reference")
    args = parser.parse_args(argv)
    facet_path = args.facet or args.download.with_name(f"{args.download.stem}.api-facet.json")
    reference = args.reference

    try:
        df = read_download(args.download)
        api_total, api = read_facet(facet_path)
        print(f"{args.download.name}: {df.height:,} rows, {df['n'].sum():,} records")

        print(f"\n1. Status values, SQL against the API's COL facet ({facet_path.name})")
        rows = check_status(df, api)
        for status, (sql_count, api_count) in rows.items():
            diff = f"{(sql_count - api_count) / api_count:+.2%}" if api_count else "-"
            print(f"   {str(status):24} sql {sql_count:>12,}   api {f'{api_count:,}' if api_count is not None else '-':>12}   {diff}")
        print(f"   api total {api_total:,}, sql total with a status {sum(s for k, (s, _) in rows.items() if k is not None):,}")

        print("\n2. Species outside the accepted species of the release")
        accepted = load_taxa(reference / "exports" / "extended" / "col.zip")["species_id"]
        species = check_species(df, accepted)
        print(f"   rows with specieskey {species.rows:,}, species {species.species:,}, none outside the {accepted.len():,} accepted")
        print(f"   records without specieskey {species.records_without_key:,}")

        print("\n3. Cell codes against GBIF's grid files")
        grids = {}
        for resolution in RESOLUTIONS:
            gpkg = reference / "grids" / f"ISEA3H-{resolution}.gpkg"
            grid = load_grid(gpkg, reference / "natural_earth" / "ne_10m_land.zip")
            if grid.resolution != resolution:
                raise CheckFailed(f"check 3, the cells: {gpkg.name} holds resolution {grid.resolution}")
            codes = check_cells(df, f"cell{resolution}", grid)
            print(f"   cell{resolution}: {codes:,} codes, all among the {grid.cells.height:,} geogrid_ids of {gpkg.name}")
            grids[resolution] = grid
    except (CheckFailed, ValueError) as error:
        raise SystemExit(f"FAILED {error}") from error

    print("\n4. Records on a part of zero area (reported, not asserted)")
    for resolution, grid in grids.items():
        for part, zero in zero_area_parts(df, f"cell{resolution}", grid).items():
            share = f"{zero.records / zero.of:.2%}" if zero.of else "-"
            print(
                f"   resolution {resolution}, {part} part of cells with no {part}: {zero.records:,} records in {zero.cells:,} cells, "
                f"{share} of the {zero.of:,} {part} records"
            )

    print("\n5. Distinct (species, cell) pairs, resolution 7 against 6 (reported, not asserted)")
    pairs6, pairs7 = pair_counts(df)
    print(f"   res 6: {pairs6:,}, res 7: {pairs7:,}, ratio {pairs7 / pairs6:.2f}")
    print("\nok")


if __name__ == "__main__":
    main()
