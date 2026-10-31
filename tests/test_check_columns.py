import importlib.util
import json
import zipfile
from pathlib import Path

import polars as pl
import pytest

from realms.grid import Grid

SCRIPT = Path(__file__).parent.parent / "scripts" / "check_columns.py"
spec = importlib.util.spec_from_file_location("check_columns", SCRIPT)
check_columns = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_columns)

HEADER = "specieskey\ttaxonomicstatus\tcell6\tcell7\tpart\tn"
# Ids that no float can tell apart, as in the real files
C6A, C6B, C6C = "654539408015826335", "654539408015826336", "654001699011250000"
C7A, C7B, C7C, C7D, C7E = (f"75781397900869151{i}" for i in range(5))
ROWS = [
    ("8VR75", "ACCEPTED", C6A, C7A, "land", 60000),
    ("8VR75", "ACCEPTED", C6A, C7B, "land", 30000),
    ("8VR75", "SYNONYM", C6A, C7A, "land", 9000),
    ("C56N", "ACCEPTED", C6B, C7C, "sea", 20000),
    ("C56N", "AMBIGUOUS_SYNONYM", C6B, C7C, "sea", 500),
    ("4K7TT", "PROVISIONALLY_ACCEPTED", C6A, C7A, "land", 1500),
    ("4K7TT", "MISAPPLIED", C6B, C7D, "land", 100),
    (None, "ACCEPTED", C6A, C7A, "land", 5000),
    ("8VR75", None, C6A, C7A, "sea", 700),
    ("C56N", "ACCEPTED", C6C, C7E, "sea", 300),
]
API = {"ACCEPTED": 115300, "SYNONYM": 9000, "PROVISIONALLY_ACCEPTED": 1500, "AMBIGUOUS_SYNONYM": 500, "MISAPPLIED": 100}
SPECIES = pl.Series(["8VR75", "C56N", "4K7TT", "UNSEEN"])


def tsv(rows, header=HEADER) -> str:
    lines = ["\t".join("" if v is None else str(v) for v in row) for row in rows]
    return "\n".join([header, *lines]) + "\n"


def write_zip(path: Path, members: dict[str, str]) -> Path:
    with zipfile.ZipFile(path, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    return path


def download(tmp_path, rows=ROWS, header=HEADER) -> Path:
    # the member is named .csv but is tab-separated, as in GBIF's SQL downloads
    return write_zip(tmp_path / "0001234-261023074918325.zip", {"0001234-261023074918325.csv": tsv(rows, header)})


def facet(counts=API) -> dict:
    """A GBIF occurrence search response with the COL status facet."""
    return {
        "offset": 0,
        "limit": 0,
        "endOfRecords": False,
        "count": sum(counts.values()),
        "results": [],
        "facets": [{"field": "TAXONOMIC_STATUS", "counts": [{"name": k, "count": v} for k, v in counts.items()]}],
    }


def frame(rows=ROWS) -> pl.DataFrame:
    return pl.DataFrame(rows, schema=check_columns.COLUMNS, orient="row")


def grid(resolution, cells) -> Grid:
    """`cells`: (cell, land_km2, sea_km2)."""
    table = pl.DataFrame(
        {
            "cell": [c for c, *_ in cells],
            "lat": 0.0,
            "lon": 0.0,
            "area_km2": [land + sea for _, land, sea in cells],
            "land_km2": [float(land) for _, land, _ in cells],
            "sea_km2": [float(sea) for _, _, sea in cells],
        }
    )
    return Grid(resolution, table, {})


GRID6 = grid(6, [(C6A, 1000, 3e-8), (C6B, 0, 1000), (C6C, 500, 500)])  # A is all land, but its sea part is float noise, not 0
GRID7 = grid(7, [(C7A, 100, 0), (C7B, 100, 0), (C7C, 0, 100), (C7D, 0, 100), (C7E, 50, 50)])


# reading the download


def test_read_download_parses_the_tab_separated_member(tmp_path):
    df = check_columns.read_download(download(tmp_path))
    assert df.columns == ["specieskey", "taxonomicstatus", "cell6", "cell7", "part", "n"]
    assert df.height == 10
    assert df["n"].dtype == pl.Int64 and df["n"].sum() == 127100
    # ids stay text, and the empty fields of the SQL's nulls are nulls
    assert df["specieskey"].dtype == pl.String and df["cell6"].dtype == pl.String
    assert df["cell6"][0] == C6A
    assert df["specieskey"].null_count() == 1 and df["taxonomicstatus"].null_count() == 1


def test_read_download_rejects_other_columns(tmp_path):
    path = download(tmp_path, header="specieskey\ttaxonomicstatus\tcell6\tcell7\tn\tpart")
    with pytest.raises(ValueError, match="columns"):
        check_columns.read_download(path)


@pytest.mark.parametrize("members", [{}, {"a.csv": tsv(ROWS), "b.csv": tsv(ROWS)}])
def test_read_download_needs_exactly_one_member(tmp_path, members):
    with pytest.raises(ValueError, match="member"):
        check_columns.read_download(write_zip(tmp_path / "x.zip", members))


def test_read_download_rejects_a_part_that_is_neither_land_nor_sea(tmp_path):
    with pytest.raises(ValueError, match="part"):
        check_columns.read_download(download(tmp_path, [("8VR75", "ACCEPTED", C6A, C7A, "Land", 1)]))


@pytest.mark.parametrize("n", ["", "many", "1.5"])
def test_read_download_rejects_a_count_that_is_not_a_whole_number(tmp_path, n):
    with pytest.raises(ValueError, match="n is not a whole number"):
        check_columns.read_download(download(tmp_path, [("8VR75", "ACCEPTED", C6A, C7A, "land", n)]))


# reading the facet


def test_read_facet(tmp_path):
    path = tmp_path / "f.json"
    path.write_text(json.dumps(facet()))
    assert check_columns.read_facet(path) == (sum(API.values()), API)


def test_read_facet_without_the_status_facet(tmp_path):
    path = tmp_path / "f.json"
    path.write_text(json.dumps({"count": 1, "facets": [{"field": "BASIS_OF_RECORD", "counts": []}]}))
    with pytest.raises(ValueError, match="TAXONOMIC_STATUS"):
        check_columns.read_facet(path)


def test_read_facet_that_is_not_a_search_response(tmp_path):
    path = tmp_path / "f.json"
    path.write_text('{"message": "not found"}')
    with pytest.raises(ValueError, match="TAXONOMIC_STATUS"):
        check_columns.read_facet(path)


# check 1: status values


def test_status_values_equal_the_apis_plus_null():
    rows = check_columns.check_status(frame(), API)
    assert rows == {
        "ACCEPTED": (115300, 115300),
        "SYNONYM": (9000, 9000),
        "PROVISIONALLY_ACCEPTED": (1500, 1500),
        "AMBIGUOUS_SYNONYM": (500, 500),
        "MISAPPLIED": (100, 100),
        None: (700, None),
    }


def test_status_value_the_api_lacks_fails():
    rows = [*ROWS, ("8VR75", "DOUBTFUL", C6A, C7A, "land", 3), ("8VR75", "HETEROTYPIC_SYNONYM", C6A, C7A, "land", 4)]
    with pytest.raises(check_columns.CheckFailed, match="DOUBTFUL.*HETEROTYPIC_SYNONYM"):
        check_columns.check_status(frame(rows), API)


def test_status_value_missing_from_sql_fails():
    with pytest.raises(check_columns.CheckFailed, match="MISAPPLIED"):
        check_columns.check_status(frame([r for r in ROWS if r[1] != "MISAPPLIED"]), API)


def test_status_without_null_fails():
    with pytest.raises(check_columns.CheckFailed, match="null"):
        check_columns.check_status(frame([r for r in ROWS if r[1] is not None]), API)


def test_status_counts_must_be_within_one_percent_of_the_api():
    drifted = {**API, "ACCEPTED": 115300 * 100 // 99}  # SQL is 1.00% below, which is allowed
    assert check_columns.check_status(frame(), drifted)["ACCEPTED"] == (115300, drifted["ACCEPTED"])
    for status, away in (("ACCEPTED", 1.011), ("SYNONYM", 0.988)):
        with pytest.raises(check_columns.CheckFailed, match=status):
            check_columns.check_status(frame(), {**API, status: round(API[status] * away)})


# check 2: species


def test_species_all_accepted():
    result = check_columns.check_species(frame(), SPECIES)
    assert (result.species, result.rows, result.records_without_key) == (3, 9, 5000)


def test_species_outside_the_release_fail():
    rows = [*ROWS, ("NOPE1", "ACCEPTED", C6A, C7A, "land", 7), ("NOPE1", "SYNONYM", C6A, C7A, "land", 8), ("NOPE2", "ACCEPTED", C6A, C7A, "land", 1)]
    with pytest.raises(check_columns.CheckFailed, match=r"2 species.*16 records.*NOPE1"):
        check_columns.check_species(frame(rows), SPECIES)


# check 3: cells


def test_cells_all_in_the_grid():
    check_columns.check_cells(frame(), "cell6", GRID6)
    check_columns.check_cells(frame(), "cell7", GRID7)


def test_cells_missing_from_the_grid_fail():
    with pytest.raises(check_columns.CheckFailed, match=rf"cell7.*{C7E}"):
        check_columns.check_cells(frame(), "cell7", grid(7, [(c, 1, 1) for c in (C7A, C7B, C7C, C7D)]))


def test_cells_are_compared_as_text():
    # the ids differ from the grid's only beyond the 15th digit, where a float would call them equal
    assert float(C6A) == float(C6B)
    with pytest.raises(check_columns.CheckFailed, match=C6B):
        check_columns.check_cells(frame(), "cell6", grid(6, [(C6A, 1, 1), (C6C, 1, 1)]))


def test_a_record_without_a_cell_fails():
    rows = [*ROWS, ("8VR75", "ACCEPTED", None, C7A, "land", 3)]
    with pytest.raises(check_columns.CheckFailed, match="null"):
        check_columns.check_cells(frame(rows), "cell6", GRID6)


# check 4: parts of zero area


def test_records_on_parts_of_zero_area():
    zero = check_columns.zero_area_parts(frame(), "cell6", GRID6)
    # A's sea part holds 700 records. Its sea_km2 of 3e-8 is float noise, so it counts as none. B has no land but 100 land records.
    assert zero["sea"] == check_columns.ZeroArea(records=700, cells=1, of=21500)
    assert zero["land"] == check_columns.ZeroArea(records=100, cells=1, of=105600)


def test_a_small_but_real_part_is_not_zero_area():
    cells = grid(6, [(C6A, 1000, 0.26), (C6B, 0.006, 1000), (C6C, 500, 500)])  # the smallest real parts on the grids
    zero = check_columns.zero_area_parts(frame(), "cell6", cells)
    assert (zero["sea"].records, zero["land"].records) == (0, 0)


def test_records_on_parts_of_zero_area_at_resolution_7():
    zero = check_columns.zero_area_parts(frame(), "cell7", GRID7)
    # C7C and C7D are all sea, and C7D holds 100 land records. C7A and C7B are all land, and C7A holds 700 sea records.
    assert zero["sea"].records == 700 and zero["sea"].cells == 1
    assert zero["land"].records == 100 and zero["land"].cells == 1


# check 5: growth


def test_growth_counts_distinct_species_cell_pairs():
    assert check_columns.pair_counts(frame()) == (5, 6)


# the whole script


@pytest.fixture
def run(tmp_path, monkeypatch):
    """Runs `main` on the fixture download, with the reference data replaced by the fixture's."""
    grids = {"ISEA3H-6.gpkg": GRID6, "ISEA3H-7.gpkg": GRID7}
    monkeypatch.setattr(check_columns, "load_taxa", lambda path: pl.DataFrame({"species_id": SPECIES}))
    monkeypatch.setattr(check_columns, "load_grid", lambda gpkg, land: grids[Path(gpkg).name])
    zip_path = download(tmp_path)
    (tmp_path / "0001234-261023074918325.api-facet.json").write_text(json.dumps(facet()))
    return lambda *args: check_columns.main([str(zip_path), *args])


def test_main_reports_every_check(run, capsys):
    run()
    out = capsys.readouterr().out
    assert "127,100" in out  # records
    assert "ACCEPTED" in out and "ok" in out.splitlines()[-1]
    assert "1.20" in out  # 6 pairs against 5


def test_main_reads_the_facet_snapshot_beside_the_zip_by_default(run, tmp_path):
    (tmp_path / "0001234-261023074918325.api-facet.json").write_text(json.dumps(facet({**API, "ACCEPTED": 50})))
    with pytest.raises(SystemExit) as failed:
        run()
    assert "ACCEPTED" in str(failed.value.code)


def test_main_takes_the_facet_snapshot_as_an_argument(run, tmp_path, capsys):
    (tmp_path / "0001234-261023074918325.api-facet.json").write_text("not json")
    other = tmp_path / "elsewhere.json"
    other.write_text(json.dumps(facet()))
    run("--facet", str(other))
    assert capsys.readouterr().out.splitlines()[-1] == "ok"


def test_main_exits_with_a_message_when_a_check_fails(run, monkeypatch):
    monkeypatch.setattr(check_columns, "load_taxa", lambda path: pl.DataFrame({"species_id": ["8VR75"]}))
    with pytest.raises(SystemExit) as failed:
        run()
    assert isinstance(failed.value.code, str) and "C56N" in failed.value.code


def test_main_exits_with_a_message_on_a_malformed_download(run, tmp_path):
    download(tmp_path, header="a\tb")
    with pytest.raises(SystemExit) as failed:
        run()
    assert "columns" in failed.value.code


def test_main_refuses_a_grid_file_of_another_resolution(run, monkeypatch):
    monkeypatch.setattr(check_columns, "load_grid", lambda gpkg, land: GRID7)
    with pytest.raises(SystemExit) as failed:
        run()
    assert "resolution" in failed.value.code
