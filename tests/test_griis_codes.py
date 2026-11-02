"""`scripts/griis_codes.py`: freezing the country code each sub-national GRIIS list applies to, from half A's counts.

The derivation itself is `realms.clean.derive_list_codes` (tested in `test_clean`, with the frozen files of the repository against the real data). Here are the script's
own rules: half A only, the resolution of the marker, every place with its code and its share in the log, and a frozen table that is neither replaced by a different one
nor rewritten by the same one.
"""

import importlib.util
import logging
from pathlib import Path

import polars as pl
import pytest

from clean_helpers import counts, row
from realms.clean import LIST_CODES_SCHEMA, read_list_codes
from realms.grid import Grid

ROOT = Path(__file__).parent.parent
CELLS = pl.DataFrame(
    [("11", "FR", "Reunion", "centre"), ("12", "US", "Texas", "centre"), ("13", "NO", "Svalbard", "overlap")],
    schema={"cell": pl.String, "country": pl.String, "subdivision": pl.String, "via": pl.String},
    orient="row",
)
COUNTS = counts(
    row("s1", "11:land", country="RE"),
    row("s1", "12:land", country="US"),
    row("s1", "13:land", country="SJ"),
    row("s2", "13:land", country="SJ"),
    row("s1", "13:sea", country="NO"),
    row("s2", "13:sea", country="RU"),
)


def load_script():
    spec = importlib.util.spec_from_file_location("griis_codes_script", ROOT / "scripts" / "griis_codes.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def world(tmp_path, monkeypatch):
    """The script with the loaders replaced by a made-up world: the grid and the admin-1 cells of `CELLS`, half A's `COUNTS`, and stubs for the files it checks.

    `loads` has the half of every call of `load_country_counts`.
    """
    script = load_script()
    for stub in ("grids/ISEA3H-7.gpkg", "natural_earth/ne_10m_land.zip", "manifest.json"):
        (tmp_path / "reference" / stub).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / "reference" / stub).touch()
    (tmp_path / "split.parquet").touch()
    (tmp_path / "counts").mkdir()
    loads = []

    def counts_of(out_dir, split, half, **kwargs):
        loads.append((half, kwargs.get("unsealed", False)))
        return COUNTS

    monkeypatch.setattr(script, "resolution_of", lambda out_dir, stated=None: 7)
    monkeypatch.setattr(script, "load_grid", lambda gpkg, land: Grid(7, pl.DataFrame({"cell": ["11", "12", "13"]}), {}))
    monkeypatch.setattr(script, "load_admin1_cells", lambda grid, reference: CELLS)
    monkeypatch.setattr(script, "load_split", lambda path: pl.DataFrame())
    monkeypatch.setattr(script, "load_country_counts", counts_of)

    def run(*extra):
        argv = ["--counts", str(tmp_path / "counts"), "--reference", str(tmp_path / "reference"), "--split", str(tmp_path / "split.parquet"), "--out", str(tmp_path / "out"), *extra]
        monkeypatch.setattr("sys.argv", ["griis_codes.py", *argv])
        script.main()

    return type("World", (), {"script": script, "run": staticmethod(run), "loads": loads, "path": tmp_path / "out" / "griis_uses_code-res7.csv", "tmp": tmp_path})


def test_it_derives_from_half_a_and_writes_the_table_sorted_for_the_resolution_of_the_marker(world):
    world.run()
    assert world.loads == [("A", False)]  # half A, sealed: half B is never asked for
    # RE is the code of Reunion's records and FR of its list, SJ that of Svalbard's (2 rows against NO's 1 and RU's 1)
    assert world.path.read_text().splitlines() == ["country,subdivision,code", "FR,Reunion,RE", "NO,Svalbard,SJ", "US,Texas,US"]
    assert read_list_codes(world.path).schema == pl.Schema(LIST_CODES_SCHEMA)


def test_it_logs_every_place_with_its_code_and_its_share_and_warns_of_a_low_one(world, caplog):
    with caplog.at_level(logging.INFO, logger="griis_codes"):
        world.run()
    text = caplog.text
    assert "NO/Svalbard SJ: 50.0% of 4 rows, next NO 1" in text and "US/Texas US: 100.0% of 1 rows" in text
    assert "FR/Reunion RE: 100.0% of 1 rows" in text
    assert "2 places have an effective code with under 80% of the rows in their cells: NO/Svalbard" not in text  # only Svalbard is under, and it is the one
    assert "1 places have an effective code with under 80% of the rows in their cells: NO/Svalbard" in text
    assert "lists applied to another code than their own: FR/Reunion as RE, NO/Svalbard as SJ" in text


def test_a_table_that_is_there_and_the_same_is_left_alone(world):
    world.run()
    before = world.path.stat().st_mtime_ns
    world.run()
    assert world.path.stat().st_mtime_ns == before


def test_the_comparison_does_not_depend_on_the_order_of_the_rows(world):
    world.run()
    unsorted = read_list_codes(world.path).sort("country", descending=True)  # the same table, its rows in another order
    unsorted.write_csv(world.path)
    before = world.path.stat().st_mtime_ns
    world.run()  # it is the same decision: nothing is refused and nothing is written
    assert world.path.stat().st_mtime_ns == before


def test_a_frozen_table_that_differs_is_refused_unless_it_is_meant_to_be_replaced(world):
    world.path.parent.mkdir()
    world.path.write_text("country,subdivision,code\nFR,Reunion,FR\nNO,Svalbard,SJ\nUS,Texas,US\n")  # frozen with another code for Reunion
    with pytest.raises(SystemExit, match=r"FAILED .*is a frozen decision, and it differs in 1 places, such as FR/Reunion.*--overwrite"):
        world.run()
    assert "FR,Reunion,FR" in world.path.read_text()  # untouched
    world.run("--overwrite")
    assert "FR,Reunion,RE" in world.path.read_text()


def test_a_table_of_the_first_form_is_refused_like_any_other_that_differs(world):
    world.path.parent.mkdir()
    world.path.write_text("country,subdivision,uses_code\nFR,Reunion,false\n")
    with pytest.raises(SystemExit, match=r"FAILED .*is a frozen decision, and it cannot be read as a table.*columns.*--overwrite"):
        world.run()
    world.run("--overwrite")
    assert read_list_codes(world.path).height == 3


def test_it_refuses_before_loading_anything_when_a_file_is_missing(world):
    (world.tmp / "reference" / "grids" / "ISEA3H-7.gpkg").unlink()
    with pytest.raises(SystemExit, match=r"FAILED these files are missing: .*ISEA3H-7\.gpkg"):
        world.run()
    assert world.loads == [] and not world.path.exists()


def test_a_grid_of_another_resolution_than_the_marker_is_refused(world, monkeypatch):
    monkeypatch.setattr(world.script, "load_grid", lambda gpkg, land: Grid(6, pl.DataFrame({"cell": ["11"]}), {}))
    with pytest.raises(SystemExit, match="holds resolution 6, not 7"):
        world.run()
    assert world.loads == [] and not world.path.exists()
