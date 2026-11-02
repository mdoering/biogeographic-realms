"""`scripts/griis_codes.py`: freezing where GBIF uses the country code of a sub-national GRIIS list, from half A's counts.

The derivation itself is `realms.clean.derive_uses_code` (tested in `test_clean`, with the frozen files of the repository against the real data). Here are the script's
own rules: half A only, the resolution of the marker, and a frozen table that is neither replaced by a different one nor rewritten by the same one.
"""

import importlib.util
from pathlib import Path

import polars as pl
import pytest

from clean_helpers import counts, row
from realms.clean import USES_CODE_SCHEMA, read_uses_code
from realms.grid import Grid

ROOT = Path(__file__).parent.parent
CELLS = pl.DataFrame(
    [("11", "FR", "Reunion", "centre"), ("12", "US", "Texas", "centre")],
    schema={"cell": pl.String, "country": pl.String, "subdivision": pl.String, "via": pl.String},
    orient="row",
)
COUNTS = counts(row("s1", "11:land", country="RE"), row("s1", "12:land", country="US"))


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
    monkeypatch.setattr(script, "load_grid", lambda gpkg, land: Grid(7, pl.DataFrame({"cell": ["11", "12"]}), {}))
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
    assert world.path.read_text().splitlines() == ["country,subdivision,uses_code", "FR,Reunion,false", "US,Texas,true"]  # RE is the code of the records, FR of the list
    assert read_uses_code(world.path).schema == pl.Schema(USES_CODE_SCHEMA)


def test_a_table_that_is_there_and_the_same_is_left_alone(world):
    world.run()
    before = world.path.stat().st_mtime_ns
    world.run()
    assert world.path.stat().st_mtime_ns == before


def test_a_frozen_table_that_differs_is_refused_unless_it_is_meant_to_be_replaced(world, capsys):
    world.path.parent.mkdir()
    world.path.write_text("country,subdivision,uses_code\nFR,Reunion,true\nUS,Texas,true\n")  # frozen the other way for Reunion
    with pytest.raises(SystemExit, match=r"FAILED .*is a frozen decision and differs.*FR/Reunion.*--overwrite"):
        world.run()
    assert "FR,Reunion,true" in world.path.read_text()  # untouched
    world.run("--overwrite")
    assert "FR,Reunion,false" in world.path.read_text()


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
