import ast
import importlib.util
import json
import re
import zipfile
from pathlib import Path

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from realms import grid as grid_module
from realms import occurrences
from realms.grid import Grid, part_id
from realms.occurrences import (
    COLUMNS,
    COMPLETE,
    COUNTS,
    IncompleteIngest,
    SealedHalf,
    ingest,
    load_country_counts,
    load_presence,
    mark_complete,
    qa,
    read_complete,
    resolution_of,
    scan_counts,
    sha256,
)

SCRIPT = Path(__file__).parent.parent / "scripts" / "ingest.py"
spec = importlib.util.spec_from_file_location("ingest_script", SCRIPT)
ingest_script = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ingest_script)

# Ids that no float can tell apart, as in the real files
C1, C2, C3 = "-658282526011250000", "654539408015826335", "654539408015826336"
UNKNOWN = "111111111111111111"
KEY = "0001234-261023074918325"


def row(species, country, cell, part, **counts):
    """One row of the download. A count that isn't given is 1, except that `n_raw` is the larger of `n` and 1."""
    values = {name: 1 for name in COUNTS}
    values["n_raw"] = counts.get("n", 1)
    values.update(counts)
    return (species, country, cell, part, *(values[name] for name in COUNTS))


def tsv(rows, header=None) -> str:
    lines = ["\t".join("" if v is None else str(v) for v in r) for r in rows]
    return "\n".join(["\t".join(header or COLUMNS), *lines]) + "\n"


def write_zip(path: Path, members: dict[str, str]) -> Path:
    with zipfile.ZipFile(path, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    return path


def download(tmp_path, rows, header=None, key=KEY) -> Path:
    # the member is named .csv but is tab-separated, as in GBIF's SQL downloads
    return write_zip(tmp_path / f"{key}.zip", {f"{key}.csv": tsv(rows, header)})


def ingested(tmp_path, rows, chunk_lines=3) -> Path:
    """The parts of an ingest of `rows`, with the marker of a complete ingest that `scripts/ingest.py` writes after its QA."""
    out = tmp_path / "out"
    zip_path = download(tmp_path, rows)
    count = ingest(zip_path, out, chunk_lines)
    mark_complete(out, zip_path, 6, count)
    return out


def grid(cells) -> Grid:
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
    return Grid(6, table, {})


# C1 is all land, with float noise for its sea part. C2 is coastal, and C3 has no land.
GRID = grid([(C1, 1000, 3e-8), (C2, 0.006, 0.26), (C3, 0, 1000)])


def taxa(*species) -> pl.DataFrame:
    return pl.DataFrame({"species_id": list(species)}, schema={"species_id": pl.String})


def split(**halves) -> pl.DataFrame:
    """`split(S1="A", S2="B")`."""
    return pl.DataFrame(
        {"species_id": list(halves), "lineage": [f"family:F{i}" for i in range(len(halves))], "half": list(halves.values())},
        schema={"species_id": pl.String, "lineage": pl.String, "half": pl.String},
    )


# Seven rows, one with every count different and one with no country
SEVEN = [
    row("S1", "DK", C1, "land", n_raw=120, n=100, n_datasets=7, n_years=11, n_nc=13, n_open=17, n_datasets_open=19, n_years_open=23, n_ebird=29, n_inat=31, n_synonym=37, n_human_text=41),
    row("S1", "SE", C1, "land", n=4),
    row("S1", "DK", C2, "sea", n=5),
    row("S2", None, C2, "land", n=6),
    row("S2", "DK", C3, "sea", n=7),
    row("8VR75", "NO", C2, "land", n=8),
    row("8VR75", "DK", C2, "land", n=9),
]


# ingest


def test_ingest_round_trips_seven_rows_in_chunks_of_three(tmp_path):
    out = tmp_path / "out"
    assert ingest(download(tmp_path, SEVEN), out, chunk_lines=3) == 7
    # three parts, and nothing else: the TSV is never written
    assert sorted(p.name for p in out.iterdir()) == ["part-00000.parquet", "part-00001.parquet", "part-00002.parquet"]
    expected = pl.DataFrame(SEVEN, schema=COLUMNS, orient="row")
    parts = [pl.read_parquet(out / f"part-0000{i}.parquet") for i in range(3)]
    assert [p.height for p in parts] == [3, 3, 1]
    assert_frame_equal(pl.concat(parts), expected)


def test_ingest_keeps_ids_as_text_and_empty_fields_as_nulls(tmp_path):
    out = ingested(tmp_path, SEVEN)
    df = pl.read_parquet(out / "part-00000.parquet")
    assert [df.schema[c] for c in ("specieskey", "countrycode", "cell", "part")] == [pl.String] * 4
    assert all(df.schema[c] == pl.Int64 for c in COUNTS)
    assert df["cell"][0] == C1 and float(C2) == float(C3)
    assert pl.read_parquet(out / "part-00001.parquet")["countrycode"].to_list() == [None, "DK", "NO"]


@pytest.mark.parametrize("chunk_lines", [1, 7, 8, 5_000_000])
def test_ingest_gives_the_same_rows_for_any_chunk_size(tmp_path, chunk_lines):
    out = ingested(tmp_path, SEVEN, chunk_lines)
    assert_frame_equal(scan_counts(out).rename({"species_id": "specieskey"}).collect(), pl.DataFrame(SEVEN, schema=COLUMNS, orient="row"))


def test_ingest_replaces_the_parts_of_an_earlier_run(tmp_path):
    out = ingested(tmp_path, SEVEN)
    assert len(list(out.glob("part-*.parquet"))) == 3
    assert ingest(download(tmp_path, SEVEN[:2]), out, chunk_lines=5) == 2
    assert [p.name for p in out.glob("part-*.parquet")] == ["part-00000.parquet"]
    assert scan_counts(out).collect().height == 2
    assert not (out / COMPLETE).exists()  # the earlier run's marker went with its parts


def test_ingest_rejects_other_columns(tmp_path):
    rows = [row("S1", "DK", C1, "land")]
    with pytest.raises(ValueError, match="columns"):
        ingest(download(tmp_path, rows, header=[*list(COLUMNS)[1:], "specieskey"]), tmp_path / "out")
    assert not list((tmp_path / "out").glob("part-*"))


@pytest.mark.parametrize("members", [{}, {"a.csv": "x\n", "b.csv": "y\n"}])
def test_ingest_needs_exactly_one_member(tmp_path, members):
    with pytest.raises(ValueError, match="member"):
        ingest(write_zip(tmp_path / "x.zip", members), tmp_path / "out")


def test_ingest_rejects_a_download_without_rows(tmp_path):
    with pytest.raises(ValueError, match="no rows"):
        ingest(download(tmp_path, []), tmp_path / "out")


def test_ingest_rejects_a_part_that_is_neither_land_nor_sea(tmp_path):
    rows = [row("S1", "DK", C1, "land"), row("S1", "DK", C1, "Land")]
    with pytest.raises(ValueError, match="part"):
        ingest(download(tmp_path, rows), tmp_path / "out")


def test_a_failure_in_a_later_chunk_leaves_no_parts_behind(tmp_path):
    rows = [row("S1", "DK", C1, "land")] * 3 + [row("S1", "DK", C1, "land", n="many")]
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="many"):
        ingest(download(tmp_path, rows), out, chunk_lines=3)
    assert not list(out.glob("part-*"))


def test_the_row_count_may_exceed_one_chunk_several_times_over(tmp_path):
    rows = [row(f"S{i}", "DK", C1, "land", n=i + 1) for i in range(25)]
    out = ingested(tmp_path, rows, chunk_lines=4)
    assert len(list(out.glob("part-*.parquet"))) == 7
    assert scan_counts(out).select(pl.col("n").sum()).collect().item() == sum(range(1, 26))


# the marker of a complete ingest


def test_the_marker_holds_the_rows_the_resolution_and_the_zips_sha256(tmp_path):
    zip_path = download(tmp_path, SEVEN)
    out = tmp_path / "out"
    rows = ingest(zip_path, out, chunk_lines=3)
    assert not (out / COMPLETE).exists()  # ingest alone does not say it is complete: the QA has not run
    written = mark_complete(out, zip_path, 6, rows)
    assert read_complete(out) == written == json.loads((out / COMPLETE).read_text())
    assert (written["rows"], written["resolution"], written["zip"], written["how"]) == (7, 6, zip_path.name, "ingest")
    assert written["zip_sha256"] == sha256(zip_path) and re.fullmatch(r"[0-9a-f]{64}", written["zip_sha256"])
    assert not list(out.glob(f"{COMPLETE}.*"))  # the temporary file is gone


def test_the_loaders_refuse_parts_without_the_marker(tmp_path):
    out = tmp_path / "out"
    ingest(download(tmp_path, SEVEN), out, chunk_lines=3)
    for load in (load_presence, load_country_counts):
        with pytest.raises(IncompleteIngest, match="_COMPLETE"):
            load(out, split(S1="A"), "A")
        with pytest.raises(IncompleteIngest):
            load(out, split(S1="A", S2="B"), "B", unsealed=True)  # nor does unsealing get past it


def test_the_loaders_refuse_parts_that_hold_other_rows_than_the_marker_says(tmp_path):
    out = ingested(tmp_path, SEVEN)
    (out / "part-00002.parquet").unlink()  # as a part lost or a run killed after its marker was copied would leave it
    with pytest.raises(IncompleteIngest, match="6 rows"):
        load_presence(out, split(S1="A"), "A")
    for path in out.glob("part-*.parquet"):
        path.unlink()
    with pytest.raises(IncompleteIngest, match="no part"):
        load_presence(out, split(S1="A"), "A")


def test_a_marker_that_lacks_a_key_is_refused(tmp_path):
    out = ingested(tmp_path, SEVEN)
    marker = json.loads((out / COMPLETE).read_text())
    del marker["resolution"]
    (out / COMPLETE).write_text(json.dumps(marker))
    with pytest.raises(IncompleteIngest, match="resolution"):
        read_complete(out)


def test_a_marker_is_not_written_for_parts_that_hold_other_rows(tmp_path):
    zip_path = download(tmp_path, SEVEN)
    out = tmp_path / "out"
    ingest(zip_path, out, chunk_lines=3)
    with pytest.raises(IncompleteIngest, match="7 rows"):
        mark_complete(out, zip_path, 6, 8)
    assert not (out / COMPLETE).exists()


def test_ingest_removes_the_marker_of_an_earlier_run_before_it_writes_anything(tmp_path):
    out = ingested(tmp_path, SEVEN)
    assert (out / COMPLETE).exists()
    bad = [row("S1", "DK", C1, "land")] * 3 + [row("S1", "DK", C1, "land", n="many")]
    with pytest.raises(ValueError, match="many"):
        ingest(download(tmp_path, bad), out, chunk_lines=3)
    assert not (out / COMPLETE).exists() and not list(out.glob("part-*"))
    # and a run that succeeds leaves no marker either, until the QA has passed
    ingest(download(tmp_path, SEVEN), out, chunk_lines=3)
    assert not (out / COMPLETE).exists()
    with pytest.raises(IncompleteIngest):
        load_presence(out, split(S1="A"), "A")


def test_the_marker_is_gone_before_the_first_part_is_written(tmp_path, monkeypatch):
    """A kill (SIGKILL, out of memory) runs no cleanup, so whatever parts it leaves must already stand without a marker: the earlier run's goes first,
    and the loaders refuse what is left."""
    out = ingested(tmp_path, SEVEN)
    seen = []
    original = pl.DataFrame.write_parquet

    def spy(df, path, *args, **kwargs):
        seen.append((out / COMPLETE).exists())
        return original(df, path, *args, **kwargs)

    monkeypatch.setattr(pl.DataFrame, "write_parquet", spy)
    ingest(download(tmp_path, SEVEN), out, chunk_lines=3)
    assert seen == [False, False, False]
    with pytest.raises(IncompleteIngest):
        load_country_counts(out, split(S1="A"), "A")


def test_the_resolution_is_the_markers_and_a_stated_one_must_agree(tmp_path):
    out = ingested(tmp_path, SEVEN)
    assert resolution_of(out) == 6 == resolution_of(out, 6)
    with pytest.raises(ValueError, match="resolution 6, not 7"):
        resolution_of(out, 7)
    with pytest.raises(IncompleteIngest):
        resolution_of(tmp_path / "nothing")


def load_clean_counts():
    spec = importlib.util.spec_from_file_location("clean_counts_script", Path(__file__).parent.parent / "scripts" / "clean_counts.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_clean_counts_takes_the_resolution_from_the_marker_before_it_loads_anything(tmp_path, monkeypatch):
    script = load_clean_counts()

    def loaded(*args, **kwargs):
        raise AssertionError("something was loaded")

    for name in ("twins_of", "load_grid", "load_references", "load_split", "load_country_counts"):
        monkeypatch.setattr(script, name, loaded)
    out = ingested(tmp_path, SEVEN)
    without = tmp_path / "without"
    ingest(download(tmp_path, SEVEN), without)
    for argv, message in (([str(out), "--resolution", "7"], "resolution 6, not 7"), ([str(without)], "_COMPLETE")):
        monkeypatch.setattr("sys.argv", ["clean_counts.py", "--counts", *argv, "--results", str(tmp_path / "results")])
        with pytest.raises(SystemExit, match=message):
            script.main()
    assert not (tmp_path / "results").exists()
    # with the marker and no resolution stated, it goes on, as far as the reference files that are not there
    monkeypatch.setattr("sys.argv", ["clean_counts.py", "--counts", str(out), "--reference", str(tmp_path / "none"), "--results", str(tmp_path / "results")])
    with pytest.raises(SystemExit, match="ISEA3H-6.gpkg"):
        script.main()


# scan_counts


def test_scan_counts_is_lazy_and_names_the_species_id(tmp_path):
    lazy = scan_counts(ingested(tmp_path, SEVEN))
    assert isinstance(lazy, pl.LazyFrame)
    assert lazy.collect_schema().names()[:4] == ["species_id", "countrycode", "cell", "part"]
    assert "specieskey" not in lazy.collect_schema().names()


def test_scan_counts_without_parts_fails(tmp_path):
    with pytest.raises(FileNotFoundError, match="part"):
        scan_counts(tmp_path)


# the seal and scan_counts

ROOT = Path(__file__).parent.parent
# `scan_counts` reads every species, half B included, without the seal. These are its only readers: `qa`, called by the script, which gives
# global numbers only, and `_select`, which checks the seal first. DECISIONS.md lists them in its seal paragraph, so a new reader needs a ruling.
SCAN_COUNTS_READERS = {"src/realms/occurrences.py", "scripts/ingest.py"}


PARQUET_READERS = {"scan_parquet", "read_parquet", "read_parquet_schema", "read_parquet_metadata", "scan_pyarrow_dataset", "read_table", "ParquetFile", "dataset"}


def _mentions_parts(node: ast.AST, tainted: set[str]) -> bool:
    """Whether `node` names the parts of an ingest: a string with `part-` in it, a name made from one, or the counts' directory (a name or attribute with
    `counts` in it, as `--counts` is the directory of the parts)."""
    for sub in ast.walk(node):
        if isinstance(sub, ast.Constant) and isinstance(sub.value, str) and "part-" in sub.value:
            return True
        if isinstance(sub, ast.Name) and (sub.id in tainted or "counts" in sub.id):
            return True
        if isinstance(sub, ast.Attribute) and "counts" in sub.attr:
            return True
    return False


def direct_part_reads(source: str) -> list[int]:
    """The lines of `source` that read the parts of an ingest with a Parquet reader of their own, past `scan_counts` and the seal."""
    tree = ast.parse(source)
    tainted: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign | ast.AnnAssign | ast.AugAssign) and node.value is not None and _mentions_parts(node.value, set()):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            tainted |= {t.id for target in targets for t in ast.walk(target) if isinstance(t, ast.Name)}
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = node.func.id if isinstance(node.func, ast.Name) else node.func.attr if isinstance(node.func, ast.Attribute) else None
            if name in PARQUET_READERS and any(_mentions_parts(arg, tainted) for arg in [*node.args, *(k.value for k in node.keywords)]):
                found.append(node.lineno)
    return found


def test_the_scan_finds_a_direct_read_of_the_parts():
    caught = [
        'pl.scan_parquet(sorted(d.glob("part-*.parquet")))',
        'parts = list(d.glob("part-*.parquet"))\npl.read_parquet(parts[0])',
        'pl.read_parquet(d / "part-00000.parquet")',
        'pl.scan_parquet(f"{d}/part-{i:05d}.parquet")',
        "pl.scan_parquet(args.counts)",
        "pl.read_parquet(counts_dir)",
        'pyarrow.parquet.read_table(source=str(d / "part-00001.parquet"))',
    ]
    for code in caught:
        assert direct_part_reads(code), code
    for code in ("pl.read_parquet(args.removed)", 'pl.read_parquet(out / "presence.parquet")', 'sum(p.stat().st_size for p in d.glob("part-*.parquet"))'):
        assert not direct_part_reads(code), code


def test_no_module_but_occurrences_reads_the_parts_directly():
    """The parts hold every species, half B included, and only `occurrences.py` reads them, behind the seal (`_select`) or for global numbers (`qa`)."""
    found = {
        path.relative_to(ROOT).as_posix(): lines
        for folder in ("src", "scripts")
        for path in (ROOT / folder).rglob("*.py")
        if path.name != "occurrences.py" and (lines := direct_part_reads(path.read_text()))
    }
    assert found == {}, "a direct read of the parts bypasses the seal on half B: read them through occurrences.py"


def test_only_the_qa_script_and_the_sealed_loaders_read_scan_counts():
    readers = {
        path.relative_to(ROOT).as_posix()
        for folder in ("src", "scripts")
        for path in (ROOT / folder).rglob("*.py")
        if "scan_counts" in path.read_text()
    }
    assert readers == SCAN_COUNTS_READERS, "scan_counts bypasses the seal on half B: add a reader only with a ruling listed in DECISIONS.md"


# qa

QA_ROWS = [
    row("S1", "DK", C1, "land", n_raw=120, n=100, n_synonym=4),
    row("S1", "SE", C1, "land", n_raw=12, n=10, n_synonym=1),
    row("S2", "DK", C2, "sea", n_raw=33, n=30, n_synonym=0),
    row("GONE", "DK", C1, "land", n_raw=55, n=50, n_synonym=5),
    row("GONE", "DK", C2, "land", n_raw=5, n=4),
    row("S1", "DK", UNKNOWN, "land", n_raw=9, n=8),
    row("S2", "DK", "222222222222222222", "sea", n_raw=3, n=2),
    row("S1", "DK", C1, "sea", n_raw=20, n=15, n_synonym=1),  # the sea part of an all-land cell: float noise only
    row("S2", "DK", C3, "land", n_raw=6, n=5),  # land in a cell with no land
    row("S2", "DK", C3, "sea", n_raw=2, n=1),
]


def test_qa_reports_the_global_numbers(tmp_path):
    report = qa(scan_counts(ingested(tmp_path, QA_ROWS)), taxa("S1", "S2"), GRID)
    assert report["rows"] == 10
    assert report["n_raw"] == 265
    assert report["n"] == 225
    assert report["n_synonym"] == 16  # a row that doesn't say has 1
    assert report["synonym_share"] == 16 / 225


def test_qa_synonym_share_is_the_sum_over_the_sum_not_a_mean_of_shares(tmp_path):
    rows = [row("S1", "DK", C1, "land", n=1000, n_synonym=10), row("S1", "SE", C1, "land", n=10, n_synonym=5)]
    assert qa(scan_counts(ingested(tmp_path, rows)), taxa("S1"), GRID)["synonym_share"] == 15 / 1010


def test_qa_shows_species_that_are_not_accepted(tmp_path):
    report = qa(scan_counts(ingested(tmp_path, QA_ROWS)), taxa("S1", "S2"), GRID)
    assert report["unaccepted"] == {"species": 1, "rows": 2, "n_raw": 60, "n": 54}
    assert qa(scan_counts(ingested(tmp_path, QA_ROWS)), taxa("S1", "S2", "GONE"), GRID)["unaccepted"] == {"species": 0, "rows": 0, "n_raw": 0, "n": 0}


def test_qa_a_row_without_a_species_is_not_accepted(tmp_path):
    rows = [row("S1", "DK", C1, "land"), row(None, "DK", C1, "land", n=3)]
    assert qa(scan_counts(ingested(tmp_path, rows)), taxa("S1"), GRID)["unaccepted"] == {"species": 1, "rows": 1, "n_raw": 3, "n": 3}


def test_qa_shows_cell_codes_missing_from_the_grid(tmp_path):
    report = qa(scan_counts(ingested(tmp_path, QA_ROWS)), taxa("S1", "S2"), GRID)
    assert report["unknown_cells"] == {"cells": 2, "rows": 2, "n_raw": 12, "n": 10, "examples": [UNKNOWN, "222222222222222222"]}


def test_qa_compares_cell_codes_as_text(tmp_path):
    # C2 and C3 are the same float, and only one of them is in this grid
    report = qa(scan_counts(ingested(tmp_path, [row("S1", "DK", C3, "sea")])), taxa("S1"), grid([(C2, 1, 1)]))
    assert report["unknown_cells"]["cells"] == 1 and report["unknown_cells"]["examples"] == [C3]


def test_qa_lists_a_few_unknown_codes_at_most(tmp_path):
    rows = [row("S1", "DK", f"9{i:017d}", "land") for i in range(25)]
    unknown = qa(scan_counts(ingested(tmp_path, rows)), taxa("S1"), GRID)["unknown_cells"]
    assert unknown["cells"] == 25 and len(unknown["examples"]) == 10


def test_qa_counts_records_on_parts_of_zero_area(tmp_path):
    report = qa(scan_counts(ingested(tmp_path, QA_ROWS)), taxa("S1", "S2"), GRID)
    # the sea part of C1 has 3e-8 km2, which is none, and its 15 clean records count. The land part of C3 has none, and holds 5.
    assert report["zero_area"]["sea"] == {"cells": 1, "rows": 1, "n_raw": 20, "n": 15}
    assert report["zero_area"]["land"] == {"cells": 1, "rows": 1, "n_raw": 6, "n": 5}


def test_qa_zero_area_is_the_one_threshold_of_the_grid():
    assert occurrences.ZERO_KM2 is grid_module.ZERO_KM2 == 1e-3


def test_qa_a_small_part_is_not_zero_area(tmp_path):
    rows = [row("S1", "DK", C2, "land"), row("S1", "DK", C2, "sea")]  # 0.006 km2 of land and 0.26 km2 of sea, the smallest on the grids
    assert qa(scan_counts(ingested(tmp_path, rows)), taxa("S1"), GRID)["zero_area"] == {
        "land": {"cells": 0, "rows": 0, "n_raw": 0, "n": 0},
        "sea": {"cells": 0, "rows": 0, "n_raw": 0, "n": 0},
    }


def test_qa_does_not_count_an_unknown_cell_as_a_zero_area_part(tmp_path):
    report = qa(scan_counts(ingested(tmp_path, [row("S1", "DK", UNKNOWN, "land")])), taxa("S1"), GRID)
    assert report["zero_area"]["land"]["rows"] == 0 and report["unknown_cells"]["rows"] == 1


def test_qa_names_no_species_and_no_lineage(tmp_path):
    rows = [row("SECRET-ONE", "DK", C1, "land"), row("SECRET-TWO", "DK", C2, "sea"), row("NOT-ACCEPTED", "DK", C1, "land")]
    report = qa(scan_counts(ingested(tmp_path, rows)), taxa("SECRET-ONE", "SECRET-TWO"), GRID)
    text = json.dumps(report)
    assert "SECRET" not in text and "NOT-ACCEPTED" not in text
    assert report["unaccepted"]["species"] == 1


def test_qa_a_download_of_no_records_has_no_synonym_share(tmp_path):
    rows = [row("S1", "DK", C1, "land", n_raw=0, n=0, n_synonym=0)]
    assert qa(scan_counts(ingested(tmp_path, rows)), taxa("S1"), GRID)["synonym_share"] is None


# load_presence


def test_two_countries_in_one_part_sum_their_records_and_take_the_larger_distinct_count(tmp_path):
    rows = [
        row("S1", "DK", C1, "land", n_raw=12, n=10, n_datasets=2, n_years=3, n_nc=1, n_open=5, n_datasets_open=1, n_years_open=3, n_ebird=2, n_inat=3, n_synonym=1, n_human_text=0),
        row("S1", "SE", C1, "land", n_raw=6, n=5, n_datasets=1, n_years=7, n_nc=2, n_open=4, n_datasets_open=3, n_years_open=2, n_ebird=1, n_inat=1, n_synonym=2, n_human_text=1),
    ]
    presence = load_presence(ingested(tmp_path, rows), split(S1="A"), "A")
    assert presence.columns == ["species_id", "part", *COUNTS]
    assert presence.height == 1
    assert presence.row(0, named=True) == {
        "species_id": "S1",
        "part": part_id(C1, "land"),
        "n_raw": 18,
        "n": 15,
        "n_datasets": 2,  # the larger, not 3
        "n_years": 7,  # the larger, not 10
        "n_nc": 3,
        "n_open": 9,
        "n_datasets_open": 3,
        "n_years_open": 3,
        "n_ebird": 3,
        "n_inat": 4,
        "n_synonym": 3,
        "n_human_text": 1,
    }


def test_the_land_and_sea_parts_of_a_cell_stay_apart(tmp_path):
    rows = [row("S1", "DK", C2, "land", n=3), row("S1", "DK", C2, "sea", n=4), row("S1", "SE", C2, "sea", n=5)]
    presence = load_presence(ingested(tmp_path, rows), split(S1="A"), "A")
    assert dict(zip(presence["part"], presence["n"])) == {f"{C2}:land": 3, f"{C2}:sea": 9}


def test_load_presence_returns_counts_as_integers_and_is_sorted(tmp_path):
    presence = load_presence(ingested(tmp_path, SEVEN), split(S1="A", S2="A", **{"8VR75": "A"}), "A")
    assert all(presence.schema[c] == pl.Int64 for c in COUNTS)
    assert presence["species_id"].to_list() == sorted(presence["species_id"].to_list())
    assert presence.sort("species_id", "part").equals(presence)


def test_a_species_outside_the_taxa_is_absent_from_presence(tmp_path):
    out = ingested(tmp_path, QA_ROWS)
    assert qa(scan_counts(out), taxa("S1", "S2"), GRID)["unaccepted"]["species"] == 1
    presence = load_presence(out, split(S1="A", S2="A"), "A")
    assert "GONE" not in presence["species_id"].to_list()
    assert set(presence["species_id"]) == {"S1", "S2"}


def test_half_a_holds_only_the_species_of_half_a(tmp_path):
    out = ingested(tmp_path, SEVEN)
    presence = load_presence(out, split(S1="A", S2="B", **{"8VR75": "B"}), "A")
    assert set(presence["species_id"]) == {"S1"}
    assert presence["n"].sum() == 100 + 4 + 5


def test_the_halves_make_up_all_and_do_not_overlap(tmp_path):
    out = ingested(tmp_path, SEVEN)
    halves = split(S1="A", S2="B", **{"8VR75": "B"})
    a = load_presence(out, halves, "A")
    b = load_presence(out, halves, "B", unsealed=True)
    everything = load_presence(out, halves, "all", unsealed=True)
    assert set(a["species_id"]).isdisjoint(b["species_id"])
    assert_frame_equal(pl.concat([a, b]).sort("species_id", "part"), everything)


@pytest.mark.parametrize("half", ["B", "all"])
def test_a_sealed_half_raises_without_unsealed(tmp_path, half):
    out = ingested(tmp_path, SEVEN)
    with pytest.raises(SealedHalf):
        load_presence(out, split(S1="A", S2="B"), half)
    with pytest.raises(SealedHalf):
        load_country_counts(out, split(S1="A", S2="B"), half)
    with pytest.raises(SealedHalf):
        load_presence(out, split(S1="A", S2="B"), half, unsealed=False)


def test_the_seal_holds_before_anything_is_read(tmp_path):
    # no parts exist, so a read would raise FileNotFoundError
    with pytest.raises(SealedHalf):
        load_presence(tmp_path / "nothing", split(S1="B"), "B")


def test_the_seal_message_names_no_taxon(tmp_path):
    with pytest.raises(SealedHalf) as sealed:
        load_presence(ingested(tmp_path, SEVEN), split(S1="B"), "B")
    assert "S1" not in str(sealed.value) and "unsealed" in str(sealed.value)


def test_sealed_halves_open_with_unsealed(tmp_path):
    out = ingested(tmp_path, SEVEN)
    assert set(load_presence(out, split(S1="A", S2="B"), "B", unsealed=True)["species_id"]) == {"S2"}
    assert set(load_presence(out, split(S1="A", S2="B"), "all", unsealed=True)["species_id"]) == {"S1", "S2"}


def test_half_a_needs_no_unsealed(tmp_path):
    assert load_presence(ingested(tmp_path, SEVEN), split(S1="A", S2="B"), "A").height == 2


def test_an_unknown_half_is_an_error(tmp_path):
    with pytest.raises(ValueError, match="half"):
        load_presence(ingested(tmp_path, SEVEN), split(S1="A"), "C", unsealed=True)


def test_a_split_that_repeats_a_species_is_an_error(tmp_path):
    twice = pl.concat([split(S1="A"), split(S1="B")])
    with pytest.raises(ValueError, match="species_id"):
        load_presence(ingested(tmp_path, SEVEN), twice, "A")


# load_country_counts


def test_country_counts_keep_the_country(tmp_path):
    rows = [
        row("S1", "DK", C1, "land", n=10, n_years=3),
        row("S1", "SE", C1, "land", n=5, n_years=7),
        row("S1", None, C1, "land", n=2, n_years=1),
        row("GONE", "DK", C1, "land"),
    ]
    counts = load_country_counts(ingested(tmp_path, rows), split(S1="A"), "A")
    assert counts.columns == ["species_id", "countrycode", "part", *COUNTS]
    assert counts.height == 3
    by_country = {c: (n, y) for c, n, y in counts.select("countrycode", "n", "n_years").iter_rows()}
    assert by_country == {"DK": (10, 3), "SE": (5, 7), None: (2, 1)}
    assert set(counts["part"]) == {part_id(C1, "land")}


def test_country_counts_are_sealed_like_presence(tmp_path):
    with pytest.raises(SealedHalf):
        load_country_counts(ingested(tmp_path, SEVEN), split(S1="B"), "B")
    assert load_country_counts(ingested(tmp_path, SEVEN), split(S1="B"), "B", unsealed=True)["species_id"].unique().to_list() == ["S1"]


# the script

REFERENCE_FILES = ("grids/ISEA3H-6.gpkg", "natural_earth/ne_10m_land.zip", "exports/extended/col.zip")
CLEAN = [
    row("S1", "DK", C1, "land", n_raw=1000, n=950, n_synonym=30),
    row("S1", "SE", C1, "land", n_raw=100, n=50, n_synonym=0),
    row("S2", "DK", C2, "sea", n_raw=10, n=10, n_synonym=0),
]


@pytest.fixture
def run(tmp_path, monkeypatch):
    """Runs the script's `main` on a fixture download, with the reference data replaced by the fixtures above."""
    for name in REFERENCE_FILES:
        (tmp_path / "reference" / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / "reference" / name).touch()
    monkeypatch.setattr(ingest_script, "load_grid", lambda gpkg, land: GRID)
    monkeypatch.setattr(ingest_script, "read_centres", lambda gpkg: GRID.cells.select("cell", "lat", "lon"))
    monkeypatch.setattr(ingest_script, "load_taxa", lambda path: taxa("S1", "S2"))

    def go(rows=CLEAN, *args):
        zip_path = download(tmp_path, rows)
        argv = [str(zip_path), "--resolution", "6", "--reference", str(tmp_path / "reference"), "--results", str(tmp_path / "results")]
        return ingest_script.main([*argv, *args])

    return go


def report_of(tmp_path) -> dict:
    return json.loads((tmp_path / "results" / "res6" / "qa.json").read_text())


def test_the_script_writes_the_parts_and_qa_json(run, tmp_path):
    run()
    assert [p.name for p in (tmp_path / KEY).glob("part-*.parquet")] == ["part-00000.parquet"]
    report = report_of(tmp_path)
    assert report["rows"] == 3 and report["n"] == 1010 and report["synonym_share"] == 30 / 1010
    assert report["unaccepted"]["species"] == 0 and report["unknown_cells"]["cells"] == 0
    assert report["ingest"]["rows"] == 3 and report["ingest"]["seconds"] >= 0 and report["ingest"]["parquet_bytes"] > 0
    assert report["ingest"]["peak_rss_mb"] > 0


def test_the_script_writes_the_marker_once_the_qa_has_passed(run, tmp_path):
    run()
    out = tmp_path / KEY
    marker = read_complete(out)
    assert (marker["rows"], marker["resolution"], marker["how"]) == (3, 6, "ingest")
    assert marker["zip_sha256"] == sha256(tmp_path / f"{KEY}.zip")
    assert load_presence(out, split(S1="A", S2="A"), "A").height == 2


def test_the_script_writes_no_marker_when_a_pass_mark_fails(run, tmp_path):
    run()
    with pytest.raises(SystemExit):
        run([*CLEAN, row("GONE", "DK", C1, "land", n=40)])
    assert list((tmp_path / KEY).glob("part-*.parquet")) and not (tmp_path / KEY / COMPLETE).exists()  # the first run's marker is gone too


def test_mark_complete_writes_the_marker_for_an_earlier_ingest_without_ingesting_again(run, tmp_path):
    run()
    out = tmp_path / KEY
    (out / COMPLETE).unlink()  # an ingest made before the marker existed
    parts = {p.name: p.stat().st_mtime_ns for p in out.glob("part-*.parquet")}
    run(CLEAN, "--mark-complete")
    assert {p.name: p.stat().st_mtime_ns for p in out.glob("part-*.parquet")} == parts  # nothing ingested
    marker = read_complete(out)
    assert (marker["rows"], marker["resolution"], marker["how"]) == (3, 6, "mark-complete")
    assert marker["zip_sha256"] == sha256(tmp_path / f"{KEY}.zip")
    with pytest.raises(SystemExit, match="already"):
        run(CLEAN, "--mark-complete")


def earlier_ingest(run, tmp_path) -> Path:
    run()
    out = tmp_path / KEY
    (out / COMPLETE).unlink()
    return out


def test_mark_complete_refuses_parts_that_hold_other_rows_than_qa_json(run, tmp_path):
    out = earlier_ingest(run, tmp_path)
    qa_path = tmp_path / "results" / "res6" / "qa.json"
    report = json.loads(qa_path.read_text())
    qa_path.write_text(json.dumps({**report, "rows": 4, "ingest": {**report["ingest"], "rows": 4}}))
    with pytest.raises(SystemExit, match="3 rows"):
        run(CLEAN, "--mark-complete")
    assert not (out / COMPLETE).exists()


def test_mark_complete_refuses_a_qa_json_that_fails_the_pass_marks(run, tmp_path):
    out = earlier_ingest(run, tmp_path)
    qa_path = tmp_path / "results" / "res6" / "qa.json"
    qa_path.write_text(json.dumps({**json.loads(qa_path.read_text()), "synonym_share": 0.5}))
    with pytest.raises(SystemExit, match="synonym"):
        run(CLEAN, "--mark-complete")
    assert not (out / COMPLETE).exists()


def test_mark_complete_refuses_parts_with_a_cell_of_another_grid(run, tmp_path, monkeypatch):
    out = earlier_ingest(run, tmp_path)
    monkeypatch.setattr(ingest_script, "load_grid", lambda gpkg, land: grid([(C1, 1000, 0)]))  # C2 is not in it, as a grid of another resolution
    with pytest.raises(SystemExit, match="not in the grid"):
        run(CLEAN, "--mark-complete")
    assert not (out / COMPLETE).exists()


def test_the_script_prints_global_numbers_only(run, capsys):
    run()
    out = capsys.readouterr().out
    assert "1,010" in out and "S1" not in out and "S2" not in out


def test_the_qa_file_names_no_species(run, tmp_path):
    run()
    text = (tmp_path / "results" / "res6" / "qa.json").read_text()
    assert '"S1"' not in text and '"S2"' not in text


def failure(run, rows):
    with pytest.raises(SystemExit) as failed:
        run(rows)
    assert isinstance(failed.value.code, str)
    return failed.value.code


def test_the_script_fails_on_an_unaccepted_species_and_still_writes_the_numbers(run, tmp_path):
    message = failure(run, [*CLEAN, row("GONE", "DK", C1, "land", n=40)])
    assert "unaccepted" in message and "S1" not in message and "GONE" not in message
    assert report_of(tmp_path)["unaccepted"]["rows"] == 1


def test_the_script_fails_on_an_unknown_cell(run, tmp_path):
    assert "unknown" in failure(run, [*CLEAN, row("S1", "DK", UNKNOWN, "land")])
    assert report_of(tmp_path)["unknown_cells"]["rows"] == 1


@pytest.mark.parametrize("synonyms", [0, 19, 51, 200])
def test_the_script_fails_on_a_synonym_share_outside_two_to_five_percent(run, synonyms):
    # 1000 records of S1 and 10 of S2, so 20 synonyms is just under 2%, and 51 just over 5%
    rows = [row("S1", "DK", C1, "land", n=990, n_synonym=synonyms), row("S2", "DK", C2, "sea", n=10, n_synonym=0)]
    assert "synonym" in failure(run, rows)


@pytest.mark.parametrize("synonyms", [20, 36, 50])
def test_the_script_passes_a_synonym_share_from_two_to_five_percent(run, synonyms):
    run([row("S1", "DK", C1, "land", n=990, n_synonym=synonyms), row("S2", "DK", C2, "sea", n=10, n_synonym=0)])


def test_the_script_reports_every_failed_mark_at_once(run):
    message = failure(run, [row("GONE", "DK", UNKNOWN, "land", n=10, n_synonym=0)])
    assert "unaccepted" in message and "unknown" in message and "synonym" in message


def test_the_script_checks_the_reference_files_before_ingesting(run, tmp_path):
    (tmp_path / "reference" / "exports" / "extended" / "col.zip").unlink()
    with pytest.raises(SystemExit) as failed:
        run()
    assert "col.zip" in failed.value.code
    assert not (tmp_path / KEY).exists()


def test_the_script_refuses_a_grid_file_of_another_resolution(run, monkeypatch):
    monkeypatch.setattr(ingest_script, "load_grid", lambda gpkg, land: Grid(7, GRID.cells, {}))
    with pytest.raises(SystemExit) as failed:
        run()
    assert "resolution" in failed.value.code


def test_the_pass_marks_are_those_of_the_plan():
    none = {"rows": 0, "n_raw": 0, "n": 0}
    some = {"rows": 1, "n_raw": 2, "n": 1}
    clean = {"unaccepted": {"species": 0, **none}, "unknown_cells": {"cells": 0, **none}, "synonym_share": 0.036}
    assert ingest_script.problems(clean) == []
    assert ingest_script.problems({**clean, "synonym_share": 0.02}) == []
    assert ingest_script.problems({**clean, "synonym_share": 0.05}) == []
    assert len(ingest_script.problems({**clean, "synonym_share": 0.0199})) == 1
    assert len(ingest_script.problems({**clean, "synonym_share": 0.0501})) == 1
    assert len(ingest_script.problems({**clean, "synonym_share": None})) == 1
    assert len(ingest_script.problems({**clean, "unaccepted": {"species": 1, **some}})) == 1
    assert len(ingest_script.problems({**clean, "unknown_cells": {"cells": 1, **some}})) == 1


def test_the_counts_and_columns_are_those_of_the_query():
    assert len(COUNTS) == 12
    assert list(COLUMNS) == ["specieskey", "countrycode", "cell", "part", *COUNTS]
    sql = (Path(__file__).parent.parent / "queries" / "occurrence_res6.sql").read_text()
    assert re.findall(r"\bAS (\w+)", sql) == list(COLUMNS)[2:]


# near-miss cell codes: GBIF's code for a point differs from the grid file's id of the same cell by a few hundred metres in a handful of cells


def cell_id(code: int, lat: float, lon: float, negative: bool = False) -> str:
    """A cell code as GBIF writes it: the 2-digit code, |lat| in millionths of a degree in 8 digits and |lon| in 9, as an integer (a code below 10 loses its zero)."""
    digits = f"{code:02d}{round(abs(lat) * 1e6):08d}{round(abs(lon) * 1e6):09d}"
    return ("-" if negative else "") + str(int(digits))


def centres_of(places) -> pl.DataFrame:
    """`places`: (cell, lat, lon) of the cells of a made-up grid, as `read_centres` gives them."""
    return pl.DataFrame(
        {"cell": [c for c, _, _ in places], "lat": [float(la) for _, la, _ in places], "lon": [float(lo) for _, _, lo in places]},
        schema={"cell": pl.String, "lat": pl.Float64, "lon": pl.Float64},
    )


def place_grid(places) -> Grid:
    """The same cells as a `Grid`, all land, for the QA. Its resolution is 6, which is what the script fixture asks for: the codes are made up."""
    table = centres_of(places).with_columns(area_km2=pl.lit(1000.0), land_km2=pl.lit(1000.0), sea_km2=pl.lit(0.0))
    return Grid(6, table, {})


# Code 07 is north and east, 29 south and east, 51 north and west and 73 south and west, as on the real grids
CELL_A, FAR_A = cell_id(7, 10.0, 20.0), cell_id(7, 10.5, 20.0)  # FAR_A is 55 km north of A
CELL_S, CELL_W, CELL_SW = cell_id(29, -35.0, 150.0), cell_id(51, 45.0, -100.0), cell_id(73, -35.0, -70.0)
PLACES = [
    (CELL_A, 10.0, 20.0),
    (FAR_A, 10.5, 20.0),
    (CELL_S, -35.0, 150.0),
    (CELL_W, 45.0, -100.0),
    (CELL_SW, -35.0, -70.0),
    # A's coordinates, with the signs of the other three quadrants: a decoder that reads the digits and not the code would find them as near as A
    (cell_id(29, -10.0, 20.0), -10.0, 20.0),
    (cell_id(51, 10.0, -20.0), 10.0, -20.0),
    (cell_id(73, -10.0, -20.0), -10.0, -20.0),
]
NEAR_A = cell_id(7, 10.003, 20.0)  # 0.33 km north of A
NEAR_S = cell_id(29, -35.0, 150.004)  # 0.36 km east of S
NEAR_W = cell_id(51, 45.0, -100.003)  # 0.24 km west of W
NEAR_SW = cell_id(73, -35.003, -70.0)  # 0.33 km south of SW
KM = 111.195  # a degree of latitude


def parts_of(out) -> pl.DataFrame:
    return pl.concat([pl.read_parquet(p) for p in sorted(out.glob("part-*.parquet"))])


def remap_of(tmp_path, code, places, chunk_lines=5):
    """The ingest of one row in `code`, with the cells of `places`: what it was remapped to, or None."""
    out = tmp_path / "out"
    done = occurrences.ingest_remapped(download(tmp_path, [row("S1", "DK", code, "land")]), out, centres_of(places), chunk_lines)
    assert done.rows == done.downloaded == 1
    cell = parts_of(out)["cell"].item()
    assert cell == code or [r.target for r in done.remaps] == [cell]
    return None if cell == code else cell


def test_a_near_miss_code_is_remapped_to_the_grid_cell_with_its_code_and_nearest_centre(tmp_path):
    rows = [
        row("S1", "DK", NEAR_A, "land"),
        row("S1", "DK", NEAR_S, "land"),
        row("S1", "DK", NEAR_W, "land"),
        row("S1", "DK", NEAR_SW, "land"),
        row("S2", "DK", CELL_A, "land"),  # a code of the grid stays as it is
    ]
    out = tmp_path / "out"
    done = occurrences.ingest_remapped(download(tmp_path, rows), out, centres_of(PLACES), chunk_lines=2)
    assert (done.rows, done.downloaded) == (5, 5)
    found = parts_of(out)
    assert set(zip(found["specieskey"], found["cell"], strict=True)) == {("S1", CELL_A), ("S1", CELL_S), ("S1", CELL_SW), ("S1", CELL_W), ("S2", CELL_A)}
    assert {r.code: r.target for r in done.remaps} == {NEAR_A: CELL_A, NEAR_S: CELL_S, NEAR_W: CELL_W, NEAR_SW: CELL_SW}
    assert set(parts_of(out)["cell"]) <= {c for c, _, _ in PLACES}


def test_the_target_is_a_cell_with_the_same_code_and_not_the_nearest_cell_of_any_code(tmp_path):
    # 0.4 km north of the equator in code 07 there is no cell, and 0.4 km south of it a cell of code 29: another cell, and not a rounding of this one
    places = [(cell_id(7, 0.5, 30.0), 0.5, 30.0), (cell_id(29, -0.001, 30.0), -0.001, 30.0)]
    assert remap_of(tmp_path, cell_id(7, 0.003, 30.0), places) is None


@pytest.mark.parametrize("negative", [False, True])
def test_the_sign_of_a_code_is_dropped_before_it_is_read(tmp_path, negative):
    # twelve cells of a grid have a negative id (here, one at 58.28 N, 11.25 E): the sign is not part of what the digits say
    vertex = cell_id(7, 58.282526, 11.25, negative=True)
    places = [(vertex, 58.282526, 11.25), (cell_id(7, 58.8, 11.25), 58.8, 11.25)]
    assert remap_of(tmp_path, cell_id(7, 58.2825, 11.25, negative=negative), places) == vertex


@pytest.mark.parametrize(("km", "remapped"), [(4.9, True), (5.1, False), (30.0, False)])
def test_a_near_miss_is_within_five_kilometres(tmp_path, km, remapped):
    places = [(CELL_A, 10.0, 20.0)]
    assert (remap_of(tmp_path, cell_id(7, 10.0 + km / KM, 20.0), places) == CELL_A) is remapped


@pytest.mark.parametrize(("second_km", "remapped"), [(2.5, False), (2.9, False), (3.2, True), (4.0, True)])
def test_the_next_cell_must_be_more_than_three_times_as_far(tmp_path, second_km, remapped):
    near, second = cell_id(7, 10.009, 20.0), cell_id(7, round(10.0 - second_km / KM, 6), 20.0)  # 1.0 km and `second_km` from (10, 20)
    places = [(near, 10.009, 20.0), (second, round(10.0 - second_km / KM, 6), 20.0)]
    assert (remap_of(tmp_path, cell_id(7, 10.0, 20.0), places) == near) is remapped


def test_a_code_in_the_grid_is_never_remapped_even_if_another_cell_is_next_to_it(tmp_path):
    places = [(CELL_A, 10.0, 20.0), (cell_id(7, 10.003, 20.0), 10.003, 20.0)]
    assert remap_of(tmp_path, CELL_A, places) is None
    assert remap_of(tmp_path, cell_id(7, 10.003, 20.0), places) is None


def test_two_cells_of_one_centre_are_ambiguous(tmp_path):
    places = [(CELL_A, 10.0, 20.0), (cell_id(7, 10.0, 20.0, negative=True), 10.0, 20.0)]
    assert remap_of(tmp_path, cell_id(7, 10.001, 20.0), places) is None


@pytest.mark.parametrize(
    "code",
    [
        "not-a-cell",
        "7" * 25,  # more than 19 digits
        "²⁰" + "0" * 17,  # digits that are not ASCII
        cell_id(7, 90.003, 0.0),  # beyond the pole, 0.33 km from the cell at it on a sphere that goes on past the pole
        cell_id(51, 10.0, 180.003),  # beyond the antimeridian, 0.3 km from the cell at -180 on a sphere that goes on past it
        cell_id(8, 10.003, 20.0),  # a code that no cell of the grid has, so that its hemisphere is not known
    ],
)
def test_a_code_that_does_not_decode_to_a_place_is_not_remapped(tmp_path, code):
    edges = [(cell_id(7, 90.0, 0.0), 90.0, 0.0), (cell_id(51, 10.0, -180.0), 10.0, -180.0)]
    assert remap_of(tmp_path, code, [*PLACES, *edges]) is None


def test_decoding_reads_the_code_the_latitude_and_the_longitude():
    assert occurrences.decode("704809617102037511") == ("07", 4.809617, 102.037511)  # 18 digits: the first zero is gone
    assert occurrences.decode("-758282526011250000") == ("07", 58.282526, 11.25)
    assert occurrences.decode("5100000000079533471") == ("51", 0.0, 79.533471)
    assert occurrences.decode("7300815160144987711") == ("73", 0.815160, 144.987711)
    assert occurrences.decode("") is None and occurrences.decode("-") is None


def test_a_code_whose_hemisphere_the_grid_does_not_show_is_not_remapped(tmp_path):
    # the cells of code 07 straddle the equator, so the hemisphere of the code is not known from the grid: nothing is remapped
    places = [(cell_id(7, 5.0, 20.0), 5.0, 20.0), (cell_id(7, -5.0, 20.0), -5.0, 20.0)]
    assert remap_of(tmp_path, cell_id(7, 5.003, 20.0), places) is None


def test_rows_without_a_cell_and_codes_that_are_not_near_a_cell_stay_as_they_are(tmp_path):
    rows = [row("S1", "DK", None, "land"), row("S1", "DK", UNKNOWN, "land"), row("S1", "DK", NEAR_A, "land")]
    out = tmp_path / "out"
    done = occurrences.ingest_remapped(download(tmp_path, rows), out, centres_of(PLACES))
    assert (done.rows, done.downloaded, [r.code for r in done.remaps]) == (3, 3, [NEAR_A])
    assert sorted(parts_of(out)["cell"].to_list(), key=str) == sorted([None, UNKNOWN, CELL_A], key=str)


def test_ingest_without_the_cells_of_a_grid_remaps_nothing(tmp_path):
    out = tmp_path / "out"
    assert ingest(download(tmp_path, [row("S1", "DK", NEAR_A, "land")]), out) == 1
    assert parts_of(out)["cell"].to_list() == [NEAR_A]


def test_a_remap_lists_its_code_target_rows_records_and_distance(tmp_path):
    rows = [
        row("S1", "DK", NEAR_A, "land", n_raw=12, n=10),
        row("S2", "SE", NEAR_A, "land", n_raw=8, n=7),
        row("S1", "DK", NEAR_S, "land", n_raw=4, n=3),
    ]
    done = occurrences.ingest_remapped(download(tmp_path, rows), tmp_path / "out", centres_of(PLACES))
    first, second = sorted(done.remaps, key=lambda r: r.code == NEAR_S)  # they come sorted by code, as text
    assert (first.code, first.target, first.rows, first.n_raw, first.n) == (NEAR_A, CELL_A, 2, 20, 17)
    assert (second.code, second.target, second.rows, second.n_raw, second.n) == (NEAR_S, CELL_S, 1, 4, 3)
    assert first.km == pytest.approx(0.334, abs=1e-3) and second.km == pytest.approx(0.364, abs=1e-3)
    assert [r.code for r in done.remaps] == sorted([NEAR_A, NEAR_S])


def counts(**kw):
    return {"n_raw": 12, "n": 10, "n_datasets": 2, "n_years": 3, "n_nc": 1, "n_open": 5, "n_datasets_open": 1, "n_years_open": 3, "n_ebird": 2, "n_inat": 3, "n_synonym": 1, "n_human_text": 0, **kw}


NEAR_A2 = cell_id(7, 10.0, 20.003)  # a second code for A, 0.33 km east of it


def test_rows_that_collide_after_a_remap_are_merged_as_the_loaders_merge_them(tmp_path):
    rows = [
        row("S1", "DK", CELL_A, "land", **counts()),  # chunk 0
        row("S1", "SE", CELL_A, "land", **counts(n=1)),
        row("S1", "DK", CELL_S, "land"),
        row("S3", None, CELL_A, "land", **counts(n=2)),  # chunk 1, with row 5 after the remap, a row with no country
        row("S1", "DK", NEAR_A, "land", **counts(n_raw=6, n=5, n_datasets=1, n_years=7, n_nc=2, n_open=4, n_datasets_open=3, n_years_open=2, n_ebird=1, n_inat=1, n_synonym=2, n_human_text=1)),
        row("S3", None, NEAR_A, "land", **counts(n=3)),
        row("S1", "DK", NEAR_A2, "land", **counts(n_raw=1, n=1, n_datasets=9, n_years=1, n_nc=1, n_open=1, n_datasets_open=1, n_years_open=1, n_ebird=1, n_inat=1, n_synonym=0, n_human_text=1)),  # chunk 2, a second code, onto the same row as row 4
        row("S4", "DK", NEAR_A, "land"),  # collides with nothing
    ]
    out = tmp_path / "out"
    done = occurrences.ingest_remapped(download(tmp_path, rows), out, centres_of(PLACES), chunk_lines=3)
    assert (done.downloaded, done.rows) == (8, 5)  # three rows went into others
    merged = parts_of(out)
    assert merged.height == 5 and merged.select("specieskey", "countrycode", "cell", "part").unique().height == 5
    assert occurrences.part_rows(out) == 5
    one = merged.filter((pl.col("specieskey") == "S1") & (pl.col("countrycode") == "DK") & (pl.col("cell") == CELL_A)).row(0, named=True)
    # the counts add up, and the four distinct counts take the largest of the three rows
    assert {k: one[k] for k in COUNTS} == {
        "n_raw": 12 + 6 + 1,
        "n": 10 + 5 + 1,
        "n_datasets": 9,
        "n_years": 7,
        "n_nc": 1 + 2 + 1,
        "n_open": 5 + 4 + 1,
        "n_datasets_open": 3,
        "n_years_open": 3,
        "n_ebird": 2 + 1 + 1,
        "n_inat": 3 + 1 + 1,
        "n_synonym": 1 + 2 + 0,
        "n_human_text": 0 + 1 + 1,
    }
    assert merged.filter(pl.col("specieskey") == "S3")["n"].to_list() == [5]
    # and it is what the loaders make of the same download with the codes of the grid in it
    replaced = [(sp, c, CELL_A if cell in (NEAR_A, NEAR_A2) else cell, *rest) for sp, c, cell, *rest in rows]
    for name, source in (("remapped", rows), ("replaced", replaced)):
        (tmp_path / name).mkdir()
        zip_path = download(tmp_path / name, source)
        if name == "remapped":
            count = occurrences.ingest_remapped(zip_path, tmp_path / name / "out", centres_of(PLACES), 3).rows
        else:
            count = ingest(zip_path, tmp_path / name / "out", 3)
        mark_complete(tmp_path / name / "out", zip_path, 7, count)
    by = split(S1="A", S3="A", S4="A")
    for loader in (load_country_counts, load_presence):
        assert_frame_equal(loader(tmp_path / "remapped" / "out", by, "A"), loader(tmp_path / "replaced" / "out", by, "A"))


def test_a_remap_ingest_with_no_collision_keeps_every_row(tmp_path):
    rows = [row("S1", "DK", CELL_A, "land"), row("S2", "DK", NEAR_A, "land"), row("S2", "DK", NEAR_A, "sea")]
    done = occurrences.ingest_remapped(download(tmp_path, rows), tmp_path / "out", centres_of(PLACES), chunk_lines=2)
    assert (done.downloaded, done.rows) == (3, 3)


def test_a_part_that_was_being_rewritten_when_an_earlier_run_was_killed_is_removed(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "part-00000.parquet.tmp").write_bytes(b"half a part")
    occurrences.ingest_remapped(download(tmp_path, [row("S1", "DK", CELL_A, "land")]), out, centres_of(PLACES))
    assert sorted(p.name for p in out.iterdir()) == ["part-00000.parquet"]


def test_more_near_misses_than_the_cap_stop_the_ingest_and_leave_no_parts(tmp_path, monkeypatch):
    monkeypatch.setattr(occurrences, "MAX_REMAPS", 2)
    out = tmp_path / "out"
    rows = [row("S1", "DK", code, "land") for code in (NEAR_A, NEAR_S, NEAR_W)]
    with pytest.raises(ValueError, match="3 cell codes"):
        occurrences.ingest_remapped(download(tmp_path, rows), out, centres_of(PLACES), chunk_lines=2)
    assert not list(out.glob("part-*.parquet"))
    done = occurrences.ingest_remapped(download(tmp_path, rows[:2]), out, centres_of(PLACES))  # two are within the cap
    assert len(done.remaps) == 2


def test_a_failed_remap_leaves_no_parts_behind(tmp_path):
    out = tmp_path / "out"
    rows = [row("S1", "DK", NEAR_A, "land"), row("S1", "DK", CELL_A, "water")]  # `water` is neither land nor sea
    with pytest.raises(ValueError, match="part is one of"):
        occurrences.ingest_remapped(download(tmp_path, rows), out, centres_of(PLACES), chunk_lines=1)
    assert not list(out.glob("part-*.parquet"))


REAL_GRIDS = Path(__file__).parent.parent / "data" / "reference" / "grids"


@pytest.mark.parametrize("resolution", [5, 6, 7])
def test_decoding_every_cell_of_a_real_grid_gives_its_centre(resolution):
    gpkg = REAL_GRIDS / f"ISEA3H-{resolution}.gpkg"
    if not gpkg.exists():
        pytest.skip(f"{gpkg} is not there")
    cells = grid_module.read_centres(gpkg)
    centres = occurrences.Centres(cells)
    wrong = [
        cell
        for cell, lat, lon in cells.iter_rows()
        if (placed := centres.place(cell)) is None or abs(placed[0] - lat) > 1e-6 or abs(placed[1] - lon) > 1e-6
    ]
    assert cells.height > 2000 and wrong == []


def test_the_seven_near_misses_of_the_resolution_7_download_find_their_cells():
    gpkg = REAL_GRIDS / "ISEA3H-7.gpkg"
    if not gpkg.exists():
        pytest.skip(f"{gpkg} is not there")
    centres = occurrences.Centres(grid_module.read_centres(gpkg))
    expected = {
        "2902405731102034479": ("2902411730102033531", 0.67),
        "5100000000079533471": ("5100000000079532484", 0.11),
        "5121434201079619262": ("5121436312079621336", 0.32),
        "704809617102037511": ("704821168102036680", 1.29),
        "714370309102070552": ("714390873102070952", 2.28),
        "7300815160144987711": ("7300817749144998464", 1.23),
        "7302405731079534479": ("7302411730079533531", 0.67),
    }
    for code, (target, km) in expected.items():
        found = centres.nearest(code)
        assert found is not None and found[0] == target and found[1] == pytest.approx(km, abs=0.011), code


# the script with near misses


@pytest.fixture
def near_run(run, monkeypatch):
    """`run` with the grid of `places` for both the centres the ingest reads and the grid the QA loads."""

    def go(rows, places=PLACES):
        monkeypatch.setattr(ingest_script, "load_grid", lambda gpkg, land: place_grid(places))
        monkeypatch.setattr(ingest_script, "read_centres", lambda gpkg: centres_of(places))
        return run(rows)

    return go


BASE = row("S1", "DK", CELL_S, "land", n_raw=1000, n=950, n_synonym=30)


def test_remaps_are_counted_in_qa_json(near_run, tmp_path, capsys):
    near_run([BASE, row("S1", "DK", NEAR_A, "land", n_raw=12, n=10, n_synonym=0), row("S2", "SE", NEAR_A, "land", n_raw=8, n=7, n_synonym=0), row("S2", "DK", NEAR_W, "land", n_raw=1, n=1, n_synonym=0)])
    report = report_of(tmp_path)
    remaps = report["remaps"]
    assert (remaps["codes"], remaps["rows"], remaps["n_raw"], remaps["n"], remaps["merged_rows"]) == (2, 3, 21, 18, 0)
    listed = {r["code"]: r for r in remaps["list"]}  # sorted by code, as text
    assert list(listed) == sorted([NEAR_A, NEAR_W])
    assert [(r["target"], r["rows"], r["n_raw"], r["n"]) for r in (listed[NEAR_A], listed[NEAR_W])] == [(CELL_A, 2, 20, 17), (CELL_W, 1, 1, 1)]
    assert listed[NEAR_A]["km"] == pytest.approx(0.334, abs=1e-3) and listed[NEAR_W]["km"] == pytest.approx(0.236, abs=1e-3)
    # the QA runs on the parts, which hold no unknown code, and the ingest passes
    assert report["unknown_cells"]["rows"] == 0 and report["rows"] == report["ingest"]["rows"] == 4
    assert read_complete(tmp_path / KEY)["rows"] == 4
    assert set(parts_of(tmp_path / KEY)["cell"]) == {CELL_S, CELL_A, CELL_W}
    printed = capsys.readouterr().out
    assert NEAR_A in printed and CELL_A in printed and "S1" not in printed and "S2" not in printed


def test_qa_json_counts_the_rows_the_remaps_merged(near_run, tmp_path):
    near_run([BASE, row("S1", "DK", CELL_A, "land", n_synonym=0), row("S1", "DK", NEAR_A, "land", n_synonym=0)])
    report = report_of(tmp_path)
    assert report["remaps"]["merged_rows"] == 1 and report["rows"] == report["ingest"]["rows"] == 2 and report["ingest"]["downloaded_rows"] == 3
    assert read_complete(tmp_path / KEY)["rows"] == 2


def test_a_download_without_near_misses_lists_no_remap(run, tmp_path):
    run()
    assert report_of(tmp_path)["remaps"] == {"codes": 0, "rows": 0, "n_raw": 0, "n": 0, "merged_rows": 0, "list": []}


def test_a_code_with_no_grid_cell_within_tolerance_still_fails_qa(near_run, tmp_path):
    far = cell_id(7, 10.0 + 6.0 / KM, 20.0)  # 6 km from A
    message = failure(near_run, [BASE, row("S1", "DK", far, "land", n_synonym=0)])
    assert "unknown cells: 1 codes in 1 rows" in message
    report = report_of(tmp_path)
    assert report["unknown_cells"]["examples"] == [far] and report["remaps"]["codes"] == 0
    assert not (tmp_path / KEY / COMPLETE).exists()


def test_an_ambiguous_near_miss_still_fails_qa(near_run, tmp_path):
    two = [(cell_id(7, 10.009, 20.0), 10.009, 20.0), (cell_id(7, round(10.0 - 2.5 / KM, 6), 20.0), round(10.0 - 2.5 / KM, 6), 20.0)]  # 1.0 and 2.5 km away
    unknown = cell_id(7, 10.0, 20.0)
    with pytest.raises(SystemExit) as failed:
        near_run([BASE, row("S1", "DK", unknown, "land", n_synonym=0)], [*two, (CELL_S, -35.0, 150.0)])
    assert "unknown cells: 1 codes in 1 rows" in failed.value.code
    assert report_of(tmp_path)["unknown_cells"]["examples"] == [unknown] and report_of(tmp_path)["remaps"]["codes"] == 0
    assert not (tmp_path / KEY / COMPLETE).exists()
