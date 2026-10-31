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
from realms.occurrences import COLUMNS, COUNTS, SealedHalf, ingest, load_country_counts, load_presence, qa, scan_counts

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
    out = tmp_path / "out"
    ingest(download(tmp_path, rows), out, chunk_lines)
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
