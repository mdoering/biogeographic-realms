import zipfile

import polars as pl
import pytest

from realms.reference import griis_datasets, latest_export, load_tax_groups, load_taxa, rank_coverage, repeated_names


def export(key, status="finished", created="2026-07-01", size=10, **request):
    return {
        "key": key,
        "status": status,
        "created": created,
        "size": size,
        "download": f"https://download.example/{key}.zip",
        "request": {"datasetKey": 1, "format": "COLDP", "synonyms": True, "tabFormat": "tsv", **request},
    }


def test_latest_export_matches_flags_exactly(api):
    api.routes["/export"] = (
        200,
        {"total": 3, "result": [export("ext", extended=True), export("plain", created="2026-06-01"), export("cls", created="2026-05-01", classification=True)]},
    )
    assert latest_export(1, classification=True, api=api.base).key == "cls"
    assert latest_export(1, api=api.base).key == "plain"
    assert latest_export(1, extended=True, api=api.base).key == "ext"
    assert "datasetKey=1" in api.requests[0][1] and "format=COLDP" in api.requests[0][1]


def test_latest_export_takes_the_newest(api):
    api.routes["/export"] = (200, {"total": 2, "result": [export("old", created="2026-01-01"), export("new", created="2026-09-01")]})
    found = latest_export(1, api=api.base)
    assert (found.key, found.url, found.size, found.created) == ("new", "https://download.example/new.zip", 10, "2026-09-01")


def test_latest_export_skips_unfinished(api):
    api.routes["/export"] = (200, {"total": 1, "result": [export("running", status="running")]})
    assert latest_export(1, api=api.base) is None


def test_latest_export_needs_synonyms(api):
    # an export without synonyms can't resolve what GBIF matched to
    row = export("nosyn")
    row["request"]["synonyms"] = False
    api.routes["/export"] = (200, {"total": 1, "result": [row]})
    assert latest_export(1, api=api.base) is None


def test_latest_export_none_when_absent_and_does_not_post(api):
    api.routes["/export"] = (200, {"total": 0, "result": []})
    assert latest_export(26853, api=api.base) is None
    assert all(method == "GET" for method, _, _ in api.requests)


def griis_listing(api, results):
    api.routes["/dataset"] = (200, {"total": len(results), "last": True, "result": results})


def test_griis_parses_country_and_subdivision(api):
    griis_listing(
        api,
        [
            {"key": 20359, "alias": "GRIIS-US-Hawaii", "title": "Global Register of Introduced and Invasive Species - Hawaii, United States"},
            {"key": 26226, "alias": "GRIIS-EC", "title": "GRIIS Checklist of Introduced and Invasive Species - Ecuador"},
        ],
    )
    found = {d.alias: d for d in griis_datasets(api.base)}
    assert (found["GRIIS-US-Hawaii"].key, found["GRIIS-US-Hawaii"].country, found["GRIIS-US-Hawaii"].subdivision) == (20359, "US", "Hawaii")
    assert (found["GRIIS-EC"].country, found["GRIIS-EC"].subdivision) == ("EC", None)


def test_griis_excludes_protected_areas_and_others(api):
    griis_listing(
        api,
        [
            {"key": 23243, "alias": "GRIIS-EC-Galapagos-1", "title": "Protected Areas - Global Register of Introduced and Invasive Species - Galapagos"},
            {"key": 26853, "alias": None, "title": "World Register of Introduced Marine Species (WRiMS)"},
            {"key": 1, "alias": "Something", "title": "A GRIIS mention"},
            {"key": 25058, "alias": "GRIIS-ZW", "title": "GRIIS Checklist of Introduced and Invasive Species - Zimbabwe"},
        ],
    )
    assert [d.alias for d in griis_datasets(api.base)] == ["GRIIS-ZW"]


COLUMNS = ["col:ID", "col:parentID", "col:status", "col:rank", "col:scientificName", "col:extinct", "col:environment", "col:genus", "col:family", "col:order", "col:class", "col:phylum", "col:kingdom"]
ROWS = [
    ["F1", "", "accepted", "family", "Fidae", "", "", "", "", "Oales", "Cea", "Phy", "Animalia"],
    ["G1", "F1", "accepted", "genus", "Gus", "", "", "", "Fidae", "Oales", "Cea", "Phy", "Animalia"],
    ["SG", "G1", "accepted", "subgenus", "Gus (Hus)", "", "", "Gus", "Fidae", "Oales", "Cea", "Phy", "Animalia"],
    ["S1", "G1", "accepted", "species", 'Gus "quoted" bus', "false", "marine", "Gus", "Fidae", "Oales", "Cea", "Phy", "Animalia"],
    ["S2", "SG", "accepted", "species", "Gus cus", "true", "", "Gus", "Fidae", "Oales", "Cea", "Phy", "Animalia"],
    ["G2", "", "accepted", "genus", "Dus", "", "", "", "", "", "", "", "Plantae"],
    ["S3", "G2", "accepted", "species", "Dus eus", "false", "", "Dus", "", "", "", "", "Plantae"],
    ["S4", "G2", "provisionally accepted", "species", "Dus fus", "false", "", "Dus", "", "", "", "", "Plantae"],
    ["SY", "S1", "synonym", "species", "Gus old", "", "", "", "", "", "", "", ""],
]
# a second accepted Gus in the same kingdom, as COL has for thousands of genera
HOMONYM = [
    ["G3", "", "accepted", "genus", "Gus", "", "", "", "", "", "", "", "Animalia"],
    ["S5", "G3", "accepted", "species", "Gus hus", "false", "", "Gus", "", "", "", "", "Animalia"],
]


def coldp(tmp_path, columns=COLUMNS, rows=ROWS):
    path = tmp_path / "coldp.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("NameUsage.tsv", "\n".join("\t".join(r) for r in [columns, *rows]) + "\n")
    return path


def test_load_taxa_keeps_accepted_species_only(tmp_path):
    taxa = load_taxa(coldp(tmp_path))
    assert sorted(taxa["species_id"]) == ["S1", "S2", "S3", "S4"]


def test_load_taxa_takes_flat_classification(tmp_path):
    row = load_taxa(coldp(tmp_path)).filter(pl.col("species_id") == "S1").row(0, named=True)
    assert row == {
        "species_id": "S1",
        "genus": "Gus",
        "genus_id": "G1",
        "family": "Fidae",
        "family_id": "F1",
        "order": "Oales",
        "class": "Cea",
        "phylum": "Phy",
        "kingdom": "Animalia",
        "environment": "marine",
        "extinct": False,
    }


def test_load_taxa_reads_empty_as_null(tmp_path):
    row = load_taxa(coldp(tmp_path)).filter(pl.col("species_id") == "S3").row(0, named=True)
    assert row["family"] is None and row["family_id"] is None and row["environment"] is None and row["extinct"] is False


def test_load_taxa_walks_up_through_ranks_between(tmp_path):
    row = load_taxa(coldp(tmp_path)).filter(pl.col("species_id") == "S2").row(0, named=True)
    assert (row["genus_id"], row["family_id"]) == ("G1", "F1")


def test_load_taxa_tells_homonym_genera_apart_by_id(tmp_path):
    taxa = load_taxa(coldp(tmp_path, rows=ROWS + HOMONYM)).filter(pl.col("genus") == "Gus")
    assert dict(zip(taxa["species_id"], taxa["genus_id"])) == {"S1": "G1", "S2": "G1", "S5": "G3"}


def test_load_taxa_species_directly_under_a_family_has_no_genus_id(tmp_path):
    rows = ROWS + [["S6", "F1", "accepted", "species", "Xus yus", "false", "", "", "Fidae", "Oales", "Cea", "Phy", "Animalia"]]
    row = load_taxa(coldp(tmp_path, rows=rows)).filter(pl.col("species_id") == "S6").row(0, named=True)
    assert (row["genus_id"], row["family_id"]) == (None, "F1")


def test_load_taxa_refuses_a_genus_that_is_not_accepted(tmp_path):
    rows = ROWS + [
        ["GS", "", "synonym", "genus", "Old", "", "", "", "", "", "", "", "Animalia"],
        ["S7", "GS", "accepted", "species", "Old bus", "false", "", "Old", "", "", "", "", "Animalia"],
    ]
    with pytest.raises(ValueError, match="S7.*synonym genus GS"):
        load_taxa(coldp(tmp_path, rows=rows))


def test_repeated_names_count_ids_holding_species_within_a_kingdom(tmp_path):
    taxa = load_taxa(coldp(tmp_path, rows=ROWS + HOMONYM))
    assert repeated_names(taxa, "genus") == {"names": 1, "species": 3, "outside_largest": 1}
    assert repeated_names(taxa, "family") == {"names": 0, "species": 0, "outside_largest": 0}


def test_load_taxa_is_not_confused_by_quotes(tmp_path):
    assert load_taxa(coldp(tmp_path)).height == 4


def test_load_taxa_names_columns_on_mismatch(tmp_path):
    path = coldp(tmp_path, [c for c in COLUMNS if c != "col:family"], [r[:8] + r[9:] for r in ROWS])
    with pytest.raises(ValueError, match=r"col:family.*col:genus"):
        load_taxa(path)


def test_load_taxa_without_name_usage_is_an_error(tmp_path):
    path = tmp_path / "empty.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("metadata.yaml", "title: x\n")
    with pytest.raises(ValueError, match="NameUsage.tsv"):
        load_taxa(path)


def test_load_tax_groups_reads_the_flat_classification(tmp_path):
    path = tmp_path / "flat.zip"
    columns = ["col:ID", "col:status", "col:rank", "col:scientificName", "clb:taxGroupFromName", "clb:taxGroup"]
    rows = [
        ["F1", "accepted", "family", "Fidae", "", "arthropods"],
        ["S1", "accepted", "species", "Gus bus", "", "arthropods"],
        ["S2", "provisionally accepted", "species", "Dus eus", "", ""],
        ["SY", "synonym", "species", "Gus old", "", "arthropods"],
    ]
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("metadata.yaml", "title: x\n")
        z.writestr("dataset-1.tsv", "\n".join("\t".join(r) for r in [columns, *rows]) + "\n")
    groups = load_tax_groups(path)
    assert dict(zip(groups["species_id"], groups["tax_group"])) == {"S1": "arthropods", "S2": None}


def test_rank_coverage_counts(tmp_path):
    coverage = rank_coverage(load_taxa(coldp(tmp_path)))
    assert coverage["species"] == 4
    assert coverage["with_genus"] == 4
    assert coverage["with_family"] == 2
    assert coverage["with_environment"] == 1
    assert coverage["extinct"] == 1
    assert coverage["per_kingdom"]["Plantae"] == {"species": 2, "with_genus": 2, "with_family": 0}
    assert coverage["per_kingdom"]["Animalia"] == {"species": 2, "with_genus": 2, "with_family": 2}
