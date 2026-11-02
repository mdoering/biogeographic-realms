import json
import logging
import zipfile
from pathlib import Path

import geopandas as gpd
import polars as pl
import pytest
from shapely import MultiPolygon, box

from realms.gbif import Match
from realms.grid import Grid
from realms.lists import NO_WILD_RANGE, NO_WILD_RANGE_DROPPED, resolve_no_wild_range
from realms.natives import (
    admin1_cells,
    admin1_places,
    gbif_fallback,
    griis_lists,
    griis_names,
    griis_name_only_pairs,
    griis_pairs,
    load_admin1_cells,
    load_references,
    marine_country,
    match_names,
    match_stats,
    griis_without_cells,
    name_index,
    read_fallback,
    read_tsv,
    tdwg3_parts,
    unknown_tdwg3,
    wcvp_regions,
    write_fallback,
    wrims_pairs,
)

COL_COLUMNS = ["col:ID", "col:parentID", "col:status", "col:rank", "col:scientificName", "col:extinct", "col:environment", "col:genus", "col:family", "col:order", "col:class", "col:phylum", "col:kingdom"]


def usage(id, parent, status, rank, name, kingdom="", extinct="", environment="", genus="", family=""):
    return [id, parent, status, rank, name, extinct, environment, genus, family, "", "", "", kingdom]


def write_zip(path, files):
    with zipfile.ZipFile(path, "w") as z:
        for name, (columns, rows) in files.items():
            z.writestr(name, "\n".join("\t".join(r) for r in [columns, *rows]) + "\n")
    return path


def col_zip(tmp_path, extra=(), name="col.zip"):
    rows = [
        usage("F1", "", "accepted", "family", "Fidae", "Animalia"),
        usage("G1", "F1", "accepted", "genus", "Aus", "Animalia"),
        usage("A1", "G1", "accepted", "species", "Aus bus", "Animalia", genus="Aus", family="Fidae"),
        usage("A2", "G1", "accepted", "species", "Aus cus", "Animalia", genus="Aus", family="Fidae"),
        usage("A2s", "A2", "accepted", "subspecies", "Aus cus dus", "Animalia"),
        usage("Y1", "A1", "synonym", "species", "Aus old"),
        usage("Y2", "A2s", "synonym", "species", "Aus sub"),
        usage("Y3", "A1", "ambiguous synonym", "species", "Aus amb"),
        usage("Y4", "A1", "misapplied", "species", "Aus mis"),
        usage("Y5", "G1", "synonym", "species", "Aus gen"),
        usage("Y7", "A2", "synonym", "species", "Aus old"),
        usage("P1", "", "accepted", "species", "Plantus eus", "Plantae"),
        usage("P2", "", "accepted", "species", "Aus bus", "Plantae"),
        usage("A3", "G1", "accepted", "species", "Aus twin", "Animalia", genus="Aus"),
        usage("A4", "G1", "accepted", "species", "Aus twin", "Animalia", genus="Aus"),
        usage("A5", "G1", "accepted", "species", "Aus prov", "Animalia", genus="Aus"),
        usage("A6", "G1", "provisionally accepted", "species", "Aus prov", "Animalia", genus="Aus"),
        usage("A7", "G1", "accepted", "species", "Aus syn", "Animalia", genus="Aus"),
        usage("Y6", "A1", "synonym", "species", "Aus syn"),
        usage("B1", "", "accepted", "species", "Bactus bus", "Bacillati"),
        usage("P3", "", "accepted", "species", "Plantus twin", "Plantae"),
        usage("P4", "", "accepted", "species", "Plantus twin", "Plantae"),
        *extra,
    ]
    return write_zip(tmp_path / name, {"NameUsage.tsv": (COL_COLUMNS, rows)})


@pytest.fixture
def index(tmp_path):
    return name_index(col_zip(tmp_path))


def ids(index, name, kingdom=None):
    rows = index.filter(pl.col("name") == name)
    if kingdom:
        rows = rows.filter(pl.col("kingdom") == kingdom)
    return sorted(rows["species_id"].to_list())


def names_frame(*rows):
    return pl.DataFrame(rows, schema={"name": pl.String, "kingdom": pl.String}, orient="row")


def one(matched, name):
    return matched.filter(pl.col("name") == name).row(0, named=True)


# ---------------------------------------------------------------- name_index


def test_index_holds_accepted_and_provisionally_accepted_species(index):
    assert index.columns == ["name", "kingdom", "species_id", "status"]
    assert ids(index, "Aus bus", "Animalia") == ["A1"]
    assert index.filter(pl.col("species_id") == "A6")["status"].to_list() == ["provisionally accepted"]
    assert "Aus cus dus" not in index["name"].to_list() and "Aus" not in index["name"].to_list()


def test_a_synonym_resolves_to_its_accepted_species(index):
    row = index.filter(pl.col("name") == "Aus old", pl.col("species_id") == "A1").row(0, named=True)
    assert row == {"name": "Aus old", "kingdom": "Animalia", "species_id": "A1", "status": "synonym"}


def test_a_synonym_of_a_subspecies_rolls_up_to_the_species(index):
    assert ids(index, "Aus sub") == ["A2"]


def test_ambiguous_synonyms_misapplied_names_and_synonyms_of_a_genus_are_left_out(index):
    assert not {"Aus amb", "Aus mis", "Aus gen"} & set(index["name"].to_list())


def test_an_accepted_name_keeps_its_row_over_the_same_species_synonym_of_it(tmp_path):
    index = name_index(col_zip(tmp_path, [usage("Y8", "A1", "synonym", "species", "Aus bus")]))
    assert index.filter(pl.col("species_id") == "A1", pl.col("name") == "Aus bus")["status"].to_list() == ["accepted"]


def test_index_names_the_columns_it_lacks(tmp_path):
    path = write_zip(tmp_path / "bad.zip", {"NameUsage.tsv": (["col:ID", "col:status"], [["A", "accepted"]])})
    with pytest.raises(ValueError, match=r"col:rank.*col:scientificName"):
        name_index(path)


def test_read_tsv_without_the_file_is_an_error(tmp_path):
    path = write_zip(tmp_path / "other.zip", {"Taxon.tsv": (["col:ID"], [["A"]])})
    with pytest.raises(ValueError, match="NameUsage.tsv"):
        read_tsv(path, "NameUsage.tsv", ["col:ID"])


def test_read_tsv_reads_a_directory_too(tmp_path):
    (tmp_path / "NameUsage.tsv").write_text("col:ID\tcol:rank\nA\tspecies\n")
    assert read_tsv(tmp_path, "NameUsage.tsv", ["col:ID"]).to_dicts() == [{"col:ID": "A"}]


# ---------------------------------------------------------------- match_names


def test_an_exact_name_matches(index):
    matched = match_names(names_frame(("Aus cus", "Animalia"), ("Plantus eus", None), ("Nothing here", None)), index)
    assert [(r["name"], r["status"], r["species_id"]) for r in matched.iter_rows(named=True)] == [
        ("Aus cus", "matched", "A2"),
        ("Plantus eus", "matched", "P1"),
        ("Nothing here", "unmatched", None),
    ]


def test_a_synonym_name_matches_its_accepted_species(index):
    assert one(match_names(names_frame(("Aus sub", "Animalia")), index), "Aus sub")["species_id"] == "A2"


def test_an_ambiguous_name_is_counted_and_never_guessed(index):
    matched = match_names(names_frame(("Aus twin", "Animalia"), ("Aus cus", "Animalia"), ("Aus old", "Animalia"), ("Nothing", "Animalia")), index)
    assert one(matched, "Aus twin") == {"name": "Aus twin", "kingdom": "Animalia", "species_id": None, "status": "ambiguous", "n_species": 2, "via": None}
    assert one(matched, "Aus old")["status"] == "ambiguous"  # a synonym of two species
    assert match_stats(matched) == {"names": 4, "matched": 1, "ambiguous": 2, "unmatched": 1, "matched_share": 0.25}


def test_the_kingdom_tells_homonyms_apart(index):
    matched = match_names(names_frame(("Aus bus", "Animalia"), ("Aus bus", "Plantae"), ("Aus bus", None)), index)
    assert [(r["kingdom"], r["status"], r["species_id"]) for r in matched.iter_rows(named=True)] == [
        ("Animalia", "matched", "A1"),
        ("Plantae", "matched", "P2"),
        (None, "ambiguous", None),
    ]


def test_a_clash_of_kingdoms_is_no_match_but_other_kingdoms_are_not_clashes(index):
    matched = match_names(names_frame(("Plantus eus", "Animalia"), ("Bactus bus", "Bacteria"), ("Plantus eus", "Chromista"), ("Aus cus", "Plantae")), index)
    assert [r["status"] for r in matched.iter_rows(named=True)] == ["unmatched", "matched", "matched", "unmatched"]


def test_an_accepted_species_beats_a_provisional_one_and_a_synonym(index):
    matched = match_names(names_frame(("Aus prov", "Animalia"), ("Aus syn", "Animalia")), index)
    assert one(matched, "Aus prov")["species_id"] == "A5"
    assert one(matched, "Aus syn")["species_id"] == "A7"  # A1 holds the same name as a synonym


def test_match_names_keeps_other_columns_and_order(index):
    names = pl.DataFrame({"name": ["Aus cus", "Zzz", "Aus bus"], "alias": ["x", "y", "z"], "kingdom": ["Animalia"] * 3})
    matched = match_names(names, index)
    assert matched["alias"].to_list() == ["x", "y", "z"] and matched["name"].to_list() == ["Aus cus", "Zzz", "Aus bus"]
    assert match_names(matched, index).equals(matched)  # matching a result again changes nothing


def test_names_without_a_kingdom_column_still_match(index):
    matched = match_names(pl.DataFrame({"name": ["Aus cus"]}), index)
    assert matched["species_id"].to_list() == ["A2"]


def test_a_name_needs_a_name_column(index):
    with pytest.raises(ValueError, match="name"):
        match_names(pl.DataFrame({"species": ["x"]}), index)


def test_the_fallback_places_names_that_stayed_unmatched_and_nothing_else(index):
    fallback = pl.DataFrame({"name": ["Nothing here", "Aus twin"], "kingdom": ["Animalia", "Animalia"], "species_id": ["A1", "A3"]})
    matched = match_names(names_frame(("Nothing here", "Animalia"), ("Aus twin", "Animalia"), ("Aus cus", "Animalia")), index, fallback)
    assert [(r["status"], r["species_id"], r["via"]) for r in matched.iter_rows(named=True)] == [("matched", "A1", "gbif"), ("ambiguous", None, None), ("matched", "A2", "name")]


# ---------------------------------------------------------------- WCVP

WCVP_NAMES = ["col:ID", "col:status", "col:rank", "col:scientificName"]
WCVP_DIST = ["col:taxonID", "col:areaID", "col:gazetteer", "col:establishmentMeans"]


def wcvp_zip(tmp_path, usages, distributions):
    return write_zip(tmp_path / "wcvp.zip", {"NameUsage.tsv": (WCVP_NAMES, usages), "Distribution.tsv": (WCVP_DIST, distributions)})


WCVP_USAGES = [
    ["w1", "accepted", "species", "Plantus eus"],
    ["w2", "accepted", "species", "Aus cus"],  # an animal's name: WCVP lists plants, so it is read as a plant and finds nothing
    ["w3", "accepted", "species", "Aus bus"],  # a homonym that only the kingdom separates
    ["w4", "accepted", "genus", "Plantus"],
    ["w5", "synonym", "species", "Plantus old"],
    ["w6", "accepted", "species", "Plantus twin"],  # two plants of that name in COL
]


def test_wcvp_empty_means_native_and_introduced_is_not(tmp_path, index):
    zip_ = wcvp_zip(tmp_path, WCVP_USAGES, [["w1", "BZE", "tdwg", ""], ["w1", "FRA", "tdwg", "introduced"], ["w3", "NZS", "tdwg", ""]])
    regions = wcvp_regions(zip_, index)
    assert regions.columns == ["species_id", "tdwg3", "native"]
    assert regions.rows() == [("P1", "BZE", True), ("P1", "FRA", False), ("P2", "NZS", True)]


def test_uncertain_regions_are_dropped(tmp_path, index):
    zip_ = wcvp_zip(tmp_path, WCVP_USAGES, [["w1", "BZE", "tdwg", ""], ["w1", "FRA", "tdwg", "uncertain"]])
    assert wcvp_regions(zip_, index)["tdwg3"].to_list() == ["BZE"]


def test_wcvp_regions_above_level_3_are_dropped(tmp_path, index):
    zip_ = wcvp_zip(tmp_path, WCVP_USAGES, [["w1", "BZE", "tdwg", ""], ["w1", "84", "tdwg", ""], ["w1", "7", "tdwg", "introduced"]])
    assert wcvp_regions(zip_, index)["tdwg3"].to_list() == ["BZE"]


def test_wcvp_ranks_other_than_species_and_synonyms_have_no_regions(tmp_path, index):
    zip_ = wcvp_zip(tmp_path, WCVP_USAGES, [["w4", "BZE", "tdwg", ""], ["w5", "BZE", "tdwg", ""]])
    assert wcvp_regions(zip_, index).height == 0


def test_a_wcvp_name_that_is_ambiguous_or_unmatched_has_no_regions(tmp_path, index):
    zip_ = wcvp_zip(tmp_path, WCVP_USAGES + [["w7", "accepted", "species", "Plantus unknownus"]], [["w6", "BZE", "tdwg", ""], ["w2", "BZE", "tdwg", ""], ["w7", "BZE", "tdwg", ""]])
    assert wcvp_regions(zip_, index).height == 0


def test_two_wcvp_names_for_one_species_are_native_where_either_says_so(tmp_path):
    usages = [["w1", "accepted", "species", "Plantus eus"], ["w8", "accepted", "species", "Plantus synonymus"]]
    synonym = usage("S1", "P1", "synonym", "species", "Plantus synonymus")
    index = name_index(col_zip(tmp_path, [synonym], "col2.zip"))
    zip_ = wcvp_zip(tmp_path, usages, [["w1", "BZE", "tdwg", "introduced"], ["w8", "BZE", "tdwg", ""], ["w8", "FRA", "tdwg", "introduced"]])
    assert wcvp_regions(zip_, index).rows() == [("P1", "BZE", True), ("P1", "FRA", False)]


def test_wcvp_values_nobody_judged_are_refused(tmp_path, index):
    with pytest.raises(ValueError, match="doubtful"):
        wcvp_regions(wcvp_zip(tmp_path, WCVP_USAGES, [["w1", "BZE", "tdwg", "doubtful"]]), index)
    with pytest.raises(ValueError, match="gazetteer"):
        wcvp_regions(wcvp_zip(tmp_path, WCVP_USAGES, [["w1", "BZE", "text", ""]]), index)


# ---------------------------------------------------------------- GRIIS

PLAIN = ["col:ID", "col:parentID", "col:status", "col:rank", "col:scientificName", "col:authorship"]
EXTENDED = ["col:ID", "col:parentID", "col:status", "col:rank", "col:scientificName", "col:kingdom"]


GRIIS_DIST = ["col:taxonID", "col:sourceID", "col:areaID", "col:area", "col:gazetteer", "col:establishmentMeans"]


def griis_zip(tmp_path, alias, rows, columns=PLAIN, means="introduced"):
    """A GRIIS export. `means` is the `establishmentMeans` of every accepted species entry ("" for none), `None` for a plain export with no
    `Distribution.tsv`, or a dict of entry id to a means (or a list of them, one distribution row each) for the entries it names, the rest `introduced`."""
    files = {"NameUsage.tsv": (columns, rows)}
    if means is not None:
        given = means if isinstance(means, dict) else {}
        default = means if isinstance(means, str) else "introduced"
        dist = []
        for entry, _, status, rank, *_ in rows:
            if rank == "species" and status == "accepted":
                values = given.get(entry, default)
                dist += [[entry, "", "ZW", "Zedland", "iso", v] for v in (values if isinstance(values, list) else [values])]
        files["Distribution.tsv"] = (GRIIS_DIST, dist)
    path = tmp_path / f"{alias}.zip"
    write_zip(path, files)
    return path


def plain_rows(*species):
    rows = [["k", "", "accepted", "kingdom", "Animalia", ""], ["c", "k", "accepted", "class", "Insecta", ""], ["f", "c", "accepted", "family", "Fidae", ""]]
    return rows + [[f"s{i}", "f", "accepted", "species", s, ""] for i, s in enumerate(species)]


def test_ranks_above_species_are_dropped(tmp_path, index):
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus") + [["g", "f", "accepted", "genus", "Aus", ""], ["ss", "s0", "accepted", "subspecies", "Aus bus x", ""]])
    names = griis_names([path])
    assert names["name"].to_list() == ["Aus bus"] and names["rank"].to_list() == ["species"]
    assert sorted(griis_names([path], species_only=False)["rank"].to_list()) == ["class", "family", "genus", "kingdom", "species", "subspecies"]
    pairs = griis_pairs([path], index)
    assert pairs.columns == ["species_id", "country", "subdivision"]
    assert pairs.rows() == [("A1", "ZW", None)]


def test_a_synonym_in_an_alien_list_is_not_an_alien_taxon(tmp_path, index):
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus") + [["sy", "s0", "synonym", "species", "Aus cus", ""]])
    assert griis_pairs([path], index)["species_id"].to_list() == ["A1"]


def test_country_and_subdivision_come_from_the_file_name(tmp_path, index):
    one = griis_zip(tmp_path, "GRIIS-US-Hawaii", plain_rows("Aus bus"))
    two = griis_zip(tmp_path, "GRIIS-NO", plain_rows("Aus bus", "Aus cus"))
    pairs = griis_pairs([one, two], index)
    assert pairs.rows() == [("A1", "NO", None), ("A1", "US", "Hawaii"), ("A2", "NO", None)]


def test_a_place_list_and_its_country_list_stay_apart_when_the_place_has_two_words(tmp_path, index):
    # `fetch_reference.py` names the list of Cayman Brac `GRIIS-KY-CaymanBrac`, as ChecklistBank gives it the alias of the country's `GRIIS-KY`
    country = griis_zip(tmp_path, "GRIIS-KY", plain_rows("Aus bus"))
    place = griis_zip(tmp_path, "GRIIS-KY-CaymanBrac", plain_rows("Aus cus"))
    assert griis_pairs([country, place], index).rows() == [("A1", "KY", None), ("A2", "KY", "CaymanBrac")]


def test_a_plain_list_takes_its_kingdom_from_the_kingdom_row_above(tmp_path, index):
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus"))
    assert griis_names([path])["kingdom"].to_list() == ["Animalia"]
    # a plant by that name is not the animal GRIIS lists
    plant = griis_zip(tmp_path, "GRIIS-PT", [["k", "", "accepted", "kingdom", "Plantae", ""], ["s", "k", "accepted", "species", "Aus bus", ""]])
    assert griis_pairs([plant], index)["species_id"].to_list() == ["P2"]


def test_an_extended_list_uses_its_own_kingdom_and_falls_back_on_the_row_above(tmp_path, index):
    rows = [["k", "", "accepted", "kingdom", "Plantae", "Plantae"], ["s1", "k", "accepted", "species", "Aus bus", "Plantae"], ["s2", "k", "accepted", "species", "Aus cus", ""]]
    names = griis_names([griis_zip(tmp_path, "GRIIS-EC", rows, EXTENDED)])
    assert names.sort("name")["kingdom"].to_list() == ["Plantae", "Plantae"]


def test_a_list_with_no_kingdom_row_matches_on_the_name_alone(tmp_path, index):
    path = griis_zip(tmp_path, "GRIIS-ZW", [["s", "", "accepted", "species", "Aus cus", ""]])
    assert griis_names([path])["kingdom"].to_list() == [None]
    assert griis_pairs([path], index)["species_id"].to_list() == ["A2"]


def test_a_file_that_is_not_a_griis_list_is_refused(tmp_path):
    with pytest.raises(ValueError, match="GRIIS"):
        griis_names([griis_zip(tmp_path, "checklist", plain_rows("Aus bus"))])


def test_no_lists_give_no_pairs(index):
    assert griis_pairs([], index).height == 0
    assert griis_name_only_pairs([], index).height == 0


def listed(lists, alias):
    return lists.filter(pl.col("alias") == alias).row(0, named=True)


def test_griis_native_entries_are_not_alien(tmp_path, index):
    # GRIIS lists some species that are native where it lists them (invasive natives): they are no aliens there
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus", "Aus cus"), means={"s0": "native"})
    assert griis_pairs([path], index).rows() == [("A2", "ZW", None)]
    # and the list says how many it left out
    counted = listed(griis_lists(griis_names([path])), "GRIIS-ZW")
    assert (counted["entries"], counted["introduced"], counted["native"], counted["uncertain"], counted["unstated"]) == (2, 1, 1, 0, 0)


def test_griis_uncertain_entries_are_not_alien(tmp_path, index):
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus", "Aus cus"), means={"s1": "uncertain"})
    assert griis_pairs([path], index).rows() == [("A1", "ZW", None)]
    counted = listed(griis_lists(griis_names([path])), "GRIIS-ZW")
    assert (counted["entries"], counted["introduced"], counted["uncertain"]) == (2, 1, 1)


def test_griis_empty_means_counts_as_introduced(tmp_path, index):
    # a GRIIS list is a list of aliens, so an entry with no value is introduced, and the list says how many had none
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus", "Aus cus"), means={"s0": ""})
    assert griis_pairs([path], index).rows() == [("A1", "ZW", None), ("A2", "ZW", None)]
    counted = listed(griis_lists(griis_names([path])), "GRIIS-ZW")
    assert (counted["introduced"], counted["unstated"], counted["native"]) == (1, 1, 0)


def test_griis_list_without_distributions_is_not_applied(tmp_path, index):
    # a plain export has no distribution table: nothing says what the entries are, so none is applied. It is reported, and its pairs kept apart
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus", "Aus cus"), means=None)
    assert griis_pairs([path], index).height == 0
    assert griis_name_only_pairs([path], index).rows() == [("A1", "ZW", None), ("A2", "ZW", None)]
    counted = listed(griis_lists(griis_names([path])), "GRIIS-ZW")
    assert counted["table"] is False and counted["entries"] == 2
    assert (counted["introduced"], counted["unstated"], counted["native"], counted["uncertain"]) == (0, 0, 0, 0)  # nothing is known of them
    # a list with a table is no name-only list
    assert griis_name_only_pairs([griis_zip(tmp_path, "GRIIS-NO", plain_rows("Aus bus"))], index).height == 0


def test_an_entry_native_in_any_distribution_row_is_native(tmp_path, index):
    # as WCVP: a species that is native somewhere in the list is not an alien of it. An uncertain row beside an introduced one leaves it introduced
    path = griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus", "Aus cus", "Aus syn"), means={"s0": ["introduced", "native"], "s1": ["introduced", "uncertain"], "s2": ["uncertain"]})
    assert griis_pairs([path], index).rows() == [("A2", "ZW", None)]  # A1 is native in one row, A7 only uncertain


def test_two_entries_of_one_species_are_native_if_either_says_so(tmp_path, index):
    # `Aus sub` is a synonym of a subspecies of A2, which the index rolls up to A2: two entries for one species
    rows = plain_rows("Aus cus", "Aus sub")
    path = griis_zip(tmp_path, "GRIIS-ZW", rows, means={"s1": "native"})
    assert griis_pairs([path], index).height == 0
    assert griis_pairs([griis_zip(tmp_path, "GRIIS-NO", rows)], index).rows() == [("A2", "NO", None)]


def test_an_establishment_means_nobody_judged_is_refused(tmp_path):
    with pytest.raises(ValueError, match="doubtful"):
        griis_names([griis_zip(tmp_path, "GRIIS-ZW", plain_rows("Aus bus"), means="doubtful")])


def test_the_counts_of_a_list_say_what_each_entry_is(tmp_path):
    rows = plain_rows("Aus bus", "Aus cus", "Aus twin", "Aus prov", "Aus syn")
    means = {"s0": "native", "s1": "native", "s2": "uncertain", "s3": ""}
    one = griis_zip(tmp_path, "GRIIS-US-Hawaii", rows, means=means)
    two = griis_zip(tmp_path, "GRIIS-NO", rows, means=None)
    lists = griis_lists(griis_names([one, two]))
    assert lists.columns == ["alias", "country", "subdivision", "table", "entries", "introduced", "unstated", "native", "uncertain"]
    assert lists.sort("alias").rows() == [
        ("GRIIS-NO", "NO", None, False, 5, 0, 0, 0, 0),
        ("GRIIS-US-Hawaii", "US", "Hawaii", True, 5, 1, 1, 2, 1),
    ]


# ---------------------------------------------------------------- WRiMS

WRIMS_USAGES = ["col:ID", "col:status", "col:rank", "col:scientificName", "col:kingdom"]
WRIMS_DIST = ["col:taxonID", "col:areaID", "col:area", "col:gazetteer", "col:establishmentMeans"]


def test_wrims_keeps_introduced_areas_of_species(tmp_path, index):
    usages = [["m1", "accepted", "species", "Aus bus", "Animalia"], ["m2", "accepted", "genus", "Aus", "Animalia"], ["m3", "synonym", "species", "Aus cus", "Animalia"]]
    dist = [
        ["m1", "25607", "Spanish part of the Balearic (Iberian Sea)", "mrgid", "introduced"],
        ["m1", "3257", "Firth of Forth", "mrgid", "introduced"],
        ["m1", "25365", "United States part of the North Atlantic Ocean", "mrgid", "native"],
        ["m1", "25611", "Italian part of the Tyrrhenian Sea", "mrgid", ""],
        ["m1", "25623", "Italian part of the Ionian Sea", "mrgid", "uncertain"],
        ["m2", "1", "Italian part of the Ionian Sea", "mrgid", "introduced"],
        ["m3", "2", "Italian part of the Ionian Sea", "mrgid", "introduced"],
    ]
    path = write_zip(tmp_path / "wrims.zip", {"NameUsage.tsv": (WRIMS_USAGES, usages), "Distribution.tsv": (WRIMS_DIST, dist)})
    pairs = wrims_pairs(path, index)
    assert pairs.columns == ["species_id", "area_id", "area", "country"]
    assert pairs.rows() == [("A1", "25607", "Spanish part of the Balearic (Iberian Sea)", "ES"), ("A1", "3257", "Firth of Forth", None)]


@pytest.mark.parametrize(
    "area, country",
    [
        ("Greek part of the Aegean Sea", "GR"),
        ("Tunisian Exclusive Economic Zone", "TN"),
        ("French Exclusive Economic Zone (Crozet Islands)", "FR"),
        ("United States part of the North Atlantic Ocean", "US"),
        ("New Zealand part of the South Pacific Ocean", "NZ"),
        ("Spain - Demarcación Marina Noratlántica", "ES"),
        ("Germany", "DE"),
        ("Mediterranean Sea", None),
        ("Overlapping claim Ukraine: Ukraine / Russia", None),
        ("Dominican part of the Caribbean Sea", None),
        (None, None),
    ],
)
def test_marine_country(area, country):
    assert marine_country(area) == country


# ---------------------------------------------------------------- the grid joins


def write_cells(path, squares):
    """A grid file of square cells given as (id, lon, lat, half-width), as GBIF's files hold them."""
    shapes = [box(lon - h, lat - h, lon + h, lat + h) if not isinstance(h, MultiPolygon) else h for _, lon, lat, h in squares]
    gpd.GeoDataFrame(
        {"lat_c": [lat for _, _, lat, _ in squares], "lon_c": [lon for _, lon, _, _ in squares], "res": [6.0] * len(squares), "geogrid_id": [i for i, *_ in squares]},
        geometry=shapes,
        crs="EPSG:4326",
    ).to_file(path, driver="GPKG")
    return path


def make_grid(squares, land):
    cells = pl.DataFrame(
        {
            "cell": [i for i, *_ in squares],
            "lat": [lat for _, _, lat, _ in squares],
            "lon": [lon for _, lon, _, _ in squares],
            "area_km2": [100.0] * len(squares),
            "land_km2": [land.get(i, 0.0) for i, *_ in squares],
            "sea_km2": [100.0 - land.get(i, 0.0) for i, *_ in squares],
        }
    )
    return Grid(6, cells, {i: set() for i, *_ in squares})


def write_level3(path, regions):
    gpd.GeoDataFrame({"LEVEL3_NAM": [c for c, _ in regions], "LEVEL3_COD": [c for c, _ in regions]}, geometry=[g for _, g in regions], crs="EPSG:4326").to_file(path, driver="GeoJSON")
    return path


def test_each_land_part_gets_the_regions_its_cell_overlaps(tmp_path):
    squares = [("c1", 5.0, 5.0, 1.0), ("c2", 10.0, 5.0, 1.0), ("c3", 7.5, 5.0, 1.0), ("c4", 50.0, 5.0, 1.0), ("c5", 5.0, 50.0, 1.0)]
    grid = make_grid(squares, {"c1": 10.0, "c2": 10.0, "c3": 10.0, "c4": 10.0})  # c5 is all sea
    level3 = write_level3(tmp_path / "level3.geojson", [("AAA", box(0, 0, 7.5, 10)), ("BBB", box(7.5, 0, 15, 10))])
    parts = tdwg3_parts(grid, level3, write_cells(tmp_path / "g.gpkg", squares))
    assert parts.columns == ["part", "tdwg3"]
    assert parts.rows() == [("c1:land", "AAA"), ("c2:land", "BBB"), ("c3:land", "AAA"), ("c3:land", "BBB"), ("c4:land", None)]


def test_a_cell_that_only_touches_a_region_does_not_overlap_it(tmp_path):
    squares = [("c1", 5.0, 5.0, 1.0)]
    level3 = write_level3(tmp_path / "level3.geojson", [("AAA", box(6, 0, 10, 10))])  # shares the edge x = 6
    parts = tdwg3_parts(make_grid(squares, {"c1": 10.0}), level3, write_cells(tmp_path / "g.gpkg", squares))
    assert parts.rows() == [("c1:land", None)]


def test_a_cell_across_the_antimeridian_meets_a_region_on_the_other_side(tmp_path):
    # the file holds the cell in two pieces at 180 and -180, and the region lies east of the line
    cell = MultiPolygon([box(176, 0, 180, 4), box(-180, 0, -176, 4)])
    squares = [("c1", 178.0, 2.0, cell)]
    level3 = write_level3(tmp_path / "level3.geojson", [("FJI", box(-180, 0, -170, 10)), ("RUS", box(150, 0, 177, 10))])
    parts = tdwg3_parts(make_grid(squares, {"c1": 10.0}), level3, write_cells(tmp_path / "g.gpkg", squares))
    assert parts.rows() == [("c1:land", "FJI"), ("c1:land", "RUS")]


def test_polygons_of_another_grid_are_refused(tmp_path):
    squares = [("c1", 5.0, 5.0, 1.0)]
    level3 = write_level3(tmp_path / "level3.geojson", [("AAA", box(0, 0, 10, 10))])
    with pytest.raises(ValueError, match="not the 2 of the grid"):
        tdwg3_parts(make_grid(squares + [("c2", 6.0, 6.0, 1.0)], {"c1": 1.0}), level3, write_cells(tmp_path / "g.gpkg", squares))


def write_admin1(path, places):
    gpd.GeoDataFrame(
        {"name": [n for n, *_ in places], "name_alt": [alt for _, _, alt, _ in places], "iso_a2": [iso for _, iso, _, _ in places]},
        geometry=[g for *_, g in places],
        crs="EPSG:4326",
    ).to_file(path, driver="GeoJSON")
    return path


ADMIN1 = [
    ("Hawaii", "US", None, box(-160, 19, -155, 22)),
    ("Alaska", "US", None, box(-170, 55, -140, 70)),
    ("Texas", "US", None, box(-105, 26, -94, 36)),
    ("Madeira", "PT", "Madeira Islands", box(-17.3, 32.6, -16.6, 32.9)),
    ("Madeira", "BR", None, box(-60, -5, -58, -3)),
    ("Saba", "NL", None, box(-63.3, 17.6, -63.2, 17.7)),
    ("Ha'apai", "TO", "Haapai|Ha`apai", box(-175, -20, -174, -19)),
]


def test_a_subdivision_is_found_by_its_name_in_any_form(tmp_path):
    places = admin1_places(write_admin1(tmp_path / "a.geojson", ADMIN1), [("US", "Hawaii"), ("TO", "Haapai"), ("PT", "Madeira"), ("XX", "Nowhere")])
    assert sorted(zip(places["country"], places["subdivision"], strict=True)) == [("PT", "Madeira"), ("TO", "Haapai"), ("US", "Hawaii")]


def test_a_subdivision_named_without_spaces_is_found_by_its_name_with_them(tmp_path):
    places = admin1_places(write_admin1(tmp_path / "a.geojson", [("Cayman Brac", "KY", None, box(-79.9, 19.68, -79.7, 19.76))]), [("KY", "CaymanBrac")])
    assert list(zip(places["country"], places["subdivision"], strict=True)) == [("KY", "CaymanBrac")]


def test_a_name_in_several_countries_takes_the_lists_own(tmp_path):
    places = admin1_places(write_admin1(tmp_path / "a.geojson", ADMIN1), [("PT", "Madeira")])
    assert places.geometry.iloc[0].bounds[0] < -16  # the Portuguese one, not the Brazilian


def test_a_name_that_only_one_polygon_has_is_taken_whatever_its_country(tmp_path):
    assert len(admin1_places(write_admin1(tmp_path / "a.geojson", ADMIN1), [("BQ", "Saba")])) == 1


def test_united_is_the_conterminous_united_states(tmp_path):
    places = admin1_places(write_admin1(tmp_path / "a.geojson", ADMIN1), [("US", "United")])
    assert len(places) == 1 and places.geometry.iloc[0].bounds[0] == -105  # Texas alone; Alaska and Hawaii are out


def test_subnational_list_covers_cells_overlapping_its_subdivision(tmp_path):
    # Hawaii is the box (-160, 19)-(-155, 22). `inside` has its centre in it, `coast` has its centre at sea but reaches 0.5 degrees into it,
    # `touching` shares the edge x = -155 and nothing more, and `far` is elsewhere
    squares = [("inside", -157.5, 20.5, 1.0), ("coast", -153.0, 20.5, 2.5), ("touching", -154.0, 20.5, 1.0), ("far", -100.0, 20.0, 1.0)]
    grid = make_grid(squares, {c: 1.0 for c, *_ in squares})
    small = [("Hawaii", "US", None, box(-160, 19, -155, 22))]
    cells = admin1_cells(grid, write_admin1(tmp_path / "a.geojson", small), [("US", "Hawaii")], write_cells(tmp_path / "g.gpkg", squares))
    assert cells.columns == ["cell", "country", "subdivision", "via"]
    # `via` says how the cell came by it: by its centre, or by overlapping only
    assert cells.rows() == [("coast", "US", "Hawaii", "overlap"), ("inside", "US", "Hawaii", "centre")]


def test_united_has_the_cells_that_overlap_the_conterminous_states(tmp_path):
    # Texas stands for the conterminous states in this fixture; `edge` overlaps it with its centre outside
    squares = [("in", -100.0, 30.0, 1.0), ("edge", -93.0, 30.0, 2.0), ("out", -80.0, 30.0, 1.0)]
    grid = make_grid(squares, {"in": 1.0, "edge": 1.0, "out": 1.0})
    cells = admin1_cells(grid, write_admin1(tmp_path / "a.geojson", ADMIN1), [("US", "United")], write_cells(tmp_path / "g.gpkg", squares))
    assert cells.rows() == [("edge", "US", "United", "overlap"), ("in", "US", "United", "centre")]


def test_a_cell_that_overlaps_two_subdivisions_belongs_to_both(tmp_path):
    # Texas and a second state share the line x = -105: a cell across it overlaps both, one on either side only the one beside it
    places = [("Texas", "US", None, box(-105, 26, -94, 36)), ("Utah", "US", None, box(-115, 26, -105, 36))]
    squares = [("across", -105.0, 30.0, 1.0), ("east", -100.0, 30.0, 1.0), ("west", -110.0, 30.0, 1.0)]
    grid = make_grid(squares, {c: 1.0 for c, *_ in squares})
    cells = admin1_cells(grid, write_admin1(tmp_path / "a.geojson", places), [("US", "Texas"), ("US", "Utah")], write_cells(tmp_path / "g.gpkg", squares))
    assert cells.rows() == [
        ("across", "US", "Texas", "overlap"),
        ("east", "US", "Texas", "centre"),
        ("across", "US", "Utah", "overlap"),
        ("west", "US", "Utah", "centre"),
    ]


def test_an_island_no_centre_falls_on_gets_the_cells_that_overlap_it(tmp_path):
    # a speck of island at (-161.5, 22.5) inside a cell centred at (-157, 20) that reaches 6 degrees each way
    squares = [("sea", -157.0, 20.0, 6.0), ("far", -100.0, 20.0, 1.0)]
    grid = make_grid(squares, {"sea": 5.0})
    small = [("Hawaii", "US", None, box(-162, 22, -161, 23))]
    cells = admin1_cells(grid, write_admin1(tmp_path / "a.geojson", small), [("US", "Hawaii")], write_cells(tmp_path / "g.gpkg", squares))
    assert cells.rows() == [("sea", "US", "Hawaii", "overlap")]


def test_a_subdivision_with_no_polygon_has_no_cells(tmp_path):
    squares = [("c", -100.0, 30.0, 1.0)]
    cells = admin1_cells(make_grid(squares, {}), write_admin1(tmp_path / "a.geojson", ADMIN1), [("IC", "Tenerife")], write_cells(tmp_path / "g.gpkg", squares))
    assert cells.height == 0 and cells.columns == ["cell", "country", "subdivision", "via"]


# ---------------------------------------------------------------- NO_WILD_RANGE


def test_the_list_and_what_was_left_out_do_not_overlap():
    assert not set(NO_WILD_RANGE) & set(NO_WILD_RANGE_DROPPED)
    assert "Capra hircus" in NO_WILD_RANGE_DROPPED  # COL files the Chiltan markhor under it
    assert all(kingdom in ("Animalia", "Plantae") and reason for kingdom, reason in NO_WILD_RANGE.values())
    assert all(NO_WILD_RANGE_DROPPED.values())


def test_entries_are_all_the_ids_of_an_accepted_name_of_their_kingdom(tmp_path, caplog):
    extra = [
        usage("H2", "", "provisionally accepted", "species", "Homo sapiens", "Animalia"),
        usage("W1", "", "accepted", "species", "Canis lupus", "Animalia"),
        usage("W2", "W1", "synonym", "species", "Gallus domesticus"),
        usage("H3", "", "accepted", "species", "Homo sapiens", "Plantae"),
    ]
    index = name_index(col_zip(tmp_path, [usage("H1", "", "accepted", "species", "Homo sapiens", "Animalia"), *extra]))
    entries = {"Homo sapiens": ("Animalia", "people"), "Gallus domesticus": ("Animalia", "chickens"), "Canis lupus": ("Animalia", "wolves"), "Unknown unknown": ("Plantae", "x")}
    with caplog.at_level(logging.WARNING):
        found, left_out = resolve_no_wild_range(index, entries)
    assert found == frozenset({"H1", "H2", "W1"})  # not the plant of that name
    assert sorted(left_out) == ["Gallus domesticus", "Unknown unknown"]
    assert left_out["Gallus domesticus"] == "a synonym of ['W1']"
    assert "Gallus domesticus dropped" in caplog.text


# ---------------------------------------------------------------- the GBIF fallback


def stub(answers):
    asked = []

    def match(name):
        asked.append(name)
        return answers.get(name, Match(None, None, "NONE"))

    return match, asked


def test_the_fallback_takes_exact_species_that_are_in_this_release(index):
    match, asked = stub(
        {
            "Aus new": Match("A2", "SPECIES", "EXACT"),
            "Aus var": Match("A1", "SPECIES", "VARIANT"),
            "Aus syn2": Match("Y9", "SUBSPECIES", "EXACT", "A2", "SPECIES"),  # a synonym of a species: its accepted species
            "Aus fuz": Match("A2", "SPECIES", "FUZZY"),
            "Aus gen": Match("G1", "GENUS", "HIGHERRANK"),
            "Aus sub2": Match("S9", "SUBSPECIES", "EXACT"),
            "Aus gone": Match("ZZZ9", "SPECIES", "EXACT"),  # not a species of this release
            "Aus syn3": Match("Y1", "SPECIES", "EXACT", "A2s", "SUBSPECIES"),  # synonym of a subspecies
            "Plantus bad": Match("A1", "SPECIES", "EXACT"),  # an animal for a plant: a clash
        }
    )
    names = names_frame(*[(n, "Plantae" if n == "Plantus bad" else "Animalia") for n in ["Aus new", "Aus var", "Aus syn2", "Aus fuz", "Aus gen", "Aus sub2", "Aus gone", "Aus syn3", "Plantus bad", "Aus none"]])
    found = gbif_fallback(names, index, match=match, sleep=lambda s: None)
    assert sorted(found["name"].to_list()) == ["Aus new", "Aus syn2", "Aus var"]
    assert dict(zip(found["name"], found["species_id"], strict=True)) == {"Aus new": "A2", "Aus var": "A1", "Aus syn2": "A2"}
    assert len(asked) == 10


def test_the_fallback_asks_once_per_name_and_at_most_five_times_a_second(index):
    match, asked = stub({})
    now = [0.0]
    waits = []

    def sleep(seconds):
        waits.append(seconds)
        now[0] += seconds

    names = names_frame(("Aus a", "Animalia"), ("Aus a", "Animalia"), ("Aus b", "Animalia"), ("Aus c", "Animalia"))
    gbif_fallback(names, index, per_second=5.0, match=match, sleep=sleep, clock=lambda: now[0])
    assert asked == ["Aus a", "Aus b", "Aus c"]
    assert waits == pytest.approx([0.2, 0.2])  # the stub takes no time, so each request waits out the 0.2 s


# ---------------------------------------------------------------- load_references


@pytest.fixture
def reference_dir(tmp_path):
    root = tmp_path / "reference"
    root.mkdir()
    squares = [("c1", 5.0, 5.0, 1.0), ("c2", 10.0, 5.0, 1.0)]
    col_rows = [
        usage("F1", "", "accepted", "family", "Fidae", "Animalia"),
        usage("G1", "F1", "accepted", "genus", "Aus", "Animalia"),
        usage("A1", "G1", "accepted", "species", "Aus bus", "Animalia", genus="Aus", family="Fidae"),
        usage("H1", "G1", "accepted", "species", "Homo sapiens", "Animalia", genus="Aus", family="Fidae"),
        usage("P1", "", "accepted", "species", "Plantus eus", "Plantae"),
    ]
    col = write_zip(root / "col.zip", {"NameUsage.tsv": (COL_COLUMNS, col_rows)})
    groups = write_zip(root / "groups.zip", {"dataset-1.tsv": (["col:ID", "col:status", "col:rank", "col:scientificName", "clb:taxGroup"], [["A1", "accepted", "species", "Aus bus", "arthropods"]])})
    wcvp = write_zip(root / "wcvp.zip", {"NameUsage.tsv": (WCVP_NAMES, [["w1", "accepted", "species", "Plantus eus"]]), "Distribution.tsv": (WCVP_DIST, [["w1", "BZE", "tdwg", ""]])})
    griis = griis_zip(root, "GRIIS-ZW", plain_rows("Aus bus"))
    sub = griis_zip(root, "GRIIS-US-Hawaii", plain_rows("Aus bus"))
    island = griis_zip(root, "GRIIS-IC-Tenerife", plain_rows("Aus bus"))  # Natural Earth has no polygon for it
    manifest = {
        "col": {"path": str(col), "release": 315557},
        "col_groups": {"path": str(groups)},
        "wcvp": {"path": str(wcvp)},
        "wrims": {"dataset": 26853, "status": "missing"},
        "griis/GRIIS-ZW": {"path": str(griis)},
        "griis/GRIIS-US-Hawaii": {"path": str(sub)},
        "griis/GRIIS-IC-Tenerife": {"path": str(island)},
        "griis/GRIIS-GQ": {"dataset": 1, "status": "missing"},
        "natural_earth_admin1": {"path": str(write_admin1(root / "admin1.geojson", [("Hawaii", "US", None, box(4.5, 4.5, 5.5, 5.5))]))},
        "wgsrpd_level3": {"path": str(write_level3(root / "level3.geojson", [("BZE", box(0, 0, 7, 10))]))},
        "isea3h_6": {"path": str(write_cells(root / "grid.gpkg", squares))},
    }
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root, make_grid(squares, {"c1": 10.0, "c2": 10.0})


def test_load_references_joins_everything_and_reports_what_is_missing(reference_dir, caplog):
    root, grid = reference_dir
    with caplog.at_level(logging.WARNING):
        refs = load_references(grid, root)
    assert refs.taxa["species_id"].sort().to_list() == ["A1", "H1", "P1"]
    assert refs.groups.to_dicts() == [{"species_id": "A1", "tax_group": "arthropods"}]
    assert refs.wcvp.rows() == [("P1", "BZE", True)]
    assert refs.griis.rows() == [("A1", "IC", "Tenerife"), ("A1", "US", "Hawaii"), ("A1", "ZW", None)]
    assert refs.griis_no_cells.rows() == [("A1", "IC", "Tenerife")]  # Hawaii has cells, ZW is no sub-national list
    assert refs.wrims.height == 0 and refs.wrims.columns == ["species_id", "area_id", "area", "country"]
    assert refs.wrims_no_country.height == 0 and refs.wrims_no_country.columns == refs.wrims.columns
    assert refs.tdwg3_parts.rows() == [("c1:land", "BZE"), ("c2:land", None)]
    assert refs.admin1_cells.rows() == [("c1", "US", "Hawaii", "centre")]
    assert refs.no_wild_range == frozenset({"H1"})
    assert refs.missing == ("griis/GRIIS-GQ", "wrims")
    assert "griis/GRIIS-GQ has no export" in caplog.text
    assert "1 GRIIS pairs of 1 sub-national lists have no admin-1 cells" in caplog.text


def test_the_admin1_cells_load_without_the_rest_of_the_references_and_are_the_ones_load_references_has(reference_dir):
    root, grid = reference_dir
    cells = load_admin1_cells(grid, root)
    assert cells.rows() == [("c1", "US", "Hawaii", "centre")]
    assert cells.equals(load_references(grid, root).admin1_cells)
    # a list the manifest has no export for names no place, and the polygons and the grid's file are needed
    manifest = json.loads((root / "manifest.json").read_text())
    manifest["griis/GRIIS-US-Hawaii"] = {"dataset": 2, "status": "missing"}
    (root / "manifest.json").write_text(json.dumps(manifest))
    assert load_admin1_cells(grid, root).height == 0
    for item in ("natural_earth_admin1", "isea3h_6"):
        broken = {k: v for k, v in manifest.items() if k != item}
        (root / "manifest.json").write_text(json.dumps(broken))
        with pytest.raises(ValueError, match=item):
            load_admin1_cells(grid, root)


def test_load_references_applies_griis_by_its_own_means_and_reports_the_rest(reference_dir, caplog):
    root, grid = reference_dir
    manifest = json.loads((root / "manifest.json").read_text())
    # ZW names Aus bus as native and Homo sapiens as introduced, Hawaii has Aus bus as introduced, and the Tenerife list is a plain export with no distributions
    manifest["griis/GRIIS-ZW"] = {"path": str(griis_zip(root, "GRIIS-ZW", plain_rows("Aus bus", "Homo sapiens"), means={"s0": "native"}))}
    manifest["griis/GRIIS-IC-Tenerife"] = {"path": str(griis_zip(root, "GRIIS-IC-Tenerife", plain_rows("Aus bus"), means=None))}
    (root / "manifest.json").write_text(json.dumps(manifest))
    with caplog.at_level(logging.WARNING):
        refs = load_references(grid, root)
    assert refs.griis.rows() == [("A1", "US", "Hawaii"), ("H1", "ZW", None)]  # Aus bus is native in ZW, and Tenerife says nothing
    assert refs.griis_name_only.rows() == [("A1", "IC", "Tenerife")]
    assert refs.griis_no_cells.height == 0  # Hawaii has cells
    assert refs.griis_name_only_no_cells.rows() == [("A1", "IC", "Tenerife")]  # Natural Earth has no polygon for it
    assert refs.griis_lists.columns == ["alias", "country", "subdivision", "table", "entries", "introduced", "unstated", "native", "uncertain"]
    assert refs.griis_lists.sort("alias").rows() == [
        ("GRIIS-IC-Tenerife", "IC", "Tenerife", False, 1, 0, 0, 0, 0),
        ("GRIIS-US-Hawaii", "US", "Hawaii", True, 1, 1, 0, 0, 0),
        ("GRIIS-ZW", "ZW", None, True, 2, 1, 0, 1, 0),
    ]
    assert "1 GRIIS pairs of 1 sub-national lists with no distributions have no admin-1 cells" in caplog.text


def test_load_references_makes_its_griis_pairs_with_the_functions_the_tests_check(reference_dir):
    # one path: what the references hold is what `griis_pairs` and `griis_name_only_pairs` give for the same lists and index
    root, grid = reference_dir
    manifest = json.loads((root / "manifest.json").read_text())
    manifest["griis/GRIIS-ZW"] = {"path": str(griis_zip(root, "GRIIS-ZW", plain_rows("Aus bus", "Homo sapiens"), means={"s0": "native"}))}
    manifest["griis/GRIIS-IC-Tenerife"] = {"path": str(griis_zip(root, "GRIIS-IC-Tenerife", plain_rows("Aus bus", "Homo sapiens"), means=None))}
    (root / "manifest.json").write_text(json.dumps(manifest))
    refs = load_references(grid, root)
    paths = [Path(v["path"]) for k, v in manifest.items() if k.startswith("griis/") and v.get("path")]
    index = name_index(Path(manifest["col"]["path"]))
    assert refs.griis.equals(griis_pairs(paths, index)) and refs.griis.height == 2
    assert refs.griis_name_only.equals(griis_name_only_pairs(paths, index)) and refs.griis_name_only.height == 2


def test_load_references_uses_wrims_when_there_is_an_export(reference_dir, caplog):
    root, grid = reference_dir
    manifest = json.loads((root / "manifest.json").read_text())
    usages = [["m1", "accepted", "species", "Aus bus", "Animalia"]]
    dist = [["m1", "25607", "Spanish part of the Balearic (Iberian Sea)", "mrgid", "introduced"], ["m1", "3257", "Firth of Forth", "mrgid", "introduced"]]
    manifest["wrims"] = {"path": str(write_zip(root / "wrims.zip", {"NameUsage.tsv": (WRIMS_USAGES, usages), "Distribution.tsv": (WRIMS_DIST, dist)}))}
    (root / "manifest.json").write_text(json.dumps(manifest))
    with caplog.at_level(logging.WARNING):
        refs = load_references(grid, root)
    assert refs.wrims.rows() == [("A1", "25607", "Spanish part of the Balearic (Iberian Sea)", "ES"), ("A1", "3257", "Firth of Forth", None)]
    assert refs.wrims_no_country.rows() == [("A1", "3257", "Firth of Forth", None)]
    assert "1 of 2 WRiMS pairs are in areas that name no country" in caplog.text
    assert refs.missing == ("griis/GRIIS-GQ",)


def with_plantus_eus_authored(root):
    """The reference directory with a WCVP whose one name carries an author, which only the GBIF fallback can place."""
    manifest = json.loads((root / "manifest.json").read_text())
    names = [["w1", "accepted", "species", "Plantus eus (L.)"]]
    manifest["wcvp"] = {"path": str(write_zip(root / "wcvp2.zip", {"NameUsage.tsv": (WCVP_NAMES, names), "Distribution.tsv": (WCVP_DIST, [["w1", "BZE", "tdwg", ""]])}))}
    (root / "manifest.json").write_text(json.dumps(manifest))


def fallback_of(*species):
    return pl.DataFrame({"name": ["Plantus eus (L.)"] * len(species), "kingdom": ["Plantae"] * len(species), "species_id": list(species)})


def test_load_references_applies_the_gbif_fallback_when_one_was_saved(reference_dir):
    root, grid = reference_dir
    with_plantus_eus_authored(root)
    assert load_references(grid, root).wcvp.height == 0
    write_fallback(fallback_of("P1"), root / "gbif_fallback.parquet", 315557)
    assert load_references(grid, root).wcvp.rows() == [("P1", "BZE", True)]


def test_a_fallback_made_for_another_release_is_refused(reference_dir):
    root, grid = reference_dir
    with_plantus_eus_authored(root)
    write_fallback(fallback_of("P1"), root / "gbif_fallback.parquet", 300000)
    with pytest.raises(ValueError, match="release 300000, not 315557"):
        load_references(grid, root)


def test_a_fallback_with_no_release_is_refused(reference_dir):
    root, grid = reference_dir
    fallback_of("P1").write_parquet(root / "gbif_fallback.parquet")
    with pytest.raises(ValueError, match="release None, not 315557"):
        load_references(grid, root)


def test_a_fallback_cannot_be_checked_when_the_manifest_has_no_release(reference_dir):
    root, grid = reference_dir
    manifest = json.loads((root / "manifest.json").read_text())
    del manifest["col"]["release"]
    (root / "manifest.json").write_text(json.dumps(manifest))
    write_fallback(fallback_of("P1"), root / "gbif_fallback.parquet", 315557)
    with pytest.raises(ValueError, match="not None"):
        load_references(grid, root)


def test_a_fallback_species_the_index_lacks_is_dropped_and_counted(reference_dir, tmp_path, caplog):
    root, grid = reference_dir
    with_plantus_eus_authored(root)
    write_fallback(fallback_of("P1", "GONE1"), root / "gbif_fallback.parquet", 315557)
    with caplog.at_level(logging.WARNING):
        refs = load_references(grid, root)
    assert refs.wcvp.rows() == [("P1", "BZE", True)]
    assert "1 of the 2 names GBIF placed" in caplog.text and "are dropped" in caplog.text
    # and the same on the function itself
    index = name_index(col_zip(tmp_path))
    kept = read_fallback(root / "gbif_fallback.parquet", "315557", index)  # the release may be text or a number
    assert kept["species_id"].to_list() == ["P1"]


def test_wcvp_codes_the_polygon_file_lacks_are_warned_about_with_examples(reference_dir, caplog):
    root, grid = reference_dir
    manifest = json.loads((root / "manifest.json").read_text())
    distributions = [["w1", "BZE", "tdwg", ""], ["w1", "ZZZ", "tdwg", ""], ["w1", "YYY", "tdwg", "introduced"]]
    manifest["wcvp"] = {"path": str(write_zip(root / "wcvp3.zip", {"NameUsage.tsv": (WCVP_NAMES, [["w1", "accepted", "species", "Plantus eus"]]), "Distribution.tsv": (WCVP_DIST, distributions)}))}
    (root / "manifest.json").write_text(json.dumps(manifest))
    with caplog.at_level(logging.WARNING):
        refs = load_references(grid, root)
    assert refs.wcvp.height == 3  # the rows stay: only the warning is new
    assert "2 WCVP level 3 codes are not in LEVEL3_COD" in caplog.text and "(2 rows)" in caplog.text and "YYY (1 rows), ZZZ (1 rows)" in caplog.text


def test_unknown_tdwg3_counts_rows_and_species(tmp_path):
    level3 = write_level3(tmp_path / "level3.geojson", [("AAA", box(0, 0, 1, 1))])
    wcvp = pl.DataFrame({"species_id": ["s1", "s1", "s2", "s2", "s3"], "tdwg3": ["AAA", "ZZZ", "ZZZ", "YYY", "AAA"], "native": [True] * 5})
    found = unknown_tdwg3(wcvp, level3)
    assert found.rows() == [("ZZZ", 2, 2), ("YYY", 1, 1)]
    assert unknown_tdwg3(wcvp.filter(pl.col("tdwg3") == "AAA"), level3).height == 0


def test_griis_pairs_of_subdivisions_with_no_cells_are_the_unplaced_ones():
    griis = pl.DataFrame({"species_id": ["a", "b", "c", "d"], "country": ["US", "US", "ZW", "IC"], "subdivision": ["Hawaii", "Alaska", None, "Tenerife"]})
    cells = pl.DataFrame({"cell": ["c1"], "country": ["US"], "subdivision": ["Hawaii"], "via": ["centre"]})
    assert griis_without_cells(griis, cells).rows() == [("b", "US", "Alaska"), ("d", "IC", "Tenerife")]
    assert griis_without_cells(griis.head(0), cells).height == 0


@pytest.mark.parametrize("item", ["col", "col_groups", "wcvp", "natural_earth_admin1", "wgsrpd_level3", "isea3h_6"])
def test_load_references_needs_the_essential_items(reference_dir, item):
    root, grid = reference_dir
    manifest = json.loads((root / "manifest.json").read_text())
    manifest[item] = {"dataset": 1, "status": "missing"}
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match=item):
        load_references(grid, root)
