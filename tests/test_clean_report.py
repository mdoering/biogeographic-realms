"""The summaries of a `Cleaned` and the report on provisionally accepted twins."""

import polars as pl

from clean_helpers import GRID, counts, make_grid, refs, row, taxon
from realms.clean import clean
from realms.clean_report import (
    griis_change,
    griis_feral_candidates,
    griis_totals,
    neighbour_countries,
    provisional_twin_records,
    provisional_twins,
    removal_summary,
)

# ---- the summary

def test_removal_summary_counts_records_and_presences_per_layer_rule_group_and_country():
    r = refs(
        [taxon("ext1", extinct=True), taxon("ext2", extinct=True), taxon("alien")],
        groups=[("ext1", "Birds"), ("ext2", "Insects"), ("alien", "Insects")],
        griis=[("alien", "AA", None)],
    )
    c = counts(
        row("ext1", "10:land", country="AA", n=5), row("ext1", "10:land", country="BB", n=3), row("ext1", "13:land", country="AA", n=2),
        row("ext2", "10:land", country="AA", n=7),
        row("alien", "10:land", country="AA", n=11),
    )
    out = clean(c, r, GRID)
    summary = removal_summary(out, r)
    assert summary.columns == ["layer", "rule", "unplaced", "tax_group", "countrycode", "rows", "presences", "records"]

    def cell(**k):
        return [(d["rows"], d["presences"], d["records"]) for d in summary.filter(*(pl.col(a) == b for a, b in k.items())).to_dicts()]

    assert cell(rule="extinct", tax_group="Birds", countrycode="AA") == [(2, 2, 7)]
    assert cell(rule="extinct", tax_group="Birds", countrycode="BB") == [(1, 1, 3)]  # 10:land is also AA's, so a country's presences don't add up
    assert cell(rule="extinct", tax_group="Insects") == [(1, 1, 7)]
    assert cell(rule="griis") == [(1, 1, 11)]
    assert summary.filter(pl.col("rule") == "griis")["layer"].to_list() == [3]
    assert set(summary["unplaced"]) == {False}
    assert summary["records"].sum() == 5 + 3 + 2 + 7 + 11
    # with no country in the key, the presences are the true ones: 10:land and 13:land of ext1
    by_group = removal_summary(out, r, by=("tax_group",))
    assert by_group.columns == ["layer", "rule", "unplaced", "tax_group", "rows", "presences", "records"]
    assert by_group.filter(pl.col("rule") == "extinct", pl.col("tax_group") == "Birds").row(0, named=True)["presences"] == 2
    by_rule = removal_summary(out, r, by=())
    assert by_rule.filter(pl.col("rule") == "extinct").row(0, named=True)["presences"] == 3  # ext1's two and ext2's one
    assert by_rule.filter(pl.col("rule") == "extinct").row(0, named=True)["records"] == 17


def test_a_presence_that_survives_in_another_country_is_not_a_presence_lost():
    r = refs([taxon("s")], griis=[("s", "AA", None)])
    c = counts(row("s", "10:land", country="AA", n=5), row("s", "10:land", country="BB", n=3))
    out = clean(c, r, GRID)
    assert out.removed["presence_lost"].to_list() == [False]
    assert removal_summary(out, r, by=()).row(0, named=True)["presences"] == 0


def test_removal_summary_of_a_clean_with_no_removals_is_empty_and_has_the_columns():
    r = refs([taxon("s")])
    summary = removal_summary(clean(counts(row("s", "10:land")), r, GRID), r)
    assert summary.height == 0 and summary.columns == ["layer", "rule", "unplaced", "tax_group", "countrycode", "rows", "presences", "records"]


def test_removal_summary_names_a_species_with_no_group():
    r = refs([taxon("s", extinct=True)], groups=[])
    summary = removal_summary(clean(counts(row("s", "10:land")), r, GRID), r)
    assert summary["tax_group"].to_list() == ["unknown"]


# ---- provisionally accepted twins

def test_provisional_twins_are_the_provisional_ids_whose_name_also_has_an_accepted_id():
    index = pl.DataFrame(
        {
            "name": ["Aus bus", "Aus bus", "Aus cus", "Aus dus", "Aus dus", "Aus eus", "Aus eus", "Aus fus", "Aus fus"],
            "kingdom": ["Animalia"] * 7 + ["Animalia", "Plantae"],
            "species_id": ["A1", "A2", "B1", "C1", "C2", "D1", "D2", "E1", "E2"],
            "status": ["accepted", "provisionally accepted", "provisionally accepted", "accepted", "synonym", "accepted", "accepted", "accepted", "provisionally accepted"],
        }
    )
    twins = provisional_twins(index)
    assert twins.columns == ["species_id", "name", "kingdom", "accepted_ids"]
    # A2 has an accepted twin; B1 has none; a synonym is not a twin; E2 is a plant, which an animal does not twin
    assert twins["species_id"].to_list() == ["A2"]
    assert twins["accepted_ids"].to_list() == [1]


def test_provisional_twin_records_count_the_records_of_those_ids():
    twins = pl.DataFrame({"species_id": ["A2", "A3"], "name": ["x", "y"], "kingdom": ["Animalia"] * 2, "accepted_ids": [1, 2]})
    c = counts(row("A2", "10:land", n=5), row("A2", "13:land", n=2), row("A1", "10:land", n=100), row("A3", "10:land", n=0, n_raw=4))
    found = provisional_twin_records(c, twins)
    assert found == {"species": 2, "rows": 3, "n": 7, "n_raw": 11, "ambiguous_species": 1}


# ---- GRIIS

# a star: cell 1 touches 2, 3 and 4, and 5 touches nothing
STAR = make_grid({c: (1000.0, 0.0) for c in "12345"}, [("1", "2"), ("1", "3"), ("1", "4")])


def test_countries_are_neighbours_when_their_rows_are_in_the_same_or_adjacent_cells():
    c = counts(
        row("s", "10:land", country="AA"), row("s", "11:land", country="BB"), row("s", "13:land", country="CC"),
        row("t", "10:land", country="DD"), row("t", "12:sea", country="EE"), row("t", "12:sea", country=None),
    )
    found = neighbour_countries(c, GRID)
    assert found.columns == ["country", "neighbour"]
    pairs = set(found.iter_rows())
    # 10 - 11 - 12 - 13: AA and DD share a cell, BB touches both of them and EE, EE touches CC
    assert pairs == {("AA", "BB"), ("AA", "DD"), ("BB", "AA"), ("BB", "DD"), ("BB", "EE"), ("DD", "AA"), ("DD", "BB"), ("EE", "BB"), ("EE", "CC"), ("CC", "EE")}
    assert ("AA", "CC") not in pairs and ("AA", "AA") not in pairs  # no cell between 10 and 13, and a country is not its own neighbour


def test_the_largest_griis_removals_of_species_clean_in_three_neighbouring_countries_are_listed():
    r = refs(
        [taxon(s) for s in ("mal", "two", "both", "far")],
        griis=[("mal", "AA", None), ("two", "AA", None), ("both", "AA", None), ("both", "BB", None), ("far", "AA", None)],
    )
    c = counts(
        # mal: listed in AA (100 records), clean and unlisted in BB, CC and DD, which touch AA's cell, and in EE, which does not
        row("mal", "1:land", country="AA", n=100), row("mal", "2:land", country="BB", n=5), row("mal", "3:land", country="CC", n=6),
        row("mal", "4:land", country="DD", n=7), row("mal", "5:land", country="EE", n=8),
        # two: only two neighbours
        row("two", "1:land", country="AA", n=500), row("two", "2:land", country="BB", n=5), row("two", "3:land", country="CC", n=6),
        # both: BB is listed too, so only CC and DD count
        row("both", "1:land", country="AA", n=300), row("both", "2:land", country="BB", n=5), row("both", "3:land", country="CC", n=6), row("both", "4:land", country="DD", n=7),
        # far: records only where AA's cell has no neighbour
        row("far", "1:land", country="AA", n=900), row("far", "5:land", country="EE", n=5),
        # and a neighbour's records that are all human text are not clean
        row("two", "4:land", country="DD", n=3, n_human_text=3),
    )
    removed = clean(c, r, STAR, frozenset({1, 3})).removed
    found = griis_feral_candidates(c, removed, STAR)
    assert found.columns == ["species_id", "country", "records", "rows", "neighbours", "neighbour_records", "neighbour_countries"]
    assert found.rows() == [("mal", "AA", 100, 1, 3, 18, ["DD:7", "CC:6", "BB:5"])]
    # lowering the bar to two neighbours lets `two` and `both` in, in order of what GRIIS took, and `n` cuts the list
    several = griis_feral_candidates(c, removed, STAR, min_neighbours=2)
    assert [(a, b) for a, b, *_ in several.rows()] == [("two", "AA"), ("both", "AA"), ("mal", "AA")]
    assert griis_feral_candidates(c, removed, STAR, n=1, min_neighbours=2).height == 1
    # a neighbour needs `min_records` of the species' own clean, unlisted records: with 6, mal keeps CC and DD only, which is too few
    assert griis_feral_candidates(c, removed, STAR, min_records=6).height == 0
    assert griis_feral_candidates(c, removed, STAR, min_neighbours=2, min_records=6)["species_id"].to_list() == ["both", "mal"]  # `two` has CC only


def test_griis_totals_count_rows_presences_and_records_once_across_the_two_rules():
    r = refs([taxon("s1"), taxon("s2")], griis=[("s1", "AA", None), ("s2", "FR", "Reunion")], admin1=[("11", "FR", "Reunion", "centre")])
    c = counts(row("s1", "10:land", n=5), row("s1", "13:land", n=7), row("s1", "13:land", country="BB", n=2), row("s2", "11:land", country="FR", n=4))
    out = clean(c, r, GRID)
    totals = griis_totals(out.removed)
    assert (totals["rows"], totals["records"], totals["species"]) == (3, 16, 2)
    assert totals["presences"] == 2  # s1's 10:land and s2's 11:land: 13:land keeps s1's BB row, so that part is not lost
    assert totals["griis"]["rows"] == 2 and totals["griis_subnational"] == {"rows": 1, "presences": 1, "records": 4, "species": 1}


def test_the_species_whose_griis_removals_changed_most_come_first():
    def removed(*rows):
        return pl.DataFrame(
            [{"layer": 3, "rule": rule, "species_id": s, "part": "10:land", "countrycode": "AA", "records": n, "to_part": None, "presence_lost": True} for s, rule, n in rows]
        )

    before = removed(("a", "griis", 100), ("b", "griis", 50), ("c", "griis_subnational", 10), ("x", "extinct", 999))
    after = removed(("a", "griis", 100), ("b", "griis", 5), ("c", "griis", 30), ("d", "griis", 40), ("x", "extinct", 1))
    found = griis_change(before, after, n=3)
    assert found.columns == ["species_id", "before", "after", "change"]
    assert found.rows() == [("b", 50, 5, -45), ("d", 0, 40, 40), ("c", 10, 30, 20)]  # `a` did not change and `x` is not GRIIS
