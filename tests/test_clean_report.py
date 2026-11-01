"""The summaries of a `Cleaned` and the report on provisionally accepted twins."""

import polars as pl

from clean_helpers import GRID, counts, refs, row, taxon
from realms.clean import clean
from realms.clean_report import provisional_twin_records, provisional_twins, removal_summary

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
