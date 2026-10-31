import polars as pl
import pytest

from helpers import presence
from realms.endemism import endemism_score, presence_by_rank, realm_endemism


def test_tolerant_counts_95_percent_inside():
    # 19 of the species' 20 cells are in realm 0
    rows = [(f"in{i}", "s") for i in range(19)] + [("out", "s"), ("out", "t")]
    labels = {f"in{i}": 0 for i in range(19)} | {"out": 1}
    assert realm_endemism(labels, presence(rows), tolerance=0.95)[0] == 1.0
    assert realm_endemism(labels, presence(rows), tolerance=1.0)[0] == 0.0


def test_below_tolerance_is_not_endemic_anywhere():
    # 3 of 4 cells inside is under 95%, and the species stays non-endemic in both realms
    rows = [("a", "s"), ("b", "s"), ("c", "s"), ("d", "s")]
    result = realm_endemism({"a": 0, "b": 0, "c": 0, "d": 1}, presence(rows))
    assert result == {0: 0.0, 1: 0.0}


def test_strict_needs_all_cells():
    rows = [("a", "s"), ("b", "s"), ("c", "t")]
    assert realm_endemism({"a": 0, "b": 0, "c": 1}, presence(rows), tolerance=1.0) == {0: 1.0, 1: 1.0}


def test_realm_endemism_share_hand_example():
    # realm 0 holds s1, s2 (endemic) and w (shared): 2 of 3. Realm 1 holds s3 (endemic) and w: 1 of 2
    rows = [("a", "s1"), ("a", "s2"), ("a", "w"), ("b", "s3"), ("b", "w")]
    assert realm_endemism({"a": 0, "b": 1}, presence(rows)) == {0: pytest.approx(2 / 3), 1: 0.5}


def test_realm_with_no_species_scores_zero():
    assert realm_endemism({"a": 0, "b": 1}, presence([("a", "s")])) == {0: 1.0, 1: 0.0}


def test_score_is_area_weighted():
    # realm 0 scores 1.0 on 3 units of area, realm 1 scores 0.5 on 1 unit: (3 + 0.5) / 4
    rows = [("a", "s1"), ("b", "s2"), ("b", "w"), ("a", "w2"), ("c", "w")]
    labels = {"a": 0, "c": 0, "b": 1}
    # realm 0 holds s1, w2 and w (2 of 3 endemic); realm 1 holds s2 and w (1 of 2)
    areas = {"a": 2.0, "c": 1.0, "b": 1.0}
    assert endemism_score(labels, presence(rows), areas) == pytest.approx((3 * (2 / 3) + 1 * 0.5) / 4)


def test_presence_by_rank_rolls_up_and_skips_null_ranks():
    p = presence([("a", "s1"), ("a", "s2"), ("b", "s3")])
    taxa = pl.DataFrame(
        {"species_id": ["s1", "s2", "s3"], "genus_id": ["G", "G", "H"], "family_id": [None, None, None]},
        schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String},
    )
    ranks = presence_by_rank(p, taxa)
    assert set(ranks) == {"species", "genus"}
    assert sorted(ranks["genus"].iter_rows()) == [("a", "G"), ("b", "H")]
    assert ranks["species"].height == 3


def test_presence_by_rank_without_taxa_is_species_only():
    assert set(presence_by_rank(presence([("a", "s1")]), None)) == {"species"}


def test_cells_outside_the_map_count_against_endemism():
    # s also occurs in a cell the map doesn't cover, so it isn't known to be unique to realm 0
    rows = [("a", "s"), ("outside", "s"), ("b", "t")]
    assert realm_endemism({"a": 0, "b": 1}, presence(rows), tolerance=1.0) == {0: 0.0, 1: 1.0}


def test_taxon_only_outside_the_map_changes_nothing():
    rows = [("a", "s"), ("b", "t"), ("outside", "u")]
    assert realm_endemism({"a": 0, "b": 1}, presence(rows), tolerance=1.0) == {0: 1.0, 1: 1.0}


def test_null_species_is_not_a_taxon():
    rows = [("a", "s1"), ("a", None), ("b", None), ("b", "s2")]
    assert realm_endemism({"a": 0, "b": 1}, presence(rows)) == {0: 1.0, 1: 1.0}


def test_presence_by_rank_drops_null_species():
    assert presence_by_rank(presence([("a", "s1"), ("b", None)]), None)["species"].rows() == [("a", "s1")]


def test_presence_by_rank_keeps_homonyms_apart_by_id():
    p = presence([("a", "s1"), ("b", "s2")])
    taxa = pl.DataFrame(
        {"species_id": ["s1", "s2"], "genus": ["Drosophila", "Drosophila"], "genus_id": ["BXWZ6", "8HFX7"], "family_id": [None, None]},
        schema={"species_id": pl.String, "genus": pl.String, "genus_id": pl.String, "family_id": pl.String},
    )
    assert sorted(presence_by_rank(p, taxa)["genus"].iter_rows()) == [("a", "BXWZ6"), ("b", "8HFX7")]


def test_presence_by_rank_counts_accepted_species_only():
    # a synonym id must not count as a species when comparing cells or scoring endemism
    p = presence([("a", "s1"), ("b", "synonym")])
    taxa = pl.DataFrame({"species_id": ["s1"], "genus_id": ["G"], "family_id": ["F"]})
    ranks = presence_by_rank(p, taxa)
    assert ranks["species"].rows() == [("a", "s1")]
    assert ranks["genus"].rows() == [("a", "G")]
