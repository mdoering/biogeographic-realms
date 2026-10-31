import numpy as np
import polars as pl
import pytest

import reference_impl as ref
from helpers import matrix, presence
from realms.endemism import area_weighted, endemism_score, realm_endemism, realm_share
from realms.matrix import RANKS, cell_taxa
from realms.synth import planted_network


def links(ct, rank):
    """The `(cell, taxon id)` pairs of one rank of the matrix."""
    coo = ct.x.tocoo()
    return sorted((ct.cells[r], ct.taxa[c].removeprefix(f"{rank}:")) for r, c in zip(coo.row, coo.col) if ct.taxa[c].startswith(f"{rank}:"))


def test_tolerant_counts_95_percent_inside():
    # 19 of the species' 20 cells are in realm 0
    rows = [(f"in{i}", "s") for i in range(19)] + [("out", "s"), ("out", "t")]
    labels = {f"in{i}": 0 for i in range(19)} | {"out": 1}
    assert realm_endemism(labels, matrix(rows), tolerance=0.95)[0] == 1.0
    assert realm_endemism(labels, matrix(rows), tolerance=1.0)[0] == 0.0


def test_below_tolerance_is_not_endemic_anywhere():
    # 3 of 4 cells inside is under 95%, and the species stays non-endemic in both realms
    rows = [("a", "s"), ("b", "s"), ("c", "s"), ("d", "s")]
    result = realm_endemism({"a": 0, "b": 0, "c": 0, "d": 1}, matrix(rows))
    assert result == {0: 0.0, 1: 0.0}


def test_strict_needs_all_cells():
    rows = [("a", "s"), ("b", "s"), ("c", "t")]
    assert realm_endemism({"a": 0, "b": 0, "c": 1}, matrix(rows), tolerance=1.0) == {0: 1.0, 1: 1.0}


def test_realm_endemism_share_hand_example():
    # realm 0 holds s1, s2 (endemic) and w (shared): 2 of 3. Realm 1 holds s3 (endemic) and w: 1 of 2
    rows = [("a", "s1"), ("a", "s2"), ("a", "w"), ("b", "s3"), ("b", "w")]
    assert realm_endemism({"a": 0, "b": 1}, matrix(rows)) == {0: pytest.approx(2 / 3), 1: 0.5}


def test_realm_with_no_species_scores_zero():
    assert realm_endemism({"a": 0, "b": 1}, matrix([("a", "s")])) == {0: 1.0, 1: 0.0}


def test_score_is_area_weighted():
    # realm 0 scores 1.0 on 3 units of area, realm 1 scores 0.5 on 1 unit: (3 + 0.5) / 4
    rows = [("a", "s1"), ("b", "s2"), ("b", "w"), ("a", "w2"), ("c", "w")]
    labels = {"a": 0, "c": 0, "b": 1}
    # realm 0 holds s1, w2 and w (2 of 3 endemic); realm 1 holds s2 and w (1 of 2)
    areas = {"a": 2.0, "c": 1.0, "b": 1.0}
    assert endemism_score(labels, matrix(rows), areas) == pytest.approx((3 * (2 / 3) + 1 * 0.5) / 4)


def test_cell_taxa_rolls_up_and_skips_null_ranks():
    p = presence([("a", "s1"), ("a", "s2"), ("b", "s3")])
    taxa = pl.DataFrame(
        {"species_id": ["s1", "s2", "s3"], "genus_id": ["G", "G", "H"], "family_id": [None, None, None]},
        schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String},
    )
    ct = cell_taxa(p, taxa)
    assert {RANKS[r] for r in ct.rank} == {"species", "genus"}
    assert links(ct, "genus") == [("a", "G"), ("b", "H")]
    assert len(links(ct, "species")) == 3


def test_cell_taxa_without_taxa_is_species_only():
    assert {RANKS[r] for r in cell_taxa(presence([("a", "s1")]), None).rank} == {"species"}


def test_cells_outside_the_map_count_against_endemism():
    # s also occurs in a cell the map doesn't cover, so it isn't known to be unique to realm 0
    rows = [("a", "s"), ("outside", "s"), ("b", "t")]
    assert realm_endemism({"a": 0, "b": 1}, matrix(rows), tolerance=1.0) == {0: 0.0, 1: 1.0}


def test_taxon_only_outside_the_map_changes_nothing():
    rows = [("a", "s"), ("b", "t"), ("outside", "u")]
    assert realm_endemism({"a": 0, "b": 1}, matrix(rows), tolerance=1.0) == {0: 1.0, 1: 1.0}


def test_null_species_is_not_a_taxon():
    rows = [("a", "s1"), ("a", None), ("b", None), ("b", "s2")]
    assert realm_endemism({"a": 0, "b": 1}, matrix(rows)) == {0: 1.0, 1: 1.0}


def test_cell_taxa_drops_null_species():
    assert links(cell_taxa(presence([("a", "s1"), ("b", None)]), None), "species") == [("a", "s1")]


def test_cell_taxa_keeps_homonyms_apart_by_id():
    p = presence([("a", "s1"), ("b", "s2")])
    taxa = pl.DataFrame(
        {"species_id": ["s1", "s2"], "genus": ["Drosophila", "Drosophila"], "genus_id": ["BXWZ6", "8HFX7"], "family_id": [None, None]},
        schema={"species_id": pl.String, "genus": pl.String, "genus_id": pl.String, "family_id": pl.String},
    )
    assert links(cell_taxa(p, taxa), "genus") == [("a", "BXWZ6"), ("b", "8HFX7")]


def test_cell_taxa_counts_accepted_species_only():
    # a synonym id must not count as a species when comparing cells or scoring endemism
    p = presence([("a", "s1"), ("b", "synonym")])
    taxa = pl.DataFrame({"species_id": ["s1"], "genus_id": ["G"], "family_id": ["F"]})
    ct = cell_taxa(p, taxa)
    assert links(ct, "species") == [("a", "s1")]
    assert links(ct, "genus") == [("a", "G")]


def noisy_labels(truth, seed, outside=0.2, wrong=0.15):
    """The true realms with some cells given a random realm and some left out of the map."""
    rng = np.random.default_rng(seed)
    n_realms = len(set(truth.values()))
    return {c: (int(rng.integers(n_realms)) if rng.random() < wrong else r) for c, r in truth.items() if rng.random() >= outside}


def reference_shares(labels, presence, taxa, tolerance):
    """The shares from the string-keyed table and `incidence` as they were, in the realm order of the labels."""
    cells = sorted(labels)
    realms = sorted(set(labels.values()))
    cell_realm = np.array([realms.index(labels[c]) for c in cells])
    table = ref.incidence(ref.rank_presence(presence, taxa, RANKS), cells)
    return cells, cell_realm, dict(zip(realms, map(float, realm_share(table, cell_realm, len(realms), tolerance))))


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("tolerance", [0.95, 1.0])
def test_realm_endemism_matches_reference(seed, tolerance):
    presence, taxa, truth = planted_network(n_realms=4, cells_per_realm=15, endemic_species=30, widespread_species=8, seed=seed)
    labels = noisy_labels(truth, seed)
    _, _, old = reference_shares(labels, presence, taxa, tolerance)
    assert len(labels) < len(truth) and len(set(old.values())) > 1  # cells outside the map count against endemism, and shares differ
    assert realm_endemism(labels, cell_taxa(presence, taxa), tolerance) == old


@pytest.mark.parametrize("seed", range(5))
def test_endemism_score_matches_reference(seed):
    presence, taxa, truth = planted_network(n_realms=4, cells_per_realm=15, endemic_species=30, widespread_species=8, seed=seed)
    labels = noisy_labels(truth, seed)
    areas = {c: float(a) for c, a in zip(labels, np.random.default_rng(seed).uniform(1, 9, size=len(labels)))}
    cells, cell_realm, shares = reference_shares(labels, presence, taxa, 0.95)
    old = area_weighted(np.array(list(shares.values())), cell_realm, np.array([areas[c] for c in cells]))
    assert endemism_score(labels, cell_taxa(presence, taxa), areas) == old
