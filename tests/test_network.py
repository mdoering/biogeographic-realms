import math

import polars as pl
import pytest

from realms.network import build_network, information_weight


def presence(rows):
    return pl.DataFrame(rows, schema={"cell": pl.String, "species_id": pl.String}, orient="row")


def taxa(rows):
    return pl.DataFrame(rows, schema={"species_id": pl.String, "genus": pl.String, "family": pl.String}, orient="row")


def links(net):
    return {(net.cells[c], net.taxa[t]): float(w) for c, t, w in zip(net.cell_idx, net.taxon_idx, net.weight)}


def test_weight_one_cell_is_one():
    assert information_weight(1, 100) == 1.0


def test_weight_every_cell_is_zero():
    assert information_weight(100, 100) == 0.0


def test_weight_halfway_on_log_scale():
    # 10 of 100 cells is half way on a log scale
    assert information_weight(10, 100) == pytest.approx(0.5)


def test_weight_single_cell_world():
    assert information_weight(1, 1) == 1.0


def test_species_links_have_weight_one():
    net = build_network(presence([("a", "s1"), ("b", "s1"), ("b", "s2")]))
    assert links(net) == {("a", "species:s1"): 1.0, ("b", "species:s1"): 1.0, ("b", "species:s2"): 1.0}


def test_ancestor_links_use_information_weight():
    # 4 cells; genus G spans two of them through two different species
    p = presence([("a", "s1"), ("b", "s2"), ("c", "s3"), ("d", "s4")])
    t = taxa([("s1", "G", "F1"), ("s2", "G", "F1"), ("s3", "H", "F2"), ("s4", "H", "F2")])
    got = links(build_network(p, t, ranks=("species", "genus")))
    assert got[("a", "genus:G")] == pytest.approx(1 - math.log(2) / math.log(4))
    assert got[("b", "genus:G")] == pytest.approx(0.5)
    assert ("c", "genus:G") not in got


def test_two_species_in_one_cell_make_one_ancestor_link():
    p = presence([("a", "s1"), ("a", "s2"), ("b", "s3")])
    t = taxa([("s1", "G", "F"), ("s2", "G", "F"), ("s3", "H", "F2")])
    net = build_network(p, t, ranks=("genus",))
    assert sorted(links(net)) == [("a", "genus:G"), ("b", "genus:H")]


def test_missing_rank_creates_no_node():
    p = presence([("a", "s1"), ("b", "s2")])
    t = taxa([("s1", None, None), ("s2", None, None)])
    net = build_network(p, t)
    assert net.taxa == ["species:s1", "species:s2"]


def test_species_missing_from_taxa_keeps_its_species_link():
    net = build_network(presence([("a", "s1"), ("b", "s2")]), taxa([("s1", "G", "F")]))
    assert ("b", "species:s2") in links(net)
    assert "genus:G" in net.taxa


def test_same_name_at_two_ranks_stays_two_nodes():
    p = presence([("a", "s1"), ("b", "s2")])
    t = taxa([("s1", "Foo", "Foo"), ("s2", "Bar", "Baz")])
    net = build_network(p, t)
    assert "genus:Foo" in net.taxa and "family:Foo" in net.taxa


def test_zero_weight_links_are_dropped():
    p = presence([("a", "s1"), ("b", "s2")])
    t = taxa([("s1", "G1", "Cosmo"), ("s2", "G2", "Cosmo")])
    net = build_network(p, t)
    assert "family:Cosmo" not in net.taxa
    assert not any(w == 0 for w in net.weight)


def test_species_only_when_taxa_is_none():
    net = build_network(presence([("a", "s1")]))
    assert net.taxa == ["species:s1"]


def test_cell_without_species_stays_in_the_network():
    net = build_network(presence([("a", "s1"), ("empty", None)]))
    assert net.cells == ["a", "empty"]
    assert list(links(net)) == [("a", "species:s1")]


def test_cell_order_does_not_depend_on_row_order():
    a = build_network(presence([("b", "s2"), ("a", "s1")]))
    b = build_network(presence([("a", "s1"), ("b", "s2")]))
    assert a.cells == b.cells == ["a", "b"]
    assert links(a) == links(b)
