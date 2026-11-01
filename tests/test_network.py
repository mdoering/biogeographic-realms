import math

import infomap
import polars as pl
import pytest

import reference_impl as ref
from realms.cluster import Partition, run_infomap
from realms.network import build_network, information_weight
from realms.synth import planted_network


def presence(rows):
    return pl.DataFrame(rows, schema={"cell": pl.String, "species_id": pl.String}, orient="row")


def taxa(rows):
    return pl.DataFrame(rows, schema={"species_id": pl.String, "genus_id": pl.String, "family_id": pl.String}, orient="row")


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


def test_species_not_listed_as_accepted_never_becomes_a_node():
    # taxa holds the accepted species, so a synonym id arriving in the presence table must not be clustered
    net = build_network(presence([("a", "s1"), ("b", "synonym")]), taxa([("s1", "G", "F")]))
    assert net.taxa == ["family:F", "genus:G", "species:s1"]
    assert net.cells == ["a", "b"]


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


def test_homonym_genera_are_separate_nodes_by_id():
    # three accepted Drosophila, and names repeat within a kingdom too: only the id says which genus a species is in
    p = presence([("a", "fly"), ("b", "fungus"), ("c", "other")])
    t = pl.DataFrame(
        {
            "species_id": ["fly", "fungus", "other"],
            "genus": ["Drosophila", "Drosophila", "Gus"],
            "genus_id": ["BXWZ6", "8HFX7", "G9"],
            "kingdom": ["Animalia", "Animalia", "Animalia"],
        }
    )
    net = build_network(p, t, ranks=("genus",))
    assert net.taxa == ["genus:8HFX7", "genus:BXWZ6", "genus:G9"]
    assert all(w == 1.0 for w in net.weight)


def test_species_of_one_genus_id_share_a_node():
    p = presence([("a", "s1"), ("b", "s2"), ("c", "s3")])
    net = build_network(p, taxa([("s1", "G", None), ("s2", "G", None), ("s3", "H", None)]), ranks=("genus",))
    assert net.taxa == ["genus:G", "genus:H"]


def planted_net(ranks=("species", "genus", "family")):
    presence, taxa, _ = planted_network(n_realms=6, cells_per_realm=10, endemic_species=30, widespread_species=5)
    return build_network(presence, taxa, ranks)


def test_bulk_links_give_the_reference_partition():
    net = planted_net()
    assert len(set(ref.run_infomap(net).leaf_labels().values())) > 1  # a real clustering, so the comparison says something
    assert run_infomap(net) == ref.run_infomap(net)
    assert run_infomap(net, seed=4, trials=2) == ref.run_infomap(net, seed=4, trials=2)


def test_bulk_links_keep_a_cell_without_links_unassigned():
    presence, taxa, _ = planted_network(n_realms=2, cells_per_realm=5, endemic_species=10, widespread_species=0)
    with_empty = pl.concat([presence, pl.DataFrame({"cell": ["empty"], "species_id": [None]}, schema=presence.schema)])
    net = build_network(with_empty)
    assert run_infomap(net) == ref.run_infomap(net)
    assert run_infomap(net).unassigned == ["empty"]


def test_a_network_with_no_links_leaves_every_cell_unassigned():
    presence, _, _ = planted_network(n_realms=2, cells_per_realm=2, endemic_species=1, widespread_species=0)
    empty = pl.DataFrame({"cell": ["b", "a"], "species_id": [None, None]}, schema=presence.schema)
    net = build_network(empty)
    assert len(net.cell_idx) == 0
    assert run_infomap(net) == ref.run_infomap(net) == Partition(paths={}, unassigned=["a", "b"])


def leaf_modules(partition):
    return len(set(partition.leaf_labels().values()))


def test_markov_time_below_one_gives_more_modules():
    net = planted_net(ranks=("species",))
    assert leaf_modules(run_infomap(net, markov_time=0.5)) > leaf_modules(run_infomap(net, markov_time=1.0))


def test_default_markov_time_leaves_the_command_line_alone(monkeypatch):
    flags = []
    real = infomap.Infomap

    def spy(options=""):
        flags.append(options)
        return real(options)

    monkeypatch.setattr(infomap, "Infomap", spy)
    net = planted_net()
    run_infomap(net)
    run_infomap(net, markov_time=1.0)
    run_infomap(net, markov_time=0.5)
    assert flags[0] == flags[1] == "--flow-model undirected --seed 1 --num-trials 10 --silent"
    assert flags[2] == flags[0] + " --markov-time 0.5"
