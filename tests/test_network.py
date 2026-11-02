import math
from dataclasses import fields, replace

import infomap
import polars as pl
import pytest

import reference_impl as ref
from realms import cluster
from realms.candidates import candidate_maps
from realms.cluster import InfomapOptions, Partition, run_infomap
from realms.matrix import cell_taxa
from realms.network import build_network, information_weight, network_from
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


DEFAULT_LINE = "--flow-model undirected --seed 1 --num-trials 10 --silent"

# what to set for each option, with the switch that its damping or strength needs, and the flags that adds to the command line,
# which are those of infomap 2.16
EACH_OPTION = {
    "markov_time": ({"markov_time": 0.5}, "--markov-time 0.5"),
    "variable_markov_time": ({"variable_markov_time": True}, "--variable-markov-time"),
    "variable_markov_damping": (
        {"variable_markov_time": True, "variable_markov_damping": 0.5},
        "--variable-markov-time --variable-markov-damping 0.5",
    ),
    "regularized": ({"regularized": True}, "--regularized"),
    "regularization_strength": ({"regularized": True, "regularization_strength": 2.0}, "--regularized --regularization-strength 2.0"),
    "entropy_corrected": ({"entropy_corrected": True}, "--entropy-corrected"),
    "entropy_correction_strength": (
        {"entropy_corrected": True, "entropy_correction_strength": 2.0},
        "--entropy-corrected --entropy-correction-strength 2.0",
    ),
    "recurse_share": ({"recurse_share": 0.5}, ""),  # the recursion is ours: Infomap is not told about it
}


@pytest.fixture
def command_lines(monkeypatch):
    """The option strings `infomap.Infomap` is constructed with, in order."""
    flags = []
    real = infomap.Infomap

    def spy(options=""):
        flags.append(options)
        return real(options)

    monkeypatch.setattr(infomap, "Infomap", spy)
    return flags


def test_default_options_give_todays_command_line(command_lines):
    assert InfomapOptions().args() == ""
    net = planted_net()
    run_infomap(net)
    run_infomap(net, options=InfomapOptions())
    run_infomap(net, options=InfomapOptions(markov_time=1.0))
    run_infomap(net, seed=4, trials=2)
    assert command_lines[:3] == [DEFAULT_LINE] * 3
    assert command_lines[3] == "--flow-model undirected --seed 4 --num-trials 2 --silent"
    assert len(command_lines) == 4


def test_every_option_has_a_case_below():
    assert set(EACH_OPTION) == {f.name for f in fields(InfomapOptions)}


@pytest.mark.parametrize("field", list(EACH_OPTION))
def test_each_option_adds_exactly_its_flag(field, command_lines):
    settings, flags = EACH_OPTION[field]
    options = InfomapOptions(**settings)
    assert options.args() == flags
    run_infomap(planted_net(), options=options)  # Infomap accepts it
    assert set(command_lines) == {DEFAULT_LINE + (f" {flags}" if flags else "")}


def test_the_options_together_keep_the_order_of_their_flags():
    options = InfomapOptions(**{name: value for settings, _ in EACH_OPTION.values() for name, value in settings.items()})
    assert options.args() == (
        "--markov-time 0.5 --variable-markov-time --variable-markov-damping 0.5 --regularized --regularization-strength 2.0 "
        "--entropy-corrected --entropy-correction-strength 2.0"
    )


@pytest.mark.parametrize(
    "settings, message",
    [
        ({"recurse_share": 0.0}, "recurse_share"),
        ({"recurse_share": 1.0}, "recurse_share"),
        ({"recurse_share": -0.1}, "recurse_share"),
        ({"recurse_share": 1.5}, "recurse_share"),
        ({"recurse_share": float("nan")}, "recurse_share"),
        ({"variable_markov_damping": 0.5}, "variable_markov_time"),
        ({"regularization_strength": 2.0}, "regularized"),
        ({"entropy_correction_strength": 2.0}, "entropy_corrected"),
        ({"variable_markov_time": True, "regularization_strength": 2.0}, "regularized"),  # another option's switch does not do
    ],
)
def test_options_that_would_silently_change_nothing_are_refused(settings, message):
    with pytest.raises(ValueError, match=message):
        InfomapOptions(**settings)


def test_a_share_just_inside_the_range_is_accepted():
    assert InfomapOptions(recurse_share=0.001).recurse_share == 0.001
    assert InfomapOptions(recurse_share=0.999).recurse_share == 0.999


def test_markov_time_below_one_gives_more_modules():
    net = planted_net(ranks=("species",))
    assert leaf_modules(run_infomap(net, options=InfomapOptions(markov_time=0.5))) > leaf_modules(run_infomap(net, options=InfomapOptions(markov_time=1.0)))


def groups(labels, cells):
    """The cells of `labels` that are in `cells`, grouped by label."""
    found = {}
    for cell in cells:
        found.setdefault(labels[cell], set()).add(cell)
    return {frozenset(g) for g in found.values()}


def in_label_order(labels, cells):
    """The groups of `cells`, by the label that they have."""
    return sorted(groups(labels, cells), key=lambda group: labels[next(iter(group))])


@pytest.fixture
def giant_module():
    """Eight planted realms of ten cells, and a cell with no species. At Markov time 3 Infomap makes two of the realms modules of their
    own and merges the other five. The cell without species counts in the weights, so this is not the partition of the network
    without it."""
    presence, taxa, _ = planted_network(n_realms=8, cells_per_realm=10, endemic_species=50, widespread_species=8)
    presence = pl.concat([presence, pl.DataFrame({"cell": ["zz-empty"], "species_id": [None]}, schema=presence.schema)])  # sorts last
    ct = cell_taxa(presence, taxa)
    options = InfomapOptions(markov_time=3.0)
    plain = run_infomap(network_from(ct), options=options)
    return ct, options, plain


def the_giant(top):
    """The cells of the largest top module."""
    biggest = max(set(top.values()), key=list(top.values()).count)
    return {c for c, t in top.items() if t == biggest}


def test_recursion_splits_a_planted_giant_module_and_leaves_the_rest(giant_module):
    ct, options, plain = giant_module
    top = plain.top_labels()
    giant = the_giant(top)
    rest = set(top) - giant
    assert (len(giant), len(top)) == (50, 80)  # 62% of the 81 units in one module, which is a single leaf
    assert plain.unassigned == ["zz-empty"]
    before = plain.leaf_labels()
    assert len(groups(before, giant)) == 1 and len(groups(before, rest)) == 3

    recursed = run_infomap(network_from(ct), options=replace(options, recurse_share=0.5))
    after = recursed.leaf_labels()
    assert len(groups(after, giant)) > 1  # the giant module is split
    assert groups(after, rest) == groups(before, rest)  # and the other three stay as they were
    assert in_label_order(before, rest) == in_label_order(after, rest)  # in the same order, though their integers may be higher
    assert {after[c] for c in giant}.isdisjoint({after[c] for c in rest})  # labels stay unique
    assert len(set(after.values())) == len(groups(after, giant)) + 3
    assert recursed.top_labels() == top  # the top modules are what they were: the sub-run goes inside them
    assert recursed.unassigned == ["zz-empty"]
    assert set(recursed.paths) == set(plain.paths)


def test_recursion_runs_a_sub_network_of_the_module_with_the_original_weights(giant_module, monkeypatch):
    ct, options, plain = giant_module
    net = network_from(ct)
    asked = []
    real = cluster.load_infomap

    def spy(network, seed=1, trials=10, options=InfomapOptions()):
        asked.append((network, seed, trials, options))
        return real(network, seed, trials, options)

    monkeypatch.setattr(cluster, "load_infomap", spy)
    recursion = replace(options, recurse_share=0.5)
    run_infomap(net, seed=5, trials=2, options=recursion)

    assert len(asked) == 2  # the whole network, then the one giant module
    assert asked[0][0] is net and asked[0][1:] == (5, 2, recursion)
    sub, seed, trials, sub_options = asked[1]
    assert (seed, trials, sub_options) == (5, 2, options)  # same seed and trials, same options with the recursion off
    giant = the_giant(run_infomap(net, seed=5, trials=2, options=options).top_labels())
    assert sub.cells == [c for c in net.cells if c in giant]  # its units, in the order of the network
    whole = links(net)
    assert links(sub) == {(c, t): w for (c, t), w in whole.items() if c in giant}  # every link of them, with the weight it had
    assert set(sub.taxa) == {t for c, t in whole if c in giant}  # and the taxa that those links lead to, no more


def test_recursion_keeps_a_module_that_the_sub_run_does_not_split():
    net = planted_net()  # six realms, which Infomap finds at the first run
    plain = run_infomap(net)
    assert leaf_modules(plain) == 6
    assert run_infomap(net, options=InfomapOptions(recurse_share=0.1)) == plain  # each holds 17%, each is run again and stays one module


def test_recursion_leaves_every_module_alone_below_the_share(giant_module):
    ct, options, plain = giant_module
    net = network_from(ct)
    assert len(net.cells) == 81  # the giant has 50 of them
    assert run_infomap(net, options=replace(options, recurse_share=50 / 81)) == plain  # exactly its share, not more than that
    assert run_infomap(net, options=replace(options, recurse_share=0.62)) == plain  # 50 of the 81 cells, which is less, and not 50 of the 80 linked
    assert run_infomap(net, options=replace(options, recurse_share=0.61)) != plain  # a little under its share, so it is split


def test_a_recursed_partition_goes_through_candidate_maps(giant_module):
    ct, options, plain = giant_module
    areas = {c: 1e6 for c in ct.cells}
    adjacency = {c: set() for c in ct.cells}
    recursed = run_infomap(network_from(ct), options=replace(options, recurse_share=0.5))
    found = candidate_maps(recursed, ct, areas, adjacency, k_min=2)
    assert len(found.maps) > 0
    assert len(set(found.units.values())) == leaf_modules(recursed)


def test_recursion_on_a_network_without_links_is_still_an_empty_partition():
    presence, _, _ = planted_network(n_realms=2, cells_per_realm=2, endemic_species=1, widespread_species=0)
    empty = pl.DataFrame({"cell": ["b", "a"], "species_id": [None, None]}, schema=presence.schema)
    assert run_infomap(build_network(empty), options=InfomapOptions(recurse_share=0.5)) == Partition(paths={}, unassigned=["a", "b"])
