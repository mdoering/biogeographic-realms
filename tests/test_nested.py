import pytest

from helpers import matrix
from realms.nested import at_cuts, bootstrap_recovery, leaf_realms, levels_of, nested_merges, realm_recovery, sea_realms


def taxa_rows(unit_taxa):
    return [(unit, t) for unit, taxa in unit_taxa.items() for t in taxa]


def test_leaf_realms_merge_while_the_closest_pair_is_below_the_cut():
    # a and b share 4 of 5 (beta-sim 0.2), c shares 2 of 5 with them (0.6)
    ct = matrix(taxa_rows({"a": ["s0", "s1", "s2", "s3", "s4"], "b": ["s1", "s2", "s3", "s4", "s5"], "c": ["s0", "s1", "t2", "t3", "t4"]}))
    labels, merges = leaf_realms({"a": 0, "b": 1, "c": 2}, ct, cut=0.5)
    assert labels == {"a": 0, "b": 0, "c": 2}
    assert merges["height"].to_list() == pytest.approx([0.2]) and merges["k"].to_list() == [3]
    assert leaf_realms({"a": 0, "b": 1, "c": 2}, ct, cut=0.1)[0] == {"a": 0, "b": 1, "c": 2}


def test_sea_realms_go_by_the_area_of_sea_parts():
    labels = {"1:sea": 0, "2:sea": 0, "2:land": 0, "3:land": 1, "4:sea": 1}
    areas = {"1:sea": 5.0, "2:sea": 1.0, "2:land": 2.0, "3:land": 3.0, "4:sea": 1.0}
    assert sea_realms(labels, areas) == {0: True, 1: False}


def test_nested_merges_join_land_with_land_and_sea_with_sea_before_land_with_sea():
    # the sea realm 0 shares most with land 1 (beta-sim 0.2), but land 1 must first join land 2 (0.6)
    unit_taxa = {
        "o:sea": ["m0", "m1", "s0", "s1", "s2", "s3"],
        "x:land": ["s0", "s1", "s2", "s3", "x4"],
        "y:land": ["s0", "y1", "y2", "x4", "y4"],
    }
    labels = {"o:sea": 0, "x:land": 1, "y:land": 2}
    merges, maps = nested_merges(labels, matrix(taxa_rows(unit_taxa)), {0: True, 1: False, 2: False})
    assert merges.select("k", "keep", "drop", "kind").rows() == [(3, 1, 2, "land"), (2, 0, 1, "both")]
    assert maps[2] == {"o:sea": 0, "x:land": 1, "y:land": 1} and set(maps[1].values()) == {0}
    assert merges["height"][0] == pytest.approx(0.6)


def test_at_cuts_takes_the_map_before_the_first_merge_at_or_above_the_cut():
    unit_taxa = {
        "o:sea": ["m0", "m1", "s0", "s1", "s2", "s3"],
        "x:land": ["s0", "s1", "s2", "s3", "x4"],
        "y:land": ["s0", "y1", "y2", "x4", "y4"],
    }
    labels = {"o:sea": 0, "x:land": 1, "y:land": 2}
    merges, maps = nested_merges(labels, matrix(taxa_rows(unit_taxa)), {0: True, 1: False, 2: False})
    found = at_cuts(merges, maps, [0.5, 0.7])
    assert found[0.5] == labels  # the first merge is at 0.6
    assert found[0.7] == maps[2]  # land merged, but never land with sea, whatever that join's height


def test_realm_recovery_is_the_best_jaccard_by_area_on_the_common_cells():
    reference = {"a": 0, "b": 0, "c": 1, "d": 1}
    replicate = {"a": 5, "b": 6, "c": 6, "e": 7}  # d is missing from the replicate
    areas = {"a": 1.0, "b": 1.0, "c": 2.0, "d": 9.0, "e": 1.0}
    got = realm_recovery(reference, replicate, areas)
    assert got == {0: pytest.approx(0.5), 1: pytest.approx(2 / 3)}  # 0: {a} of {a, b}; 1: {c} against {b, c}


def test_levels_of_refuses_a_cut_at_or_below_the_leaf_cut():
    with pytest.raises(ValueError, match="above the leaf cut"):
        levels_of({"a": 0}, matrix([("a", "s")]), {"a": 1.0}, 0.5, [0.5])


def test_bootstrap_recovery_finds_a_planted_map_again():
    import numpy as np

    from helpers import grid
    from realms.cluster import InfomapOptions
    from realms.matrix import cell_taxa
    from realms.synth import planted_network

    presence, taxa, truth = planted_network(n_realms=3, cells_per_realm=10, endemic_species=40, widespread_species=3)
    positions, adjacency = grid(6, 5)
    ordered = sorted(truth, key=lambda c: (truth[c], int(c.split("-c")[1])))
    name = dict(zip(positions, ordered))
    adjacency = {name[p]: {name[n] for n in adjacency[p]} for p in positions}
    areas = {c: 1.0 for c in ordered}
    ct = cell_taxa(presence, taxa)
    reference = {0.5: truth}
    table = bootstrap_recovery(ct, areas, adjacency, reference, 0.5, [], InfomapOptions(), n_boot=3, trials=1)
    assert table.height == 3 * 3 and table["jaccard"].min() > 0.9
    assert table.equals(bootstrap_recovery(ct, areas, adjacency, reference, 0.5, [], InfomapOptions(), n_boot=3, trials=1, workers=2))
    assert np.isclose(table["jaccard"].max(), 1.0)


# the cut rules


def test_replay_heights_gives_the_reference_heights_on_the_same_matrix():
    unit_taxa = {
        "o:sea": ["m0", "m1", "s0", "s1", "s2", "s3"],
        "x:land": ["s0", "s1", "s2", "s3", "x4"],
        "y:land": ["s0", "y1", "y2", "x4", "y4"],
    }
    from realms.nested import replay_heights

    labels = {"o:sea": 0, "x:land": 1, "y:land": 2}
    ct = matrix(taxa_rows(unit_taxa))
    merges, _ = nested_merges(labels, ct, {0: True, 1: False, 2: False})
    assert replay_heights(labels, merges, ct) == pytest.approx(merges["height"].to_list())


def test_robust_levels_need_one_cut_to_fit_in_most_replicates():
    import numpy as np
    import polars as pl

    from realms.nested import MERGES, robust_levels

    merges = pl.DataFrame([{"k": 4, "height": 0.2, "keep": 0, "drop": 1, "kind": "land"}, {"k": 3, "height": 0.3, "keep": 0, "drop": 2, "kind": "land"},
                           {"k": 2, "height": 0.8, "keep": 0, "drop": 3, "kind": "land"}], schema=MERGES)
    # the first gap (0.2 to 0.3) is crossed in half of the replicates, the second (0.3 to 0.8) never
    replicates = np.array([[0.2, 0.3, 0.8], [0.35, 0.25, 0.8], [0.2, 0.3, 0.75], [0.31, 0.3, 0.8]])
    table = robust_levels(merges, replicates)
    assert table["realms"].to_list() == [3, 2]
    assert table["separated"].to_list() == pytest.approx([0.5, 1.0])
    assert table["robust"].to_list() == [False, True]
    assert table["cut"].to_list() == pytest.approx([0.25, 0.55])


def test_robust_levels_have_no_cut_where_a_later_merge_is_lower():
    import numpy as np
    import polars as pl

    from realms.nested import MERGES, robust_levels

    merges = pl.DataFrame([{"k": 3, "height": 0.5, "keep": 0, "drop": 1, "kind": "land"}, {"k": 2, "height": 0.4, "keep": 0, "drop": 2, "kind": "land"}], schema=MERGES)
    table = robust_levels(merges, np.array([[0.5, 0.4]]))
    assert table.row(0, named=True)["cut"] is None and not table["robust"][0]


def test_realm_excess_finds_a_planted_realm_and_not_a_random_one():
    import numpy as np

    from helpers import grid
    from realms.matrix import cell_taxa
    from realms.nested import realm_excess
    from realms.synth import planted_network

    presence, taxa, truth = planted_network(n_realms=3, cells_per_realm=10, endemic_species=40, widespread_species=3)
    positions, adjacency = grid(6, 5)
    ordered = sorted(truth, key=lambda c: (truth[c], int(c.split("-c")[1])))
    name = dict(zip(positions, ordered))
    adjacency = {name[p]: {name[n] for n in adjacency[p]} for p in positions}
    areas = {c: 1.0 for c in ordered}
    ct = cell_taxa(presence, taxa)
    table = realm_excess(truth, ct, adjacency, areas, n_null=19, seed=1)
    # on a 6 x 5 grid a null region can be the planted realm itself, so one realm's p is not at its floor of 1/20
    assert table.height == 3 and (table["excess"] > 0.5).all() and (table["p"] <= 0.05).sum() >= 2
    rng = np.random.default_rng(3)
    shuffled = dict(zip(ordered, rng.permutation([truth[c] for c in ordered]).tolist()))
    random = realm_excess(shuffled, ct, adjacency, areas, n_null=19, seed=1)
    assert (random["excess"].abs() < 0.05).all() and (random["p"] > 0.05).all()


@pytest.mark.parametrize("measure", ["concordance", "endemic"])
def test_class_signal_keeps_a_class_that_follows_the_realms_and_drops_one_that_does_not(measure):
    import numpy as np
    import polars as pl

    from helpers import grid
    from realms.matrix import cell_taxa
    from realms.nested import class_signal

    positions, adjacency = grid(6, 6)
    cells = [f"c{i}" for i in range(36)]
    name = dict(zip(positions, cells))
    adjacency = {name[p]: {name[n] for n in adjacency[p]} for p in positions}
    labels = {c: (0 if i % 6 < 3 else 1) for i, c in enumerate(cells)}  # two realms, the left and the right half
    rows = []
    for s in range(40):  # endemic: each species in 3 cells of one half
        half = s % 2
        rows += [(c, f"E{s}") for c in [cells[r * 6 + half * 3 + (s % 3)] for r in ((s // 3) % 4, (s // 3) % 4 + 1, (s // 3) % 4 + 2)]]
    rng = np.random.default_rng(1)
    for s in range(40):  # scattered: each species in 3 random cells anywhere
        rows += [(c, f"W{s}") for c in rng.choice(cells, 3, replace=False)]
    for s in range(40):  # single cells: endemic to any map, never scored
        rows.append((cells[s % 36], f"S{s}"))
    presence = pl.DataFrame(rows, schema=["cell", "species_id"], orient="row")
    taxa = pl.DataFrame({"species_id": sorted(set(presence["species_id"]))}).with_columns(pl.lit(None, pl.String).alias("genus_id"), pl.lit(None, pl.String).alias("family_id"))
    ct = cell_taxa(presence, taxa)
    classes = {f"E{s}": "Endemica" for s in range(40)} | {f"W{s}": "Vagantia" for s in range(40)} | {f"S{s}": "Solitaria" for s in range(40)}
    table = class_signal(labels, ct, adjacency, {c: 1.0 for c in cells}, classes, n_null=19, measure=measure)
    got = {r["class"]: r for r in table.iter_rows(named=True)}
    assert got["Endemica"]["observed"] == 1.0 and not got["Endemica"]["drop"]
    assert got["Vagantia"]["drop"] and got["Vagantia"]["testable"]
    assert got["Solitaria"]["scored"] == 0 and not got["Solitaria"]["testable"] and not got["Solitaria"]["drop"]



@pytest.mark.parametrize(("realms", "reported"), [
    ([5, 6, 7, 8, 9, 14], [5, 9, 14]),  # the eukaryotes headline: a run to 9, then five realms at once
    ([13, 5, 10, 6, 9, 7, 8], [5, 10, 13]),  # unsorted: native_only's levels
    ([5, 6, 7], [5, 7]),  # the run reaches the finest
    ([5, 9, 14], [5, 14]),  # no run past the coarsest
    ([5], [5]),
    ([], []),
])
def test_reported_levels_give_the_coarsest_the_end_of_its_run_and_the_finest(realms, reported):
    from realms.nested import reported_levels

    assert reported_levels(realms) == reported
