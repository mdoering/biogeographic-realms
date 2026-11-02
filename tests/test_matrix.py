from collections import Counter

import numpy as np
import polars as pl
import pytest

import reference_impl as ref
from realms.matrix import RANKS, cell_taxa, duplicate_names
from realms.network import build_network, network_from
from realms.synth import grid_presence, planted_network


def frame(rows, columns):
    return pl.DataFrame(rows, schema={c: pl.String for c in columns}, orient="row")


def hand_made_case():
    """Four cells and every awkward row: a duplicate, a null species, a species missing from taxa, one with no genus.

    Genus G0 is in every cell (so its weight is 0), family F2 too, and `s9` is listed in taxa but never seen.
    """
    presence = frame(
        [
            ("a", "s1"), ("a", "s1"), ("a", "s2"), ("a", "s5"),
            ("b", "s2"), ("b", "s4"), ("b", "s5"), ("b", "ghost"),
            ("c", "s3"), ("c", "s4"), ("c", None),
            ("d", "s3"), ("d", "s5"), ("d", "s1"), ("d", "ghost"),
        ],
        ["cell", "species_id"],
    )  # fmt: skip
    taxa = frame(
        [("s1", "G0", "F1"), ("s2", "G0", "F1"), ("s3", "G0", "F2"), ("s4", None, "F2"), ("s5", "G1", "F2"), ("s9", "G9", "F9")],
        ["species_id", "genus_id", "family_id"],
    )
    return presence, taxa


CASES = {"planted": planted_network()[:2], "hand made": hand_made_case(), "grid": grid_presence(300, side=20)[:2]}


def assert_same_network(new, old):
    assert new.cells == old.cells and new.taxa == old.taxa
    assert np.array_equal(new.cell_idx, old.cell_idx) and np.array_equal(new.taxon_idx, old.taxon_idx)
    assert np.array_equal(new.weight, old.weight)  # bit for bit


@pytest.mark.parametrize("presence, taxa", CASES.values(), ids=CASES)
def test_network_from_matches_reference(presence, taxa):
    assert_same_network(network_from(cell_taxa(presence, taxa)), ref.build_network(presence, taxa))


@pytest.mark.parametrize("ranks", [("species",), ("genus",), ("family",), ("species", "genus"), ("genus", "family")])
@pytest.mark.parametrize("presence, taxa", CASES.values(), ids=CASES)
def test_build_network_matches_reference_for_any_choice_of_ranks(presence, taxa, ranks):
    assert_same_network(build_network(presence, taxa, ranks), ref.build_network(presence, taxa, ranks))


@pytest.mark.parametrize("presence, taxa", CASES.values(), ids=CASES)
def test_without_taxa_only_species_become_nodes(presence, taxa):
    new = build_network(presence)
    assert_same_network(new, ref.build_network(presence))
    assert set(t.split(":")[0] for t in new.taxa) == {"species"}


def test_empty_presence_gives_an_empty_network():
    empty = frame([], ["cell", "species_id"])
    new = build_network(empty)
    assert new.cells == [] and new.taxa == [] and len(new.weight) == 0


@pytest.mark.parametrize("presence, taxa", CASES.values(), ids=CASES)
def test_matrix_layout(presence, taxa):
    ct = cell_taxa(presence, taxa)
    assert ct.cells == sorted(ct.cells) and ct.taxa == sorted(ct.taxa)
    assert ct.x.shape == (len(ct.cells), len(ct.taxa)) and ct.x.dtype == bool and ct.x.has_sorted_indices
    assert ct.rank.dtype == np.int8 and [RANKS[r] for r in ct.rank] == [t.split(":")[0] for t in ct.taxa]
    for column in (ct.genus_col, ct.family_col):
        assert column.dtype == np.int64 and len(column) == len(ct.taxa)


def test_parent_columns_point_at_the_species_own_genus_and_family():
    presence, taxa = hand_made_case()
    ct = cell_taxa(presence, taxa)
    wanted = {r["species_id"]: r for r in taxa.iter_rows(named=True)}
    seen = 0
    for j, name in enumerate(ct.taxa):
        if not name.startswith("species:"):
            assert ct.genus_col[j] == ct.family_col[j] == -1
            continue
        seen += 1
        genus, family = wanted[name.removeprefix("species:")]["genus_id"], wanted[name.removeprefix("species:")]["family_id"]
        assert (ct.taxa[ct.genus_col[j]] if ct.genus_col[j] >= 0 else None) == (f"genus:{genus}" if genus else None)
        assert (ct.taxa[ct.family_col[j]] if ct.family_col[j] >= 0 else None) == (f"family:{family}" if family else None)
    assert seen == 5  # s1 to s5, and not s9, which has no cell


def test_cells_keep_their_place_when_their_species_are_null_or_not_accepted():
    presence = frame([("a", "s1"), ("b", "ghost"), ("c", None)], ["cell", "species_id"])
    ct = cell_taxa(presence, frame([("s1", "G", "F")], ["species_id", "genus_id", "family_id"]))
    assert ct.cells == ["a", "b", "c"]
    assert ct.taxa == ["family:F", "genus:G", "species:s1"]
    net = network_from(ct)
    assert net.cells == ["a", "b", "c"]
    assert set(net.cell_idx.tolist()) == {0}  # b and c are rows without links
    assert ct.rows(["b", "c"]).x.nnz == 0


@pytest.mark.parametrize("presence, taxa", CASES.values(), ids=CASES)
def test_rows_match_reference_incidence(presence, taxa):
    ct = cell_taxa(presence, taxa)
    ranked = ref.rank_presence(presence, taxa, RANKS)
    names = sorted(set(ranked["species_id"].to_list()))  # the reference's columns, in its order
    assert ct.taxa == names

    order = np.random.default_rng(5).permutation(len(ct.cells))
    asked = [ct.cells[i] for i in order[: max(2, len(order) // 2)]] + ["nowhere"]
    new, old = ct.rows(asked), ref.incidence(ranked, asked)

    def pairs(table):
        coo = table.x.tocoo()
        return {(asked[r], names[c]) for r, c in zip(coo.row, coo.col)}

    assert new.x.shape == old.x.shape == (len(asked), len(names))
    assert pairs(new) == pairs(old) and new.x.nnz == old.x.nnz
    assert set(new.x.data.tolist()) <= {1}
    assert dict(zip(names, new.totals.tolist())) == dict(zip(names, old.totals.tolist()))
    assert all(new.x[len(asked) - 1].toarray().ravel() == 0)  # the unknown cell


def test_rows_totals_count_cells_outside_the_rows_asked_for():
    presence, taxa = hand_made_case()
    ct = cell_taxa(presence, taxa)
    assert np.array_equal(ct.rows(["a"]).totals, ct.rows(ct.cells).totals)
    assert np.array_equal(ct.rows(["a"]).totals, np.asarray(ct.x.sum(axis=0)).ravel())


def test_rows_keep_the_order_asked_for_and_repeat_a_cell_asked_twice():
    presence, taxa = hand_made_case()
    ct = cell_taxa(presence, taxa)
    everything = ct.rows(ct.cells).x.toarray()
    got = ct.rows(["d", "a", "d"]).x.toarray()
    assert np.array_equal(got, everything[[3, 0, 3]])


def test_incidence_transpose_is_the_taxa_by_cells_matrix_built_once():
    presence, taxa = hand_made_case()
    table = cell_taxa(presence, taxa).rows(["d", "a", "b"])
    assert table.xt.shape == table.x.T.shape and table.xt.format == "csr"
    assert (table.xt != table.x.T).nnz == 0
    assert table.xt is table.xt


def test_species_ids_follow_column_order():
    presence, taxa = hand_made_case()
    assert cell_taxa(presence, taxa).species_ids() == ["s1", "s2", "s3", "s4", "s5"]
    assert cell_taxa(presence, taxa).select(["genus"]).species_ids() == []


def test_select_species_only():
    presence, taxa = hand_made_case()
    ct = cell_taxa(presence, taxa)
    only = ct.select(["species"])
    species = [j for j, t in enumerate(ct.taxa) if t.startswith("species:")]
    assert only.cells == ct.cells and only.taxa == [ct.taxa[j] for j in species]
    assert (only.x != ct.x[:, species]).nnz == 0 and only.x.has_sorted_indices
    assert (only.rank == 0).all() and (only.genus_col == -1).all() and (only.family_col == -1).all()


def test_select_remaps_the_parent_columns_that_remain():
    presence, taxa = hand_made_case()
    kept = cell_taxa(presence, taxa).select(["species", "genus"])
    assert kept.taxa == ["genus:G0", "genus:G1", "species:s1", "species:s2", "species:s3", "species:s4", "species:s5"]
    assert [kept.taxa[j] if j >= 0 else None for j in kept.genus_col] == [None, None, "genus:G0", "genus:G0", "genus:G0", None, "genus:G1"]
    assert (kept.family_col == -1).all()


def test_cell_taxa_ranks_argument_picks_the_same_columns_as_select():
    presence, taxa = hand_made_case()
    full = cell_taxa(presence, taxa)
    for ranks in (("species",), ("genus",), ("family",), ("genus", "family")):
        direct, picked = cell_taxa(presence, taxa, ranks), full.select(ranks)
        assert direct.taxa == picked.taxa and (direct.x != picked.x).nnz == 0
        assert np.array_equal(direct.genus_col, picked.genus_col) and np.array_equal(direct.family_col, picked.family_col)


def test_unknown_rank_is_refused():
    presence, taxa = hand_made_case()
    with pytest.raises(ValueError, match="order"):
        cell_taxa(presence, taxa, ("species", "order"))
    with pytest.raises(ValueError, match="order"):
        cell_taxa(presence, taxa).select(["order"])


def test_a_species_listed_twice_in_taxa_is_refused():
    presence, _ = hand_made_case()
    twice = frame([("s1", "G0", "F1"), ("s1", "G1", "F1")], ["species_id", "genus_id", "family_id"])
    with pytest.raises(ValueError, match="s1"):
        cell_taxa(presence, twice)


@pytest.mark.parametrize("presence, taxa", [CASES["planted"], CASES["grid"]], ids=["planted", "grid"])
def test_resample_matches_reference(presence, taxa):
    ct = cell_taxa(presence, taxa)
    n = len(ct.species_ids())
    draws = np.random.default_rng(3).integers(n, size=n)
    new = ct.resample(draws)
    assert new.taxa == sorted(new.taxa) and new.x.has_sorted_indices
    assert_same_network(network_from(new), ref.build_network(*ref._resample(presence, taxa, np.random.default_rng(3))))


def test_resample_drops_cells_that_no_drawn_species_occupies_like_the_reference():
    presence, taxa = CASES["grid"]
    ct = cell_taxa(presence, taxa)
    draws = np.random.default_rng(3).integers(len(ct.species_ids()), size=len(ct.species_ids()))
    new = ct.resample(draws)
    old = ref.build_network(*ref._resample(presence, taxa, np.random.default_rng(3)))
    assert len(new.cells) == len(old.cells) < len(ct.cells)  # the weights count cells, so this is a result and not a detail
    assert new.cells == old.cells


def test_resample_by_hand():
    # s1 twice and s2 once. G1, F2 and the species without a genus are not drawn, so they have no column. Cell c holds only
    # s3 and s4, which are not drawn either, so it is not in the resample
    presence, taxa = hand_made_case()
    drawn = cell_taxa(presence, taxa).resample(np.array([0, 0, 1]))
    assert drawn.cells == ["a", "b", "d"]
    assert drawn.taxa == ["family:F1", "genus:G0", "species:s1#0", "species:s1#1", "species:s2#2"]
    assert drawn.x.toarray().astype(int).tolist() == [[1, 1, 1, 1, 1], [1, 1, 0, 0, 1], [1, 1, 1, 1, 0]]
    assert [drawn.taxa[j] if j >= 0 else None for j in drawn.genus_col] == [None, None, "genus:G0", "genus:G0", "genus:G0"]
    assert [drawn.taxa[j] if j >= 0 else None for j in drawn.family_col] == [None, None, "family:F1", "family:F1", "family:F1"]


def test_resample_sorts_draw_numbers_as_text():
    presence, taxa = hand_made_case()
    drawn = cell_taxa(presence, taxa).resample(np.zeros(11, dtype=int))
    assert drawn.taxa[2:] == sorted(f"species:s1#{i}" for i in range(11))
    assert drawn.taxa[2:5] == ["species:s1#0", "species:s1#1", "species:s1#10"]


def test_resample_of_species_only_matrix_has_no_parents():
    presence, taxa = hand_made_case()
    drawn = cell_taxa(presence, taxa).select(["species"]).resample(np.array([4, 4, 0]))
    assert drawn.taxa == ["species:s1#2", "species:s5#0", "species:s5#1"]
    assert (drawn.genus_col == -1).all() and (drawn.family_col == -1).all()


def test_grid_presence_reproduces_the_baseline_generator():
    # figures from the scratch generator the baseline was measured with, at 1,000 species on 90 x 90 cells
    presence, _, _ = grid_presence(1000)
    cells = presence["cell"].str.slice(1).cast(pl.Int64)
    species = presence["species_id"].str.slice(1).cast(pl.Int64)
    assert presence.height == 34751
    assert cells.sum() == 131798633
    assert (cells * (species + 1)).sum() % 1000003 == 941599


def test_grid_presence_is_repeatable_and_seeded():
    assert grid_presence(60)[0].equals(grid_presence(60)[0])
    assert not grid_presence(60)[0].equals(grid_presence(60, seed=2)[0])


def test_grid_taxa_and_adjacency():
    presence, taxa, adjacency = grid_presence(250, side=10)
    assert taxa.columns == ["species_id", "genus_id", "family_id"] and taxa.height == 250
    assert taxa.row(137) == ("s137", "g27", "f1")
    assert set(presence["species_id"]) == set(taxa["species_id"])  # a species is always somewhere
    assert set(presence["cell"]) <= set(adjacency) and len(adjacency) == 100
    assert adjacency["c11"] == {"c01", "c10", "c12", "c21"}
    assert all(a in adjacency[b] for a, neighbours in adjacency.items() for b in neighbours)
    assert Counter(len(n) for n in adjacency.values()) == {2: 4, 3: 32, 4: 64}


NAMES = ["species_id", "genus_id", "family_id", "genus", "family", "order", "kingdom"]


def duplicates_case():
    """Genus G1 and G2 are both Abax in Carabidae, family F1 and F2 are both Carabidae in Coleoptera. s5 is another kingdom's Abax."""
    presence = frame([("a", "s1"), ("a", "s3"), ("b", "s2"), ("b", "s4"), ("c", "s5"), ("c", "s1")], ["cell", "species_id"])
    taxa = frame(
        [
            ("s1", "G1", "F1", "Abax", "Carabidae", "Coleoptera", "Animalia"),
            ("s2", "G1", "F1", "Abax", "Carabidae", "Coleoptera", "Animalia"),
            ("s3", "G2", "F2", "Abax", "Carabidae", "Coleoptera", "Animalia"),
            ("s4", "G2", "F2", "Abax", "Carabidae", "Coleoptera", "Animalia"),
            ("s5", "G3", "F3", "Abax", "Carabidae", "Coleoptera", "Plantae"),
        ],
        NAMES,
    )
    return presence, taxa


def test_two_ids_with_one_name_merge():
    presence, taxa = duplicates_case()
    found = duplicate_names(taxa)
    assert found.columns == ["rank", "kingdom", "name", "parent", "ids"]
    assert found.sort("rank").rows() == [
        ("family", "Animalia", "Carabidae", "Coleoptera", ["F1", "F2"]),
        ("genus", "Animalia", "Abax", "Carabidae", ["G1", "G2"]),
    ]
    merged = cell_taxa(presence, taxa, merge_duplicate_names=True)
    assert merged.taxa == ["family:F1", "family:F3", "genus:G1", "genus:G3", "species:s1", "species:s2", "species:s3", "species:s4", "species:s5"]
    j = {t: i for i, t in enumerate(merged.taxa)}
    for s, genus, family in [("s1", "G1", "F1"), ("s2", "G1", "F1"), ("s3", "G1", "F1"), ("s4", "G1", "F1"), ("s5", "G3", "F3")]:
        assert merged.genus_col[j[f"species:{s}"]] == j[f"genus:{genus}"]
        assert merged.family_col[j[f"species:{s}"]] == j[f"family:{family}"]
    row = {c: i for i, c in enumerate(merged.cells)}
    assert merged.x[row["a"], j["genus:G1"]] and merged.x[row["b"], j["genus:G1"]]  # G2's species s3 and s4 count for G1
    assert merged.x[row["a"], j["family:F1"]] and not merged.x[row["a"], j["family:F3"]]


def test_one_name_in_two_kingdoms_does_not_merge():
    presence, taxa = duplicates_case()
    kingdoms = taxa.filter(pl.col("species_id").is_in(["s1", "s5"]))
    assert duplicate_names(kingdoms).height == 0
    merged = cell_taxa(presence, kingdoms, merge_duplicate_names=True)
    assert merged.taxa == cell_taxa(presence, kingdoms).taxa
    assert "genus:G3" in merged.taxa and "genus:G1" in merged.taxa


def test_one_genus_name_in_two_families_does_not_merge():
    presence, taxa = duplicates_case()
    homonyms = taxa.with_columns(
        pl.when(pl.col("genus_id") == "G2").then(pl.lit("Pterostichidae")).otherwise(pl.col("family")).alias("family"),
        pl.when(pl.col("genus_id") == "G2").then(pl.lit("F9")).otherwise(pl.col("family_id")).alias("family_id"),
    ).filter(pl.col("kingdom") == "Animalia")
    assert duplicate_names(homonyms, ranks=("genus",)).height == 0
    assert "genus:G2" in cell_taxa(presence, homonyms, merge_duplicate_names=True).taxa
    nameless = taxa.with_columns(pl.when(pl.col("genus_id") == "G2").then(None).otherwise(pl.col("family")).alias("family"))
    assert duplicate_names(nameless, ranks=("genus",)).height == 0  # a lineage that isn't known never merges


def test_the_default_is_unchanged():
    presence, taxa = duplicates_case()
    plain, explicit = cell_taxa(presence, taxa), cell_taxa(presence, taxa, merge_duplicate_names=False)
    assert plain.taxa == explicit.taxa == [
        "family:F1", "family:F2", "family:F3", "genus:G1", "genus:G2", "genus:G3", "species:s1", "species:s2", "species:s3", "species:s4", "species:s5",
    ]  # fmt: skip
    assert (plain.x != explicit.x).nnz == 0
    assert np.array_equal(plain.genus_col, explicit.genus_col) and np.array_equal(plain.family_col, explicit.family_col)
