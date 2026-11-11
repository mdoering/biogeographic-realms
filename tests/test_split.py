import hashlib

import polars as pl
import pytest

from realms.split import assign_halves, carry_over, crossed, id_changes, load_split, split_sha256


def taxa(rows):
    return pl.DataFrame(
        rows,
        schema={"species_id": pl.String, "group": pl.String, "family": pl.String, "order": pl.String, "class": pl.String},
        orient="row",
    )


def families(sizes, group="g"):
    rows = []
    for f, size in enumerate(sizes):
        rows += [(f"s{f}-{i}", group, f"F{f}", "O", "C") for i in range(size)]
    return taxa(rows)


def species_half(df):
    return dict(zip(df["species_id"], df["half"]))


def test_whole_lineage_stays_in_one_half():
    result = assign_halves(families([3, 5, 2, 7, 4, 1]))
    per_family = result.group_by("family").agg(pl.col("half").n_unique())
    assert per_family["half"].to_list() == [1] * 6
    assert set(result["half"]) == {"A", "B"}


def test_same_result_for_any_input_order():
    df = families(list(range(1, 30)))
    shuffled = df.sample(fraction=1.0, shuffle=True, seed=99)
    assert species_half(assign_halves(df)) == species_half(assign_halves(shuffled))


def test_different_seed_changes_assignment():
    df = families([5] * 40)
    assert species_half(assign_halves(df, seed=1)) != species_half(assign_halves(df, seed=2))


def test_balanced_within_group():
    sizes_one = list(range(1, 51)) * 2
    sizes_two = [3] * 200
    df = pl.concat([families(sizes_one, "plants"), families(sizes_two, "insects")])
    # family names repeat between the two groups in this fixture, so keep lineages separate
    df = df.with_columns((pl.col("group") + "-" + pl.col("family")).alias("family"))
    result = assign_halves(df)
    for group, largest in [("plants", 50), ("insects", 3)]:
        counts = result.filter(pl.col("group") == group)["half"].value_counts()
        by_half = dict(zip(counts["half"], counts["count"]))
        assert abs(by_half["A"] - by_half["B"]) <= largest


def test_missing_family_uses_next_rank():
    result = assign_halves(taxa([("s1", "g", None, "X", "C"), ("s2", "g", None, None, "Y")]))
    lineage = dict(zip(result["species_id"], result["lineage"]))
    assert lineage == {"s1": "order:X", "s2": "class:Y"}


def test_family_lineage_is_named_by_rank():
    result = assign_halves(taxa([("s1", "g", "F", "X", "C")]))
    assert result["lineage"].to_list() == ["family:F"]


def test_no_ranks_at_all_is_its_own_lineage():
    result = assign_halves(taxa([("s1", "g", None, None, None), ("s2", "g", None, None, None)]))
    lineage = dict(zip(result["species_id"], result["lineage"]))
    assert lineage == {"s1": "species:s1", "s2": "species:s2"}


def test_family_with_two_groups_stays_in_one_half():
    # as in the real taxGroup: Asteraceae is mostly angiosperms, with some species grouped as eukaryotes
    rows = [(f"a{i}", "angiosperms", "Asteraceae", "O", "C") for i in range(9)] + [("b0", "eukaryotes", "Asteraceae", "O", "C")]
    for seed in range(20):
        assert assign_halves(taxa(rows), seed=seed)["half"].n_unique() == 1


def test_null_group_still_gets_a_half():
    result = assign_halves(taxa([("s1", None, "F1", "O", "C"), ("s2", None, "F2", "O", "C")]))
    assert result["half"].null_count() == 0


def test_null_group_follows_the_rest_of_its_family():
    rows = [("s1", "g", "F", "O", "C"), ("s2", None, "F", "O", "C")]
    assert assign_halves(taxa(rows))["half"].n_unique() == 1


def test_same_family_name_in_two_kingdoms_is_two_lineages():
    df = taxa([("s1", "g", "Foo", "O", "C"), ("s2", "g", "Foo", "O", "C")]).with_columns(pl.Series("kingdom", ["Animalia", "Plantae"]))
    assert sorted(assign_halves(df)["lineage"]) == ["family:Animalia:Foo", "family:Plantae:Foo"]


def write_split(path, df):
    df.select("species_id", "lineage", "half").sort("species_id").write_parquet(path)


def test_split_sha256_is_the_hash_of_the_sorted_rows():
    df = pl.DataFrame(
        {"species_id": ["s2", "é1", "s10"], "lineage": ["order:O", "species:é1", "family:F"], "half": ["A", "A", "B"]},
    )
    content = "s10\tfamily:F\tB\ns2\torder:O\tA\né1\tspecies:é1\tA\n".encode()
    assert split_sha256(df) == hashlib.sha256(content).hexdigest()
    assert split_sha256(df.reverse()) == split_sha256(df)
    assert split_sha256(df.select("half", "lineage", "species_id", pl.lit(1).alias("extra"))) == split_sha256(df)


def test_split_sha256_refuses_a_missing_value():
    df = pl.DataFrame({"species_id": ["s1", "s2"], "lineage": ["family:F", None], "half": ["A", "B"]})
    with pytest.raises(ValueError, match="null"):
        split_sha256(df)


def test_load_split_returns_the_file_when_its_hash_matches(tmp_path):
    df = assign_halves(families([3, 5, 2, 7]))
    path = tmp_path / "halves.parquet"
    write_split(path, df)
    loaded = load_split(path, sha256=split_sha256(df))
    assert loaded.columns == ["species_id", "lineage", "half"]
    assert species_half(loaded) == species_half(df)


def test_load_split_refuses_a_changed_file(tmp_path):
    df = assign_halves(families([3, 5, 2, 7]))
    expected = split_sha256(df)
    other = {"A": "B", "B": "A"}[df["half"][0]]
    flipped = df.with_columns(pl.when(pl.col("species_id") == df["species_id"][0]).then(pl.lit(other)).otherwise(pl.col("half")).alias("half"))
    assert species_half(flipped) != species_half(df)
    path = tmp_path / "halves.parquet"
    write_split(path, flipped)
    with pytest.raises(ValueError, match="halves.parquet"):
        load_split(path, sha256=expected)


def previous(rows):
    """A split of an earlier release: (species_id, lineage, half) rows."""
    return pl.DataFrame(rows, schema={"species_id": pl.String, "lineage": pl.String, "half": pl.String}, orient="row")


def test_carry_over_keeps_each_known_family_in_its_half_whatever_its_species():
    old = assign_halves(families([3, 5, 2, 7, 4, 1]))
    renamed = families([3, 5, 2, 7, 4, 1]).with_columns(("new-" + pl.col("species_id")).alias("species_id"))
    result = carry_over(renamed, old.select("species_id", "lineage", "half"))
    old_half = dict(zip(old["lineage"], old["half"]))
    assert all(old_half[lineage] == half for lineage, half in zip(result["lineage"], result["half"]))


def test_carry_over_of_an_unchanged_release_is_the_same_split():
    df = families(list(range(1, 30)))
    old = assign_halves(df)
    assert species_half(carry_over(df, old.select("species_id", "lineage", "half"))) == species_half(old)


def test_carry_over_moves_a_species_with_its_new_family():
    old = previous([("s1", "family:F1", "A"), ("s2", "family:F2", "B")])
    result = carry_over(taxa([("s1", "g", "F2", "O", "C"), ("s2", "g", "F2", "O", "C")]), old)
    assert species_half(result) == {"s1": "B", "s2": "B"}


def test_carry_over_puts_a_new_family_in_the_lighter_half_of_its_group():
    old = previous([(f"a{i}", "family:Big", "A") for i in range(10)] + [("b0", "family:Small", "B")])
    rows = [(f"a{i}", "g", "Big", "O", "C") for i in range(10)] + [("b0", "g", "Small", "O", "C")]
    rows += [(f"n{i}", "g", "New", "O", "C") for i in range(3)]
    result = carry_over(taxa(rows), old)
    assert set(result.filter(pl.col("family") == "New")["half"]) == {"B"}


def test_carry_over_balances_new_families_within_each_group():
    old = previous([(f"p{i}", "family:P", "A") for i in range(6)] + [(f"q{i}", "family:Q", "B") for i in range(6)])
    rows = [(f"p{i}", "plants", "P", "O", "C") for i in range(6)] + [(f"q{i}", "insects", "Q", "O", "C") for i in range(6)]
    rows += [("n1", "plants", "NewPlant", "O", "C"), ("n2", "insects", "NewInsect", "O", "C")]
    result = species_half(carry_over(taxa(rows), old))
    assert (result["n1"], result["n2"]) == ("B", "A")


def test_carry_over_same_result_for_any_input_order():
    old = assign_halves(families([4, 4, 4, 4]))
    df = families([4, 4, 4, 4, 3, 2, 2, 1, 1, 1])
    shuffled = df.sample(fraction=1.0, shuffle=True, seed=99)
    keep = old.select("species_id", "lineage", "half")
    assert species_half(carry_over(df, keep)) == species_half(carry_over(shuffled, keep))


def test_carry_over_refuses_a_family_in_both_halves():
    old = previous([("s1", "family:F", "A"), ("s2", "family:F", "B")])
    with pytest.raises(ValueError, match="family:F"):
        carry_over(taxa([("s1", "g", "F", "O", "C")]), old)


def test_crossed_lists_the_species_of_both_splits_whose_half_changed():
    old = previous([("s1", "family:F1", "A"), ("s2", "family:F1", "A"), ("s3", "family:F2", "B")])
    new = previous([("s1", "family:F2", "B"), ("s2", "family:F1", "A"), ("s4", "family:F3", "B")])
    result = crossed(old, new)
    assert result.select("species_id", "lineage_before", "half_before", "lineage", "half").rows() == [("s1", "family:F1", "A", "family:F2", "B")]


def names(rows):
    return pl.DataFrame(rows, schema={"species_id": pl.String, "name": pl.String}, orient="row")


def test_id_changes_pairs_the_ids_of_a_name_found_once_in_each_release():
    before = names([("a", "Aus bus L."), ("b", "Cus dus L."), ("c", "Eus fus L."), ("d", "Eus fus L."), ("e", "Gus hus L.")])
    after = names([("a", "Aus bus L."), ("x", "Cus dus L."), ("y", "Eus fus L."), ("z", "Gus hus L."), ("w", "Gus hus L.")])
    # a kept its id; b became x; Eus fus is twice before and Gus hus twice after, so neither can be paired
    assert id_changes(before, after).rows() == [("b", "x")]


def test_carry_over_puts_a_new_family_where_most_of_its_known_species_were():
    # Margarodidae species move to the new family Monophlebidae; it follows them into B although A is lighter
    old = previous([(f"m{i}", "family:Margarodidae", "B") for i in range(3)] + [("x0", "family:Other", "A"), ("b0", "family:Big", "B")])
    old = pl.concat([old, previous([(f"b{i}", "family:Big", "B") for i in range(1, 10)])])
    rows = [(f"m{i}", "g", "Monophlebidae", "O", "C") for i in range(3)] + [("n0", "g", "Monophlebidae", "O", "C"), ("x0", "g", "Other", "O", "C")]
    rows += [(f"b{i}", "g", "Big", "O", "C") for i in range(10)]
    assert set(carry_over(taxa(rows), old).filter(pl.col("family") == "Monophlebidae")["half"]) == {"B"}


def test_carry_over_places_a_new_family_whose_known_species_tie_by_balance():
    old = previous([("p0", "family:P", "A"), ("q0", "family:Q", "B")] + [(f"a{i}", "family:Heavy", "A") for i in range(10)])
    rows = [("p0", "g", "New", "O", "C"), ("q0", "g", "New", "O", "C")] + [(f"a{i}", "g", "Heavy", "O", "C") for i in range(10)]
    assert set(carry_over(taxa(rows), old).filter(pl.col("family") == "New")["half"]) == {"B"}
