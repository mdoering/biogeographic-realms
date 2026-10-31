import polars as pl

from realms.split import assign_halves


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
