import polars as pl
import pytest

from realms.seeds import across_seeds, choose


def table(rows):
    return pl.DataFrame(rows, schema={"k": pl.Int64, "stability": pl.Float64, "excess": pl.Float64}, orient="row")


def test_the_median_over_seeds_counts_a_missing_map_as_zero_stability():
    tables = {
        1: table([(7, 0.81, 0.15), (8, 0.80, 0.14), (9, 0.82, 0.16)]),
        2: table([(7, 0.75, 0.16), (8, 0.83, 0.12), (9, 0.85, 0.10)]),
        3: table([(7, 0.90, 0.14), (8, 0.84, 0.13)]),  # no k=9: stability 0 there
    }
    got = across_seeds(tables)
    assert got["k"].to_list() == [7, 8, 9]
    assert got["stability"].to_list() == pytest.approx([0.81, 0.83, 0.82])
    assert got["excess"].to_list() == pytest.approx([0.15, 0.13, 0.13])  # over the seeds that have the map
    assert got["seeds"].to_list() == [3, 3, 2]
    assert got["stability_min"].to_list() == pytest.approx([0.75, 0.80, 0.0])
    assert got["stability_max"].to_list() == pytest.approx([0.90, 0.84, 0.85])


def test_a_map_most_seeds_lack_cannot_be_stable():
    tables = {s: table([(7, 0.9, 0.1)]) for s in (1, 2, 3)}
    tables[1] = table([(7, 0.9, 0.1), (10, 0.95, 0.2)])
    got = across_seeds(tables)
    assert got.filter(pl.col("k") == 10)["stability"].item() == 0.0


def test_across_seeds_wants_some_seed():
    with pytest.raises(ValueError, match="no seed"):
        across_seeds({})


def summary(rows):
    return pl.DataFrame(rows, schema={"k": pl.Int64, "stability": pl.Float64, "excess": pl.Float64}, orient="row")


def test_choose_takes_the_highest_excess_among_stable_maps_of_every_configuration():
    got = choose({
        "a": summary([(6, 0.79, 0.30), (7, 0.81, 0.15)]),
        "b": summary([(8, 0.85, 0.17), (9, 0.82, 0.18)]),
    })
    assert (got.config, got.k, got.fell_short) == ("b", 8, False)  # 9 is within 0.01 of 0.18, so the coarser 8
    assert "0.180" in got.reason and "3 stable maps" in got.reason


def test_choose_breaks_a_tie_in_k_by_the_excess():
    got = choose({"a": summary([(7, 0.81, 0.150)]), "b": summary([(7, 0.90, 0.155)])})
    assert (got.config, got.k) == ("b", 7)


def test_choose_without_a_stable_map_reports_the_most_stable_one_as_falling_short():
    got = choose({"a": summary([(6, 0.70, 0.3), (7, 0.76, 0.1)]), "b": summary([(5, 0.76, 0.2)])})
    assert (got.config, got.k, got.fell_short) == ("b", 5, True)  # a tie in stability goes to the coarser map
    assert "no map reached" in got.reason


def selected(k=6, fell_short=True, ks=(5, 6, 7)):
    from pathlib import Path

    from realms.explore import Partition, Selected

    parts = pl.DataFrame({"part": [f"p{k_}" for k_ in ks], "k": list(ks), "realm": [0] * len(ks), "selected": [k_ == k for k_ in ks]})
    return Selected(directory=Path("run/mt090/select"), stamp={}, k=k, fell_short=fell_short, parts=parts, candidates=pl.DataFrame(),
                    selection=pl.DataFrame(), partition=Partition(paths={}, unassigned=[]), merged=False)


def test_apply_choice_puts_the_k_chosen_over_seeds_in_place_of_seed_one():
    from realms.seeds import apply_choice

    got = apply_choice(selected(), {"config": "mt090", "k": 7, "fell_short": False}, "mt090")
    assert (got.k, got.fell_short) == (7, False)
    assert got.parts.filter(pl.col("selected"))["k"].unique().to_list() == [7]


def test_apply_choice_leaves_another_configuration_or_no_choice_alone():
    from realms.seeds import apply_choice

    found = selected()
    assert apply_choice(found, {"config": "mt080", "k": 7, "fell_short": False}, "mt090") is found
    assert apply_choice(found, None, "mt090") is found


def test_apply_choice_refuses_a_k_the_result_has_no_map_of():
    from realms.seeds import apply_choice

    with pytest.raises(ValueError, match="k=9"):
        apply_choice(selected(), {"config": "mt090", "k": 9, "fell_short": False}, "mt090")
