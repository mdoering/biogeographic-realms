import polars as pl
import pytest

from realms.crosscheck import area_agreement, pairwise_v, unit_labels

SCHEMA = {"part": pl.String, "k": pl.Int64, "realm": pl.Int64, "unit": pl.String, "imputed": pl.Boolean, "selected": pl.Boolean}


def parts(rows):
    return pl.DataFrame(rows, schema=SCHEMA, orient="row")


def test_unit_labels_come_from_the_parts_that_are_not_imputed():
    table = parts([
        ("a:land", 7, 0, "u1", False, True),
        ("a:sea", 7, 0, "u1", False, True),  # two parts of one unit
        ("b:land", 7, 1, "u2", False, True),
        ("c:land", 7, 4, None, True, True),  # imputed: no unit
        ("a:land", 5, 3, "u1", False, False),
        ("b:land", 5, 3, "u2", False, False),
    ])
    assert unit_labels(table, 7) == {"u1": 0, "u2": 1}
    assert unit_labels(table, 5) == {"u1": 3, "u2": 3}
    with pytest.raises(ValueError, match="k=6"):
        unit_labels(table, 6)
    split = parts([("a:land", 7, 0, "u1", False, True), ("a:sea", 7, 1, "u1", False, True)])
    with pytest.raises(ValueError, match="u1"):
        unit_labels(split, 7)
    lost = parts([("a:land", 7, 0, None, False, True)])
    with pytest.raises(ValueError, match="no unit"):
        unit_labels(lost, 7)


def test_area_agreement_matches_realms_one_to_one_by_area():
    a = {"u1": 0, "u2": 0, "u3": 1, "u4": 2}
    b = {"u1": 5, "u2": 6, "u3": 6, "u4": 7}
    areas = {"u1": 10.0, "u2": 1.0, "u3": 4.0, "u4": 5.0}
    # 0 -> 5 (10), 1 -> 6 (4), 2 -> 7 (5): 19 of 20; 0 -> 6 would match only 1 + 4 of it
    assert area_agreement(a, b, areas) == pytest.approx(19 / 20)
    assert area_agreement(a, {u: 9 - r for u, r in a.items()}, areas) == pytest.approx(1.0)
    # one realm against three: only one of them can be matched, the largest
    assert area_agreement(a, {u: 0 for u in a}, areas) == pytest.approx(11 / 20)
    # only the units both maps hold count
    assert area_agreement({**a, "u9": 0}, b, {**areas, "u9": 100.0}) == pytest.approx(19 / 20)
    with pytest.raises(ValueError, match="area"):
        area_agreement(a, b, {"u1": 1.0})


def test_pairwise_v_scores_every_pair_once():
    maps = {1: {"x": 0, "y": 0, "z": 1}, 2: {"x": 5, "y": 5, "z": 6}, 3: {"x": 0, "y": 1, "z": 1}}
    found = pairwise_v(maps)
    assert [(r["a"], r["b"]) for r in found] == [(1, 2), (1, 3), (2, 3)]
    assert found[0]["v_measure"] == pytest.approx(1.0) and found[1]["v_measure"] < 1
