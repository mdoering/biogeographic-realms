"""The projection and the module summary of `scripts/pilot.py`, on made-up measurements. The pilot itself is run on the real counts."""

import importlib.util
import math
from pathlib import Path

import polars as pl
import pytest

from realms.cluster import Partition

SCRIPT = Path(__file__).parent.parent / "scripts" / "pilot.py"


@pytest.fixture(scope="module")
def pilot():
    spec = importlib.util.spec_from_file_location("pilot_script", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def measured(**changes) -> dict:
    data = {
        "links": 10_000_000,
        "ct_gb": 0.1,
        "chosen_markov_time": 1.0,
        "infomap": {"1.0": {"seconds": 20.0, "peak_gb": 4.5}, "0.75": {"seconds": 30.0, "peak_gb": 4.0}},
        "stages": {
            "infomap_10_trials": {"child_seconds": 100.0, "child_peak_rss_gb": 5.0},
            "candidate_maps_1.0": {"child_seconds": 10.0, "child_peak_rss_gb": 2.0},
            "bootstrap": {"child_peak_rss_gb": 3.0, "workers": [{"seconds": 60.0, "peak_gb": 4.0}, {"seconds": 50.0, "peak_gb": 3.5}]},
            "excess_table": {"tasks_seconds_sum": 96.0, "maps": 28, "child_peak_rss_gb": 3.0, "workers": [{"peak_gb": 1.0}, {"peak_gb": 0.8}]},
        },
        "disk": {"zip_bytes": 1_000_000_000, "parquet_bytes": 400_000_000, "free_gb": 50.0},
    }
    data.update(changes)
    return data


def test_projection_scales_each_stage_from_its_measure(pilot):
    result = pilot.project(measured(), 2.5)
    f, c = 5.0, pilot.CELLS_RATIO
    merges = c * c**2  # c times the merges, each over a dense c squared array, which beats the 5 times the links
    assert result["link_factor"] == f
    assert result["links"] == 50_000_000
    # the worst child peak, 5.0 GB, over 5 times the links, alone and with the matrix beside it
    assert result["a_infomap_peak_gb"] == {"alone": 25.0, "with_matrix": 25.5, "limit": 40, "ok": True}
    candidate = 10.0 * merges
    one_bootstrap = f * (60.0 - 10.0) + candidate  # the worker's time, split into what is linear in links and the merges
    # one worker fits: the parent holds 5 * 3.0 and a worker 5 * 4.0, so (40 - 15) / 20 = 1.25
    assert result["b_selection_hours"]["bootstrap_workers"] == 1
    assert result["b_selection_hours"]["candidate_maps"] == pytest.approx(candidate / 3600, abs=0.01)
    assert result["b_selection_hours"]["bootstraps"] == pytest.approx(20 * one_bootstrap / 3600, abs=0.01)
    assert result["b_selection_hours"]["infomap_10_trials"] == pytest.approx(f * 200.0 / 3600, abs=0.01)  # ten of the 20 s trial beat the 100 s measured
    assert result["b_selection_hours"]["excess_workers"] == 5  # (40 - 5 * 3.0) / (5 * 1.0)
    excess = 96.0 * f * (56 / 28) / 5  # the 96 s of work, for 56 maps and 5 times the links, shared by five workers
    assert result["b_selection_hours"]["excess_table"] == pytest.approx(excess / 3600, abs=0.01)
    total = 1000.0 + candidate + 20 * one_bootstrap + excess
    assert result["b_selection_hours"]["total"] == pytest.approx(total / 3600, abs=0.01)
    assert result["b_selection_hours"]["without_candidate_maps"] == pytest.approx((1000.0 + 20 * f * 50.0 + excess) / 3600, abs=0.01)
    assert result["memory_gb"]["candidate_maps"] == pytest.approx(2.0 * c**2, abs=0.1)  # dense arrays: c squared beats the 5 times the links
    assert result["memory_gb"]["bootstrap_worker"] == 20.0
    assert result["c_disk_left_gb"] == {"after_download": 46.5, "download_gb": 3.5, "minimum": 10, "ok": True}
    assert result["b_selection_hours"]["ok"] is True  # 3.5 hours
    assert result["go"] is True


def test_the_smaller_growth_costs_less(pilot):
    large, small = pilot.project(measured(), 2.5), pilot.project(measured(), 1.38)
    assert small["link_factor"] == pytest.approx(2.76)
    assert small["a_infomap_peak_gb"]["alone"] < large["a_infomap_peak_gb"]["alone"]
    assert small["b_selection_hours"]["total"] < large["b_selection_hours"]["total"]
    assert small["c_disk_left_gb"]["after_download"] > large["c_disk_left_gb"]["after_download"]


def test_each_condition_alone_can_say_stay(pilot):
    assert pilot.project(measured(infomap={"1.0": {"seconds": 20.0, "peak_gb": 9.0}}), 2.5)["go"] is False  # 45 GB
    slow = measured()
    slow["stages"]["infomap_10_trials"] = {"child_seconds": 100_000.0, "child_peak_rss_gb": 5.0}
    assert pilot.project(slow, 2.5)["b_selection_hours"]["ok"] is False
    full = measured(disk={"zip_bytes": 1_000_000_000, "parquet_bytes": 400_000_000, "free_gb": 12.0})
    assert pilot.project(full, 2.5)["c_disk_left_gb"]["ok"] is False


def test_a_worker_that_does_not_fit_is_a_failure_whatever_the_time(pilot):
    big = measured()
    big["stages"]["bootstrap"]["workers"] = [{"seconds": 60.0, "peak_gb": 9.0}]  # 5 * 9 = 45 GB for one worker
    result = pilot.project(big, 2.5)
    assert result["b_selection_hours"]["bootstrap_workers"] == 1
    assert result["b_selection_hours"]["ok"] is False
    assert result["go"] is False


def test_workers_are_at_most_one_fewer_than_the_cores(pilot):
    tiny = measured()
    tiny["stages"]["excess_table"]["workers"] = [{"peak_gb": 0.001}]
    assert pilot.project(tiny, 2.5)["b_selection_hours"]["excess_workers"] == pilot.CORES - 1
    assert math.isfinite(pilot.project(tiny, 2.5)["b_selection_hours"]["total"])


def test_the_excess_without_zero_area_parts_stands_in_for_a_failed_one(pilot):
    assert pilot.project(measured(), 2.5)["excess_measured_on"] == "excess_table"
    data = measured()
    data["stages"]["excess_table_without_zero_area_parts"] = data["stages"]["excess_table"]
    data["stages"]["excess_table"] = {"seconds": 3.0, "error": "ZeroDivisionError: float division by zero", "child_peak_rss_gb": 3.0}
    result = pilot.project(data, 2.5)
    assert result["excess_measured_on"] == "excess_table_without_zero_area_parts"
    assert result["b_selection_hours"]["excess_table"] == pytest.approx(pilot.project(measured(), 2.5)["b_selection_hours"]["excess_table"])


def test_structure_counts_modules_and_shares(pilot):
    paths = {"a:land": (0, 1), "b:sea": (0, 2), "c:land": (1, 1), "d:sea": (1, 1), "e:land": (1, 1)}
    info = pl.DataFrame(
        {
            "part": ["a:land", "b:sea", "c:land", "d:sea", "e:land", "f:sea"],
            "kind": ["land", "sea", "land", "sea", "land", "sea"],
            "area_km2": [2_000_000.0, 500_000.0, 400_000.0, 300_000.0, 100_000.0, 1.0],
            "lat": [10.0, -10.0, 50.0, 52.0, 54.0, 0.0],
            "lon": [0.0, 1.0, 100.0, 101.0, 102.0, 5.0],
        }
    )
    found = pilot.structure(Partition(paths, ["f:sea"]), info)
    assert found["parts"] == 5
    assert found["unassigned"] == 1
    assert found["levels"] == 2
    assert found["top_modules"] == 2
    assert found["leaf_modules"] == 3
    assert found["largest_top_share_of_parts"] == 0.6  # the second module holds 3 of 5
    assert found["largest_leaf_share_of_parts"] == 0.6
    assert found["largest_top_share_of_area"] == round(2_500_000 / 3_300_000, 4)
    assert found["leaf_modules_of_one_part"] == 2
    assert found["leaf_modules_under_floor"] == 2  # the leaves (0, 2) and (1, 1) are under 1,000,000 km2, (0, 1) is not
    first, second = found["top_level_modules"]
    assert (first["module"], first["parts"], first["sea_share"], first["leaf_modules"]) == (1, 3, 0.333, 1)
    assert (second["module"], second["parts"], second["area_km2"], second["sea_share"]) == (0, 2, 2_500_000.0, 0.5)
    assert second["median_lat"] == 0.0
