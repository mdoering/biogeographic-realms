"""The projection and the module summary of `scripts/pilot.py`, on made-up measurements. The pilot itself is run on the real counts."""

import importlib.util
import json
import math
import sys
import time
import types
from pathlib import Path

import polars as pl
import pytest

from realms import select
from realms.cluster import Partition
from realms.matrix import cell_taxa
from realms.synth import planted_network

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
            "bootstrap": {
                "child_peak_rss_gb": 3.0,
                "workers": [
                    {"seconds": 60.0, "candidate_maps_seconds": 4.0, "peak_gb": 4.0},
                    {"seconds": 50.0, "candidate_maps_seconds": 3.0, "peak_gb": 3.5},
                ],
            },
            "excess_table": {"tasks_seconds_sum": 96.0, "maps": 28, "child_peak_rss_gb": 3.0, "workers": [{"peak_gb": 1.0}, {"peak_gb": 0.8}]},
        },
        "disk": {"zip_bytes": 1_000_000_000, "parquet_bytes": 400_000_000, "free_gb": 50.0},
    }
    data.update(changes)
    return data


def test_projection_scales_each_stage_from_its_measure(pilot):
    result = pilot.project(measured(), 2.5)
    f, c = 5.0, pilot.CELLS_RATIO
    merges = c * f  # c times the merges, each as slow as the 5 times the links, which beat the c times the units of a merge's vectors
    assert result["link_factor"] == f
    assert result["links"] == 50_000_000
    # the worst child peak, 5.0 GB, over 5 times the links, alone and with the matrix beside it
    assert result["a_infomap_peak_gb"] == {"alone": 25.0, "with_matrix": 25.5, "limit": 40, "ok": True}
    candidate = 10.0 * merges
    # a worker's time is split by its own timing: 56 s of the first is not candidate maps, and 4 s are
    one_bootstrap = f * 56.0 + 4.0 * merges
    # two workers fit: the parent holds 5 * 3.0 and a worker 5 * 4.0, so (40 - 15) / 20 = 1.25 is one, and the candidate maps' 2 * 5 = 10 GB is less
    assert result["b_selection_hours"]["bootstrap_workers"] == 1
    assert result["b_selection_hours"]["candidate_maps"] == pytest.approx(candidate / 3600, abs=0.01)
    assert result["b_selection_hours"]["bootstraps"] == pytest.approx(20 * one_bootstrap / 3600, abs=0.01)
    assert result["b_selection_hours"]["one_bootstrap_hours"] == pytest.approx(one_bootstrap / 3600, abs=0.01)
    assert result["b_selection_hours"]["of_which_candidate_maps"] == pytest.approx(4.0 * merges / 3600, abs=0.01)
    assert result["b_selection_hours"]["infomap_10_trials"] == pytest.approx(f * 200.0 / 3600, abs=0.01)  # ten of the 20 s trial beat the 100 s measured
    assert result["b_selection_hours"]["excess_workers"] == 5  # (40 - 5 * 3.0) / (5 * 1.0)
    excess = 96.0 * f * (56 / 28) / 5  # the 96 s of work, for 56 maps and 5 times the links, shared by five workers
    assert result["b_selection_hours"]["excess_table"] == pytest.approx(excess / 3600, abs=0.01)
    total = 1000.0 + candidate + 20 * one_bootstrap + excess
    assert result["b_selection_hours"]["total"] == pytest.approx(total / 3600, abs=0.01)
    assert result["b_selection_hours"]["without_candidate_maps"] == pytest.approx((1000.0 + 20 * f * 56.0 + excess) / 3600, abs=0.01)
    assert result["memory_gb"]["candidate_maps"] == pytest.approx(2.0 * f, abs=0.1)  # the links, with no array of units squared
    assert result["memory_gb"]["bootstrap_worker"] == 20.0
    assert result["c_disk_left_gb"] == {"after_download": 46.5, "download_gb": 3.5, "minimum": 10, "ok": True}
    assert result["unmeasured"] == []
    assert result["b_selection_hours"]["ok"] is True  # 3.5 hours
    assert result["go"] is True


def test_a_bootstrap_shorter_than_the_headline_candidate_maps_is_still_costed(pilot):
    # the real data's path: the headline's candidate_maps (500 s) took longer than a whole bootstrap worker (300 s), which holds a
    # 3-trial Infomap, the resample and the network besides its own 50 s of candidate maps, and that part must not be costed at nothing
    data = measured()
    data["stages"]["candidate_maps_1.0"]["child_seconds"] = 500.0
    data["stages"]["bootstrap"]["workers"] = [{"seconds": 300.0, "candidate_maps_seconds": 50.0, "peak_gb": 4.0}]
    result = pilot.project(data, 2.5)["b_selection_hours"]
    f, merges = 5.0, pilot.CELLS_RATIO * 5.0
    assert result["one_bootstrap_hours"] == pytest.approx((f * 250.0 + 50.0 * merges) / 3600, abs=0.01)
    assert result["without_candidate_maps"] > 20 * f * 250.0 / 3600  # the bootstraps' own work is in it, and more than the 10-trial Infomap and the excess


def test_the_candidate_maps_of_a_worker_scale_apart_from_the_rest_of_it(pilot):
    data = measured()
    data["stages"]["bootstrap"]["workers"] = [{"seconds": 60.0, "candidate_maps_seconds": 0.0, "peak_gb": 4.0}]
    plain = pilot.project(data, 2.5)["b_selection_hours"]["one_bootstrap_hours"]
    data["stages"]["bootstrap"]["workers"] = [{"seconds": 60.0, "candidate_maps_seconds": 30.0, "peak_gb": 4.0}]
    split = pilot.project(data, 2.5)["b_selection_hours"]["one_bootstrap_hours"]
    c = pilot.CELLS_RATIO
    assert plain == pytest.approx(5.0 * 60.0 / 3600, abs=0.01)
    assert split == pytest.approx((5.0 * 30.0 + 30.0 * c * 5.0) / 3600, abs=0.01)  # 30 s are scaled by the merges instead of the links


def test_a_merge_is_as_slow_as_the_slower_of_the_cells_and_the_links(pilot):
    # at growth 1.38 the headline has 2.76 times the links of half A and the grid 3.0 times the cells, so the cells are the slower
    result = pilot.project(measured(), 1.38)
    c = pilot.CELLS_RATIO
    assert result["b_selection_hours"]["candidate_maps"] == pytest.approx(10.0 * c * c / 3600, abs=0.01)


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
    big["stages"]["bootstrap"]["workers"] = [{"seconds": 60.0, "candidate_maps_seconds": 4.0, "peak_gb": 9.0}]  # 5 * 9 = 45 GB for one worker
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


def test_neither_excess_stage_measured_leaves_it_out_and_is_not_a_go(pilot):
    data = measured()
    failed = {"seconds": 3.0, "error": "ZeroDivisionError: float division by zero", "child_peak_rss_gb": 3.0}
    data["stages"]["excess_table"] = failed
    data["stages"]["excess_table_without_zero_area_parts"] = dict(failed)
    result = pilot.project(data, 2.5)  # no KeyError for the workers that never ran
    assert result["unmeasured"] == ["excess_table"]
    assert result["excess_measured_on"] is None
    assert result["b_selection_hours"]["excess_table"] is None
    c, f = pilot.CELLS_RATIO, 5.0
    assert result["b_selection_hours"]["total"] == pytest.approx((1000.0 + 10.0 * c * f + 20 * (f * 56.0 + 4.0 * c * f)) / 3600, abs=0.01)  # the rest, as a lower bound
    assert result["go"] is False  # everything else passes, and an unmeasured stage is not a pass
    pilot.add_projection(data)
    assert data["decision"].startswith("undecided")
    assert "excess_table" in data["decision"]


def test_a_stage_that_was_not_run_at_all_is_not_a_crash_either(pilot):
    data = measured()
    del data["stages"]["excess_table"]
    assert pilot.project(data, 2.5)["unmeasured"] == ["excess_table"]


def test_project_only_needs_only_the_results(pilot, tmp_path, capsys):
    (tmp_path / "res6").mkdir()
    (tmp_path / "res6" / "measure.json").write_text(json.dumps(measured()))
    pilot.main(["--project-only", "--results", str(tmp_path)])
    saved = json.loads((tmp_path / "res6" / "measure.json").read_text())
    assert set(saved["projection"]) == {"2.5", "1.38"}
    assert saved["decision"] == "go to resolution 7"
    assert "decision at growth 2.5: go to resolution 7" in capsys.readouterr().out


def test_project_only_without_a_measure_file_says_so(pilot, tmp_path):
    with pytest.raises(SystemExit, match="measure.json"):
        pilot.main(["--project-only", "--results", str(tmp_path)])


def test_measuring_needs_the_counts_and_the_workers(pilot, tmp_path):
    with pytest.raises(SystemExit) as error:
        pilot.main(["--results", str(tmp_path)])
    assert error.value.code == 2


def test_a_worker_times_the_function_it_is_asked_to_watch(pilot, monkeypatch):
    stage = types.ModuleType("watched_stage")
    stage.slow = lambda: time.sleep(0.06)
    monkeypatch.setitem(sys.modules, "watched_stage", stage)

    def task(shared, item):
        stage.slow()  # looked up when it is called, so the replacement in the worker is what runs
        stage.slow()
        time.sleep(0.1)
        return item * 2

    result, info = pilot._with_peak(task, ("watched_stage", "slow"), None, 21)
    assert result == 42
    assert 0.12 <= info["slow_seconds"] < info["seconds"]  # two calls of 0.06 s, and the 0.1 s outside them is not counted
    assert info["seconds"] >= 0.2
    assert info["peak_gb"] > 0
    assert stage.slow is not None and stage.slow.__name__ == "<lambda>"  # put back afterwards
    plain = pilot._with_peak(task, None, None, 1)[1]
    assert "slow_seconds" not in plain


def test_the_bootstrap_reports_each_workers_candidate_maps(pilot, monkeypatch):
    # in this process, which `map_tasks` does with one worker: the same function the workers run
    presence, taxa, truth = planted_network(n_realms=4, cells_per_realm=12, endemic_species=30, widespread_species=6, seed=1)
    ct = cell_taxa(presence, taxa)
    areas = {c: 4e5 for c in truth}
    adjacency = {c: set() for c in truth}
    labels = dict(truth)
    real = select.candidate_maps

    def slow(*args, **kwargs):
        time.sleep(0.2)
        return real(*args, **kwargs)

    monkeypatch.setattr(pilot, "BOOT_PROBE", 1)
    monkeypatch.setattr(select, "candidate_maps", slow)
    result = pilot._bootstrap(ct, {4: labels}, areas, adjacency, 1.0)
    [worker] = result["workers"]
    assert set(result["scores"]) == {4}
    assert 0.2 <= worker["candidate_maps_seconds"] <= worker["seconds"] + 0.05  # the seconds are rounded to 0.1, the Infomap trials are the rest
    assert worker["peak_gb"] > 0
    assert select.candidate_maps is slow  # put back afterwards


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
