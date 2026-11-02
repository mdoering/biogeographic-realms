"""The exploration pipeline on a small planted world: four realms of 30 cells, each with its own species, and three cells too thin to enter the network.

The world has the cells' land parts only (the sea parts have no area). It goes through the real `clean`, `network_units`, Infomap, `candidate_maps`, the bootstrap and the null
maps, with fewer trials and nulls than the analysis, so that a run takes seconds.
"""

import inspect
import json
from dataclasses import dataclass, replace
from pathlib import Path

import geopandas as gpd
import polars as pl
import pytest
from shapely import box

from clean_helpers import counts, make_grid, refs, row
from realms import explore
from realms.clean import SWITCHABLE, Settings
from realms.cluster import InfomapOptions
from realms.explore import Params
from realms.grid import Grid
from realms.score import v_measure
from realms.synth import planted_network
from realms.units import UnitSettings, Units

LAYER_5 = ("detached_satellite", "under_reported", "thin_long_range")
SETTINGS = Settings(rules_off=frozenset(LAYER_5))  # the planted species are spread evenly, and layer 5 has nothing to find in them: keep it out of the way
UNITS = UnitSettings(min_species=30, min_group_species=5, min_groups=2, impute_k=3)
FAST = Params(k_min=2, k_max=10, trials_screen=2, n_null_screen=5, trials_select=2, n_boot=4, boot_trials=1, n_null_select=9)
THIN_CELLS = ("9001", "9002", "9003")


@dataclass
class World:
    counts: pl.DataFrame
    refs: object
    grid: Grid
    coarse: Grid
    coarse_gpkg: Path
    regions: pl.DataFrame
    truth: dict[str, int]  # land part -> planted realm
    list_codes: pl.DataFrame
    marker: dict


@pytest.fixture(scope="module")
def world(tmp_path_factory) -> World:
    presence, taxa, truth_by_cell = planted_network()
    cells = {c: str(1000 + i) for i, c in enumerate(truth_by_cell)}
    # a realm is a block of 5 x 6 cells, side by side along the equator: longitude is the column and latitude the row
    centres, edges = {}, []
    for realm in range(4):
        for i, old in enumerate(c for c, r in truth_by_cell.items() if r == realm):
            centres[cells[old]] = (float(i // 5), float(5 * realm + i % 5))
    for thin, lon in zip(THIN_CELLS, (21.0, 22.0, 23.0), strict=True):
        centres[thin] = (0.0, lon)
    by_place = {centre: name for name, centre in centres.items()}
    for (lat, lon), name in by_place.items():
        for other in ((lat + 1, lon), (lat, lon + 1)):
            if other in by_place:
                edges.append((name, by_place[other]))
    grid = make_grid({name: (100_000.0, 0.0) for name in centres}, edges, centres)
    grid = replace(grid, resolution=7)

    rows = [row(s, f"{cells[c]}:land") for c, s in presence.select("cell", "species_id").iter_rows()]
    thin_species = [f"e0-{s}" for s in range(5)]  # five species in each thin cell, far under the 30 that the network needs
    rows += [row(s, f"{cell}:land") for cell in THIN_CELLS for s in thin_species]
    species = sorted(taxa["species_id"])
    all_taxa = pl.DataFrame(
        {
            "species_id": species,
            "kingdom": "Animalia",
            "phylum": "Arthropoda",
            "class": "Insecta",
            "environment": None,
            "extinct": False,
        },
        schema={"species_id": pl.String, "kingdom": pl.String, "phylum": pl.String, "class": pl.String, "environment": pl.String, "extinct": pl.Boolean},
    ).join(taxa.select("species_id", "genus_id", "family_id"), on="species_id")
    groups = [(s, "angiosperms" if i % 2 else "coleoptera") for i, s in enumerate(species)]
    references = replace(refs(groups=groups), taxa=all_taxa)

    tmp = tmp_path_factory.mktemp("coarse")
    gpkg = tmp / "ISEA3H-5.gpkg"
    gpd.GeoDataFrame(
        {"lat_c": [3.0], "lon_c": [12.0], "res": [5.0], "geogrid_id": ["C1"]}, geometry=[box(-1, -1, 30, 10)], crs="EPSG:4326"
    ).to_file(gpkg, driver="GPKG")
    coarse = Grid(
        5,
        pl.DataFrame({"cell": ["C1"], "lat": [3.0], "lon": [12.0], "area_km2": [1e6], "land_km2": [5e5], "sea_km2": [5e5]}),
        {"C1": set()},
    )
    parts = [f"{name}:land" for name in centres]
    regions = pl.DataFrame({"part": parts, "region": "Africa"}, schema={"part": pl.String, "region": pl.String})
    truth = {f"{cells[c]}:land": r for c, r in truth_by_cell.items()}
    codes = pl.DataFrame(schema={"country": pl.String, "subdivision": pl.String, "code": pl.String})
    return World(counts(*rows), references, grid, coarse, gpkg, regions, truth, codes, {"rows": len(rows), "resolution": 7})


def run_chain(out: Path, w: World, params: Params = FAST, *, settings: Settings = SETTINGS, stage: str = "select", force: bool = False, options=InfomapOptions()) -> dict:
    """Every stage in order up to `stage`, as `run_pipeline` does, on the loaded world. Returns each stage's result."""
    done: dict = {}
    inputs = explore.clean_inputs(settings, explore.LAYERS, w.list_codes, w.marker)
    done["clean"] = explore.clean_stage(out, w.counts, w.refs, w.grid, settings, explore.LAYERS, w.list_codes, inputs, force=force)
    done["units"] = explore.units_stage(out, w.grid, w.coarse, w.coarse_gpkg, UNITS, w.regions, done["clean"]["hash"], force=force)
    if stage in ("screen", "select"):
        units, light, ct = explore.read_units(out), explore.read_light_grid(out), explore.network_matrix(out)
        done["screen"] = explore.screen_stage(out, "default", units, ct, light, options, params, done["units"]["hash"], force=force)
        if stage == "select":
            done["select"] = explore.select_stage(out, "default", units, ct, light, options, params, done["units"]["hash"], force=force)
    return done


@pytest.fixture(scope="module")
def planted(world, tmp_path_factory):
    out = tmp_path_factory.mktemp("explore") / "res7" / "native_only"
    return out, run_chain(out, world)


def test_planted_realms_are_recovered_through_every_stage(world, planted):
    out, done = planted
    assert [k for k, v in done.items() if v.get("skipped")] == []
    chosen = json.loads((out / "default" / "selected.json").read_text())
    assert chosen["k"] == 4 and not chosen["fell_short"]

    candidates = pl.read_parquet(out / "default" / "candidates.parquet")
    assert set(candidates.columns) == {"unit", "k", "realm", "sub_realm", "sub_realm_of"}
    found = dict(candidates.filter(pl.col("k") == 4).select("unit", "realm").iter_rows())
    assert len(found) == 120 and v_measure(found, {p: world.truth[p] for p in found}) >= 0.95

    # the map by part: the selected k is marked, and the thin parts have the realm of their neighbours
    parts = pl.read_parquet(out / "default" / "parts.parquet")
    picked = parts.filter(pl.col("selected"))
    assert picked["k"].unique().to_list() == [4] and picked.height == 123
    assert v_measure(dict(picked.select("part", "realm").iter_rows()), {**world.truth, **{f"{c}:land": 3 for c in THIN_CELLS}}) >= 0.95

    selection = pl.read_csv(out / "default" / "selection.csv")
    assert {"k", "stability", "excess"} <= set(selection.columns)
    assert selection.filter(pl.col("k") == 4)["stability"][0] >= 0.8
    screen = pl.read_csv(out / "default" / "screen.csv")
    assert screen["k"].to_list() == [2, 3, 4] and screen["excess"].is_finite().all()

    structure = json.loads((out / "default" / "structure.json").read_text())
    assert structure["degenerate"] is False and structure["modules"] == 4 and structure["single_unit_modules"] == 0
    assert 0.25 <= structure["largest_module_share_of_units"] <= 0.3 and structure["levels"] >= 1
    partition = pl.read_parquet(out / "default" / "partition.parquet")
    assert partition.columns == ["unit", "leaf", "top"] and partition.height == 120

    for name in ("presence.parquet", "summary.csv", "major.parquet", "taxa.parquet"):
        assert (out / "cleaned" / name).exists()
    for name in ("members", "presence", "areas", "adjacency", "excluded", "stats"):
        assert (out / "units" / f"{name}.parquet").exists()
    assert (out / "units" / "removed_by_region.csv").exists()
    log = json.loads((out / "default" / "log.json").read_text())
    assert {"screen", "select"} <= set(log) and log["select"]["seconds"] > 0 and log["select"]["peak_rss_gb"] > 0


def test_a_stage_is_reused_when_its_inputs_are_unchanged(world, tmp_path, monkeypatch):
    out = tmp_path / "res7" / "native_only"
    calls = {"clean": 0, "network_units": 0, "run_infomap": 0}
    for name, module in (("clean", explore), ("network_units", explore), ("run_infomap", explore)):
        real = getattr(module, name)

        def counted(*args, _real=real, _name=name, **kwargs):
            calls[_name] += 1
            return _real(*args, **kwargs)

        monkeypatch.setattr(module, name, counted)

    first = run_chain(out, world, stage="screen")
    assert calls == {"clean": 1, "network_units": 1, "run_infomap": 1} and not any(v.get("skipped") for v in first.values())
    again = run_chain(out, world, stage="screen")
    assert calls == {"clean": 1, "network_units": 1, "run_infomap": 1} and all(v["skipped"] for v in again.values())
    assert again["screen"]["hash"] == first["screen"]["hash"]

    # a change of the configuration reruns the clustering and nothing before it
    other = run_chain(out, world, stage="screen", options=InfomapOptions(markov_time=0.9))
    assert calls == {"clean": 1, "network_units": 1, "run_infomap": 2}
    assert [k for k, v in other.items() if not v["skipped"]] == ["screen"] and other["screen"]["hash"] != first["screen"]["hash"]
    # a change of the cleaning reruns everything after it, through the hashes
    changed = run_chain(out, world, stage="screen", settings=replace(SETTINGS, rules_off=frozenset(LAYER_5) | {"griis"}))
    assert calls == {"clean": 2, "network_units": 2, "run_infomap": 3} and not any(v["skipped"] for v in changed.values())
    # --force reruns the stage it is given, and the ones before it stay
    forced = explore.clean_stage(out, world.counts, world.refs, world.grid, SETTINGS, explore.LAYERS, world.list_codes, explore.clean_inputs(SETTINGS, explore.LAYERS, world.list_codes, world.marker), force=True)
    assert calls["clean"] == 3 and not forced["skipped"]
    # a screen that is run again takes the select result it overwrote out of date
    run_chain(out, world, stage="select")
    assert (out / "default" / "selected.json").exists()
    explore.screen_stage(out, "default", explore.read_units(out), explore.network_matrix(out), explore.read_light_grid(out), InfomapOptions(), FAST, explore.read_stamp(out / "units")["hash"], force=True)
    assert not (out / "default" / "selected.json").exists() and not (out / "default" / "selection.csv").exists()


def test_a_degenerate_configuration_is_recorded_not_raised(world, tmp_path):
    out = tmp_path / "res7" / "native_only"
    done = run_chain(out, world, replace(FAST, floor_km2=1e12))  # no realm is as large as that, so they all merge into one
    structure = json.loads((out / "default" / "structure.json").read_text())
    assert structure["degenerate"] is True and "above the size floor" in structure["degenerate_reason"]
    assert structure["modules"] == 4
    chosen = json.loads((out / "default" / "selected.json").read_text())
    assert chosen["degenerate"] is True and "k" not in chosen
    assert pl.read_parquet(out / "default" / "candidates.parquet").height == 0
    assert pl.read_csv(out / "default" / "screen.csv").height == 0 and pl.read_csv(out / "default" / "selection.csv").height == 0
    assert done["select"]["degenerate"] and done["screen"]["degenerate"]


def test_a_network_without_links_is_degenerate_too(world, tmp_path):
    out = tmp_path / "res7" / "native_only"
    run_chain(out, world, stage="units")
    units = explore.read_units(out)
    nothing = pl.DataFrame(schema={"unit": pl.String, "species_id": pl.String, "records": pl.Int64})
    empty = replace(units, presence=nothing)
    ct = explore.network_matrix(out, empty)
    explore.screen_stage(out, "default", empty, ct, explore.read_light_grid(out), InfomapOptions(), FAST, "x", force=True)
    structure = json.loads((out / "default" / "structure.json").read_text())
    assert structure["degenerate"] is True and structure["modules"] == 0


def test_excluded_parts_are_never_in_the_network_or_the_scores(world, planted, monkeypatch):
    out, _ = planted
    units = explore.read_units(out)
    left_out = {f"{c}:land" for c in THIN_CELLS}
    assert set(units.excluded["part"]) == left_out and set(units.excluded["reason"]) == {"sparse"}
    assert not left_out & set(units.members) and not left_out & {p for u in units.adjacency.values() for p in u}
    # the matrix the network is built from has the units of the network and nothing else
    ct = explore.network_matrix(out)
    assert set(ct.cells) == set(units.areas) and len(ct.cells) == 120
    assert set(units.presence["unit"]) == set(units.areas)
    # no map that is scored has them, and no unit of a candidate or the partition is one of them
    candidates = pl.read_parquet(out / "default" / "candidates.parquet")
    assert set(candidates["unit"]) == set(units.areas)
    assert set(pl.read_parquet(out / "default" / "partition.parquet")["unit"]) <= set(units.areas)
    seen = []
    real = explore.excess_table
    monkeypatch.setattr(explore, "excess_table", lambda maps, ct, adjacency, areas, **kw: seen.append((maps, adjacency, areas)) or real(maps, ct, adjacency, areas, **kw))
    explore.select_stage(out, "default", units, ct, explore.read_light_grid(out), InfomapOptions(), FAST, explore.read_stamp(out / "units")["hash"], force=True)
    maps, adjacency, areas = seen[0]
    assert all(set(labels) == set(units.areas) for labels in maps.values()) and set(areas) == set(units.areas) and set(adjacency) == set(units.areas)
    # in the outputs only, they have a realm, and the parts file says that it was imputed
    parts = pl.read_parquet(out / "default" / "parts.parquet").filter(pl.col("k") == 4)
    imputed = parts.filter(pl.col("imputed"))
    assert set(imputed["part"]) == left_out and imputed["unit"].is_null().all()
    assert not parts.filter(~pl.col("imputed"))["part"].is_in(list(left_out)).any()


def test_units_without_a_module_are_left_out_like_the_thin_ones(world, planted):
    out, _ = planted
    units = explore.read_units(out)
    unit = sorted(units.areas)[0]
    kept = explore.labelled_units(units, set(units.areas) - {unit})
    assert unit not in kept.areas and unit not in kept.adjacency and all(unit not in near for near in kept.adjacency.values())
    assert unit not in kept.members.values()
    lost = kept.excluded.filter(pl.col("reason") == "unlinked")
    assert lost["part"].to_list() == [unit] and lost["tried_unit"].to_list() == [unit]
    assert explore.labelled_units(units, set(units.areas)) is units


def test_only_half_a_is_loaded(world, tmp_path, monkeypatch):
    asked = []

    def spy(out_dir, split, half, *, unsealed=False):
        asked.append((half, unsealed))
        return world.counts

    monkeypatch.setattr(explore, "load_country_counts", spy)
    monkeypatch.setattr(explore, "resolution_of", lambda path, stated=None: 7)
    monkeypatch.setattr(explore, "frozen_list_codes", lambda resolution, directory: world.list_codes)
    monkeypatch.setattr(explore, "load_grid", lambda gpkg, land: world.grid)
    monkeypatch.setattr(explore, "load_references", lambda grid, directory: world.refs)
    monkeypatch.setattr(explore, "load_split", lambda path: pl.DataFrame({"species_id": ["x"], "half": ["A"]}))
    spec = explore.CleanSpec(
        counts=tmp_path / "counts", reference=tmp_path / "reference", split=tmp_path / "split.parquet", frozen_dir=tmp_path / "frozen", out=tmp_path / "out",
        settings=SETTINGS, layers=explore.LAYERS, inputs=explore.clean_inputs(SETTINGS, explore.LAYERS, world.list_codes, world.marker), force=False,
    )
    explore.clean_job(spec)
    assert asked == [("A", False)]
    # a second job finds its outputs and loads nothing at all
    explore.clean_job(spec)
    assert asked == [("A", False)]
    # and nothing in the module asks for the other half
    source = inspect.getsource(explore)
    assert "unsealed" not in source and "score_half_b" not in source and '"B"' not in source and '"all"' not in source


def test_the_tracked_configs_are_the_default_options_and_the_screen_figures(tmp_path):
    configs, params = explore.load_configs(Path(__file__).parent.parent / "explore" / "configs.toml")
    assert configs == {"default": InfomapOptions()}
    assert (params.trials_screen, params.n_null_screen) == (3, 19)
    assert (params.trials_select, params.n_boot, params.boot_trials, params.n_null_select) == (10, 20, 3, 99)

    path = tmp_path / "configs.toml"
    path.write_text('trials_screen = 2\n[configs.fine]\nmarkov_time = 0.8\n[configs.recursive]\nrecurse_share = 0.5\n')
    configs, params = explore.load_configs(path)
    assert configs == {"fine": InfomapOptions(markov_time=0.8), "recursive": InfomapOptions(recurse_share=0.5)} and params.trials_screen == 2
    for text, message in (
        ("[configs.x]\nmarkov = 1\n", "markov"),
        ("[configs.x]\nregularization_strength = 2.0\n", "regularization_strength"),
        ("trials = 3\n[configs.x]\n", "trials"),
        ("trials_screen = 3\n", "no configuration"),
    ):
        path.write_text(text)
        with pytest.raises(ValueError, match=message):
            explore.load_configs(path)


def test_the_run_is_named_for_the_mode_and_the_rules_switched_off():
    assert explore.run_name("native_only", []) == "native_only"
    assert explore.run_name("introduced", ["wrims", "griis"]) == "introduced-no-griis-no-wrims"
    assert explore.run_name("off", ["griis"]) == "off-no-griis"
    with pytest.raises(ValueError, match="nonsense"):
        explore.run_name("native_only", ["nonsense"])
    assert set(SWITCHABLE) >= {"griis"}


def test_workers_that_fit_the_budget():
    assert explore.fit_workers(12, parent_gb=5.0, worker_gb=3.0, budget_gb=40.0) == 9  # 35 GB left over 3 GB * 1.25
    assert explore.fit_workers(4, parent_gb=5.0, worker_gb=3.0, budget_gb=40.0) == 4
    assert explore.fit_workers(12, parent_gb=39.0, worker_gb=3.0, budget_gb=40.0) == 1  # at least one
    assert explore.fit_workers(1, parent_gb=1.0, worker_gb=1.0, budget_gb=40.0) == 1


def test_a_second_worker_pays_its_way_in_the_log(world, planted):
    out, _ = planted
    log = json.loads((out / "default" / "log.json").read_text())
    assert log["select"]["workers_requested"] == 1 and log["select"]["workers"] == []
    units = explore.read_units(out)
    ct = explore.network_matrix(out)
    explore.screen_stage(out, "default", units, ct, explore.read_light_grid(out), InfomapOptions(), FAST, explore.read_stamp(out / "units")["hash"], workers=2, force=True)
    log = json.loads((out / "default" / "log.json").read_text())
    assert log["screen"]["workers_requested"] == 2 and 1 <= log["screen"]["workers_used"] <= 2
    assert len(log["screen"]["workers"]) == 3 and all(w["peak_gb"] > 0 and w["pid"] for w in log["screen"]["workers"])  # a task for each of k = 2, 3, 4


def test_units_round_trip_through_their_files(planted):
    out, _ = planted
    units = explore.read_units(out)
    assert isinstance(units, Units) and units.settings == UNITS
    assert list(units.areas) == sorted(units.areas) and all(units.adjacency[u] <= set(units.areas) for u in units.areas)
    assert all(u in units.adjacency for u in units.areas)
    assert explore.read_light_grid(out).cells.height == 123


def test_the_pipeline_runs_the_stages_up_to_the_one_asked_for_and_skips_the_current(world, tmp_path, monkeypatch):
    loaded = []
    monkeypatch.setattr(explore, "load_country_counts", lambda out_dir, split, half, **kw: loaded.append(half) or world.counts)
    monkeypatch.setattr(explore, "resolution_of", lambda path, stated=None: 7)
    monkeypatch.setattr(explore, "read_complete", lambda path: world.marker)
    monkeypatch.setattr(explore, "frozen_list_codes", lambda resolution, directory: world.list_codes)
    monkeypatch.setattr(explore, "load_grid", lambda gpkg, land: world.coarse if "ISEA3H-5" in str(gpkg) else world.grid)
    monkeypatch.setattr(explore, "_manifest_path", lambda reference, item: world.coarse_gpkg if item.endswith("5") else Path(f"{item}.gpkg"))
    monkeypatch.setattr(explore, "part_regions", lambda grid, path: world.regions)
    monkeypatch.setattr(explore, "load_references", lambda grid, directory: world.refs)
    monkeypatch.setattr(explore, "load_split", lambda path: pl.DataFrame({"species_id": ["x"], "half": ["A"]}))
    monkeypatch.setattr(explore, "load_configs", lambda path: ({"default": InfomapOptions(), "fine": InfomapOptions(markov_time=0.9)}, FAST))
    plan = explore.Plan(
        counts=tmp_path / "counts", reference=tmp_path / "reference", split=tmp_path / "split", frozen_dir=tmp_path / "frozen", results=tmp_path / "results",
        configs=tmp_path / "configs.toml", mode="native_only", rules_off=tuple(LAYER_5), stage="units", isolate=False,
    )
    first = explore.run_pipeline(replace(plan, rules_off=()), UNITS)
    out = first["directory"]
    assert out == tmp_path / "results" / "res7" / "native_only" and not first["clean"]["skipped"] and not first["units"]["skipped"]
    assert not (out / "default").exists()  # the stage given is the last that runs
    log = json.loads((out / "log.json").read_text())
    assert log["clean"]["child_peak_rss_gb"] > 0 and log["units"]["parts_excluded"] >= 0
    # a later stage runs the earlier ones from their stamps, and --only picks the configuration
    second = explore.run_pipeline(replace(plan, rules_off=(), stage="select", only=("fine",)), UNITS)
    assert second["clean"]["skipped"] and second["units"]["skipped"] and loaded == ["A"]
    assert set(second) == {"directory", "clean", "units", "fine"} and not second["fine"]["screen"]["skipped"] and not second["fine"]["select"]["skipped"]
    assert (out / "fine" / "selected.json").exists() and not (out / "default").exists()
    third = explore.run_pipeline(replace(plan, rules_off=(), stage="select"), UNITS)
    assert third["fine"] == {"screen": {"skipped": True}, "select": {"skipped": True}} and not third["default"]["screen"]["skipped"]
    with pytest.raises(ValueError, match="nonsense"):
        explore.run_pipeline(replace(plan, only=("nonsense",)), UNITS)
    # the rules switched off make a run of their own
    other = explore.run_pipeline(plan, UNITS)
    assert other["directory"].name == "native_only-no-detached_satellite-no-thin_long_range-no-under_reported" and not other["clean"]["skipped"]
