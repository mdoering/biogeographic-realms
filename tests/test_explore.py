"""The exploration pipeline on a small planted world: four realms of 30 cells, each with its own species, and three cells too thin to enter the network.

The world has the cells' land parts only (the sea parts have no area). It goes through the real `clean`, `network_units`, Infomap, `candidate_maps`, the bootstrap and the null
maps, with fewer trials and nulls than the analysis, so that a run takes seconds.
"""

import inspect
import json
import shutil
from dataclasses import asdict, dataclass, replace
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
    chosen = json.loads((out / "default" / "select" / "selected.json").read_text())
    assert chosen["k"] == 4 and not chosen["fell_short"]

    candidates = pl.read_parquet(out / "default" / "select" / "candidates.parquet")
    assert set(candidates.columns) == {"unit", "k", "realm", "sub_realm", "sub_realm_of"}
    found = dict(candidates.filter(pl.col("k") == 4).select("unit", "realm").iter_rows())
    assert len(found) == 120 and v_measure(found, {p: world.truth[p] for p in found}) >= 0.95

    # the map by part: the selected k is marked, and the thin parts have the realm of their neighbours
    parts = pl.read_parquet(out / "default" / "select" / "parts.parquet")
    picked = parts.filter(pl.col("selected"))
    assert picked["k"].unique().to_list() == [4] and picked.height == 123
    assert v_measure(dict(picked.select("part", "realm").iter_rows()), {**world.truth, **{f"{c}:land": 3 for c in THIN_CELLS}}) >= 0.95

    selection = pl.read_csv(out / "default" / "select" / "selection.csv")
    assert {"k", "stability", "excess"} <= set(selection.columns)
    assert selection.filter(pl.col("k") == 4)["stability"][0] >= 0.8
    screen = pl.read_csv(out / "default" / "screen.csv")
    assert screen["k"].to_list() == [2, 3, 4] and screen["excess"].is_finite().all()

    structure = json.loads((out / "default" / "select" / "structure.json").read_text())
    assert structure["degenerate"] is False and structure["modules"] == 4 and structure["single_unit_modules"] == 0
    assert 0.25 <= structure["largest_module_share_of_units"] <= 0.3 and structure["levels"] >= 1
    partition = pl.read_parquet(out / "default" / "select" / "partition.parquet")
    assert partition.columns == ["unit", "leaf", "top"] and partition.height == 120

    for name in ("presence.parquet", "summary.csv", "major.parquet", "taxa.parquet"):
        assert (out / "cleaned" / name).exists()
    for name in ("members", "presence", "areas", "adjacency", "excluded", "stats"):
        assert (out / "units" / f"{name}.parquet").exists()
    assert (out / "units" / "removed_by_region.csv").exists()
    log = json.loads((out / "default" / "select" / "log.json").read_text())
    assert set(log) == {"select"} and log["select"]["seconds"] > 0 and log["select"]["peak_rss_gb"] > 0
    assert set(json.loads((out / "default" / "log.json").read_text())) == {"screen"}
    # the screen's own files are in the configuration's directory, the select's in select/, and neither holds the other's
    assert json.loads((out / "default" / "structure.json").read_text())["stage"] == "screen" and structure.get("stage") == "select"
    assert {p.name for p in (out / "default").iterdir()} == {*explore.SCREEN_FILES, "log.json", "stage-screen.json", "select"}
    assert {p.name for p in (out / "default" / "select").iterdir()} == {*explore.SELECT_FILES, "log.json", "stage-select.json"}
    assert chosen["units"] == 120  # the units scored, all of them assigned here


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
    # and the stages after it are no longer current, although their hashes are the same: cleaning deleted their stamps, so they run again
    assert not (out / "units" / "stage.json").exists() and not (out / "default" / "stage-screen.json").exists()
    rerun = run_chain(out, world, stage="screen")
    assert rerun["clean"]["skipped"] and not rerun["units"]["skipped"] and not rerun["screen"]["skipped"]
    assert rerun["units"]["hash"] == first["units"]["hash"] and rerun["screen"]["hash"] == first["screen"]["hash"]
    # a screen that is run again leaves the select result alone, which has a directory of its own
    run_chain(out, world, stage="select")
    selected = (out / "default" / "select" / "selected.json").read_bytes()
    explore.screen_stage(out, "default", explore.read_units(out), explore.network_matrix(out), explore.read_light_grid(out), InfomapOptions(), FAST, explore.read_stamp(out / "units")["hash"], force=True)
    assert (out / "default" / "select" / "selected.json").read_bytes() == selected
    assert all(v["skipped"] for v in run_chain(out, world, stage="select").values())


def test_a_degenerate_configuration_is_recorded_not_raised(world, tmp_path):
    out = tmp_path / "res7" / "native_only"
    done = run_chain(out, world, replace(FAST, floor_km2=1e12))  # no realm is as large as that, so they all merge into one
    structure = json.loads((out / "default" / "select" / "structure.json").read_text())
    assert structure["degenerate"] is True and "above the size floor" in structure["degenerate_reason"]
    assert structure["modules"] == 4
    chosen = json.loads((out / "default" / "select" / "selected.json").read_text())
    assert chosen["degenerate"] is True and "k" not in chosen
    assert pl.read_parquet(out / "default" / "select" / "candidates.parquet").height == 0
    assert pl.read_csv(out / "default" / "screen.csv").height == 0 and pl.read_csv(out / "default" / "select" / "selection.csv").height == 0
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
    candidates = pl.read_parquet(out / "default" / "select" / "candidates.parquet")
    assert set(candidates["unit"]) == set(units.areas)
    assert set(pl.read_parquet(out / "default" / "select" / "partition.parquet")["unit"]) <= set(units.areas)
    seen = []
    real = explore.excess_table
    monkeypatch.setattr(explore, "excess_table", lambda maps, ct, adjacency, areas, **kw: seen.append((maps, adjacency, areas)) or real(maps, ct, adjacency, areas, **kw))
    explore.select_stage(out, "default", units, ct, explore.read_light_grid(out), InfomapOptions(), FAST, explore.read_stamp(out / "units")["hash"], force=True)
    maps, adjacency, areas = seen[0]
    assert all(set(labels) == set(units.areas) for labels in maps.values()) and set(areas) == set(units.areas) and set(adjacency) == set(units.areas)
    # in the outputs only, they have a realm, and the parts file says that it was imputed
    parts = pl.read_parquet(out / "default" / "select" / "parts.parquet").filter(pl.col("k") == 4)
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


def unlinked_world(out: Path, n: int):
    """The planted units, with the first `n` of them holding only a species that is not accepted: they stay rows of the matrix without a link."""
    units = explore.read_units(out)
    lost = sorted(units.areas)[:n]
    presence = units.presence.with_columns(
        pl.when(pl.col("unit").is_in(lost)).then(pl.lit("zz-not-accepted")).otherwise(pl.col("species_id")).alias("species_id")
    ).unique(["unit", "species_id"])
    units = replace(units, presence=presence)
    ct = explore.network_matrix(out, units)
    assert set(lost) <= set(ct.cells)
    return units, ct, lost


@pytest.mark.parametrize("n, warns", [(1, False), (3, True)])
def test_units_without_a_module_leave_the_scores_but_not_the_map(world, planted, monkeypatch, tmp_path, n, warns):
    out = tmp_path / "res7" / "native_only"
    shutil.copytree(planted[0], out)
    units, ct, lost = unlinked_world(out, n)
    kept = set(units.areas) - set(lost)
    seen = {"boot": [], "excess": []}
    real_boot, real_excess = explore.bootstrap_stability, explore.excess_table
    monkeypatch.setattr(
        explore, "bootstrap_stability", lambda ct, maps, areas, adjacency, **kw: seen["boot"].append((ct, maps, areas, adjacency)) or real_boot(ct, maps, areas, adjacency, **kw)
    )
    monkeypatch.setattr(
        explore, "excess_table", lambda maps, ct, adjacency, areas, **kw: seen["excess"].append((ct, maps, areas, adjacency)) or real_excess(maps, ct, adjacency, areas, **kw)
    )
    light = explore.read_light_grid(out)
    parent = explore.read_stamp(out / "units")["hash"]
    explore.screen_stage(out, "default", units, ct, light, InfomapOptions(), FAST, parent, force=True)
    explore.select_stage(out, "default", units, ct, light, InfomapOptions(), FAST, parent, force=True)
    assert len(seen["excess"]) == 2 and len(seen["boot"]) == 1  # screen and select score once each, and only select bootstraps
    for scored_ct, maps, areas, adjacency in seen["excess"] + seen["boot"]:
        assert set(scored_ct.cells) == kept and set(areas) == kept and set(adjacency) == kept
        assert all(near <= kept for near in adjacency.values()) and all(set(labels) == kept for labels in maps.values())
        assert not set(lost) & {c for c in scored_ct.cells}
    assert len(units.areas) == 120  # the units that were passed in are untouched

    structure = json.loads((out / "default" / "select" / "structure.json").read_text())
    assert structure["unassigned"] == n and structure["unassigned_units"] == lost and structure["units"] == 120
    assert ("warning" in structure) is warns
    if warns:
        assert "more than 1%" in structure["warning"] and str(n) in structure["warning"]
    # the map has no label for them, and the parts file imputes them, for every k
    partition = pl.read_parquet(out / "default" / "select" / "partition.parquet")
    assert partition.height == 120 and partition.filter(pl.col("unit").is_in(lost))["leaf"].is_null().all() and partition["leaf"].null_count() == n
    assert set(pl.read_parquet(out / "default" / "select" / "candidates.parquet")["unit"]) == kept
    parts = pl.read_parquet(out / "default" / "select" / "parts.parquet")
    for k, rows in parts.group_by("k"):
        gone = rows.filter(pl.col("part").is_in(lost))
        assert sorted(gone["part"]) == lost and gone["imputed"].all() and gone["unit"].is_null().all() and gone["realm"].is_not_null().all()
        assert rows.height == 123
    assert json.loads((out / "default" / "select" / "selected.json").read_text())["k"] == 4


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
    assert list(configs) == ["default", "mt095", "mt090", "mt085", "mt080", "mt075", "variable", "regularized", "entropy", "recurse50"]
    assert configs == {
        "default": InfomapOptions(),
        "mt095": InfomapOptions(markov_time=0.95),
        "mt090": InfomapOptions(markov_time=0.90),
        "mt085": InfomapOptions(markov_time=0.85),
        "mt080": InfomapOptions(markov_time=0.80),
        "mt075": InfomapOptions(markov_time=0.75),
        "variable": InfomapOptions(variable_markov_time=True),
        "regularized": InfomapOptions(regularized=True),
        "entropy": InfomapOptions(entropy_corrected=True),
        "recurse50": InfomapOptions(recurse_share=0.5),
    }
    assert configs["default"].args() == "" and configs["variable"].variable_markov_damping is None  # Infomap's own damping, and nothing else set
    assert all(configs[n].markov_time == 1.0 for n in ("variable", "regularized", "entropy", "recurse50"))
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
    log = json.loads((out / "default" / "select" / "log.json").read_text())
    assert log["select"]["workers_requested"] == 1 and log["select"]["workers"] == []
    one_worker = (out / "default" / "screen.csv").read_bytes()  # the planted screen ran with one worker
    assert json.loads((out / "default" / "log.json").read_text())["screen"]["workers_requested"] == 1
    units = explore.read_units(out)
    ct = explore.network_matrix(out)
    explore.screen_stage(out, "default", units, ct, explore.read_light_grid(out), InfomapOptions(), FAST, explore.read_stamp(out / "units")["hash"], workers=2, force=True)
    log = json.loads((out / "default" / "log.json").read_text())
    assert log["screen"]["workers_requested"] == 2 and 1 <= log["screen"]["workers_used"] <= 2
    assert len(log["screen"]["workers"]) == 3 and all(w["peak_gb"] > 0 and w["pid"] for w in log["screen"]["workers"])  # a task for each of k = 2, 3, 4
    # the two pools (two workers to measure, then as many as fit) score what one worker does
    assert (out / "default" / "screen.csv").read_bytes() == one_worker


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
    monkeypatch.setattr(explore, "manifest_path", lambda reference, item: world.coarse_gpkg if item.endswith("5") else Path(f"{item}.gpkg"))
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
    assert (out / "fine" / "select" / "selected.json").exists() and not (out / "default").exists()
    third = explore.run_pipeline(replace(plan, rules_off=(), stage="select"), UNITS)
    assert third["fine"] == {"screen": {"skipped": True}, "select": {"skipped": True}} and not third["default"]["screen"]["skipped"]
    with pytest.raises(ValueError, match="nonsense"):
        explore.run_pipeline(replace(plan, only=("nonsense",)), UNITS)
    # the rules switched off make a run of their own
    other = explore.run_pipeline(plan, UNITS)
    assert other["directory"].name == "native_only-no-detached_satellite-no-thin_long_range-no-under_reported" and not other["clean"]["skipped"]


# ------------------------------------------------------------------ duplicate genus and family names merged as a run option

CODES = pl.DataFrame(schema={"country": pl.String, "subdivision": pl.String, "code": pl.String})
# The hashes of these inputs, as they were before the option existed and as every run of M1 has them. Change them consciously: a new value reruns every stage of
# every run, so a change of the hashed inputs that is not meant to change the results (a new option off by default) must keep them.
# The cleaning hash changed on 2026-11-05 (was 4a21d1f2…): layer 4 takes two thirds of the records and the land classes (DECISIONS, "Terrestrial unless COL says marine"),
# and WoRMS's environment where COL has none ("WoRMS's environment").
DEFAULT_CLEAN_HASH = "0ea97c1fb4e3d104125f8b39cc5be4e18924b9efd428a36651fcdcb229ca9cfa"
DEFAULT_SELECT_HASH = "92587114b47877aecdd749bb5dd9bb72b26a2a45b77d8ea082f22fec3ea75e96"


def test_merging_duplicate_names_enters_the_hashes_only_when_it_is_on():
    settings, marker = Settings(wcvp_mode="native_only"), {"rows": 1, "resolution": 7}
    plain = explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m")
    assert "merge_duplicate_names" not in plain and explore.digest(plain) == DEFAULT_CLEAN_HASH
    assert explore.digest(explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m", merge_duplicate_names=False)) == DEFAULT_CLEAN_HASH
    merged = explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m", merge_duplicate_names=True)
    assert merged["merge_duplicate_names"] is True and explore.digest(merged) != DEFAULT_CLEAN_HASH

    options = InfomapOptions(markov_time=0.8)
    assert explore.digest(explore.config_inputs("parent", options, Params(), "select")) == DEFAULT_SELECT_HASH
    assert explore.digest(explore.config_inputs("parent", options, Params(), "select", merge_duplicate_names=False)) == DEFAULT_SELECT_HASH
    assert explore.digest(explore.config_inputs("parent", options, Params(), "select", merge_duplicate_names=True)) != DEFAULT_SELECT_HASH


def test_a_run_that_merges_duplicate_names_is_named_for_it():
    assert explore.run_name("native_only", [], merge_duplicate_names=True) == "native_only-merged-names"
    assert explore.run_name("native_only", ["griis"], merge_duplicate_names=True) == "native_only-no-griis-merged-names"
    assert explore.run_name("native_only", ["griis"], merge_duplicate_names=False) == "native_only-no-griis"


def named(world: World) -> World:
    """The world, with the names of each species' genus, family, order and kingdom: genus ids G<r>-0 and G<r>-1 of each realm share one name in one family, and
    so do G<r>-2 and G<r>-3, which merging takes for duplicates."""
    taxa = world.refs.taxa.with_columns(
        pl.col("genus_id").str.replace(r"-(\d+)$", "").alias("_stem"),
        (pl.col("genus_id").str.extract(r"-(\d+)$").cast(pl.Int64) // 2).cast(pl.String).alias("_pair"),
    ).with_columns(
        (pl.lit("Genus") + pl.col("_stem") + pl.lit("-") + pl.col("_pair")).alias("genus"),
        (pl.lit("Family") + pl.col("family_id")).alias("family"),
        pl.lit("Coleoptera").alias("order"),
    ).drop("_stem", "_pair")
    return replace(world, refs=replace(world.refs, taxa=taxa))


def patch_loaders(monkeypatch, world: World) -> None:
    monkeypatch.setattr(explore, "load_country_counts", lambda out_dir, split, half, **kw: world.counts)
    monkeypatch.setattr(explore, "resolution_of", lambda path, stated=None: 7)
    monkeypatch.setattr(explore, "read_complete", lambda path: world.marker)
    monkeypatch.setattr(explore, "frozen_list_codes", lambda resolution, directory: world.list_codes)
    monkeypatch.setattr(explore, "load_grid", lambda gpkg, land: world.coarse if "ISEA3H-5" in str(gpkg) else world.grid)
    monkeypatch.setattr(explore, "manifest_path", lambda reference, item: world.coarse_gpkg if item.endswith("5") else Path(f"{item}.gpkg"))
    monkeypatch.setattr(explore, "part_regions", lambda grid, path: world.regions)
    monkeypatch.setattr(explore, "load_references", lambda grid, directory: world.refs)
    monkeypatch.setattr(explore, "load_split", lambda path: pl.DataFrame({"species_id": ["x"], "half": ["A"]}))
    monkeypatch.setattr(explore, "load_configs", lambda path: ({"default": InfomapOptions()}, FAST))


def test_merge_duplicate_names_reaches_the_matrix_and_makes_a_run_of_its_own(world, tmp_path, monkeypatch):
    world = named(world)
    patch_loaders(monkeypatch, world)
    asked = []
    real = explore.cell_taxa
    monkeypatch.setattr(explore, "cell_taxa", lambda presence, taxa, *args, **kw: asked.append(kw.get("merge_duplicate_names", False)) or real(presence, taxa, *args, **kw))
    plan = explore.Plan(
        counts=tmp_path / "counts", reference=tmp_path / "reference", split=tmp_path / "split", frozen_dir=tmp_path / "frozen", results=tmp_path / "results",
        configs=tmp_path / "configs.toml", mode="native_only", rules_off=tuple(LAYER_5), stage="select", isolate=False,
    )
    plain = explore.run_pipeline(plan, UNITS)
    assert asked == [False]  # the screen and the select of a run share one matrix
    merged = explore.run_pipeline(replace(plan, merge_duplicate_names=True), UNITS)
    assert asked == [False, True]
    assert merged["directory"].name == plain["directory"].name + "-merged-names"
    for stage in ("clean", "units"):
        assert not merged[stage]["skipped"] and merged[stage]["hash"] != plain[stage]["hash"]
    for stage in ("screen", "select"):
        assert not merged["default"][stage]["skipped"] and merged["default"][stage]["hash"] != plain["default"][stage]["hash"]
    stamp = explore.read_stamp(merged["directory"] / "default" / "select", "stage-select")["inputs"]
    assert stamp["merge_duplicate_names"] is True and "merge_duplicate_names" not in explore.read_stamp(plain["directory"] / "default" / "select", "stage-select")["inputs"]

    # the cleaned taxa carry the names that merging needs only in the run that merges, and the presence is the same
    assert pl.read_parquet(plain["directory"] / "cleaned" / "taxa.parquet").columns == ["species_id", "genus_id", "family_id"]
    names = pl.read_parquet(merged["directory"] / "cleaned" / "taxa.parquet")
    assert {"kingdom", "order", "family", "genus"} <= set(names.columns)
    assert pl.read_parquet(merged["directory"] / "cleaned" / "presence.parquet").equals(pl.read_parquet(plain["directory"] / "cleaned" / "presence.parquet"))
    genera = {d: sum(t.startswith("genus:") for t in explore.network_matrix(d, merge_duplicate_names=m).taxa) for d, m in ((plain["directory"], False), (merged["directory"], True))}
    assert genera[merged["directory"]] < genera[plain["directory"]]
    # a matrix that merges needs the names, and says so
    with pytest.raises(ValueError, match="--merge-duplicate-names"):
        explore.network_matrix(plain["directory"], merge_duplicate_names=True)

    # the run on ids is untouched: run again, every stage is current
    again = explore.run_pipeline(plan, UNITS)
    assert again["clean"]["skipped"] and again["units"]["skipped"] and again["default"] == {"screen": {"skipped": True}, "select": {"skipped": True}}


# ------------------------------------------------------------------ method versions, a rewritten stage's dependants, the seed, the half, reading a result back


def test_the_method_versions_are_one_and_enter_the_hashes_only_above_it(monkeypatch):
    assert explore.METHOD_VERSION == {"clean": 1, "units": 1, "screen": 1, "select": 1}  # raise one when its stage's behaviour changes
    settings, marker = Settings(wcvp_mode="native_only"), {"rows": 1, "resolution": 7}
    options = InfomapOptions(markov_time=0.8)
    units_before = explore.digest({"cleaned": "parent", "unit_settings": asdict(UnitSettings())})  # the units' inputs as M1 hashed them
    assert explore.digest(explore.units_inputs("parent", UnitSettings())) == units_before
    hashes = {
        "clean": lambda: explore.digest(explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m")),
        "units": lambda: explore.digest(explore.units_inputs("parent", UnitSettings())),
        "screen": lambda: explore.digest(explore.config_inputs("parent", options, Params(), "screen")),
        "select": lambda: explore.digest(explore.config_inputs("parent", options, Params(), "select")),
    }
    today = {stage: hash_() for stage, hash_ in hashes.items()}
    assert today["clean"] == DEFAULT_CLEAN_HASH and today["select"] == DEFAULT_SELECT_HASH
    for stage in explore.METHOD_VERSION:
        monkeypatch.setitem(explore.METHOD_VERSION, stage, 2)
        assert {s: h() != today[s] for s, h in hashes.items()} == {s: s == stage for s in hashes}  # that stage's hash only; the stages after follow through it
        monkeypatch.setitem(explore.METHOD_VERSION, stage, 1)


def test_a_stage_written_again_deletes_the_stamps_of_the_stages_that_depend_on_it(world, planted, tmp_path):
    out = tmp_path / "res7" / "native_only"
    shutil.copytree(planted[0], out)
    (out / "fine-s2" / "select").mkdir(parents=True)
    for stamp in ("fine-s2/stage-screen.json", "fine-s2/select/stage-select.json"):
        (out / stamp).write_text("{}")
    screens = sorted(out.glob("*/stage-screen.json"))
    selects = sorted(out.glob("*/select/stage-select.json"))
    assert len(screens) == len(selects) == 2
    # a screen or a select takes nothing with it
    assert explore.invalidate_after(out, "screen") == [] and explore.invalidate_after(out, "select") == []
    # the units take every configuration's screen and select, and leave the cleaning's stamp
    parent = explore.read_stamp(out / "cleaned")["hash"]
    explore.units_stage(out, world.grid, world.coarse, world.coarse_gpkg, UNITS, world.regions, parent, force=True)
    assert not any(p.exists() for p in screens + selects) and (out / "cleaned" / "stage.json").exists() and (out / "units" / "stage.json").exists()
    again = run_chain(out, world)
    assert again["clean"]["skipped"] and again["units"]["skipped"] and not again["screen"]["skipped"] and not again["select"]["skipped"]
    # cleaning takes the units too
    inputs = explore.clean_inputs(SETTINGS, explore.LAYERS, world.list_codes, world.marker)
    explore.clean_stage(out, world.counts, world.refs, world.grid, SETTINGS, explore.LAYERS, world.list_codes, inputs, force=True)
    assert not (out / "units" / "stage.json").exists() and not (out / "default" / "stage-screen.json").exists()
    assert not (out / "default" / "select" / "stage-select.json").exists()


def test_a_second_seed_writes_a_directory_of_its_own_and_leaves_seed_one_as_it_was(world, tmp_path, monkeypatch):
    load_configs = explore.load_configs
    patch_loaders(monkeypatch, world)
    budgets = []
    real = explore.measured_workers
    monkeypatch.setattr(explore, "measured_workers", lambda requested, seen, used, budget_gb=explore.BUDGET_GB: budgets.append(budget_gb) or real(requested, seen, used, budget_gb))
    plan = explore.Plan(
        counts=tmp_path / "counts", reference=tmp_path / "reference", split=tmp_path / "split", frozen_dir=tmp_path / "frozen", results=tmp_path / "results",
        configs=tmp_path / "configs.toml", mode="native_only", rules_off=tuple(LAYER_5), stage="select", isolate=False,
    )
    first = explore.run_pipeline(plan, UNITS)
    out = first["directory"]
    assert budgets and set(budgets) == {explore.BUDGET_GB}
    seed_one = {p: p.read_bytes() for p in sorted((out / "default").rglob("*")) if p.is_file() and p.name != "log.json"}
    stamp = explore.read_stamp(out / "default" / "select", "stage-select")
    assert stamp["inputs"]["params"]["seed"] == 1 and stamp["hash"] == explore.digest(explore.config_inputs(first["units"]["hash"], InfomapOptions(), FAST, "select"))

    budgets.clear()
    second = explore.run_pipeline(replace(plan, seed=2, budget_gb=12.5), UNITS)
    assert second["clean"]["skipped"] and second["units"]["skipped"] and set(budgets) == {12.5}
    assert not second["default"]["screen"]["skipped"] and not second["default"]["select"]["skipped"]
    assert second["default"]["select"]["hash"] != first["default"]["select"]["hash"]
    assert explore.read_stamp(out / "default-s2" / "select", "stage-select")["inputs"]["params"]["seed"] == 2
    assert {p: p.read_bytes() for p in seed_one} == seed_one  # seed 1's files are untouched
    third = explore.run_pipeline(plan, UNITS)
    assert third["default"] == {"screen": {"skipped": True}, "select": {"skipped": True}}

    assert explore.config_dir("mt080") == "mt080" and explore.config_dir("mt080", 1) == "mt080" and explore.config_dir("mt080", 2) == "mt080-s2"
    with pytest.raises(ValueError, match="seed"):
        replace(plan, seed=0)
    path = tmp_path / "clash.toml"
    for name in ("mt080-s2", "select"):
        path.write_text(f'[configs."{name}"]\n')
        with pytest.raises(ValueError, match="like a seed"):
            load_configs(path)


def test_the_half_is_a_only_and_another_would_never_take_its_stamps_or_directory(tmp_path):
    with pytest.raises(ValueError, match="sealed until M3"):
        explore.Plan(
            counts=tmp_path, reference=tmp_path, split=tmp_path, frozen_dir=tmp_path, results=tmp_path, configs=tmp_path, mode="native_only", rules_off=(),
            stage="clean", half="B",
        )
    # the name and the hash, at the function level: half A adds nothing to either, any other half changes both
    assert explore.run_name("native_only", ["griis"]) == explore.run_name("native_only", ["griis"], half="A") == "native_only-no-griis"
    assert explore.run_name("native_only", ["griis"], half="B") == "native_only-no-griis-half-B"
    settings, marker = Settings(wcvp_mode="native_only"), {"rows": 1, "resolution": 7}
    a = explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m", half="A", split_sha256="ignored for A")
    assert explore.digest(a) == DEFAULT_CLEAN_HASH and "half" not in a
    b = explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m", half="B", split_sha256="f" * 64)
    assert b["half"] == "B" and b["split"] == "f" * 64 and explore.digest(b) != DEFAULT_CLEAN_HASH
    assert explore.digest(explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m", half="B", split_sha256="e" * 64)) != explore.digest(b)
    with pytest.raises(ValueError, match="sha256"):
        explore.clean_inputs(settings, explore.LAYERS, CODES, marker, "m", half="B")
    # a stamp of half A is not current for half B's inputs, even in the same directory
    directory = tmp_path / "cleaned"
    directory.mkdir()
    for name in explore.CLEAN_FILES:
        (directory / name).write_text("")
    explore.write_stamp(directory, "stage", explore.digest(a), a, {})
    assert explore.is_current(directory, "stage", explore.digest(a), explore.CLEAN_FILES)
    assert not explore.is_current(directory, "stage", explore.digest(b), explore.CLEAN_FILES)
    # and the units, screen and select of B hash B's cleaning, so none of them can be A's either
    assert explore.digest(explore.units_inputs(explore.digest(b), UnitSettings())) != explore.digest(explore.units_inputs(explore.digest(a), UnitSettings()))


def test_a_select_result_is_read_back_only_when_it_is_current(world, planted, tmp_path, monkeypatch):
    out = tmp_path / "res7" / "native_only"
    shutil.copytree(planted[0], out)
    monkeypatch.setattr(explore, "load_configs", lambda path: ({"default": InfomapOptions(), "fine": InfomapOptions(markov_time=0.9)}, FAST))
    found = explore.read_selected(out, "default", tmp_path / "configs.toml")
    assert found.directory == out / "default" / "select" and found.k == 4 and not found.fell_short and not found.merged
    assert found.parts.equals(pl.read_parquet(out / "default" / "select" / "parts.parquet")) and found.selection["k"].to_list() == [2, 3, 4]
    assert found.partition.unassigned == [] and len(found.partition.paths) == 120 and found.candidates.height > 0
    for kwargs, message in (
        ({"seed": 2}, "no select stamp"),  # no result at that seed
        ({"merge_duplicate_names": True}, "not the current"),  # another hash
    ):
        with pytest.raises(ValueError, match=message):
            explore.read_selected(out, "default", tmp_path / "configs.toml", **kwargs)
    with pytest.raises(ValueError, match="no configuration"):
        explore.read_selected(out, "nonsense", tmp_path / "configs.toml")
    with pytest.raises(ValueError, match="no select stamp"):
        explore.read_selected(out, "fine", tmp_path / "configs.toml")
    # other parameters than those of the result: not current
    monkeypatch.setattr(explore, "load_configs", lambda path: ({"default": InfomapOptions()}, replace(FAST, n_null_select=19)))
    with pytest.raises(ValueError, match="not the current"):
        explore.read_selected(out, "default", tmp_path / "configs.toml")
    monkeypatch.setattr(explore, "load_configs", lambda path: ({"default": InfomapOptions()}, FAST))
    # the units of another cleaning
    stamp = json.loads((out / "units" / "stage.json").read_text())
    (out / "units" / "stage.json").write_text(json.dumps(stamp | {"inputs": stamp["inputs"] | {"cleaned": "other"}}))
    with pytest.raises(ValueError, match="cleaning other"):
        explore.read_selected(out, "default", tmp_path / "configs.toml")
    (out / "units" / "stage.json").write_text(json.dumps(stamp))
    # the units written again: the select stamp is gone with them
    explore.units_stage(out, world.grid, world.coarse, world.coarse_gpkg, UNITS, world.regions, explore.read_stamp(out / "cleaned")["hash"], force=True)
    with pytest.raises(ValueError, match="no select stamp"):
        explore.read_selected(out, "default", tmp_path / "configs.toml")


def old_layout(planted: Path, out: Path) -> None:
    """A copy of the planted run as the pipeline wrote it before select had a directory of its own: select's files over the screen's, both stamps, one log."""
    shutil.copytree(planted, out)
    config = out / "default"
    logs = json.loads((config / "log.json").read_text()) | json.loads((config / "select" / "log.json").read_text())
    for path in (config / "select").iterdir():
        path.replace(config / path.name)
    (config / "select").rmdir()
    (config / "log.json").write_text(json.dumps(logs))
    shutil.copytree(config, out / "entropy")  # a configuration that only had a screen
    (out / "entropy" / "stage-select.json").unlink()


def test_migrate_select_moves_the_select_result_and_drops_the_overwritten_screen_stamp(world, planted, tmp_path, monkeypatch):
    root = tmp_path / "explore"
    out = root / "res7" / "native_only"
    old_layout(planted[0], out)
    selected = (out / "default" / "selected.json").read_bytes()
    entropy = {p.name: p.read_bytes() for p in (out / "entropy").iterdir()}
    done = explore.migrate_select(root)
    assert [d["directory"] for d in done] == [str(out / "default")] and done[0]["screen_stamp_removed"]
    config = out / "default"
    assert {p.name for p in config.iterdir()} == {"screen.csv", "log.json", "MIGRATED.txt", "select"}
    assert {p.name for p in (config / "select").iterdir()} == {*explore.SELECT_FILES, "stage-select.json", "log.json"}
    assert (config / "select" / "selected.json").read_bytes() == selected
    assert set(json.loads((config / "log.json").read_text())) == {"screen"} and set(json.loads((config / "select" / "log.json").read_text())) == {"select"}
    assert "screen.csv" in (config / "MIGRATED.txt").read_text()
    assert {p.name: p.read_bytes() for p in (out / "entropy").iterdir()} == entropy  # a screen without a select is left alone
    # the moved result reads back as current, and running it again changes nothing
    monkeypatch.setattr(explore, "load_configs", lambda path: ({"default": InfomapOptions()}, FAST))
    assert explore.read_selected(out, "default", tmp_path / "configs.toml").k == 4
    before = {p: p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}
    assert explore.migrate_select(root) == []
    assert {p: p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()} == before
    # a select stamp in both places is not guessed at
    (config / "stage-select.json").write_text("{}")
    with pytest.raises(ValueError, match="both here and in select"):
        explore.migrate_select(root)


# visits as the evidence (DECISIONS, "Species-days", 2026-11-06)


def test_run_name_ends_in_visits_for_visits_and_refuses_other_evidence():
    assert explore.run_name("native_only", ["griis"], evidence="visits") == "native_only-no-griis-visits"
    assert explore.run_name("native_only", [], merge_duplicate_names=True, evidence="visits") == "native_only-merged-names-visits"
    assert explore.run_name("native_only", [], evidence="records") == "native_only"
    with pytest.raises(ValueError, match="evidence"):
        explore.run_name("native_only", [], evidence="photos")


def test_evidence_counts_puts_visits_in_place_of_records_and_scales_the_record_subsets():
    from realms.occurrences import COUNTS

    row = {c: 0 for c in COUNTS} | {"species_id": "S1", "countrycode": "DK", "part": "c:land", "n_raw": 120, "n": 100, "n_open": 80, "n_human_text": 10,
                                   "n_ebird": 50, "n_years": 7, "n_visits": 40, "n_visits_open": 30}
    counts = pl.DataFrame([row])
    visits = explore.evidence_counts(counts, "visits")
    assert "n_visits" not in visits.columns and "n_visits_open" not in visits.columns
    got = visits.row(0, named=True)
    assert (got["n"], got["n_open"], got["n_human_text"], got["n_ebird"], got["n_raw"], got["n_years"]) == (40, 30, 4, 20, 120, 7)
    records = explore.evidence_counts(counts, "records").row(0, named=True)
    assert (records["n"], records["n_human_text"]) == (100, 10) and "n_visits" not in records
    with pytest.raises(ValueError, match="without visits"):
        explore.evidence_counts(counts.drop("n_visits", "n_visits_open"), "visits")


def test_records_as_evidence_leave_the_clean_hash_as_it_was():
    from realms.clean import Settings

    marker = {"rows": 1, "resolution": 7, "zip": "x.zip", "zip_sha256": "0"}
    codes = pl.DataFrame({"code": ["DK"]})
    plain = explore.clean_inputs(Settings(), [1, 2], codes, marker)
    assert explore.clean_inputs(Settings(), [1, 2], codes, marker, evidence="records") == plain
    assert explore.clean_inputs(Settings(), [1, 2], codes, marker, evidence="visits") == plain | {"evidence": "visits"}


# eukaryotes only (DECISIONS, "Eukaryotes only", 2026-11-07)


def test_a_eukaryotes_only_run_has_a_name_and_a_hash_of_its_own():
    from realms.clean import Settings

    assert explore.run_name("native_only", [], eukaryotes_only=True) == "native_only-eukaryotes"
    assert explore.run_name("native_only", [], evidence="visits", eukaryotes_only=True) == "native_only-visits-eukaryotes"
    marker = {"rows": 1, "resolution": 7, "zip": "x.zip", "zip_sha256": "0"}
    codes = pl.DataFrame({"code": ["DK"]})
    plain = explore.clean_inputs(Settings(), [1, 2], codes, marker)
    assert explore.clean_inputs(Settings(), [1, 2], codes, marker, eukaryotes_only=False) == plain
    assert explore.clean_inputs(Settings(), [1, 2], codes, marker, eukaryotes_only=True) == plain | {"kingdoms": list(explore.EUKARYOTES)}


def test_eukaryote_counts_keep_the_eukaryotes_and_drop_the_rest():
    counts = pl.DataFrame({"species_id": ["A", "B", "V", "N", "X"], "n": [1, 2, 3, 4, 5]})
    taxa = pl.DataFrame({"species_id": ["A", "B", "V", "N"], "kingdom": ["Animalia", "Fungi", "Orthornavirae", None]})
    assert explore.eukaryote_counts(counts, taxa)["species_id"].to_list() == ["A", "B"]
