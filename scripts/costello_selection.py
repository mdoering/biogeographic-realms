"""Try the clustering and selection rules on Costello et al.'s matrix.

This checks the method on data others have analysed. There is no target realm count, so nothing here is tuned to reach one.
"""

import argparse
import resource
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import polars as pl
from shapely.geometry import Point

from realms.candidates import DegenerateClustering, candidate_maps
from realms.cluster import run_infomap
from realms.costello import csquare_adjacency, csquare_area_km2, csquare_bounds, load_matrix
from realms.endemism import endemism_score, realm_endemism
from realms.fetch import fetch
from realms.maps import plot_cells
from realms.matrix import RANKS as ALL_RANKS
from realms.matrix import cell_taxa
from realms.network import network_from
from realms.nullmodel import excess_table
from realms.score import v_measure
from realms.select import bootstrap_stability, select_level

MATRIX = Path("data/costello/obis_2009_csq.csv")
OUT = Path("results/costello")
GAZETTEERS = Path("data/reference/gazetteers")
RAW = "https://raw.githubusercontent.com/CatalogueOfLife/col-gazetteers/main"
RANKS = ("species", "genus")


def stamp(message: str, started: float) -> None:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    print(f"[{time.time() - started:7.0f}s, peak {peak:4.1f} GB] {message}", flush=True)


def costello_realms(cells: list[str]) -> dict[str, int]:
    """Costello's realm for each cell centre. Cells whose centre falls in no polygon are left out."""
    paths = [fetch(f"{RAW}/costello/features/{i}.geojson", GAZETTEERS / "costello" / f"{i}.geojson") for i in range(1, 31)]
    realms = pl.concat(
        [pl.DataFrame({"path": [str(p)], "realm": [int(p.stem)]}) for p in paths]
    )
    polygons = gpd.GeoDataFrame(
        {"realm": realms["realm"].to_list()},
        geometry=[gpd.read_file(p).geometry.union_all() for p in realms["path"]],
        crs="EPSG:4326",
    )
    centres = []
    for cell in cells:
        b = csquare_bounds(cell)
        centres.append(Point((b[2] + b[3]) / 2, (b[0] + b[1]) / 2))
    points = gpd.GeoDataFrame({"cell": cells}, geometry=centres, crs="EPSG:4326")
    joined = gpd.sjoin(points, polygons, predicate="within", how="inner").drop_duplicates("cell")
    return dict(zip(joined["cell"], joined["realm"].astype(int)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--boot", type=int, default=20)
    parser.add_argument("--null", type=int, default=99)
    parser.add_argument("--workers", type=int, default=1, help="processes for the bootstraps and the null maps, which give the same results as one")
    args = parser.parse_args()

    started = time.time()
    matrix = load_matrix(MATRIX)
    richness = np.asarray(matrix.X.T.tocsr().getnnz(axis=1))
    cells = [c for c, r in zip(matrix.cells, richness) if r >= 5]
    presence = matrix.presence().filter(pl.col("cell").is_in(cells))
    taxa = matrix.genera()
    areas = {c: csquare_area_km2(c) for c in cells}
    adjacency = csquare_adjacency(cells)
    ct = cell_taxa(presence, taxa, RANKS)
    species_only = ct.select(["species"])
    held = [rank for i, rank in enumerate(ALL_RANKS) if (ct.rank == i).any()]
    stamp(f"{presence.height} presences, {len(cells)} cells, ranks {held}", started)

    net = network_from(ct)
    stamp(f"network: {len(net.taxa)} taxa, {len(net.weight)} links", started)
    partition = run_infomap(net, seed=1, trials=10)
    leaves = partition.leaf_labels()
    sizes = sorted(np.unique(list(leaves.values()), return_counts=True)[1].tolist(), reverse=True)
    stamp(f"infomap: {len(set(partition.top_labels().values()))} top modules, {len(sizes)} leaf modules, largest {sizes[:12]}, "
          f"{len(partition.unassigned)} unassigned", started)

    OUT.mkdir(parents=True, exist_ok=True)
    plot_cells(list(leaves), list(leaves.values()), csquare_bounds, OUT / "infomap_leaf_modules.png", "Infomap leaf modules, species and genus")

    try:
        result = candidate_maps(partition, ct, areas, adjacency)
    except DegenerateClustering as e:
        stamp(f"DEGENERATE: {e}", started)
        return
    maps = result.maps
    stamp(f"candidates: k from {min(maps)} to {max(maps)}, {len(result.sub_realms)} sub-realms: {result.sub_realms}", started)

    stability = bootstrap_stability(ct, maps, areas, adjacency, n_boot=args.boot, workers=args.workers)
    stamp("stability done", started)
    excess = excess_table(maps, ct, adjacency, areas, n_null=args.null, workers=args.workers)
    stamp("excess done", started)

    table = pl.DataFrame(
        {"k": sorted(maps), "stability": [stability[k] for k in sorted(maps)], "excess": [excess[k] for k in sorted(maps)]},
        schema={"k": pl.Int64, "stability": pl.Float64, "excess": pl.Float64},
    )
    table.write_csv(OUT / "selection.csv")
    with pl.Config(tbl_rows=70):
        print(table)
    selection = select_level(table)
    print(f"SELECTED k={selection.k}, fell_short={selection.fell_short}: {selection.reason}")

    costello = costello_realms(cells)
    print(f"Costello realm found for {len(costello)} of {len(cells)} cells")
    for k in sorted({selection.k, 7, 30} & set(maps)):
        labels = maps[k]
        plot_cells(list(labels), list(labels.values()), csquare_bounds, OUT / f"selection_k{k}.png", f"Infomap and selection rules, k={k}")
        shares = realm_endemism(labels, species_only)
        print(
            f"k={k}: V-measure with Costello's 30 = {v_measure(labels, costello):.3f}, "
            f"species endemism score = {endemism_score(labels, species_only, areas):.3f}, "
            f"strict = {endemism_score(labels, species_only, areas, tolerance=1.0):.3f}, "
            f"mean realm share = {np.mean(list(shares.values())):.3f}"
        )
    strict_shares = realm_endemism(costello, species_only, tolerance=1.0)
    print(
        f"Costello's 30 on the same cells: species endemism score = {endemism_score(costello, species_only, areas):.3f}, "
        f"strict = {endemism_score(costello, species_only, areas, tolerance=1.0):.3f}, "
        f"strict mean realm share = {np.mean(list(strict_shares.values())):.3f} (published mean: 0.42), "
        f"range {min(strict_shares.values()):.2f} to {max(strict_shares.values()):.2f} (published: 0.17 to 0.84)"
    )
    stamp("finished", started)


if __name__ == "__main__":
    main()
