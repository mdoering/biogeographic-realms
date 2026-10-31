"""Time the pipeline on a synthetic grid, to see what won't scale to the pilot.

At 1,000,000 species the grid holds about 89 million links with genera and families, as many as the pilot is expected to. Each
stage prints its own time and the peak memory of the whole process so far, so the stages can be held against their targets. Infomap
itself is not run: the grid has no modules to find. A fixed partition into blocks of cells stands in for its answer wherever the
stages after it need one.
"""

import argparse
import resource
import sys
import time
from unittest import mock

import numpy as np

from realms import select
from realms.candidates import candidate_maps
from realms.cluster import Partition, load_infomap
from realms.endemism import realm_share
from realms.matrix import cell_taxa
from realms.network import network_from
from realms.nullmodel import excess_table
from realms.select import bootstrap_stability
from realms.synth import grid_presence

CELL_KM2 = 100_000  # most blocks of cells then clear the 1,000,000 km2 size floor, and a few small ones don't
UNITS_PER_SIDE = 23  # 529 blocks stand in for the leaf modules Infomap would find
REALMS = (6, 10)  # the 60 realms of the scoring stages, as blocks


def peak_gb(who: int = resource.RUSAGE_SELF) -> float:
    peak = resource.getrusage(who).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)  # bytes on macOS, KiB on Linux


def stage(name, build, describe):
    started = time.time()
    result = build()
    print(f"{name:14s} {time.time() - started:7.1f} s   peak {peak_gb():5.1f} GB   {describe(result)}", flush=True)
    return result


def matrix_note(ct) -> str:
    return f"{len(ct.cells):,} cells, {len(ct.taxa):,} taxa, {ct.x.nnz:,} links"


def blocks(side: int, cells: list[str], down: int, across: int) -> np.ndarray:
    """The block of a `down` x `across` layout that each cell of the square grid falls in. Cells are named `c<position>`."""
    row, col = np.divmod(np.array([int(c[1:]) for c in cells]), side)
    return (row * down // side) * across + col * across // side


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=int, default=100_000)
    parser.add_argument("--side", type=int, default=90, help="cells along one side of the grid")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--null", type=int, default=99, help="null maps per map")
    parser.add_argument("--maps", type=int, default=0, help="score this many of the candidate maps, spread across them, then scale up. 0 is all")
    args = parser.parse_args()

    presence, taxa, adjacency = stage(
        "grid_presence",
        lambda: grid_presence(args.species, args.side, args.seed),
        lambda grid: f"{grid[0].height:,} presences, {grid[1].height:,} species",
    )
    ct = stage("cell_taxa", lambda: cell_taxa(presence, taxa), matrix_note)
    net = stage("network_from", lambda: network_from(ct), lambda net: f"{len(net.taxa):,} nodes, {len(net.weight):,} links")
    stage("infomap input", lambda: load_infomap(net), lambda im: f"{len(net.weight):,} links added in one call")
    table = stage("rows", lambda: ct.rows(ct.cells), lambda table: f"{table.x.shape[0]:,} rows, {table.x.nnz:,} links")
    n_species = len(ct.species_ids())
    draws = np.random.default_rng(args.seed).integers(n_species, size=n_species)
    stage("resample", lambda: ct.resample(draws), matrix_note)

    areas = {c: float(CELL_KM2) for c in ct.cells}
    unit_blocks = blocks(args.side, ct.cells, UNITS_PER_SIDE, UNITS_PER_SIDE)
    units = Partition(paths={c: (int(b),) for c, b in zip(ct.cells, unit_blocks)}, unassigned=[])
    result = stage(
        "candidate_maps",
        lambda: candidate_maps(units, ct, areas, adjacency),
        lambda result: f"{len(set(unit_blocks.tolist()))} blocks, {len(result.maps)} maps, k {min(result.maps)} to {max(result.maps)}",
    )
    maps = result.maps

    with mock.patch.object(select, "run_infomap", lambda net, **options: units):
        stage(
            "bootstrap",
            lambda: bootstrap_stability(ct, maps, areas, adjacency, n_boot=1),
            lambda scores: f"one draw, Infomap replaced by the blocks, {len(scores)} maps scored",
        )

    realm_blocks = blocks(args.side, ct.cells, *REALMS)
    scattered = np.random.default_rng(args.seed).integers(REALMS[0] * REALMS[1], size=len(ct.cells))
    n_realms = REALMS[0] * REALMS[1]
    stage("xt", lambda: table.xt, lambda xt: f"{xt.nnz:,} links, built once per table")
    for name, labels in (("contiguous", realm_blocks), ("scattered", scattered)):
        stage(
            "realm_share",
            lambda: realm_share(table, labels, n_realms, 0.95),
            lambda shares: f"{n_realms} {name} realms, {table.x.nnz:,} links, one call",
        )

    keys = sorted(maps)
    picked = keys if not args.maps else [keys[i] for i in sorted(set(np.linspace(0, len(keys) - 1, args.maps).round().astype(int).tolist()))]
    del presence, net
    started = time.time()
    stage(
        "excess_table",
        lambda: excess_table({k: maps[k] for k in picked}, ct, adjacency, areas, n_null=args.null, workers=args.workers),
        lambda excess: f"{len(excess)} of {len(keys)} maps x {args.null} null maps x 3 ranks, {args.workers} workers, "
        f"worker peak {peak_gb(resource.RUSAGE_CHILDREN):.1f} GB",
    )
    elapsed = time.time() - started
    if len(picked) < len(keys):
        print(f"{'':14s} scaled to all {len(keys)} maps: {elapsed * len(keys) / len(picked) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
