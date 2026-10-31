"""Time building and slicing the cell-taxon matrix on a synthetic grid, to see what won't scale to the pilot.

At 1,000,000 species the grid holds about 89 million links with genera and families, as many as the pilot is expected to. Each
stage prints its own time and the peak memory of the whole process so far, so the stages can be held against their targets.
"""

import argparse
import resource
import sys
import time

import numpy as np

from realms.matrix import cell_taxa
from realms.network import network_from
from realms.synth import grid_presence


def peak_gb() -> float:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1e9 if sys.platform == "darwin" else 1e6)  # bytes on macOS, KiB on Linux


def stage(name, build, describe):
    started = time.time()
    result = build()
    print(f"{name:14s} {time.time() - started:7.1f} s   peak {peak_gb():5.1f} GB   {describe(result)}", flush=True)
    return result


def matrix_note(ct) -> str:
    return f"{len(ct.cells):,} cells, {len(ct.taxa):,} taxa, {ct.x.nnz:,} links"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=int, default=100_000)
    parser.add_argument("--side", type=int, default=90, help="cells along one side of the grid")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    presence, taxa, _ = stage(
        "grid_presence",
        lambda: grid_presence(args.species, args.side, args.seed),
        lambda grid: f"{grid[0].height:,} presences, {grid[1].height:,} species",
    )
    ct = stage("cell_taxa", lambda: cell_taxa(presence, taxa), matrix_note)
    stage("network_from", lambda: network_from(ct), lambda net: f"{len(net.taxa):,} nodes, {len(net.weight):,} links")
    stage("rows", lambda: ct.rows(ct.cells), lambda table: f"{table.x.shape[0]:,} rows, {table.x.nnz:,} links")
    n_species = len(ct.species_ids())
    draws = np.random.default_rng(args.seed).integers(n_species, size=n_species)
    stage("resample", lambda: ct.resample(draws), matrix_note)


if __name__ == "__main__":
    main()
