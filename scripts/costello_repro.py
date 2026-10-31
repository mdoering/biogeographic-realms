"""Reproduce the beta-sim regionalisation of Castro-Insua et al. 2018 from Costello et al.'s matrix.

Follows the critique's own script: cells with at least 5 species, Simpson dissimilarity, R's hclust with ward.D2 and average
linkage, ANOSIM R for k from 2 to 50. Clustering runs in R because the tree isn't unique: over half of all cell pairs have a
dissimilarity of exactly 1, and the tie-breaking changes the map. The second half of this script measures by how much. The
checks aren't tuned, so a failure is a finding. The published values were read by eye from the critique's Figure 1.
"""

import csv
import math
import subprocess
from pathlib import Path

import numpy as np
from scipy.spatial.distance import squareform

from realms.betasim import anosim_r, beta_sim, hierarchical_labels
from realms.costello import csquare_bounds, load_matrix
from realms.maps import plot_cells
from realms.score import v_measure

MATRIX = Path("data/costello/obis_2009_csq.csv")
OUT = Path("results/costello")
HCLUST = Path(__file__).with_name("costello_hclust.R")
FIGURE_1_WARD = {2: 0.24, 3: 0.255, 4: 0.21, 5: 0.37, 6: 0.35, 7: 0.44, 8: 0.40, 9: 0.41, 10: 0.455, 20: 0.56, 30: 0.67, 50: 0.755}
TOLERANCE = 0.015
PLACES = {
    "North Pacific": (30, -170),
    "South Pacific": (-30, -130),
    "North Atlantic": (40, -40),
    "South Atlantic": (-30, -10),
    "Indian": (-20, 80),
    "West Antarctic": (-65, -60),
    "East Antarctic": (-65, 100),
}


def nearest_cell(cells: list[str], lat: float, lon: float) -> int:
    def distance(cell: str) -> float:
        b = csquare_bounds(cell)
        clat, clon = math.radians((b[0] + b[1]) / 2), math.radians((b[2] + b[3]) / 2)
        la, lo = math.radians(lat), math.radians(lon)
        return math.acos(min(1.0, math.sin(la) * math.sin(clat) + math.cos(la) * math.cos(clat) * math.cos(lo - clon)))

    return min(range(len(cells)), key=lambda i: distance(cells[i]))


def report(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'  ' + detail if detail else ''}")


def named(cells: list[str], labels) -> dict[str, int]:
    return {c: int(label) for c, label in zip(cells, labels)}


def main() -> None:
    matrix = load_matrix(MATRIX)
    by_cell = matrix.X.T.tocsr()
    richness = np.asarray(by_cell.getnnz(axis=1))
    keep = np.flatnonzero(richness >= 5)
    cells = [matrix.cells[i] for i in keep]
    print(f"{len(matrix.species)} species, {len(matrix.cells)} cells, {int((richness == 0).sum())} empty, {len(cells)} with >=5 species")

    dist = beta_sim(by_cell[keep])
    print(f"cell pairs with a dissimilarity of exactly 1: {(dist == 1).mean():.1%}")
    OUT.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT / "dist.txt", dist, fmt="%.10f")
    subprocess.run(["Rscript", str(HCLUST), str(OUT)], check=True)
    ward = np.loadtxt(OUT / "labels_ward.D2.txt", dtype=int)
    average = np.loadtxt(OUT / "labels_average.txt", dtype=int)
    (OUT / "dist.txt").unlink()

    curve = {k: (anosim_r(dist, ward[:, k - 2]), anosim_r(dist, average[:, k - 2])) for k in range(2, 51)}
    with (OUT / "anosim.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["k", "ward_r", "average_r"])
        writer.writerows((k, f"{w:.4f}", f"{a:.4f}") for k, (w, a) in curve.items())

    print("\nReproducing the critique, with R's clustering")
    off = {k: round(curve[k][0] - v, 3) for k, v in FIGURE_1_WARD.items() if abs(curve[k][0] - v) > TOLERANCE}
    report(f"Ward R matches Figure 1 within {TOLERANCE} at {len(FIGURE_1_WARD)} points", not off, f"off: {off}" if off else "")
    report("Ward R peaks at k=7 before falling back", curve[7][0] > curve[6][0] and curve[7][0] > curve[8][0], f"R: {curve[6][0]:.3f}, {curve[7][0]:.3f}, {curve[8][0]:.3f}")
    losing = [k for k in range(7, 51) if curve[k][0] <= curve[k][1]]
    report("Ward beats average for k above 6", not losing, f"exceptions: {losing}" if losing else "")

    k7 = ward[:, 5]
    at = {name: k7[nearest_cell(cells, *place)] for name, place in PLACES.items()}
    report("North and South Pacific differ", at["North Pacific"] != at["South Pacific"])
    report("North and South Atlantic match", at["North Atlantic"] == at["South Atlantic"])
    report("Indian differs from Atlantic", at["Indian"] != at["North Atlantic"])
    report("West and East Antarctic differ", at["West Antarctic"] != at["East Antarctic"])
    for k in (7, 30):
        plot_cells(cells, ward[:, k - 2].tolist(), csquare_bounds, OUT / f"beta_sim_ward_k{k}.png", f"Beta-sim, Ward (R ward.D2), k={k}: Costello's matrix, cells with at least 5 species")

    print("\nHow much the same method depends on tie-breaking, at k=7")
    scipy_k7 = hierarchical_labels(dist, "ward", 7)
    print(f"SciPy Ward against R: V-measure {v_measure(named(cells, k7), named(cells, scipy_k7)):.3f}, ANOSIM R {anosim_r(dist, scipy_k7):.3f} against {curve[7][0]:.3f}")
    square = squareform(dist)
    rng = np.random.default_rng(1)
    shuffled = []
    for _ in range(20):
        order = rng.permutation(len(cells))
        labels = np.empty(len(cells), dtype=int)
        labels[order] = hierarchical_labels(squareform(square[np.ix_(order, order)], checks=False), "ward", 7)
        shuffled.append(labels)
    against_r = [v_measure(named(cells, k7), named(cells, labels)) for labels in shuffled]
    between = [v_measure(named(cells, a), named(cells, b)) for i, a in enumerate(shuffled) for b in shuffled[i + 1 :]]
    rs = [anosim_r(dist, labels) for labels in shuffled]
    print(f"20 shuffled cell orders: V-measure against R {min(against_r):.3f} to {max(against_r):.3f} (median {np.median(against_r):.3f})")
    print(f"between shuffles: {min(between):.3f} to {max(between):.3f} (median {np.median(between):.3f}); ANOSIM R {min(rs):.3f} to {max(rs):.3f}")


if __name__ == "__main__":
    main()
