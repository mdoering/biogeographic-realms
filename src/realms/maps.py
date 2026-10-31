"""Working maps of realms on latitude-longitude cells.

A realm map can hold more realms than any palette can keep apart, so colour only separates neighbouring realms and the number
printed on each realm is what identifies it.
"""

from collections import defaultdict
from collections.abc import Callable, Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PatchCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.patheffects import withStroke  # noqa: E402

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def _neighbouring_realms(bounds_by_cell: dict[str, tuple], label_of: dict[str, int]) -> dict[int, set[int]]:
    by_edge: dict[tuple, set[int]] = defaultdict(set)
    for cell, (lat_min, lat_max, lon_min, lon_max) in bounds_by_cell.items():
        for edge in (("lat", lat_min, lon_min), ("lat", lat_max, lon_min), ("lon", lon_min, lat_min), ("lon", lon_max, lat_min)):
            by_edge[edge].add(label_of[cell])
    neighbours: dict[int, set[int]] = defaultdict(set)
    for realms in by_edge.values():
        for realm in realms:
            neighbours[realm] |= realms - {realm}
    return neighbours


def _separating_colours(realms: list[int], neighbours: dict[int, set[int]]) -> dict[int, str]:
    colour: dict[int, int] = {}
    for realm in sorted(realms, key=lambda r: (-len(neighbours[r]), r)):
        used = {colour[n] for n in neighbours[realm] if n in colour}
        colour[realm] = next((i for i in range(len(HUES)) if i not in used), 0)
    return {realm: HUES[i] for realm, i in colour.items()}


def plot_cells(cells: list[str], labels: Sequence[int], bounds: Callable, out: Path, title: str) -> None:
    label_of = {c: int(label) for c, label in zip(cells, labels)}
    bounds_by_cell = {c: bounds(c) for c in cells}
    realms = sorted(set(label_of.values()))
    colours = _separating_colours(realms, _neighbouring_realms(bounds_by_cell, label_of))

    fig, ax = plt.subplots(figsize=(14, 7.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    patches = [Rectangle((b[2], b[0]), b[3] - b[2], b[1] - b[0]) for b in bounds_by_cell.values()]
    ax.add_collection(
        PatchCollection(patches, facecolors=[colours[label_of[c]] for c in cells], edgecolors=SURFACE, linewidths=0.4)
    )

    for realm in realms:
        members = [c for c in cells if label_of[c] == realm]
        centres = np.array([[(bounds_by_cell[c][0] + bounds_by_cell[c][1]) / 2, (bounds_by_cell[c][2] + bounds_by_cell[c][3]) / 2] for c in members])
        anchor = centres[np.argmin(((centres - np.median(centres, axis=0)) ** 2).sum(axis=1))]
        ax.text(
            anchor[1], anchor[0], str(realm), ha="center", va="center", fontsize=10, fontweight="bold", color=TEXT,
            path_effects=[withStroke(linewidth=3, foreground=SURFACE)],
        )

    ax.set(xlim=(-180, 180), ylim=(-90, 90), aspect="equal")
    ax.set_title(title, loc="left", color=TEXT, fontsize=13)
    ax.tick_params(colors="#52514e", labelsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)
