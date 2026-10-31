"""Community detection on the cell-taxon network."""

from dataclasses import dataclass

import infomap
import numpy as np

from realms.network import Network


@dataclass(frozen=True)
class Partition:
    """Module path per cell, outermost module first, plus the cells that had no links to cluster on."""

    paths: dict[str, tuple[int, ...]]
    unassigned: list[str]

    def leaf_labels(self) -> dict[str, int]:
        """The finest modules, numbered in sorted path order so the numbering is deterministic."""
        numbers = {path: i for i, path in enumerate(sorted(set(self.paths.values())))}
        return {cell: numbers[path] for cell, path in self.paths.items()}

    def top_labels(self) -> dict[str, int]:
        """The outermost modules."""
        return {cell: path[0] for cell, path in self.paths.items()}


def load_infomap(net: Network, seed: int = 1, trials: int = 10, markov_time: float = 1.0) -> infomap.Infomap:
    """An Infomap holding every link of the network in one call, ready to run.

    Cells are nodes `0..n_cells-1` and taxa follow from `n_cells`, so the two sides of the bipartite network never share an id. A
    markov time other than 1 is passed on, and the default leaves the command line as it was.
    """
    options = f"--flow-model undirected --seed {seed} --num-trials {trials} --silent"
    if markov_time != 1.0:
        options += f" --markov-time {markov_time}"
    n_cells = len(net.cells)
    links = np.empty((len(net.weight), 3), dtype=np.float64)  # source, target, weight: all one float dtype, node ids are exact
    links[:, 0] = net.cell_idx
    links[:, 1] = net.taxon_idx
    links[:, 1] += n_cells
    links[:, 2] = net.weight
    im = infomap.Infomap(options)
    im.bipartite_start_id = n_cells
    im.add_links(links)
    return im


def run_infomap(net: Network, seed: int = 1, trials: int = 10, markov_time: float = 1.0) -> Partition:
    """Cluster the network. A `markov_time` below 1 gives more and smaller modules, above 1 fewer and larger ones."""
    n_cells = len(net.cells)
    linked = set(net.cell_idx.tolist())
    unassigned = sorted(net.cells[i] for i in range(n_cells) if i not in linked)
    if not linked:
        return Partition(paths={}, unassigned=unassigned)

    modules = load_infomap(net, seed, trials, markov_time).run().multilevel_modules()
    paths = {net.cells[node]: tuple(path) for node, path in modules.items() if node < n_cells}
    return Partition(paths=paths, unassigned=unassigned)
