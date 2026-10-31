"""Community detection on the cell-taxon network."""

from dataclasses import dataclass

import infomap

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


def run_infomap(net: Network, seed: int = 1, trials: int = 10) -> Partition:
    n_cells = len(net.cells)
    linked = set(net.cell_idx.tolist())
    unassigned = sorted(net.cells[i] for i in range(n_cells) if i not in linked)
    if not linked:
        return Partition(paths={}, unassigned=unassigned)

    im = infomap.Infomap(f"--flow-model undirected --seed {seed} --num-trials {trials} --silent")
    im.bipartite_start_id = n_cells
    for cell, taxon, weight in zip(net.cell_idx.tolist(), net.taxon_idx.tolist(), net.weight.tolist()):
        im.add_link(cell, n_cells + taxon, weight)
    modules = im.run().multilevel_modules()

    paths = {net.cells[node]: tuple(path) for node, path in modules.items() if node < n_cells}
    return Partition(paths=paths, unassigned=unassigned)
