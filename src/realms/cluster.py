"""Community detection on the cell-taxon network."""

from collections import Counter
from dataclasses import dataclass, replace

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


@dataclass(frozen=True)
class InfomapOptions:
    """What differs from Infomap's defaults when clustering. The defaults add nothing to the command line.

    `markov_time` below 1 gives more and smaller modules, above 1 fewer and larger ones. `variable_markov_time` varies it locally,
    to reduce overpartitioning in sparse areas while keeping the resolution in dense ones, and `variable_markov_damping` sets the
    damping between the local effective degree (0) and the local entropy (1). `regularized` adds a Bayesian prior network against
    overfitting to missing links, and `entropy_corrected` corrects for the entropy bias in small samples, both with a strength to
    scale them by. A damping or strength is only read by Infomap together with its switch.

    `recurse_share` is not Infomap's: every top module holding more than this share of the cells is clustered again on its own,
    as `run_infomap` says.
    """

    markov_time: float = 1.0
    variable_markov_time: bool = False
    variable_markov_damping: float | None = None
    regularized: bool = False
    regularization_strength: float | None = None
    entropy_corrected: bool = False
    entropy_correction_strength: float | None = None
    recurse_share: float | None = None

    def args(self) -> str:
        """The flags that differ from the defaults, for the end of Infomap's command line. Nothing at all for the defaults."""
        flags = []
        if self.markov_time != 1.0:
            flags.append(f"--markov-time {self.markov_time}")
        if self.variable_markov_time:
            flags.append("--variable-markov-time")
        if self.variable_markov_damping is not None:
            flags.append(f"--variable-markov-damping {self.variable_markov_damping}")
        if self.regularized:
            flags.append("--regularized")
        if self.regularization_strength is not None:
            flags.append(f"--regularization-strength {self.regularization_strength}")
        if self.entropy_corrected:
            flags.append("--entropy-corrected")
        if self.entropy_correction_strength is not None:
            flags.append(f"--entropy-correction-strength {self.entropy_correction_strength}")
        return " ".join(flags)


def load_infomap(net: Network, seed: int = 1, trials: int = 10, options: InfomapOptions = InfomapOptions()) -> infomap.Infomap:
    """An Infomap holding every link of the network in one call, ready to run.

    Cells are nodes `0..n_cells-1` and taxa follow from `n_cells`, so the two sides of the bipartite network never share an id. The
    flags of `options` that differ from the defaults come last, and the defaults leave the command line as it was.
    """
    line = f"--flow-model undirected --seed {seed} --num-trials {trials} --silent"
    if flags := options.args():
        line += f" {flags}"
    n_cells = len(net.cells)
    links = np.empty((len(net.weight), 3), dtype=np.float64)  # source, target, weight: all one float dtype, node ids are exact
    links[:, 0] = net.cell_idx
    links[:, 1] = net.taxon_idx
    links[:, 1] += n_cells
    links[:, 2] = net.weight
    im = infomap.Infomap(line)
    im.bipartite_start_id = n_cells
    im.add_links(links)
    return im


def _cluster(net: Network, seed: int, trials: int, options: InfomapOptions) -> Partition:
    """One Infomap run on the network, with no recursion."""
    n_cells = len(net.cells)
    linked = np.bincount(net.cell_idx, minlength=n_cells) > 0  # a cell with any link: a pass over the links, with no Python set of them
    unassigned = sorted(net.cells[i] for i in np.flatnonzero(~linked))
    if not linked.any():
        return Partition(paths={}, unassigned=unassigned)

    modules = load_infomap(net, seed, trials, options).run().multilevel_modules()
    paths = {net.cells[node]: tuple(path) for node, path in modules.items() if node < n_cells}
    return Partition(paths=paths, unassigned=unassigned)


def _subnetwork(net: Network, keep: np.ndarray) -> Network:
    """The cells flagged in `keep`, the taxa that they link to and every link between the two, each with the weight it had.

    The cells and taxa stay in the order of `net`.
    """
    renumber = np.cumsum(keep) - 1
    linked = keep[net.cell_idx]
    used, taxon_idx = np.unique(net.taxon_idx[linked], return_inverse=True)
    return Network(
        cells=[net.cells[i] for i in np.flatnonzero(keep)],
        taxa=[net.taxa[i] for i in used] if net.taxa else [],  # the pilot's child builds its network without the names of the taxa
        cell_idx=renumber[net.cell_idx[linked]],
        taxon_idx=taxon_idx,
        weight=net.weight[linked],
    )


def _recurse(net: Network, partition: Partition, seed: int, trials: int, options: InfomapOptions) -> Partition:
    """Cluster each top module with more than `options.recurse_share` of the network's cells again, on its own.

    The sub-network keeps the original link weights. The run has the seed, trials and options of the first one, with the recursion
    off. If it finds several modules, they replace the module's, below its top label, so that the top modules stay what they were
    and every leaf has a path of its own, and so a label of its own. A module it leaves whole keeps its paths.
    """
    top = partition.top_labels()
    big = sorted(module for module, n in Counter(top.values()).items() if n / len(net.cells) > options.recurse_share)
    if not big:
        return partition
    index = {cell: i for i, cell in enumerate(net.cells)}
    module_of = np.full(len(net.cells), -1, dtype=np.int64)  # a cell without a module, because it has no links, is -1
    for cell, module in top.items():
        module_of[index[cell]] = module
    inner = replace(options, recurse_share=None)
    paths = dict(partition.paths)
    for module in big:
        found = _cluster(_subnetwork(net, module_of == module), seed, trials, inner)
        if len(set(found.paths.values())) > 1:
            paths.update({cell: (module, *path) for cell, path in found.paths.items()})
    return Partition(paths=paths, unassigned=partition.unassigned)


def run_infomap(net: Network, seed: int = 1, trials: int = 10, options: InfomapOptions = InfomapOptions()) -> Partition:
    """Cluster the network with the given options.

    With `options.recurse_share`, each top module holding more of the network's cells than that is clustered again: on its cells,
    the taxa linked to them and all the links between them, with their original weights. That is the same run, with the same seed
    and trials, and the recursion off. Its modules replace the big one, which is what `Partition.leaf_labels` then numbers.
    """
    partition = _cluster(net, seed, trials, options)
    if options.recurse_share is None:
        return partition
    return _recurse(net, partition, seed, trials, options)
