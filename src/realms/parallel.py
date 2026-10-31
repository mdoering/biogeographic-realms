"""One function over many tasks, in this process or in a pool of workers that each hold the same tables."""

import multiprocessing
import os
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import Any

_shared: Any = None  # what this worker was started with


def _install(shared: Any) -> None:
    global _shared
    _shared = shared
    # Infomap starts a thread per core by default and gains nothing from them, so several workers would crowd each other out
    os.environ.setdefault("OMP_NUM_THREADS", "1")


def _call(function: Callable[[Any, Any], Any], task: Any) -> Any:
    return function(_shared, task)


def map_tasks(function: Callable[[Any, Any], Any], tasks: Sequence[Any], shared: Any, workers: int = 1) -> list[Any]:
    """`function(shared, task)` for each task, results in the order of the tasks.

    With more than one worker the tasks go to a pool of spawned processes, which is also how macOS starts them, so nothing is
    inherited. `shared` is sent to each worker once, as it starts, and not with every task. A worker runs OpenMP, and so Infomap,
    on one thread unless `OMP_NUM_THREADS` says otherwise. `function` has to be a top-level function that the worker can import. It
    must not read anything but its two arguments, so that the results don't depend on the number of workers.
    """
    if workers <= 1 or len(tasks) < 2:
        return [function(shared, task) for task in tasks]
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(min(workers, len(tasks)), mp_context=context, initializer=_install, initargs=(shared,)) as pool:
        return list(pool.map(partial(_call, function), tasks))
