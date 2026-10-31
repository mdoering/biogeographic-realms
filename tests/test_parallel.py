import os

from realms.parallel import map_tasks


def scaled(shared, task):
    return shared * task


def where(shared, task):
    return shared, os.getpid(), os.environ.get("OMP_NUM_THREADS")


def test_results_come_back_in_the_order_of_the_tasks():
    tasks = list(range(10))
    assert map_tasks(scaled, tasks, 3, workers=3) == [3 * t for t in tasks]
    assert map_tasks(scaled, tasks, 3) == map_tasks(scaled, tasks, 3, workers=3)


def test_workers_each_hold_the_shared_tables_on_one_openmp_thread():
    shared = {"table": [1, 2, 3]}
    answers = map_tasks(where, list(range(6)), shared, workers=2)
    assert all(table == shared for table, _, _ in answers)
    assert os.getpid() not in {pid for _, pid, _ in answers}
    assert {threads for _, _, threads in answers} == {os.environ.get("OMP_NUM_THREADS", "1")}  # the caller's choice, else one


def test_one_worker_stays_in_this_process():
    assert {pid for _, pid, _ in map_tasks(where, [1, 2], None)} == {os.getpid()}


def test_no_tasks_start_no_pool():
    assert map_tasks(scaled, [], 3, workers=4) == []
