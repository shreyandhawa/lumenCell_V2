import queue
from collections.abc import Iterator
from pathlib import Path

import pytest

from lumencell.bench.synth import generate
from lumencell.data.loader import load
from lumencell.engine.store import RunStore
from lumencell.jobs.manager import JobManager
from lumencell.prep.cache import PrepCache, dataset_dir
from lumencell.prep.stages import STAGES, PrepParams


@pytest.fixture
def jobs() -> Iterator[JobManager]:
    manager = JobManager().start()
    yield manager
    manager.shutdown()


def _drain(
    q: "queue.Queue[dict[str, object]]", until: str = "done", timeout: float = 120
) -> list[dict[str, object]]:
    events = []
    while True:
        event = q.get(timeout=timeout)
        events.append(event)
        if event["state"] in (until, "failed", "cancelled"):
            return events


def test_sleep_job_reports_progress_in_order(jobs: JobManager) -> None:
    job = jobs.submit("sleep", "s1", {"steps": 4, "seconds": 0.05})
    events = _drain(jobs.subscribe(job.job_id))
    assert events[-1]["state"] == "done"
    steps = [e["i"] for e in events if e["state"] == "running" and e["stage"] == "sleep"]
    assert steps == sorted(steps) and steps[-1] == 3
    assert jobs.get(job.job_id).fraction == 1.0  # type: ignore[union-attr]


def test_identical_active_jobs_are_deduplicated(jobs: JobManager) -> None:
    a = jobs.submit("sleep", "same", {"steps": 3, "seconds": 0.1})
    b = jobs.submit("sleep", "same", {"steps": 3, "seconds": 0.1})
    assert a.job_id == b.job_id
    jobs.wait(a.job_id)
    c = jobs.submit("sleep", "same", {"steps": 1, "seconds": 0.01})
    assert c.job_id != a.job_id, "a finished job is not reused"


def test_failure_carries_the_reason(jobs: JobManager) -> None:
    job = jobs.submit("sleep", "f", {"steps": 1, "seconds": 0.01, "fail": True})
    done = jobs.wait(job.job_id)
    assert done.state == "failed"
    assert done.error is not None and "asked to fail" in done.error


def test_cancel_queued_job_never_runs(jobs: JobManager) -> None:
    first = jobs.submit("sleep", "long", {"steps": 20, "seconds": 0.1})
    second = jobs.submit("sleep", "queued", {"steps": 1, "seconds": 0.01})
    assert jobs.cancel(second.job_id).state == "cancelled"  # type: ignore[union-attr]
    jobs.cancel(first.job_id)
    assert jobs.get(second.job_id).started is None  # type: ignore[union-attr]


def test_cancel_running_job_restarts_the_worker(jobs: JobManager) -> None:
    job = jobs.submit("sleep", "long", {"steps": 100, "seconds": 0.1})
    q = jobs.subscribe(job.job_id)
    while q.get(timeout=30)["state"] != "running":
        pass
    assert jobs.cancel(job.job_id).state == "cancelled"  # type: ignore[union-attr]
    after = jobs.submit("sleep", "after", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(after.job_id, timeout=60).state == "done", "a fresh worker takes the next job"


def test_a_worker_that_dies_mid_job_fails_the_job_and_is_replaced(jobs: JobManager) -> None:
    job = jobs.submit("sleep", "long", {"steps": 100, "seconds": 0.1})
    q = jobs.subscribe(job.job_id)
    while q.get(timeout=30)["state"] != "running":
        pass
    jobs._worker.kill()
    done = jobs.wait(job.job_id, timeout=30)
    assert done.state == "failed"
    assert done.error is not None and "stopped unexpectedly" in done.error
    after = jobs.submit("sleep", "after", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(after.job_id, timeout=60).state == "done"


def test_a_worker_that_died_while_idle_is_replaced_before_the_next_job(jobs: JobManager) -> None:
    warm = jobs.submit("sleep", "warm", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(warm.job_id, timeout=60).state == "done"
    jobs._worker.kill()
    jobs._worker.join(10)
    after = jobs.submit("sleep", "after", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(after.job_id, timeout=60).state == "done", "the first job after must not be lost"


def test_prepare_then_infer_through_the_worker(jobs: JobManager, tmp_path: Path) -> None:
    path = tmp_path / "d.h5ad"
    generate("bifurcation", "realistic", n_cells=500, n_genes=300, seed=0).write_h5ad(path)
    cache_root = tmp_path / "cache"
    params = PrepParams()
    base = {"path": str(path), "prep_params": params, "cache_root": str(cache_root)}

    prep = jobs.submit("prepare", "prep", base)
    events = _drain(jobs.subscribe(prep.job_id))
    assert events[-1]["state"] == "done", events[-1]
    seen = [e["stage"] for e in events if e["state"] == "running" and e["stage"]]
    assert [s for s in dict.fromkeys(seen) if s != "done"] == [s.name for s in STAGES]

    ds = load(path)
    assert set(PrepCache(ds, cache_root).status(params).values()) == {"done"}

    infer = jobs.submit(
        "infer",
        "run",
        {**base, "method": "paga_dpt", "root": "cluster:true_label=trunk_early", "run_id": "r1"},
    )
    done = jobs.wait(infer.job_id, timeout=120)
    assert done.state == "done", done.error
    assert done.result_id == "r1"
    store = RunStore(dataset_dir(ds.fingerprint, cache_root))
    assert store.exists("r1")
    store.load("r1").validate()
