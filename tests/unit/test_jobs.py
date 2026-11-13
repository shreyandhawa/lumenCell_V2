import queue
import signal
import threading
from collections.abc import Iterator
from pathlib import Path

import pytest

from lumencell.bench.synth import generate
from lumencell.data.loader import load
from lumencell.engine.store import RunStore
from lumencell.jobs import manager as jobs_mod
from lumencell.jobs.manager import JobManager
from lumencell.jobs.worker import HANDLERS
from lumencell.prep.cache import PrepCache, dataset_dir, prep_key
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


def test_jobs_can_be_listed_by_kind(jobs: JobManager) -> None:
    first = jobs.submit("sleep", "one", {"steps": 1, "seconds": 0.01})
    second = jobs.submit("sleep", "two", {"steps": 1, "seconds": 0.01})
    assert [j.job_id for j in jobs.of_kind("sleep")] == [first.job_id, second.job_id]
    assert jobs.of_kind("infer") == []


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
        {
            **base,
            "method": "paga_dpt",
            "root": "cluster:true_label=trunk_early",
            "run_id": "r1",
            "clonal_columns": [],
        },
    )
    done = jobs.wait(infer.job_id, timeout=120)
    assert done.state == "done", done.error
    assert done.result_id == "r1"
    store = RunStore(dataset_dir(ds.fingerprint, cache_root))
    assert store.exists("r1")
    store.load("r1").validate()
    assert store.graph("r1")["prep_key"] == prep_key(ds, params)


def test_an_infer_job_must_say_which_columns_are_scoring_only(jobs: JobManager, tmp_path: Path) -> None:
    # No default: a payload that forgot the key would otherwise let methods see clone columns.
    payload = {
        "path": str(tmp_path / "never-read.h5ad"),
        "prep_params": PrepParams(),
        "method": "paga_dpt",
        "root": "cluster:true_label=trunk_early",
        "run_id": "r2",
    }
    done = jobs.wait(jobs.submit("infer", "no-columns", payload).job_id, timeout=120)
    assert done.state == "failed"
    assert done.error is not None and "clonal_columns" in done.error
    assert not done.error.startswith("KeyError"), done.error


def test_clonal_columns_must_be_a_list_of_names() -> None:
    # A bare string would be read as one column per character, leaving the real column visible.
    payload = {"path": "never-read.h5ad", "prep_params": PrepParams(), "clonal_columns": "donor"}
    with pytest.raises(ValueError, match="must be a list of column names"):
        HANDLERS["infer"](payload, lambda *_: None)


# ---- time limits and payloads -------------------------------------------------------------


class IgnoresSigterm:
    """Unpickled in the worker, it makes the worker ignore SIGTERM (like a library with its own handler)."""

    def __reduce__(self) -> tuple[object, tuple[object, ...]]:
        return (signal.signal, (signal.SIGTERM, signal.SIG_IGN))


def _running_worker(manager: JobManager, job_id: str) -> object:
    q = manager.subscribe(job_id)
    while q.get(timeout=30)["state"] != "running":
        pass
    return manager._worker


def test_a_job_over_its_time_limit_fails_and_the_next_one_runs(capsys: pytest.CaptureFixture[str]) -> None:
    manager = JobManager(time_limits={"sleep": 2.0}).start()
    try:
        job = manager.submit("sleep", "slow", {"steps": 1000, "seconds": 0.1})
        old = _running_worker(manager, job.job_id)
        done = manager.wait(job.job_id, timeout=30)
        assert done.state == "failed"
        assert done.error is not None and "timed out" in done.error and "time limit of 2 s" in done.error
        assert not old.is_alive(), "the worker running it was stopped"  # type: ignore[attr-defined]
        assert manager._worker.pid != old.pid  # type: ignore[attr-defined]
        assert "timed out" in capsys.readouterr().err, "the terminal says so too"
        after = manager.submit("sleep", "after", {"steps": 1, "seconds": 0.01})
        assert manager.wait(after.job_id, timeout=60).state == "done", "a fresh worker takes the next job"
    finally:
        manager.shutdown()


@pytest.mark.parametrize("stop", ["time limit", "cancel"])
def test_a_worker_that_ignores_sigterm_is_still_stopped(stop: str) -> None:
    manager = JobManager(time_limits={"sleep": 2.0 if stop == "time limit" else 600}).start()
    try:
        job = manager.submit("sleep", "stubborn", {"steps": 1000, "seconds": 0.1, "x": IgnoresSigterm()})
        old = _running_worker(manager, job.job_id)
        if stop == "cancel":
            manager.cancel(job.job_id)
        manager.wait(job.job_id, timeout=60)
        assert not old.is_alive(), "no second worker runs beside a stuck one"  # type: ignore[attr-defined]
    finally:
        manager.shutdown()


def test_an_unsendable_payload_fails_at_submit_and_the_queue_keeps_moving(jobs: JobManager) -> None:
    with pytest.raises(jobs_mod.PayloadError, match="cannot be sent"):
        jobs.submit("sleep", "bad", {"steps": 1, "lock": threading.Lock()})
    assert jobs.of_kind("sleep") == [], "nothing was queued"
    after = jobs.submit("sleep", "after", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(after.job_id, timeout=30).state == "done"


def test_every_job_kind_has_a_generous_default_time_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("LUMENCELL_JOB_TIME_LIMITS", raising=False)
    manager = JobManager()
    assert manager.time_limit("prepare") == manager.time_limit("infer") == 3 * 3600
    assert manager.time_limit("sleep") == 10 * 60
    assert all(manager.time_limit(kind) >= 600 for kind in HANDLERS)


def test_time_limits_can_be_set_from_the_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_JOB_TIME_LIMITS", "prepare=60, infer=90")
    assert JobManager().time_limit("prepare") == 60
    assert JobManager().time_limit("infer") == 90
    assert JobManager(time_limits={"prepare": 5}).time_limit("prepare") == 5, "the argument wins"
    for bad in ("prepare=soon", "prepare", "prepare=-1", "prepar=60", "prepare=nan"):
        monkeypatch.setenv("LUMENCELL_JOB_TIME_LIMITS", bad)
        with pytest.raises(ValueError, match="LUMENCELL_JOB_TIME_LIMITS"):
            JobManager()


@pytest.mark.parametrize(
    "limits", [{"prepare": float("nan")}, {"prepare": "60"}, {"prepare": 0}, {"prepare": True}, {"nope": 5}]
)
def test_time_limits_given_in_code_are_checked(limits: dict[str, object]) -> None:
    with pytest.raises(ValueError, match="time limit"):
        JobManager(time_limits=limits)  # type: ignore[arg-type]
