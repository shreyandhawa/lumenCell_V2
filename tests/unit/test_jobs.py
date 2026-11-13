import queue
import signal
import threading
from collections.abc import Iterator
from pathlib import Path

import pytest

from lumencell.bench.synth import generate
from lumencell.data.loader import load
from lumencell.engine.roots import RootError
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


def test_the_named_clonal_columns_reach_the_run(tmp_path: Path) -> None:
    # In-process: a named column is scoring-only in the job, so it cannot be the root.
    path = tmp_path / "d.h5ad"
    generate("bifurcation", "realistic", n_cells=400, n_genes=200, seed=0).write_h5ad(path)
    payload = {
        "path": str(path),
        "prep_params": PrepParams(),
        "cache_root": str(tmp_path / "cache"),
        "method": "paga_dpt",
        "root": "cluster:true_label=trunk_early",
        "run_id": "r3",
        "clonal_columns": ["true_label"],
    }
    with pytest.raises(RootError, match="'true_label' is set as a clone or donor column"):
        HANDLERS["infer"](payload, lambda *_: None)


def test_foreground_jobs_overtake_queued_background_jobs(jobs: JobManager) -> None:
    block = jobs.submit("sleep", "block", {"steps": 5, "seconds": 0.1})
    bg = jobs.submit("score_clonal", "bg", {})
    fg = jobs.submit("sleep", "fg", {"steps": 1, "seconds": 0.01})
    for job in (block, bg, fg):
        jobs.wait(job.job_id, timeout=60)  # score_clonal has no handler yet, so it fails
    bg, fg = jobs.get(bg.job_id), jobs.get(fg.job_id)  # type: ignore[assignment]
    assert fg.started is not None and bg.started is not None and fg.started < bg.started  # type: ignore[union-attr]


def test_take_next_order() -> None:
    manager = JobManager()  # not started: nothing is dispatched
    ids = {
        name: manager.submit(kind, name, {}).job_id
        for name, kind in [
            ("bg1", "score_clonal"),
            ("fg1", "sleep"),
            ("bg2", "stability"),
            ("fg2", "prepare"),
        ]
    }
    with manager._lock:
        order = [manager._take_next() for _ in range(5)]
    assert order == [ids["fg1"], ids["fg2"], ids["bg1"], ids["bg2"], None]


def test_a_queued_background_job_can_be_cancelled(jobs: JobManager) -> None:
    jobs.submit("sleep", "block", {"steps": 5, "seconds": 0.1})
    bg = jobs.submit("score_clonal", "bg", {})
    assert jobs.cancel(bg.job_id).state == "cancelled"  # type: ignore[union-attr]
    with jobs._lock:
        assert bg.job_id not in jobs._background
    assert jobs.get(bg.job_id).started is None  # type: ignore[union-attr]


def test_done_callbacks_run_outside_the_lock(jobs: JobManager) -> None:
    seen: dict[str, object] = {}
    called = threading.Event()

    def callback(job: object) -> None:
        job_id = job.job_id  # type: ignore[attr-defined]
        fetched: list[object] = []
        # The thread needs the manager lock; it can only finish if the listener has let go of it.
        t = threading.Thread(target=lambda: fetched.append(jobs.get(job_id)))
        t.start()
        t.join(2)
        seen["finished"] = not t.is_alive()
        seen["state"] = job.state  # type: ignore[attr-defined]
        called.set()

    jobs.add_done_callback(callback)
    job = jobs.submit("sleep", "cb", {"steps": 1, "seconds": 0.01})
    assert called.wait(60)
    assert seen == {"finished": True, "state": "done"}
    assert jobs.get(job.job_id).state == "done"  # type: ignore[union-attr]


def test_callbacks_see_only_done_jobs(jobs: JobManager) -> None:
    seen: list[str] = []
    jobs.add_done_callback(lambda job: seen.append(job.key))
    bad = jobs.submit("sleep", "bad", {"steps": 1, "seconds": 0.01, "fail": True})
    assert jobs.wait(bad.job_id, timeout=60).state == "failed"
    good = jobs.submit("sleep", "good", {"steps": 1, "seconds": 0.01})
    jobs.wait(good.job_id, timeout=60)
    deadline = threading.Event()
    for _ in range(100):  # the callback runs just after the state changes
        if seen:
            break
        deadline.wait(0.05)
    assert seen == ["good"]


def test_a_failing_callback_is_reported_and_the_listener_keeps_going(
    jobs: JobManager, capsys: pytest.CaptureFixture[str]
) -> None:
    def boom(job: object) -> None:
        raise RuntimeError("callback exploded")

    finished = threading.Event()
    jobs.add_done_callback(boom)
    jobs.add_done_callback(lambda job: finished.set() if job.key == "two" else None)  # runs after boom
    first = jobs.submit("sleep", "one", {"steps": 1, "seconds": 0.01})
    jobs.wait(first.job_id, timeout=10)
    second = jobs.submit("sleep", "two", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(second.job_id, timeout=10).state == "done"
    assert finished.wait(10), "the listener went on to the later callbacks and the next job"
    err = capsys.readouterr().err
    assert f"callback after job {first.job_id}" in err and "callback exploded" in err


def test_scoring_has_its_own_time_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("LUMENCELL_JOB_TIME_LIMITS", raising=False)
    assert JobManager().time_limit("score_clonal") == 3600


def test_a_queued_foreground_job_can_be_cancelled(jobs: JobManager) -> None:
    jobs.submit("sleep", "block", {"steps": 5, "seconds": 0.1})
    fg = jobs.submit("sleep", "fg", {"steps": 1, "seconds": 0.01})
    assert jobs.cancel(fg.job_id).state == "cancelled"  # type: ignore[union-attr]
    with jobs._lock:
        assert fg.job_id not in jobs._foreground
    assert jobs.get(fg.job_id).started is None  # type: ignore[union-attr]


@pytest.mark.parametrize("how", ["cancel queued", "cancel running", "time limit", "worker death"])
def test_callbacks_do_not_run_for_jobs_that_did_not_finish(how: str) -> None:
    limits = {"sleep": 1.0} if how == "time limit" else None
    manager = JobManager(time_limits=limits).start()
    try:
        seen: list[str] = []
        sentinel_seen = threading.Event()

        def callback(job: object) -> None:
            seen.append(job.key)  # type: ignore[attr-defined]
            if job.key == "ok":  # type: ignore[attr-defined]
                sentinel_seen.set()

        manager.add_done_callback(callback)
        job = manager.submit("sleep", "victim", {"steps": 100, "seconds": 0.1})
        if how == "cancel queued":
            blocker = manager.submit("sleep", "blocker", {"steps": 100, "seconds": 0.1})
            q = manager.subscribe(job.job_id)
            while q.get(timeout=30)["state"] != "running":
                pass
            queued = manager.submit("sleep", "queued", {"steps": 1, "seconds": 0.01})
            manager.cancel(queued.job_id)
            manager.cancel(job.job_id)
            manager.cancel(blocker.job_id)
            victim = queued
        else:
            _running_worker(manager, job.job_id)
            if how == "cancel running":
                manager.cancel(job.job_id)
            elif how == "worker death":
                manager._worker.kill()
            victim = job
        final = manager.wait(victim.job_id, timeout=30)
        assert final.state == ("cancelled" if "cancel" in how else "failed")
        ok = manager.submit("sleep", "ok", {"steps": 1, "seconds": 0.01})
        assert manager.wait(ok.job_id, timeout=60).state == "done"
        assert sentinel_seen.wait(10)  # callbacks run in job order, so earlier ones have had their chance
        assert seen == ["ok"]
    finally:
        manager.shutdown()


def test_a_callback_that_raises_systemexit_does_not_kill_the_listener(
    jobs: JobManager, capsys: pytest.CaptureFixture[str]
) -> None:
    after = threading.Event()

    def quit_now(job: object) -> None:
        raise SystemExit(3)

    jobs.add_done_callback(quit_now)
    jobs.add_done_callback(lambda job: after.set() if job.key == "two" else None)
    first = jobs.submit("sleep", "one", {"steps": 1, "seconds": 0.01})
    jobs.wait(first.job_id, timeout=10)
    second = jobs.submit("sleep", "two", {"steps": 1, "seconds": 0.01})
    assert jobs.wait(second.job_id, timeout=10).state == "done"
    assert after.wait(10)
    err = capsys.readouterr().err
    assert f"callback after job {first.job_id}" in err and "SystemExit" in err
