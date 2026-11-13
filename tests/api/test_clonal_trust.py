"""Every run's GET carries its clonal agreement, scored in the background and healed on GET.

Each test builds its own job manager and session. Tests that check what a GET queues keep the
manager unstarted until those checks are done, so no job can finish early and race them. A test
waits for the score_clonal job itself: a waiter on the infer job is told it is done before the
session's done-callback has queued the scoring.
"""

import shutil
import threading
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient

from lumencell.api.app import create_app
from lumencell.api.models import ClonalSettings
from lumencell.api.session import Session
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
from lumencell.data.loader import load
from lumencell.jobs.manager import Job, JobManager
from lumencell.jobs.worker import HANDLERS
from lumencell.prep.cache import PrepCache
from lumencell.prep.stages import PrepParams
from lumencell.trust.clonal import settings_columns, settings_from_state, trust_hash, trust_name
from lumencell.trust.constants import N_PERMUTATIONS

ROOT = {"kind": "cluster", "key": "true_label", "value": "trunk_early"}
A = {"clone_key": "clone_id", "donor_key": "donor"}
B = {"clone_key": "clone_id", "donor_key": "donor", "strata_keys": ["true_segment"]}
SCORED = ("ok", "too_few_clones")
WAIT = 300.0


@pytest.fixture(scope="module")
def prepared(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    base = tmp_path_factory.mktemp("trust")
    path = base / "d.h5ad"
    adata = generate("bifurcation", "realistic", n_cells=1000, n_genes=300, seed=0)
    add_clones(adata, DEFAULT_LAYER, seed=0).write_h5ad(path)
    PrepCache(load(path), base / "cache").run(PrepParams())
    return path, base / "cache"


@pytest.fixture
def cache_root(prepared: tuple[Path, Path], tmp_path: Path) -> Path:
    """A cache of this test's own: the shared prep stages, and a run store no other test writes."""
    _, shared = prepared
    own = tmp_path / "cache"
    shutil.copytree(shared, own, ignore=shutil.ignore_patterns("runs"))
    return own


@pytest.fixture
def jobs() -> Iterator[JobManager]:
    manager = JobManager()  # not started: jobs wait in the queue until the test starts it
    yield manager
    manager.shutdown()


def _session(prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path) -> Session:
    path, _ = prepared
    s = Session(path, jobs, cache_root=cache_root, state_dir=tmp_path / "state")
    assert s.prep_ready()
    return s


def _client(session: Session) -> TestClient:
    return TestClient(create_app(session), raise_server_exceptions=False)


def _save_run(session: Session, rid: str, params: dict[str, Any] | None = None) -> None:
    """Compute and save a run in this process, as `lumencell run` does; no job is involved."""
    payload = {
        **session._job_payload(),
        "method": "paga_dpt",
        "root": "cluster:true_label=trunk_early",
        "params": params or {},
        "run_id": rid,
        "label_key": session.view.label_key,
        "clonal_columns": list(session.view.clonal_columns),
    }
    assert HANDLERS["infer"](payload, lambda *_: None) == rid


def _scoring(session: Session) -> list[Job]:
    return session.jobs.of_kind("score_clonal")


def _clonal(client: TestClient, rid: str) -> dict[str, Any]:
    r = client.get(f"/api/v1/trajectories/{rid}")
    assert r.status_code == 200, r.text
    clonal: dict[str, Any] = r.json()["trust"]["clonal"]
    return clonal


def _put(client: TestClient, settings: dict[str, Any]) -> None:
    r = client.put("/api/v1/dataset/clonal", json=settings)
    assert r.status_code == 200, r.text


def _wait_all(jobs: JobManager, which: list[Job]) -> None:
    for job in which:
        done = jobs.wait(job.job_id, timeout=WAIT)
        assert done.state not in ("queued", "running"), f"job {job.job_id} did not end in {WAIT} s"


def _after_done(jobs: JobManager, kind: str) -> threading.Event:
    """Set once a job of this kind is done and every callback registered before this one ran."""
    event = threading.Event()

    def note(job: Job) -> None:
        if job.kind == kind:
            event.set()

    jobs.add_done_callback(note)
    return event


def test_a_finished_infer_job_queues_scoring(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    inferred = _after_done(jobs, "infer")  # registered after the session's own callback
    client = _client(session)
    _put(client, A)
    jobs.start()
    created = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": ROOT}).json()
    rid = created["run_id"]
    assert jobs.wait(created["job_id"], timeout=WAIT).state == "done"
    assert inferred.wait(WAIT)
    (job,) = _scoring(session)
    assert job.key == f"score_clonal:{rid}:{trust_hash(ClonalSettings(**A))}"
    assert job.payload["trust_name"] == trust_name(ClonalSettings(**A))
    assert jobs.wait(job.job_id, timeout=WAIT).state == "done", job.error
    clonal = _clonal(client, rid)
    assert clonal["status"] in SCORED
    assert clonal["primary"]["n_permutations"] == N_PERMUTATIONS == 1000
    assert clonal["settings"] == ClonalSettings(**A).model_dump(mode="json")


def test_without_settings_the_status_is_not_configured(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    client = _client(session)
    _save_run(session, "r1")
    clonal = _clonal(client, "r1")
    assert clonal["status"] == "not_configured"
    assert clonal["settings"] is None and clonal["primary"] is None
    assert _scoring(session) == []
    _put(client, A)
    assert _clonal(client, "r1")["status"] == "pending"
    queued = len(_scoring(session))
    # Clearing the settings submits nothing; the next GET says so.
    assert client.delete("/api/v1/dataset/clonal").status_code == 204
    assert _clonal(client, "r1")["status"] == "not_configured"
    assert len(_scoring(session)) == queued == 1


def test_changing_settings_rescores_every_complete_run(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    client = _client(session)
    _save_run(session, "r1")
    _save_run(session, "r2", {"n_dcs": 8})
    _put(client, A)
    _put(client, B)
    hb = trust_hash(ClonalSettings(**B))
    keys = {j.key for j in _scoring(session) if j.key.endswith(f":{hb}")}
    assert keys == {f"score_clonal:r1:{hb}", f"score_clonal:r2:{hb}"}
    assert all(j.state == "queued" for j in _scoring(session))
    # The same settings again queue nothing new: the jobs are already waiting.
    _put(client, B)
    assert len(_scoring(session)) == 4


def test_a_restart_while_scoring_is_pending_heals_itself(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    first = _session(prepared, cache_root, JobManager(), tmp_path)  # never started; its jobs are lost
    _put(_client(first), A)
    _save_run(first, "r1")
    # The server restarts: a new session on a new manager, with the settings in the state file.
    session = _session(prepared, cache_root, jobs, tmp_path)
    client = _client(session)
    assert session.clonal_settings() == ClonalSettings(**A)
    pending = _clonal(client, "r1")
    assert pending["status"] == "pending"
    assert pending["settings"] == ClonalSettings(**A).model_dump(mode="json")
    assert _clonal(client, "r1")["status"] == "pending"
    (job,) = _scoring(session)
    jobs.start()
    assert jobs.wait(job.job_id, timeout=WAIT).state == "done", job.error
    assert _clonal(client, "r1")["status"] in SCORED
    assert _clonal(client, "r1")["status"] in SCORED
    assert len(_scoring(session)) == 1
    # The file goes after its job is done: the next GET queues it again.
    session.runs.store.trust_path("r1", trust_name(ClonalSettings(**A))).unlink()
    assert _clonal(client, "r1")["status"] == "pending"
    again = _scoring(session)
    assert len(again) == 2
    _wait_all(jobs, again)


def test_concurrent_gets_queue_one_job(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    session.set_clonal(ClonalSettings(**A))
    _save_run(session, "r1")
    client = _client(session)
    with ThreadPoolExecutor(8) as pool:
        statuses = list(pool.map(lambda _: _clonal(client, "r1")["status"], range(16)))
    assert statuses == ["pending"] * 16
    assert len(_scoring(session)) == 1


def test_a_failed_score_is_reported_and_not_resubmitted_on_every_get(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    session.set_clonal(ClonalSettings(**A))
    _save_run(session, "r1")
    (session.runs.store.path("r1") / "cells.lcb").unlink()
    client = _client(session)
    assert _clonal(client, "r1")["status"] == "pending"
    (job,) = _scoring(session)
    jobs.start()
    assert jobs.wait(job.job_id, timeout=WAIT).state == "failed"
    for _ in range(2):
        clonal = _clonal(client, "r1")
        assert clonal["status"] == "failed"
        assert clonal["primary"] is None
        assert clonal["notes"] == [job.error]
        assert "run r1 was replaced or is incomplete" in clonal["notes"][0]
    assert len(_scoring(session)) == 1


def test_a_cancelled_score_is_reported_as_failed(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    session.set_clonal(ClonalSettings(**A))
    _save_run(session, "r1")
    client = _client(session)
    assert _clonal(client, "r1")["status"] == "pending"
    (job,) = _scoring(session)
    jobs.cancel(job.job_id)
    clonal = _clonal(client, "r1")
    assert clonal["status"] == "failed" and clonal["notes"] == ["cancelled"]
    assert len(_scoring(session)) == 1
    # A settings PUT retries it.
    _put(client, A)
    assert _clonal(client, "r1")["status"] == "pending"
    assert len(_scoring(session)) == 2


def test_a_result_for_old_settings_is_never_served(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    client = _client(session)
    _put(client, A)
    _put(client, B)
    want = ClonalSettings(**B).model_dump(mode="json")
    pending = _clonal(client, "r1")
    assert pending["status"] == "pending" and pending["settings"] == want
    both = _scoring(session)
    assert len(both) == 2
    jobs.start()
    _wait_all(jobs, both)
    assert all(j.state == "done" for j in both), [j.error for j in both]
    clonal = _clonal(client, "r1")
    assert clonal["status"] in SCORED and clonal["settings"] == want
    trust = session.runs.store.path("r1") / "trust"
    assert (trust / f"clonal-{trust_hash(ClonalSettings(**B))}.json").exists()
    assert (trust / f"clonal-{trust_hash(ClonalSettings(**A))}.json").exists()  # stored, never served


def test_a_trust_file_that_does_not_validate_counts_as_missing(
    prepared: tuple[Path, Path],
    cache_root: Path,
    jobs: JobManager,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    settings = ClonalSettings(**A)
    _save_run(session, "r1")
    store = session.runs.store
    path = str(store.trust_path("r1", trust_name(settings)))
    # An older store version's shape, or a hand edit: valid JSON, not a ClonalAgreement.
    store.write_trust("r1", trust_name(settings), {"status": "excellent"})
    capsys.readouterr()
    session.set_clonal(settings)  # the settings PUT counts it as missing
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and path in err[0]
    assert len(_scoring(session)) == 1
    client = _client(session)
    assert _clonal(client, "r1")["status"] == "pending"  # and so does a GET
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and path in err[0]
    assert len(_scoring(session)) == 1  # the queued job is reused


def test_legacy_runs_are_never_scored(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    store = session.runs.store
    store.update_graph("r1", {**store.graph("r1"), "store_version": 1})
    session.set_clonal(ClonalSettings(**A))
    clonal = session.runs.clonal_agreement("r1", store.graph("r1"))
    assert clonal.status == "refused"
    assert clonal.notes == ["made before inference isolation; run it again"]
    assert _scoring(session) == []


def test_job_info_reports_kinds_by_name(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    client = _client(session)
    for kind in ("score_clonal", "sleep", "prepare"):
        job = jobs.submit(kind, f"k-{kind}", {})
        body = client.get(f"/api/v1/jobs/{job.job_id}").json()
        assert body["kind"] == kind


def test_an_infer_job_for_another_file_queues_nothing(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    session.set_clonal(ClonalSettings(**A))
    _save_run(session, "r1")
    mine = {"path": str(session.path), "cache_root": str(cache_root)}
    base = Job(job_id="x", kind="infer", key="k", payload=mine, result_id="r1")
    session._on_job_done(replace(base, payload={**mine, "path": "/elsewhere.h5ad"}))
    session._on_job_done(replace(base, payload={**mine, "cache_root": str(tmp_path / "other")}))
    session._on_job_done(replace(base, payload={"path": str(session.path)}))  # the default cache
    for kind in ("score_clonal", "prepare", "heatmap", "sleep"):
        session._on_job_done(replace(base, kind=kind))
    session._on_job_done(replace(base, result_id=None))
    assert _scoring(session) == []
    session._on_job_done(base)
    (job,) = _scoring(session)
    assert job.key == f"score_clonal:r1:{trust_hash(ClonalSettings(**A))}"


def test_submit_trust_jobs_queues_clonal_only_with_settings(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    assert session.runs.submit_trust_jobs("r1") == []
    assert _scoring(session) == []
    session.set_clonal(ClonalSettings(**A))
    (queued,) = _scoring(session)
    assert session.runs.submit_trust_jobs("r1") == [queued]  # the queued job, not a second one
    assert queued.key == f"score_clonal:r1:{trust_hash(ClonalSettings(**A))}"
    assert len(_scoring(session)) == 1


def test_submit_clonal_all_counts_the_runs_it_queues(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    _save_run(session, "r2", {"n_dcs": 8})
    assert session.runs.submit_clonal_all() == 0
    assert _scoring(session) == []
    session._clonal = ClonalSettings(**A)  # settings in place without set_clonal's own submit
    assert session.runs.submit_clonal_all() == 2
    hb = trust_hash(ClonalSettings(**A))
    assert {j.key for j in _scoring(session)} == {f"score_clonal:r1:{hb}", f"score_clonal:r2:{hb}"}


@pytest.mark.parametrize("second", ["put", "delete"])
def test_concurrent_settings_writes_leave_memory_disk_and_view_in_agreement(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path, second: str
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    client = _client(session)
    inside = entries = 0
    overlapped = False
    count = threading.Lock()
    both_entered = threading.Event()
    update = session.state.update

    def held_update(change: Any) -> None:
        # The first writer waits inside its state write until the second enters one (up to 1 s):
        # without serialising, the two writes overlap; with it, the wait runs out alone.
        nonlocal inside, entries, overlapped
        with count:
            inside += 1
            entries += 1
            overlapped = overlapped or inside > 1
            if entries == 2:
                both_entered.set()
        try:
            update(change)
            both_entered.wait(1.0)
        finally:
            with count:
                inside -= 1

    session.state.update = held_update  # type: ignore[method-assign]
    barrier = threading.Barrier(2, timeout=30)

    def write(which: str) -> int:
        barrier.wait()
        if which == "A":
            return client.put("/api/v1/dataset/clonal", json=A).status_code
        if second == "put":
            return client.put("/api/v1/dataset/clonal", json=B).status_code
        return client.delete("/api/v1/dataset/clonal").status_code

    with ThreadPoolExecutor(2) as pool:
        codes = list(pool.map(write, ["A", "B"]))
    assert codes == [200, 200 if second == "put" else 204]
    assert not overlapped, "two settings writes ran at once"
    memory = session.clonal_settings()
    on_disk = settings_from_state(session.state.read())
    assert memory == on_disk
    assert session.view.clonal_columns == (settings_columns(memory) if memory else ())
    clonal = _clonal(client, "r1")
    if memory is None:
        assert clonal["status"] == "not_configured"
    else:
        assert clonal["status"] == "pending" and clonal["settings"] == memory.model_dump(mode="json")
        key = f"score_clonal:r1:{trust_hash(memory)}"
        assert [j.state for j in _scoring(session) if j.key == key] == ["queued"]
