"""Every run's GET carries its stability and its tree-support verdict, healed on GET.

Each test builds its own job manager and session, as tests/api/test_cache_healing.py does. Tests
that check what a GET queues keep the manager unstarted until those checks are done.

The recorded verification failed (``VERDICT_VERIFIED = False``), so every served verdict that would
be supported or weak is capped at provisional. Tests that need the uncapped level patch
``lumencell.engine.tree_support.VERDICT_VERIFIED`` to True and say so. The calibration set no tau
(C2, C3 and C4 are excluded), so tests of the C2 and stability paths patch the thresholds
``tree_support_for_run`` uses and say so too.
"""

import shutil
import threading
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from fastapi.testclient import TestClient
from scipy import sparse

from lumencell.api.app import create_app
from lumencell.api.models import ClonalSettings
from lumencell.api.session import Session
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
from lumencell.data.loader import load
from lumencell.engine.tree_support import Thresholds
from lumencell.jobs.manager import BACKGROUND_KINDS, Job, JobManager
from lumencell.jobs.worker import HANDLERS
from lumencell.prep.cache import PrepCache
from lumencell.prep.stages import PrepParams
from lumencell.trust.stability import stability_hash, stability_name

ROOT = {"kind": "cluster", "key": "true_label", "value": "trunk_early"}
CLONES = {"clone_key": "clone_id", "donor_key": "donor"}
WAIT = 300.0
UNVERIFIED = "failed its verification"


@pytest.fixture(scope="module")
def prepared(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    base = tmp_path_factory.mktemp("robust")
    path = base / "d.h5ad"
    adata = generate("bifurcation", "realistic", n_cells=400, n_genes=300, seed=0)
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


@pytest.fixture
def verified(monkeypatch: pytest.MonkeyPatch) -> None:
    """Lift the failed-verification cap, so a test sees the level the rule itself gives."""
    monkeypatch.setattr("lumencell.engine.tree_support.VERDICT_VERIFIED", True)


def _thresholds(monkeypatch: pytest.MonkeyPatch, **taus: float | None) -> None:
    """Serve the verdict with these taus (the calibrated ones are all None, which excludes C2-C4)."""
    t = Thresholds(tau_tree=taus.get("tau_tree"), tau_lin=taus.get("tau_lin"), tau_pt=taus.get("tau_pt"))
    monkeypatch.setattr("lumencell.trust.verdict.calibrated_thresholds", lambda: t)


def _session(prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path) -> Session:
    path, _ = prepared
    s = Session(path, jobs, cache_root=cache_root, state_dir=tmp_path / "state")
    assert s.prep_ready()
    return s


def _client(session: Session) -> TestClient:
    return TestClient(create_app(session), raise_server_exceptions=False)


def _save_run(session: Session, rid: str) -> None:
    """Compute and save a run in this process, as `lumencell run` does; no job is involved."""
    payload = {
        **session._job_payload(),
        "method": "paga_dpt",
        "root": "cluster:true_label=trunk_early",
        "params": {},
        "run_id": rid,
        "label_key": session.view.label_key,
        "clonal_columns": list(session.view.clonal_columns),
    }
    assert HANDLERS["infer"](payload, lambda *_: None) == rid


def _measure(session: Session, rid: str) -> None:
    """Run the stability job's handler in this process and store its trust file."""
    payload = {**session._job_payload(), "run_id": rid, "trust_name": stability_name()}
    assert HANDLERS["stability"](payload, lambda *_: None) == rid


def _trust(client: TestClient, rid: str) -> dict[str, Any]:
    r = client.get(f"/api/v1/trajectories/{rid}")
    assert r.status_code == 200, r.text
    trust: dict[str, Any] = r.json()["trust"]
    return trust


def _check(trust: dict[str, Any], cid: str) -> dict[str, Any]:
    (check,) = [c for c in trust["tree_support"]["checks"] if c["id"] == cid]
    return check


def _stability_jobs(session: Session) -> list[Job]:
    return session.jobs.of_kind("stability")


def _after_done(jobs: JobManager, kind: str) -> threading.Event:
    """Set once a job of this kind is done and every callback registered before this one ran."""
    event = threading.Event()

    def note(job: Job) -> None:
        if job.kind == kind:
            event.set()

    jobs.add_done_callback(note)
    return event


def test_a_finished_infer_job_queues_stability_at_background_priority(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path, verified: None
) -> None:
    """With the cap lifted (``verified``): the level is then one of supported, weak or not_supported."""
    session = _session(prepared, cache_root, jobs, tmp_path)
    inferred = _after_done(jobs, "infer")  # registered after the session's own callback
    client = _client(session)
    jobs.start()
    created = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": ROOT}).json()
    rid = created["run_id"]
    assert jobs.wait(created["job_id"], timeout=WAIT).state == "done"
    assert inferred.wait(WAIT)
    (job,) = _stability_jobs(session)
    assert job.key == f"stability:{rid}:{stability_hash()}"
    assert job.payload["trust_name"] == stability_name()
    assert "stability" in BACKGROUND_KINDS
    assert jobs.wait(job.job_id, timeout=WAIT).state == "done", job.error
    trust = _trust(client, rid)
    stability = trust["stability"]
    assert stability["status"] == "ok"
    assert stability["reseeds"]["r"] == 2
    assert stability["subsamples"]["r"] == 3
    assert stability["subsamples"]["fraction"] == 0.8
    assert trust["tree_support"]["level"] in ("supported", "weak", "not_supported")
    assert _check(trust, "C5")["state"] == "not_applicable"  # no clonal settings
    assert len(_stability_jobs(session)) == 1


def test_the_served_verdict_is_capped_while_the_verification_has_failed(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    """Unpatched: the recorded verification failed, so nothing is served above provisional."""
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    _measure(session, "r1")
    trust = _trust(_client(session), "r1")
    assert trust["stability"]["status"] == "ok"
    support = trust["tree_support"]
    assert support["level"] in ("provisional", "not_supported")
    if support["level"] == "provisional":
        assert any(UNVERIFIED in r for r in support["reasons"]), support["reasons"]
    assert _stability_jobs(session) == []


def test_foreground_jobs_overtake_a_queued_stability_job(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    (stability,) = session.runs.submit_trust_jobs("r1")
    assert stability.kind == "stability"
    prepare = session.start_prep()
    with jobs._lock:
        assert jobs._take_next() == prepare.job_id
        assert jobs._take_next() == stability.job_id


def test_a_restart_heals_a_missing_stability_job(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    first = _session(prepared, cache_root, JobManager(), tmp_path)  # never started; its jobs are lost
    _save_run(first, "r1")
    session = _session(prepared, cache_root, jobs, tmp_path)
    client = _client(session)
    trust = _trust(client, "r1")
    assert trust["stability"]["status"] == "pending"
    assert trust["tree_support"]["level"] != "supported"
    (job,) = _stability_jobs(session)
    assert job.payload == {**session._job_payload(), "run_id": "r1", "trust_name": stability_name()}
    assert _trust(client, "r1")["stability"]["status"] == "pending"
    assert len(_stability_jobs(session)) == 1
    jobs.start()
    assert jobs.wait(job.job_id, timeout=WAIT).state == "done", job.error
    assert _trust(client, "r1")["stability"]["status"] == "ok"
    assert _trust(client, "r1")["stability"]["status"] == "ok"
    assert len(_stability_jobs(session)) == 1


def test_a_cancelled_stability_keeps_the_verdict_provisional_and_is_not_resubmitted(
    prepared: tuple[Path, Path],
    cache_root: Path,
    jobs: JobManager,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    verified: None,
) -> None:
    """Review Focus 4. With the cap lifted and C3/C4 given taus, so stability decides the level."""
    _thresholds(monkeypatch, tau_lin=0.5, tau_pt=0.5)
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    client = _client(session)
    assert _trust(client, "r1")["stability"]["status"] == "pending"
    (job,) = _stability_jobs(session)
    r = client.delete(f"/api/v1/jobs/{job.job_id}")
    assert r.status_code == 200, r.text
    assert r.json()["kind"] == "stability" and r.json()["state"] == "cancelled"
    for _ in range(2):
        trust = _trust(client, "r1")
        assert trust["stability"]["status"] == "failed"
        assert trust["stability"]["notes"] == ["cancelled"]
        support = trust["tree_support"]
        assert support["level"] == "provisional"
        assert any("stability" in reason for reason in support["reasons"]), support["reasons"]
        assert _check(trust, "C3")["state"] == "unavailable"
    assert len(jobs.of_kind("stability")) == 1


def test_a_cancelled_stability_is_never_served_as_supported(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    """Review Focus 4 on the served defaults: no tau and the cap still in place."""
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    client = _client(session)
    _trust(client, "r1")
    (job,) = _stability_jobs(session)
    jobs.cancel(job.job_id)
    trust = _trust(client, "r1")
    assert trust["stability"]["status"] == "failed"
    assert trust["tree_support"]["level"] in ("provisional", "not_supported")
    assert len(_stability_jobs(session)) == 1


def test_a_failed_stability_job_is_reported_with_its_error_and_not_resubmitted(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    client = _client(session)
    _trust(client, "r1")
    (job,) = _stability_jobs(session)
    (session.runs.store.path("r1") / ".complete").unlink()  # the handler then refuses the run
    jobs.start()
    assert jobs.wait(job.job_id, timeout=WAIT).state == "failed"
    session.runs.store.path("r1").joinpath(".complete").touch()
    for _ in range(2):
        stability = _trust(client, "r1")["stability"]
        assert stability["status"] == "failed"
        assert stability["notes"] == [job.error]
        assert "run r1 was replaced or is incomplete" in stability["notes"][0]
    assert len(_stability_jobs(session)) == 1


def test_a_stability_file_that_does_not_validate_counts_as_missing(
    prepared: tuple[Path, Path],
    cache_root: Path,
    jobs: JobManager,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    store = session.runs.store
    path = str(store.trust_path("r1", stability_name()))
    store.write_trust("r1", stability_name(), {"status": "excellent"})
    capsys.readouterr()
    assert _trust(_client(session), "r1")["stability"]["status"] == "pending"
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and path in err[0]
    assert len(_stability_jobs(session)) == 1


def test_legacy_runs_get_no_stability_job(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    store = session.runs.store
    store.update_graph("r1", {**store.graph("r1"), "store_version": 1})
    # The store does not serve such a run at all; the lookup itself still refuses it.
    stability = session.runs.stability("r1", store.graph("r1"))
    assert stability.status == "not_available"
    assert stability.notes == ["made before inference isolation; run it again"]
    assert _stability_jobs(session) == []


def test_dataset_info_carries_tree_likeness(
    prepared: tuple[Path, Path], tmp_path: Path, jobs: JobManager
) -> None:
    path, _ = prepared
    empty = Session(path, jobs, cache_root=tmp_path / "empty", state_dir=tmp_path / "s1")
    assert _client(empty).get("/api/v1/dataset").json()["tree_likeness"] is None
    assert empty.view.tree_likeness() is None
    own = tmp_path / "cache"
    shutil.copytree(prepared[1], own, ignore=shutil.ignore_patterns("runs"))
    ready = _session(prepared, own, jobs, tmp_path)
    tree = _client(ready).get("/api/v1/dataset").json()["tree_likeness"]
    assert 0 < tree <= 1
    assert ready.view.tree_likeness() == tree


def test_clonal_support_enters_the_verdict(
    prepared: tuple[Path, Path], cache_root: Path, jobs: JobManager, tmp_path: Path
) -> None:
    session = _session(prepared, cache_root, jobs, tmp_path)
    settings = ClonalSettings(**CLONES)
    session.set_clonal(settings)
    _save_run(session, "r1")
    _measure(session, "r1")
    assert HANDLERS["score_clonal"](session.runs._clonal_payload("r1", settings), lambda *_: None) == "r1"
    trust = _trust(_client(session), "r1")
    assert trust["stability"]["status"] == "ok"
    assert trust["clonal"]["status"] not in ("pending", "not_configured")
    assert _check(trust, "C5")["state"] not in ("not_applicable", "pending")


def test_a_dataset_without_paga_edges_is_not_supported(
    prepared: tuple[Path, Path],
    cache_root: Path,
    jobs: JobManager,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Review Focus 5. C2 is excluded under the calibrated (None) tau, so it is given one here."""
    session = _session(prepared, cache_root, jobs, tmp_path)
    _save_run(session, "r1")
    prep = session.view.loaded_prep()
    assert prep is not None
    prep.artifacts["paga_connectivities"] = sparse.csr_matrix((1, 1))
    prep.artifacts["paga_tree"] = sparse.csr_matrix((1, 1))
    client = _client(session)
    c2 = _check(_trust(client, "r1"), "C2")
    assert c2["state"] == "excluded" and c2["value"] is None  # the served defaults
    _thresholds(monkeypatch, tau_tree=0.5)
    trust = _trust(client, "r1")
    assert trust["tree_support"]["level"] == "not_supported"
    c2 = _check(trust, "C2")
    assert c2["state"] == "fail" and c2["value"] is None and c2["reason"]
    assert client.get("/api/v1/dataset").json()["tree_likeness"] is None
