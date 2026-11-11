"""A damaged stage in the prep cache goes back to the worker instead of failing every route.

Each test's job manager is started only after the assertions about the damaged state,
so the heal job cannot finish early and race them.
"""

from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any
from unittest import mock

import pytest
from fastapi.testclient import TestClient

from lumencell.api.app import create_app
from lumencell.api.session import Session
from lumencell.bench.synth import generate
from lumencell.data.loader import load
from lumencell.jobs.manager import Job, JobManager
from lumencell.prep import cache as cache_mod
from lumencell.prep.cache import PrepCache
from lumencell.prep.stages import PrepParams

ROOT = {"kind": "cluster", "key": "true_label", "value": "trunk_early"}


@pytest.fixture(scope="module")
def prepared(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    base = tmp_path_factory.mktemp("heal")
    path = base / "d.h5ad"
    generate("bifurcation", "realistic", n_cells=400, n_genes=300, seed=0).write_h5ad(path)
    PrepCache(load(path), base / "cache").run(PrepParams())
    return path, base / "cache"


@pytest.fixture
def jobs() -> Iterator[JobManager]:
    manager = JobManager()  # not started: jobs wait in the queue until the test starts it
    yield manager
    manager.shutdown()


@pytest.fixture
def session(prepared: tuple[Path, Path], jobs: JobManager, tmp_path: Path) -> Session:
    path, cache_root = prepared
    s = Session(path, jobs, cache_root=cache_root, state_dir=tmp_path / "state")
    # `lumencell launch` asks this at startup; from then on the answer is kept.
    assert s.prep_ready()
    return s


def _client(session: Session) -> TestClient:
    return TestClient(create_app(session), raise_server_exceptions=False)


def _damage(session: Session, how: str) -> None:
    key = session.cache.keys(session.params)["pca"]
    artifact = session.cache.stages_dir / f"pca-{key}" / "pca.npy"
    if how == "truncate":
        artifact.write_bytes(artifact.read_bytes()[: artifact.stat().st_size // 2])
    else:
        artifact.unlink()


def _prepare_jobs(session: Session) -> list[Job]:
    return session.jobs.of_kind("prepare")


def _heal(session: Session) -> None:
    (job,) = _prepare_jobs(session)
    session.jobs.start()
    assert session.jobs.wait(job.job_id, timeout=120).state == "done"


@pytest.mark.parametrize("how", ["truncate", "delete"])
def test_a_damaged_stage_is_sent_back_to_the_worker(
    session: Session, how: str, capsys: pytest.CaptureFixture[str]
) -> None:
    client = _client(session)
    _damage(session, how)

    res = client.get("/api/v1/root/suggestions")
    assert res.status_code == 409, res.text
    assert res.json()["error"]["code"] == "prep_not_ready"
    assert len(_prepare_jobs(session)) == 1, "the server starts the prepare job itself"
    err = capsys.readouterr().err
    assert "pca" in err and "preparing the dataset again" in err, "the server says what it did"

    _heal(session)
    assert client.get("/api/v1/root/suggestions").status_code == 200


def test_every_data_route_waits_while_a_damaged_stage_heals(session: Session) -> None:
    client = _client(session)
    _damage(session, "truncate")
    requests: list[tuple[str, str, dict[str, Any]]] = [
        ("GET", "/api/v1/obs/lumencell_leiden", {}),
        ("GET", "/api/v1/expression", {"params": {"genes": "G00001"}}),
        ("GET", "/api/v1/root/suggestions", {}),
        ("POST", "/api/v1/root/resolve", {"json": {"root": ROOT}}),
        ("POST", "/api/v1/trajectories", {"json": {"method": "paga_dpt", "root": ROOT}}),
    ]
    for method, url, kw in requests:
        res = client.request(method, url, **kw)
        assert res.status_code == 409, (url, res.text)
        assert res.json()["error"]["code"] == "prep_not_ready", url
    # As on a dataset that was never prepared: the computed embeddings are not there yet.
    assert client.get("/api/v1/embeddings/dc12").status_code == 404
    assert len(_prepare_jobs(session)) == 1

    _heal(session)
    for method, url, kw in requests[:4]:
        assert client.request(method, url, **kw).status_code == 200, url


def test_the_dataset_route_reports_a_damaged_stage_as_preparing(session: Session) -> None:
    client = _client(session)
    _damage(session, "truncate")

    res = client.get("/api/v1/dataset")
    assert res.status_code == 200, res.text
    assert res.json()["prep"]["state"] == "running"
    assert res.json()["cluster_key"] is None

    _heal(session)
    assert client.get("/api/v1/dataset").json()["prep"]["state"] == "ready"


def test_concurrent_requests_on_a_damaged_stage_start_one_job(session: Session) -> None:
    client = _client(session)
    _damage(session, "truncate")

    with ThreadPoolExecutor(4) as pool:
        codes = list(pool.map(lambda _: client.get("/api/v1/root/suggestions").status_code, range(4)))
    assert codes == [409] * 4
    assert len(_prepare_jobs(session)) == 1

    _heal(session)
    assert client.get("/api/v1/root/suggestions").status_code == 200


def test_a_stage_that_stays_unreadable_is_not_prepared_again_and_again(
    session: Session, capsys: pytest.CaptureFixture[str]
) -> None:
    client = _client(session)
    _damage(session, "truncate")
    assert client.get("/api/v1/root/suggestions").status_code == 409
    _heal(session)

    # The stage was recomputed, and reading it still fails (say the disk keeps cutting files short).
    real = cache_mod._load_artifact

    def unreadable_pca(directory: Path, name: str) -> object:
        if name == "pca":
            raise ValueError("cut short again")
        return real(directory, name)

    capsys.readouterr()
    with mock.patch.object(cache_mod, "_load_artifact", unreadable_pca):
        res = client.get("/api/v1/root/suggestions")
    assert res.status_code == 409 and res.json()["error"]["code"] == "prep_not_ready"
    assert len(_prepare_jobs(session)) == 1, "no second automatic prepare job"
    assert "not trying again" in capsys.readouterr().err
    assert client.get("/api/v1/prepare").json()["state"] == "missing", "the app offers to prepare it"
