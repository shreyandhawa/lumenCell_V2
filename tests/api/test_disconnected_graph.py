"""A run on a disconnected kNN graph fails as a job with the clear message, not an internal error."""

from collections.abc import Iterator

import pytest
from fastapi.testclient import TestClient

from lumencell.api.app import create_app
from lumencell.api.session import Session
from lumencell.bench.synth_controls import control_root_cell, generate_control
from lumencell.jobs.manager import JobManager

pytest.importorskip("scFates")


@pytest.fixture(scope="module")
def split(tmp_path_factory: pytest.TempPathFactory) -> Iterator[tuple[TestClient, Session, str]]:
    base = tmp_path_factory.mktemp("split")
    adata = generate_control("two_groups", n_cells=1000, n_genes=500, seed=0)
    adata.write_h5ad(base / "two_groups.h5ad")
    jobs = JobManager().start()
    session = Session(base / "two_groups.h5ad", jobs, cache_root=base / "cache", state_dir=base / "state")
    client = TestClient(create_app(session))
    prep = client.post("/api/v1/prepare").json()
    jobs.wait(prep["job_id"], timeout=300)
    yield client, session, control_root_cell(adata)
    jobs.shutdown()


def test_scfates_run_fails_with_the_disconnected_graph_message(
    split: tuple[TestClient, Session, str],
) -> None:
    client, session, root = split
    res = client.post(
        "/api/v1/trajectories", json={"method": "scfates_ppt", "root": {"kind": "cell", "obs_name": root}}
    )
    assert res.status_code == 202, res.text
    done = session.jobs.wait(res.json()["job_id"], timeout=300)
    assert done.state == "failed"
    assert done.error is not None
    assert done.error.startswith("lumencell.engine.graph.DisconnectedGraphError:")
    assert "2 disconnected parts" in done.error and "Subset" in done.error and root not in done.error
    listed = client.get("/api/v1/trajectories").json()["runs"]
    assert [r["error"] for r in listed if r["run_id"] == res.json()["run_id"]] == [done.error]
