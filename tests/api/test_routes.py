"""The API against a real prepared dataset, worker process included."""

import json
from collections.abc import Iterator

import numpy as np
import pytest
from fastapi.testclient import TestClient

from lumencell.api import binary
from lumencell.api.app import create_app
from lumencell.api.session import Session
from lumencell.bench.synth import generate
from lumencell.jobs.manager import JobManager

ROOT = {"kind": "cluster", "key": "true_label", "value": "trunk_early"}


@pytest.fixture(scope="module")
def session(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Session]:
    base = tmp_path_factory.mktemp("api")
    path = base / "synthetic.h5ad"
    generate("bifurcation", "realistic", n_cells=800, n_genes=400, seed=0, store_root=True).write_h5ad(path)
    jobs = JobManager().start()
    s = Session(path, jobs, cache_root=base / "cache", state_dir=base / "state")
    yield s
    jobs.shutdown()


@pytest.fixture(scope="module")
def client(session: Session) -> TestClient:
    return TestClient(create_app(session))


def _lcb(res: object) -> binary.Decoded:
    assert res.status_code == 200, res.text  # type: ignore[attr-defined]
    assert res.headers["content-type"] == binary.MEDIA_TYPE  # type: ignore[attr-defined]
    return binary.decode(res.content)  # type: ignore[attr-defined]


# ---- before prep (tests run in file order) ---------------------------------------------


def test_dataset_before_prep(client: TestClient) -> None:
    info = client.get("/api/v1/dataset").json()
    assert info["n_obs"] == 800 and info["x_kind"] == "counts"
    assert info["prep"]["state"] in ("missing", "running")
    assert info["embeddings"] == []
    assert info["root"] is None


def test_root_routes_wait_for_prep(client: TestClient) -> None:
    res = client.get("/api/v1/root/suggestions")
    assert res.status_code == 409 and res.json()["error"]["code"] == "prep_not_ready"


def test_prepare_job_runs_and_reports_ready(client: TestClient, session: Session) -> None:
    job = client.post("/api/v1/prepare").json()
    assert job["kind"] == "prepare"
    session.jobs.wait(job["job_id"], timeout=300)
    assert client.get("/api/v1/prepare").json()["state"] == "ready"
    again = client.post("/api/v1/prepare").json()
    assert again["job_id"] != job["job_id"] or again["state"] == "done"


# ---- dataset ---------------------------------------------------------------------------


def test_dataset_after_prep(client: TestClient) -> None:
    info = client.get("/api/v1/dataset").json()
    names = {e["name"]: e["source"] for e in info["embeddings"]}
    assert names == {"umap": "computed", "dc12": "diffusion"}
    assert info["default_embedding"] == "umap"
    cols = {c["name"]: c for c in info["obs"]}
    assert cols["true_label"]["kind"] == "categorical"
    assert cols["true_pseudotime"]["kind"] == "continuous"
    assert cols["lumencell_leiden"]["categories"]
    assert info["cluster_key"] == "lumencell_leiden"


def test_embedding_payload(client: TestClient) -> None:
    coords = _lcb(client.get("/api/v1/embeddings/umap")).arrays["coords"]
    assert coords.shape == (800, 2) and coords.dtype == np.float32
    missing = client.get("/api/v1/embeddings/tsne")
    assert missing.status_code == 404 and missing.json()["error"]["code"] == "not_found"


def test_obs_payloads(client: TestClient) -> None:
    cat = _lcb(client.get("/api/v1/obs/true_label"))
    assert cat.arrays["codes"].dtype == np.int16
    assert cat.meta["kind"] == "categorical" and "trunk_early" in cat.meta["categories"]
    cont = _lcb(client.get("/api/v1/obs/true_pseudotime"))
    assert cont.arrays["values"].dtype == np.float32 and cont.meta["kind"] == "continuous"
    leiden = _lcb(client.get("/api/v1/obs/lumencell_leiden"))
    assert leiden.arrays["codes"].max() == len(leiden.meta["categories"]) - 1


def test_gene_search(client: TestClient) -> None:
    res = client.get("/api/v1/genes", params={"q": "G0001", "limit": 3}).json()
    assert res["genes"][0].startswith("G0001")
    assert len(res["genes"]) == 3 and res["truncated"] is True


def test_expression(client: TestClient, session: Session) -> None:
    decoded = _lcb(client.get("/api/v1/expression", params={"genes": "G00001,g00002"}))
    expr = decoded.arrays["expression"]
    assert expr.shape == (2, 800)
    assert decoded.meta["genes"] == ["G00001", "G00002"]
    lognorm = session.prep().artifacts["lognorm"]
    np.testing.assert_allclose(expr[0], lognorm[:, 1].toarray().ravel(), rtol=1e-6)
    unknown = client.get("/api/v1/expression", params={"genes": "NOPE"})
    assert unknown.status_code == 404 and unknown.json()["error"]["code"] == "unknown_gene"


def test_methods(client: TestClient) -> None:
    methods = {m["name"]: m for m in client.get("/api/v1/methods").json()["methods"]}
    m = methods["paga_dpt"]
    assert m["available"] is True and m["experimental"] is True
    assert "n_dcs" in m["params_schema"]["properties"]
    assert "palantir" in methods


# ---- root and runs ---------------------------------------------------------------------


def test_running_without_a_root_is_refused(client: TestClient) -> None:
    res = client.post("/api/v1/trajectories", json={"method": "paga_dpt"})
    assert res.status_code == 409
    assert res.json()["error"]["code"] == "root_required"


def test_root_suggestions_are_offered_not_applied(client: TestClient) -> None:
    suggestions = client.get("/api/v1/root/suggestions").json()["suggestions"]
    assert suggestions[0]["source"] == "file"
    assert client.get("/api/v1/root").json() == {"root": None}


def test_resolve_and_set_root(client: TestClient) -> None:
    preview = client.post("/api/v1/root/resolve", json={"root": ROOT}).json()
    assert preview["spec"] == ROOT and preview["cluster_label"] == "trunk_early"
    assert client.get("/api/v1/root").json() == {"root": None}, "resolving is only a preview"
    saved = client.put("/api/v1/root", json={"root": ROOT}).json()
    assert client.get("/api/v1/root").json()["root"] == saved
    bad = client.put("/api/v1/root", json={"root": {"kind": "cluster", "key": "true_label", "value": "nope"}})
    assert bad.status_code == 422 and bad.json()["error"]["code"] == "invalid_root"


@pytest.fixture(scope="module")
def run(client: TestClient, session: Session) -> dict[str, object]:
    created = client.post("/api/v1/trajectories", json={"method": "paga_dpt"})
    assert created.status_code == 202, created.text
    body = created.json()
    done = session.jobs.wait(body["job_id"], timeout=300)
    assert done.state == "done", done.error
    return body


def test_run_uses_the_saved_root(client: TestClient, run: dict[str, object]) -> None:
    traj = client.get(f"/api/v1/trajectories/{run['run_id']}").json()
    assert traj["state"] == "done" and traj["topology_kind"] == "tree"
    assert traj["provenance"]["root"]["spec"] == ROOT
    assert len(traj["lineages"]) >= 2


def test_repeating_a_run_returns_the_existing_one(client: TestClient, run: dict[str, object]) -> None:
    again = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": ROOT})
    assert again.status_code == 200
    assert again.json() == {"run_id": run["run_id"], "job_id": None, "state": "done"}
    listed = client.get("/api/v1/trajectories").json()["runs"]
    assert run["run_id"] in [r["run_id"] for r in listed]


def test_cells_payload(client: TestClient, run: dict[str, object]) -> None:
    cells = _lcb(client.get(f"/api/v1/trajectories/{run['run_id']}/cells")).arrays
    assert cells["pseudotime"].shape == (800,)
    assert cells["membership"].shape[1] == 800


def test_curves(client: TestClient, run: dict[str, object]) -> None:
    curves = client.get(f"/api/v1/trajectories/{run['run_id']}/curves", params={"embedding": "umap"}).json()
    assert curves["curve_source"] == "binned"
    assert all(len(c["points"]) > 10 for c in curves["curves"])


def test_rename_lineage(client: TestClient, run: dict[str, object]) -> None:
    url = f"/api/v1/trajectories/{run['run_id']}"
    renamed = client.patch(f"{url}/lineages/0", json={"name": "Erythroid"}).json()
    assert renamed["name"] == "Erythroid"
    assert client.get(url).json()["lineages"][0]["name"] == "Erythroid"
    clash = client.patch(f"{url}/lineages/1", json={"name": "Erythroid"})
    assert clash.status_code == 409 and clash.json()["error"]["code"] == "name_taken"


def test_job_routes(client: TestClient, run: dict[str, object]) -> None:
    job = client.get(f"/api/v1/jobs/{run['job_id']}").json()
    assert job["state"] == "done" and job["result_id"] == run["run_id"]
    with client.stream("GET", f"/api/v1/jobs/{run['job_id']}/events") as res:
        assert res.headers["content-type"].startswith("text/event-stream")
        events = [json.loads(line[6:]) for line in res.iter_lines() if line.startswith("data: ")]
    assert events[-1]["state"] == "done"
    assert client.delete(f"/api/v1/jobs/{run['job_id']}").json()["state"] == "done"
    assert client.get("/api/v1/jobs/nope").status_code == 404


def test_unknown_run_and_method(client: TestClient) -> None:
    assert client.get("/api/v1/trajectories/nope").status_code == 404
    res = client.post("/api/v1/trajectories", json={"method": "slingshot", "root": ROOT})
    assert res.status_code == 404 and res.json()["error"]["code"] == "unknown_method"


def test_validation_errors_use_the_error_body(client: TestClient) -> None:
    res = client.put("/api/v1/root", json={"root": {"kind": "cluster"}})
    assert res.status_code == 422
    body = res.json()
    assert body["error"]["code"] == "invalid_request"
    assert body["error"]["detail"]["errors"]


def test_trends(client: TestClient, run: dict[str, object]) -> None:
    url = f"/api/v1/trajectories/{run['run_id']}/trends"
    res = client.get(url, params={"lineage": 0, "genes": "G00001,G00002"}).json()
    assert [g["gene"] for g in res["genes"]] == ["G00001", "G00002"]
    assert len(res["x"]) == len(res["counts"]) == len(res["genes"][0]["y"]) == 50
    assert all(0.0 <= g["r2"] <= 1.0 for g in res["genes"])
    assert client.get(url, params={"lineage": 9, "genes": "G00001"}).status_code == 404


def test_heatmap(client: TestClient, run: dict[str, object]) -> None:
    url = f"/api/v1/trajectories/{run['run_id']}/heatmap"
    decoded = _lcb(client.post(url, json={"lineage": 0, "top_n": 12}))
    assert decoded.arrays["heatmap"].shape == (12, 50)
    assert len(decoded.meta["genes"]) == 12 and len(decoded.meta["bins"]) == 50
    again = client.post(url, json={"lineage": 0, "top_n": 12})
    assert again.content == client.post(url, json={"lineage": 0, "top_n": 12}).content


def test_compare_two_methods(client: TestClient, session: Session, run: dict[str, object]) -> None:
    other = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "params": {"n_dcs": 8}}).json()
    session.jobs.wait(other["job_id"], timeout=300)
    res = client.get("/api/v1/compare", params={"a": run["run_id"], "b": other["run_id"]}).json()
    assert res["spearman"] > 0.9 and res["n_common"] > 700
    assert res["lineage_matches"]
    cells = _lcb(client.get("/api/v1/compare/cells", params={"a": run["run_id"], "b": other["run_id"]}))
    assert cells.arrays["rank_difference"].shape == (800,)


def test_no_session_means_503() -> None:
    res = TestClient(create_app()).get("/api/v1/dataset")
    assert res.status_code == 503 and res.json()["error"]["code"] == "no_dataset"
