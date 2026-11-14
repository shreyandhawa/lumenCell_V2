"""The API against a real prepared dataset, worker process included."""

import json
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
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
    adata = generate("bifurcation", "realistic", n_cells=800, n_genes=400, seed=0, store_root=True)
    adata.obs["clone_id"] = adata.obs["true_lineage"].copy()  # scoring-only: no run may use it
    adata.write_h5ad(path)
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
    lognorm = session.view.prep().artifacts["lognorm"]
    np.testing.assert_allclose(expr[0], lognorm[:, 1].toarray().ravel(), rtol=1e-6)
    unknown = client.get("/api/v1/expression", params={"genes": "NOPE"})
    assert unknown.status_code == 404 and unknown.json()["error"]["code"] == "unknown_gene"


def _matrices_held_by(root: object, shape: tuple[int, int]) -> dict[int, object]:
    """Every sparse matrix of this shape reachable from root, keyed by its data buffer."""
    import gc
    import types

    from scipy import sparse

    found: dict[int, object] = {}
    seen: set[int] = set()
    stack = [root]
    while stack:
        obj = stack.pop()
        if id(obj) in seen:
            continue
        seen.add(id(obj))
        if sparse.issparse(obj):
            if obj.shape == shape:  # type: ignore[attr-defined]
                found[obj.data.__array_interface__["data"][0]] = obj  # type: ignore[attr-defined]
            continue
        skip = (np.ndarray, str, bytes, type, types.ModuleType, types.FunctionType, types.BuiltinFunctionType)
        if not isinstance(obj, skip):
            stack.extend(gc.get_referents(obj))
    return found


def test_the_server_holds_expression_once(client: TestClient, session: Session) -> None:
    assert client.get("/api/v1/expression", params={"genes": "G00001"}).status_code == 200
    assert client.get("/api/v1/root/suggestions").status_code == 200
    held = _matrices_held_by(session, (session.ds.n_obs, session.ds.n_vars))
    assert len(held) == 1, f"{len(held)} copies of the expression matrix"
    (matrix,) = held.values()
    assert matrix.format == "csc" and matrix.dtype == np.float32  # type: ignore[attr-defined]


def test_the_server_still_offers_the_gene_count_root(client: TestClient) -> None:
    # It loads the file without its matrix, keeping only the genes detected per cell.
    sources = [s["source"] for s in client.get("/api/v1/root/suggestions").json()["suggestions"]]
    assert "gene_count" in sources


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
    blank = client.patch(f"{url}/lineages/1", json={"name": "   "})
    assert blank.status_code == 422 and blank.json()["error"]["code"] == "invalid_request"


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
    # Run ids name folders in the cache; a path is never looked up as one.
    for a in ("../stages", "..", "a/b"):
        res = client.get("/api/v1/compare", params={"a": a, "b": a})
        assert res.status_code == 404 and res.json()["error"]["code"] == "not_found", a
    res = client.post("/api/v1/trajectories", json={"method": "slingshot", "root": ROOT})
    assert res.status_code == 404 and res.json()["error"]["code"] == "unknown_method"


def test_validation_errors_use_the_error_body(client: TestClient) -> None:
    res = client.put("/api/v1/root", json={"root": {"kind": "cluster"}})
    assert res.status_code == 422
    body = res.json()
    assert body["error"]["code"] == "invalid_request"
    assert body["error"]["detail"]["errors"]


def test_invalid_method_params_are_a_422_not_a_500(client: TestClient) -> None:
    for params in ({"n_dcs": "many"}, {"no_such_param": 1}):
        res = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": ROOT, "params": params})
        assert res.status_code == 422, res.text
        body = res.json()
        assert body["error"]["code"] == "invalid_params"
        assert body["error"]["detail"]["errors"]


def test_a_root_that_cannot_be_written_as_text_is_a_422(client: TestClient) -> None:
    for root in ({"kind": "cell", "obs_name": ""}, {"kind": "cluster", "key": "", "value": "x"}):
        res = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": root})
        assert res.status_code == 422, res.text
        assert res.json()["error"]["code"] == "invalid_root"


def test_a_root_is_never_changed_on_its_way_to_the_worker(client: TestClient) -> None:
    # The trailing space used to be stripped on the way, so this started a run from "trunk_early".
    root = {"kind": "cluster", "key": "true_label", "value": "trunk_early "}
    res = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": root})
    assert res.status_code == 422, res.text
    assert res.json()["error"]["code"] == "invalid_root"
    assert "no cell has true_label = 'trunk_early '" in res.json()["error"]["message"], "refused at resolve"


def _runs_without_a_server(
    columns: list[str], clonal_columns: tuple[str, ...] = (), label_key: str = "l"
) -> tuple[object, list[dict[str, object]], list[object]]:
    """Runs over a stand-in dataset view: (runs, submitted job payloads, roots resolved)."""
    from types import SimpleNamespace

    from lumencell.api.runs import Runs
    from lumencell.clonal.isolation import scoring_only_columns

    submitted: list[dict[str, object]] = []
    resolved: list[object] = []

    def submit(kind: str, key: str, payload: dict[str, object]) -> SimpleNamespace:
        submitted.append(payload)
        return SimpleNamespace(job_id="j", state="queued")

    view = SimpleNamespace(
        prep_key="p",
        label_key=label_key,
        ds=SimpleNamespace(obs=pd.DataFrame(columns=columns)),
        clonal_columns=clonal_columns,
        ctx=lambda: SimpleNamespace(scoring_only=scoring_only_columns(columns, clonal_columns)),
    )
    runs = Runs(
        view,  # type: ignore[arg-type]
        SimpleNamespace(exists=lambda rid: False),  # type: ignore[arg-type]
        SimpleNamespace(submit=submit),  # type: ignore[arg-type]
        SimpleNamespace(read=lambda: {}),  # type: ignore[arg-type]
        dict,
        resolved.append,  # type: ignore[arg-type]
    )
    return runs, submitted, resolved


def test_a_root_the_plain_text_cannot_carry_reaches_the_worker_unchanged() -> None:
    from lumencell.api import models
    from lumencell.engine.roots import ClusterRoot, parse_root

    runs, submitted, resolved = _runs_without_a_server(["batch=1"])
    root = models.ClusterRoot(kind="cluster", key="batch=1", value="A ")
    runs.create_run(models.TrajectoryRequest(method="paga_dpt", root=root))  # type: ignore[attr-defined]
    assert resolved == [ClusterRoot("batch=1", "A ")]
    assert parse_root(str(submitted[0]["root"])) == ClusterRoot("batch=1", "A ")
    assert submitted[0]["clonal_columns"] == []

    named, sent, _ = _runs_without_a_server(["batch=1", "donor"], ("donor",))
    named.create_run(models.TrajectoryRequest(method="paga_dpt", root=root))  # type: ignore[attr-defined]
    assert sent[0]["clonal_columns"] == ["donor"], "the view's clonal columns reach the worker"


# ---- clone columns are scoring-only ------------------------------------------------------


def test_create_run_refuses_a_scoring_only_root_or_cluster_key(
    client: TestClient, session: Session, run: dict[str, object]
) -> None:
    assert "clone_id" in session.view.ds.obs.columns
    for body in (
        {"method": "paga_dpt", "root": {"kind": "cluster", "key": "clone_id", "value": "A"}},
        {"method": "paga_dpt", "root": ROOT, "params": {"cluster_key": "clone_id"}},
    ):
        res = client.post("/api/v1/trajectories", json=body)
        assert res.status_code == 422, res.text
        error = res.json()["error"]
        assert error["code"] == "scoring_only_column"
        assert error["detail"]["columns"] == ["clone_id"]
        assert "'clone_id' holds clone or receptor data" in error["message"]
    # The run made with clone_id in the file was steered by other columns only, and says so.
    graph = session.runs.store.graph(str(run["run_id"]))
    assert graph["store_version"] == 2 and graph["inference_isolation"] == 1
    assert graph["prep_key"] == session.prep_key
    assert "clone_id" not in graph["obs_keys"] and "true_label" in graph["obs_keys"]


def test_every_scoring_only_column_is_named_once_and_nothing_starts() -> None:
    from lumencell.api import models
    from lumencell.api.errors import ApiError

    runs, submitted, resolved = _runs_without_a_server(["clone_id", "donor", "cell_type"], ("donor",))
    root = models.ClusterRoot(kind="cluster", key="clone_id", value="A")
    with pytest.raises(ApiError) as caught:
        runs.create_run(  # type: ignore[attr-defined]
            models.TrajectoryRequest(method="paga_dpt", root=root, params={"cluster_key": "donor"})
        )
    assert caught.value.status == 422 and caught.value.code == "scoring_only_column"
    assert caught.value.detail == {"columns": ["clone_id", "donor"]}
    assert caught.value.message.endswith("it cannot be the root column")
    assert submitted == [] and resolved == [], "refused before the root is resolved"

    with pytest.raises(ApiError) as caught:
        runs.create_run(  # type: ignore[attr-defined]
            models.TrajectoryRequest(method="paga_dpt", root=root, params={"cluster_key": "clone_id"})
        )
    assert caught.value.detail == {"columns": ["clone_id"]}, "each column is named once"

    donor_root = models.ClusterRoot(kind="cluster", key="donor", value="d1")
    with pytest.raises(ApiError) as caught:
        runs.create_run(  # type: ignore[attr-defined]
            models.TrajectoryRequest(method="paga_dpt", root=donor_root, params={"cluster_key": "donor"})
        )
    assert caught.value.detail == {"columns": ["donor"]}
    assert caught.value.message == (
        "'donor' is set as a clone or donor column, which lumenCell uses only to score runs; "
        "it cannot be the root column"
    )
    with pytest.raises(ApiError, match="it cannot be the cluster column"):
        runs.create_run(  # type: ignore[attr-defined]
            models.TrajectoryRequest(
                method="paga_dpt",
                root=models.ClusterRoot(kind="cluster", key="cell_type", value="t"),
                params={"cluster_key": "donor"},
            )
        )
    labelled, _, _ = _runs_without_a_server(["clone_id", "cell_type"], label_key="clone_id")
    with pytest.raises(ApiError, match="it cannot be the label column"):
        labelled.create_run(  # type: ignore[attr-defined]
            models.TrajectoryRequest(
                method="paga_dpt", root=models.ClusterRoot(kind="cluster", key="cell_type", value="t")
            )
        )


def test_unexpected_errors_use_the_error_body() -> None:
    class Broken:
        def dataset_info(self) -> None:
            raise RuntimeError("secret detail")

    client = TestClient(create_app(Broken()), raise_server_exceptions=False)  # type: ignore[arg-type]
    res = client.get("/api/v1/dataset")
    assert res.status_code == 500
    assert res.json()["error"]["code"] == "internal_error"
    assert "secret detail" not in res.text


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


def test_state_changes_from_many_threads_are_all_kept(session: Session) -> None:
    def add(i: int) -> None:
        session.state.update(lambda state: state.setdefault("probe", {}).update({str(i): i}))

    with ThreadPoolExecutor(8) as pool:
        list(pool.map(add, range(40)))
    assert len(session.state.read()["probe"]) == 40
    session.state.update(lambda state: state.pop("probe"))
    leftovers = [p.name for p in session.state.path.parent.iterdir() if p.suffix == ".tmp"]
    assert leftovers == [], "no temporary files are left behind"


def test_a_corrupt_state_file_is_set_aside(client: TestClient, session: Session) -> None:
    saved = session.state.path.read_text(encoding="utf-8")
    try:
        session.state.path.write_text("{broken", encoding="utf-8")
        assert client.get("/api/v1/root").json() == {"root": None}
        assert client.get("/api/v1/dataset").status_code == 200
        assert session.state.path.with_suffix(".corrupt").read_text(encoding="utf-8") == "{broken"
    finally:
        session.state.path.write_text(saved, encoding="utf-8")


def test_nullable_boolean_columns_are_listed_and_served_alike(tmp_path: Path) -> None:
    # pandas' nullable boolean is numeric to pandas; it used to be listed as
    # continuous by /dataset while /obs served it as category codes.
    adata = generate("linear", "easy", n_cells=120, n_genes=60, seed=0)
    flags = pd.array(np.arange(adata.n_obs) % 2 == 0, dtype="boolean")
    flags[0] = pd.NA
    adata.obs["flag"] = flags
    path = tmp_path / "flags.h5ad"
    adata.write_h5ad(path)
    session = Session(path, JobManager(), cache_root=tmp_path / "cache", state_dir=tmp_path / "state")
    listed = {c.name: c for c in session.dataset_info().obs}["flag"]
    served = binary.decode(session.view.obs_payload("flag"))
    assert listed.kind == served.meta["kind"] == "boolean"
    assert listed.categories == served.meta["categories"] == ["False", "True"]
    assert served.arrays["codes"][:3].tolist() == [-1, 0, 1]


def test_no_session_means_503() -> None:
    res = TestClient(create_app()).get("/api/v1/dataset")
    assert res.status_code == 503 and res.json()["error"]["code"] == "no_dataset"


def test_a_finished_run_is_read_from_disk_once(
    session: Session, run: dict[str, object], monkeypatch: pytest.MonkeyPatch
) -> None:
    rid = str(run["run_id"])
    embedding = session.view.default_embedding()
    assert embedding is not None
    session.runs._loaded.clear()
    reads: list[str] = []
    real = session.runs.store.load
    monkeypatch.setattr(session.runs.store, "load", lambda r: (reads.append(r), real(r))[1])
    session.runs.run_curves(rid, embedding)
    session.runs.trends(rid, 0, ["G00001"])
    session.runs.compare_cells(rid, rid)
    assert reads == [rid]


def test_only_the_most_recent_runs_are_kept_in_memory(
    session: Session, monkeypatch: pytest.MonkeyPatch
) -> None:
    from lumencell.api.runs import MAX_RUNS_IN_MEMORY

    monkeypatch.setattr(session.runs.store, "load", lambda r: object())
    session.runs._loaded.clear()
    try:
        ids = [f"r{i}" for i in range(MAX_RUNS_IN_MEMORY + 3)]
        for rid in ids:
            session.runs._load(rid)
        session.runs._load(ids[-2])  # used again, so it is now the newest
        assert list(session.runs._loaded) == [*ids[-MAX_RUNS_IN_MEMORY:-2], ids[-1], ids[-2]]
    finally:
        session.runs._loaded.clear()


def test_what_depends_only_on_the_file_is_worked_out_once(
    client: TestClient, session: Session, monkeypatch: pytest.MonkeyPatch
) -> None:
    assert session.view.prep_ready()
    first, again = session.view.obs_columns(), session.view.obs_columns()
    # Every column from the file is the same object; only the Leiden column is added per call.
    assert len(first) > 1 and all(a is b for a, b in zip(first[:-1], again[:-1], strict=True))
    assert session.view._gene_lookup is session.view._gene_lookup

    # Once prep is ready the answer is kept, so later requests do not go back to the cache on disk.
    def no_more_checks(_params: object) -> dict[str, str]:
        raise AssertionError("the cache was checked again")

    monkeypatch.setattr(session.view.cache, "status", no_more_checks)
    assert session.view.prep_ready()
    assert client.get("/api/v1/embeddings/dc12").status_code == 200


def test_the_view_keeps_named_columns_out_of_label_detection(tmp_path: Path) -> None:
    from lumencell.api.dataset_view import DatasetView
    from lumencell.data.loader import load
    from lumencell.prep.cache import PrepCache
    from lumencell.prep.stages import PrepParams

    adata = generate("linear", "easy", n_cells=120, n_genes=60, seed=0)
    adata.obs["cell_type"] = adata.obs["true_label"].copy()
    path = tmp_path / "ct.h5ad"
    adata.write_h5ad(path)
    ds = load(path, keep_matrix=False)
    cache = PrepCache(ds, tmp_path / "cache")
    assert DatasetView(path, ds, cache, PrepParams()).label_key == "cell_type"
    named = DatasetView(path, ds, cache, PrepParams(), clonal_columns=("cell_type",))
    assert named.label_key != "cell_type"


def test_the_views_context_holds_its_named_columns(session: Session) -> None:
    from lumencell.api.dataset_view import DatasetView

    view = DatasetView(
        session.path, session.ds, session.cache, session.params, clonal_columns=("true_segment",)
    )
    assert {"true_segment", "clone_id"} <= view.ctx().scoring_only
    assert "true_segment" not in view.ctx().obs.columns
