"""Dataset-level clonal settings: the PUT/DELETE routes, what they change, and what they refuse.

The job manager is never started: a run is only ever created, never executed.
"""

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from lumencell.api import binary
from lumencell.api.app import create_app
from lumencell.api.models import ClonalSettings
from lumencell.api.session import Session
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
from lumencell.data.loader import load
from lumencell.engine.store import run_id
from lumencell.jobs.manager import JobManager
from lumencell.prep.cache import PrepCache
from lumencell.prep.stages import PrepParams
from lumencell.trust.clonal import check_settings
from lumencell.trust.constants import MAX_DONORS

N_CELLS = 1200  # more than the 1,001 distinct values the high-cardinality columns need
SETTINGS = {"clone_key": "clone_id", "donor_key": "donor"}
TRUNK = {"kind": "cluster", "key": "true_segment", "value": "trunk"}


@pytest.fixture(scope="module")
def prepared(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    base = tmp_path_factory.mktemp("clonal")
    path = base / "d.h5ad"
    adata = add_clones(
        generate("bifurcation", "realistic", n_cells=N_CELLS, n_genes=300, seed=0), DEFAULT_LAYER, 0
    )
    # Plain string columns, so the listing and serving rules meet real high cardinality.
    adata.obs["barcode2"] = [f"bc{i}" for i in range(MAX_DONORS + 1)] + ["bc0"] * (N_CELLS - MAX_DONORS - 1)
    adata.obs["many_donors"] = [f"m{i}" for i in range(MAX_DONORS + 1)] + ["m0"] * (N_CELLS - MAX_DONORS - 1)
    # Categoricals with categories that no cell uses.
    adata.obs["cat_clone"] = pd.Categorical(
        adata.obs["clone_id"].astype(object), categories=[f"c{i}" for i in range(500)] + ["unused_clone"]
    )
    adata.obs["cat_donor"] = pd.Categorical(
        adata.obs["donor"].astype(object), categories=["d0", "d1", "d2", "unused_a", "unused_b"]
    )
    adata.write_h5ad(path)
    PrepCache(load(path), base / "cache").run(PrepParams())
    return path, base / "cache"


@pytest.fixture
def jobs() -> Iterator[JobManager]:
    manager = JobManager()
    yield manager
    manager.shutdown()


def _session(prepared: tuple[Path, Path], jobs: JobManager, state: Path) -> Session:
    path, cache_root = prepared
    s = Session(path, jobs, cache_root=cache_root, state_dir=state)
    assert s.prep_ready()
    return s


@pytest.fixture
def session(prepared: tuple[Path, Path], jobs: JobManager, tmp_path: Path) -> Session:
    return _session(prepared, jobs, tmp_path / "state")


@pytest.fixture
def client(session: Session) -> TestClient:
    return TestClient(create_app(session), raise_server_exceptions=False)


def _put(client: TestClient, body: dict[str, Any]) -> Any:
    return client.put("/api/v1/dataset/clonal", json=body)


def _obs(client: TestClient, column: str) -> dict[str, Any]:
    res = client.get(f"/api/v1/obs/{column}")
    assert res.status_code == 200
    return binary.decode(res.content).meta


def test_put_returns_the_summary(client: TestClient, session: Session) -> None:
    res = _put(client, SETTINGS)
    assert res.status_code == 200
    body = res.json()
    obs = session.ds.obs
    labelled = obs["clone_id"].notna()
    assert body["settings"]["clone_key"] == "clone_id"
    assert body["n_donors"] == 3
    assert body["n_cells_labelled"] == int(labelled.sum())
    pairs = obs.loc[labelled, ["donor", "clone_id"]].astype(str).drop_duplicates()
    assert body["n_clones"] == len(pairs)
    assert body["size_max"] >= body["size_p99"] >= body["size_p90"] >= body["size_median"] >= 1


@pytest.mark.parametrize(
    "body",
    [
        {"clone_key": "clone_id", "donor_key": "donor", "single_donor": True},
        {"clone_key": "clone_id"},
    ],
)
def test_donor_or_single_donor_exactly_one(client: TestClient, body: dict[str, Any]) -> None:
    res = _put(client, body)
    assert res.status_code == 422 and res.json()["error"]["code"] == "invalid_request"


@pytest.mark.parametrize(
    "body",
    [
        {"clone_key": "nope", "donor_key": "donor"},
        {"clone_key": "clone_id", "donor_key": "clone_id"},
        {"clone_key": "true_pseudotime", "donor_key": "donor"},
        {"clone_key": "clone_id", "donor_key": "many_donors"},
        {"clone_key": "clone_id", "donor_key": "donor", "strata_keys": ["nope"]},
        {"clone_key": "clone_id", "donor_key": "donor", "strata_keys": ["clone_id"]},
    ],
)
def test_invalid_settings_answer_422_with_the_error_body(client: TestClient, body: dict[str, Any]) -> None:
    res = _put(client, body)
    assert res.status_code == 422
    error = res.json()["error"]
    assert error["code"] == "invalid_clonal_settings"
    assert error["detail"]["problems"]


def test_dataset_info_lists_candidates_and_never_selects_one(client: TestClient) -> None:
    info = client.get("/api/v1/dataset").json()
    assert info["clonal"] is None
    assert info["clonal_candidates"] == ["clone_id"]
    assert _put(client, SETTINGS).status_code == 200
    assert client.get("/api/v1/dataset").json()["clonal"] == {
        **SETTINGS,
        "single_donor": False,
        "strata_keys": [],
    }
    assert client.delete("/api/v1/dataset/clonal").status_code == 204
    assert client.get("/api/v1/dataset").json()["clonal"] is None


def test_settings_survive_a_new_session(
    prepared: tuple[Path, Path], jobs: JobManager, tmp_path: Path
) -> None:
    first = _session(prepared, jobs, tmp_path / "state")
    assert _put(TestClient(create_app(first)), SETTINGS).status_code == 200
    second = _session(prepared, jobs, tmp_path / "state")
    assert second.clonal_settings() == ClonalSettings(**SETTINGS)
    assert second.view.clonal_columns == ("clone_id", "donor")
    assert second.dataset_info().clonal == ClonalSettings(**SETTINGS)


def test_put_rebuilds_the_inference_context(client: TestClient, session: Session) -> None:
    assert "donor" in session.view.ctx().obs  # not scoring-only until it is named
    assert _put(client, SETTINGS).status_code == 200
    ctx = session.view.ctx()
    assert ctx.scoring_only >= {"clone_id", "donor"}
    assert "donor" not in ctx.obs


def test_settings_that_capture_the_saved_root_are_refused(client: TestClient, session: Session) -> None:
    assert client.put("/api/v1/root", json={"root": TRUNK}).status_code == 200
    res = _put(client, {"clone_key": "clone_id", "donor_key": "true_segment"})
    assert res.status_code == 422 and res.json()["error"]["code"] == "scoring_only_column"
    assert session.clonal_settings() is None
    assert session.view.clonal_columns == ()


def test_strata_columns_stay_visible_to_methods(client: TestClient, session: Session) -> None:
    res = _put(client, {**SETTINGS, "strata_keys": ["true_segment"]})
    assert res.status_code == 200
    ctx = session.view.ctx()
    assert "true_segment" in ctx.obs and "true_segment" not in ctx.scoring_only
    assert client.put("/api/v1/root", json={"root": TRUNK}).status_code == 200
    assert client.post("/api/v1/trajectories", json={"method": "paga_dpt"}).status_code == 202


def test_create_run_refuses_a_settings_column(client: TestClient) -> None:
    assert _put(client, SETTINGS).status_code == 200
    root = {"kind": "cluster", "key": "donor", "value": "d0"}
    res = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": root})
    assert res.status_code == 422 and res.json()["error"]["code"] == "scoring_only_column"


def test_listing_and_serving_agree_on_every_column_kind(client: TestClient) -> None:
    listed = {c["name"]: c["kind"] for c in client.get("/api/v1/dataset").json()["obs"]}
    assert listed["barcode2"] == "string"
    for name, kind in listed.items():
        assert _obs(client, name)["kind"] == kind, name
    assert _obs(client, "barcode2")["categories"], "codes and categories are still sent for a string column"


def test_unused_categories_are_not_counted(client: TestClient, session: Session) -> None:
    body = _put(client, {"clone_key": "cat_clone", "donor_key": "cat_donor"}).json()
    plain = _put(client, SETTINGS).json()
    assert body["n_donors"] == 3 == plain["n_donors"]
    assert body["n_clones"] == plain["n_clones"]
    # More donor categories than the limit, three in use: the limit counts the three.
    obs = session.ds.obs.copy()
    obs["cat_donor"] = pd.Categorical(
        obs["donor"].astype(object), categories=[f"d{i}" for i in range(MAX_DONORS + 50)]
    )
    assert check_settings(obs, ClonalSettings(clone_key="cat_clone", donor_key="cat_donor")) == []


def test_run_id_depends_on_saved_settings(client: TestClient, session: Session) -> None:
    def created() -> str:
        res = client.post("/api/v1/trajectories", json={"method": "paga_dpt", "root": TRUNK})
        assert res.status_code == 202, res.text
        return str(res.json()["run_id"])

    before = created()
    assert created() == before, "no settings: the id is stable"
    assert _put(client, SETTINGS).status_code == 200
    with_settings = created()
    assert with_settings != before
    assert created() == with_settings
    assert client.delete("/api/v1/dataset/clonal").status_code == 204
    assert created() == before, "clearing the settings gives the old id back"


def test_run_id_unchanged_without_columns() -> None:
    args = ("paga_dpt", {"a": 1}, "cluster:k=v", None, "prepkey")
    assert run_id(*args) == run_id(*args, clonal_columns=())
    assert run_id(*args, clonal_columns=("b", "a")) == run_id(*args, clonal_columns=("a", "b"))
    assert run_id(*args, clonal_columns=("a",)) != run_id(*args)
