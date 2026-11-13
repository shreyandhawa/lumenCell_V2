"""Dataset-level clonal settings: the PUT/DELETE routes, what they change, and what they refuse.

The job manager is never started: a run is only ever created, never executed.
"""

import json
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
    # Label-key detection looks for these names; the second is the next candidate once the first is
    # a donor column.
    adata.obs["cell_type"] = adata.obs["donor"].astype(str).astype("category")
    adata.obs["celltype"] = adata.obs["true_segment"].astype(str).astype("category")
    # A donor only where there is no clone label: every labelled cell lacks a donor.
    adata.obs["donor_of_unlabelled"] = pd.Categorical(
        adata.obs["donor"].astype(object).where(adata.obs["clone_id"].isna())
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


# ---- fix round 1 -----------------------------------------------------------------------


def _state(session: Session) -> dict[str, Any]:
    return session.state.read()


@pytest.mark.parametrize(
    "body, word",
    [
        ({**SETTINGS, "strata_keys": ["true_segment", "true_segment"]}, "distinct"),
        ({**SETTINGS, "strata_keys": ["donor"]}, "strata"),
        ({"clone_key": "clone_id", "donor_key": "donor_of_unlabelled"}, "no cell has both"),
    ],
)
def test_more_invalid_settings(client: TestClient, body: dict[str, Any], word: str) -> None:
    res = _put(client, body)
    assert res.status_code == 422
    error = res.json()["error"]
    assert error["code"] == "invalid_clonal_settings"
    assert any(word in p for p in error["detail"]["problems"]), error["detail"]


def test_single_donor_through_the_api(client: TestClient, session: Session) -> None:
    res = _put(client, {"clone_key": "clone_id", "single_donor": True})
    assert res.status_code == 200
    body = res.json()
    assert body["n_donors"] == 1
    labelled = session.ds.obs["clone_id"].dropna().astype(str)
    assert body["n_clones"] == labelled.nunique()
    assert body["n_cells_labelled"] == len(labelled)
    assert session.view.clonal_columns == ("clone_id",)
    assert "donor" in session.view.ctx().obs
    assert client.get("/api/v1/dataset").json()["clonal"]["single_donor"] is True


@pytest.mark.parametrize("key", ["cat_clone", "donor"])
def test_a_saved_root_on_the_clone_or_donor_column_is_refused(
    client: TestClient, session: Session, key: str
) -> None:
    value = "c0" if key != "donor" else "d0"
    root = {"kind": "cluster", "key": key, "value": value}
    assert client.put("/api/v1/root", json={"root": root}).status_code == 200
    body = {"clone_key": "cat_clone", "donor_key": "donor"}
    res = _put(client, body)
    assert res.status_code == 422 and res.json()["error"]["code"] == "scoring_only_column"
    assert client.get("/api/v1/dataset").json()["clonal"] is None
    assert "clonal" not in _state(session)
    assert session.view.clonal_columns == ()


def test_a_root_on_an_unrelated_or_strata_column_is_accepted(client: TestClient) -> None:
    assert client.put("/api/v1/root", json={"root": TRUNK}).status_code == 200
    assert _put(client, SETTINGS).status_code == 200
    assert _put(client, {**SETTINGS, "strata_keys": ["true_segment"]}).status_code == 200


def test_refusals_leave_the_saved_state_alone(
    client: TestClient, session: Session, prepared: tuple[Path, Path], jobs: JobManager, tmp_path: Path
) -> None:
    assert _put(client, {"clone_key": "nope", "donor_key": "donor"}).status_code == 422
    assert "clonal" not in _state(session)
    again = _session(prepared, jobs, tmp_path / "state")
    assert again.clonal_settings() is None


def test_an_unparseable_saved_root_does_not_block_the_settings(client: TestClient, session: Session) -> None:
    session.state.update(lambda s: s.update(root="garbage"))
    assert _put(client, SETTINGS).status_code == 200


def test_label_key_is_detected_again_after_set_and_clear(client: TestClient, session: Session) -> None:
    assert session.view.label_key == "cell_type"
    assert _put(client, {"clone_key": "clone_id", "donor_key": "cell_type"}).status_code == 200
    assert session.view.label_key == "celltype"
    assert client.get("/api/v1/dataset").json()["label_key"] == "celltype"
    assert client.delete("/api/v1/dataset/clonal").status_code == 204
    assert session.view.label_key == "cell_type"


def test_delete_puts_the_donor_column_back_in_the_context(client: TestClient, session: Session) -> None:
    assert _put(client, SETTINGS).status_code == 200
    assert "donor" not in session.view.ctx().obs
    assert client.delete("/api/v1/dataset/clonal").status_code == 204
    ctx = session.view.ctx()
    assert "donor" in ctx.obs and "donor" not in ctx.scoring_only


def test_delete_is_remembered_by_a_new_session(
    client: TestClient, session: Session, prepared: tuple[Path, Path], jobs: JobManager, tmp_path: Path
) -> None:
    assert _put(client, SETTINGS).status_code == 200
    assert client.delete("/api/v1/dataset/clonal").status_code == 204
    assert "clonal" not in _state(session)
    again = _session(prepared, jobs, tmp_path / "state")
    assert again.clonal_settings() is None and again.view.clonal_columns == ()


def _write_state(session: Session, state: dict[str, Any]) -> None:
    session.state.path.parent.mkdir(parents=True, exist_ok=True)
    session.state.path.write_text(json.dumps(state), encoding="utf-8")


@pytest.mark.parametrize(
    "saved",
    [
        {"clonal": {"clone_key": "clone_id"}},  # neither donor_key nor single_donor
        {"clonal": "x"},
    ],
)
def test_an_invalid_saved_value_is_ignored_with_one_line(
    session: Session,
    prepared: tuple[Path, Path],
    jobs: JobManager,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    saved: dict[str, Any],
) -> None:
    _write_state(session, saved)
    capsys.readouterr()
    again = _session(prepared, jobs, tmp_path / "state")
    assert again.clonal_settings() is None and again.view.clonal_columns == ()
    lines = [x for x in capsys.readouterr().err.splitlines() if "ignoring the saved clonal settings" in x]
    assert len(lines) == 1


@pytest.mark.parametrize("saved", [{}, {"clonal": None}])
def test_a_missing_saved_value_is_silent(
    session: Session,
    prepared: tuple[Path, Path],
    jobs: JobManager,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    saved: dict[str, Any],
) -> None:
    _write_state(session, saved)
    capsys.readouterr()
    again = _session(prepared, jobs, tmp_path / "state")
    assert again.clonal_settings() is None
    assert capsys.readouterr().err == ""


def test_saved_settings_that_fail_the_check_are_ignored(
    session: Session,
    prepared: tuple[Path, Path],
    jobs: JobManager,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _write_state(session, {"clonal": {"clone_key": "not_a_column", "donor_key": "donor"}})
    capsys.readouterr()
    again = _session(prepared, jobs, tmp_path / "state")
    assert again.clonal_settings() is None and again.view.clonal_columns == ()
    assert "ignoring the saved clonal settings" in capsys.readouterr().err
    assert again.dataset_info().clonal is None


def test_run_id_hashes_every_column_in_sorted_order() -> None:
    args = ("paga_dpt", {"a": 1}, "cluster:k=v", None, "prepkey")
    assert run_id(*args, ("a", "b")) != run_id(*args, ("a", "c"))
    assert run_id(*args, ("a", "b")) != run_id(*args, ("a",))
    assert run_id(*args, ("b", "a")) == run_id(*args, ("a", "b"))


def test_run_id_without_columns_is_the_id_from_before_this_task() -> None:
    # Computed at 4a8cd8a, before run_id took clonal columns. It also depends on the store
    # version and the package version: if either changes, this literal changes with it.
    args = ("paga_dpt", {"a": 1}, "cluster:k=v", None, "prepkey")
    assert run_id(*args) == "443cbd44c784478d"
    assert run_id(*args, ()) == "443cbd44c784478d"
