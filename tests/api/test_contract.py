import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from lumencell import __version__
from lumencell.api.app import create_app, openapi_document

CONTRACT = Path(__file__).resolve().parents[2] / "contracts" / "openapi.json"

client = TestClient(create_app())

# The M1 contract: (method, path) pairs that must exist.
M1_ROUTES = {
    ("get", "/api/v1/healthz"),
    ("get", "/api/v1/dataset"),
    ("get", "/api/v1/embeddings/{name}"),
    ("get", "/api/v1/obs/{column}"),
    ("get", "/api/v1/genes"),
    ("get", "/api/v1/expression"),
    ("get", "/api/v1/prepare"),
    ("post", "/api/v1/prepare"),
    ("get", "/api/v1/methods"),
    ("get", "/api/v1/root"),
    ("put", "/api/v1/root"),
    ("post", "/api/v1/root/resolve"),
    ("get", "/api/v1/root/suggestions"),
    ("post", "/api/v1/trajectories"),
    ("get", "/api/v1/trajectories"),
    ("get", "/api/v1/trajectories/{run_id}"),
    ("get", "/api/v1/trajectories/{run_id}/cells"),
    ("get", "/api/v1/trajectories/{run_id}/curves"),
    ("patch", "/api/v1/trajectories/{run_id}/lineages/{lineage_id}"),
    ("get", "/api/v1/trajectories/{run_id}/trends"),
    ("post", "/api/v1/trajectories/{run_id}/heatmap"),
    ("get", "/api/v1/compare"),
    ("get", "/api/v1/compare/cells"),
    ("get", "/api/v1/jobs/{job_id}"),
    ("get", "/api/v1/jobs/{job_id}/events"),
    ("delete", "/api/v1/jobs/{job_id}"),
}


def test_committed_contract_matches_the_models() -> None:
    assert CONTRACT.exists(), "run `make contract` to create contracts/openapi.json"
    committed = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert committed == openapi_document(), "API models changed: run `make contract` and commit the result"


def test_every_m1_route_is_declared() -> None:
    paths = openapi_document()["paths"]
    declared = {(method, path) for path, ops in paths.items() for method in ops}
    assert declared >= M1_ROUTES, sorted(M1_ROUTES - declared)


def test_lcb1_routes_advertise_the_lcb1_media_type() -> None:
    paths = openapi_document()["paths"]
    for method, path in [
        ("get", "/api/v1/embeddings/{name}"),
        ("get", "/api/v1/obs/{column}"),
        ("get", "/api/v1/expression"),
        ("get", "/api/v1/trajectories/{run_id}/cells"),
        ("post", "/api/v1/trajectories/{run_id}/heatmap"),
    ]:
        content = paths[path][method]["responses"]["200"]["content"]
        assert "application/vnd.lumencell.lcb1" in content, path


def test_trajectory_creation_documents_root_required() -> None:
    op = openapi_document()["paths"]["/api/v1/trajectories"]["post"]
    assert "409" in op["responses"]


def test_healthz() -> None:
    res = client.get("/api/v1/healthz")
    assert res.status_code == 200
    assert res.json() == {"status": "ok", "version": __version__}


@pytest.mark.parametrize("path", ["/api/v1/dataset", "/api/v1/methods", "/api/v1/trajectories"])
def test_dataset_routes_without_a_dataset_answer_503(path: str) -> None:
    res = client.get(path)
    assert res.status_code == 503
    assert res.json()["error"]["code"] == "no_dataset"


def test_unknown_route_uses_the_error_body() -> None:
    res = client.get("/api/v1/nope")
    assert res.status_code == 404
    assert res.json()["error"]["code"] == "not_found"


def test_root_spec_is_a_discriminated_union() -> None:
    schema = openapi_document()["components"]["schemas"]["RootRequest"]["properties"]["root"]
    assert schema["discriminator"]["propertyName"] == "kind"
