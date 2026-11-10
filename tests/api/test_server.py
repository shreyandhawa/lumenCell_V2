"""The assembled server: static web app, host checking, and `lumencell launch` end to end."""

import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections.abc import Iterator
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from lumencell.bench.synth import generate
from lumencell.server import build_app, free_port, web_dist

LUMENCELL = Path(sys.executable).parent / "lumencell"


def test_unknown_hosts_are_rejected() -> None:
    client = TestClient(
        build_app(None, trusted_hosts=["127.0.0.1", "localhost"]), base_url="http://127.0.0.1"
    )
    assert client.get("/api/v1/healthz").status_code == 200
    assert client.get("/api/v1/healthz", headers={"host": "evil.example"}).status_code == 400


def test_root_serves_the_web_app_or_explains_how_to_build_it() -> None:
    res = TestClient(build_app(None)).get("/")
    assert res.status_code == 200
    if web_dist() is None:
        assert "make web-dist" in res.text
    else:
        assert '<div id="root">' in res.text


def test_free_port_skips_a_busy_one() -> None:
    import socket

    with socket.socket() as busy:
        busy.bind(("127.0.0.1", 0))
        taken = busy.getsockname()[1]
        assert free_port("127.0.0.1", taken) != taken


def _get(url: str) -> dict[str, object]:
    with urllib.request.urlopen(url, timeout=5) as res:
        data: dict[str, object] = json.loads(res.read())
        return data


@pytest.fixture
def launched(tmp_path: Path, request: pytest.FixtureRequest) -> Iterator[str]:
    data = tmp_path / "d.h5ad"
    generate("bifurcation", "realistic", n_cells=500, n_genes=300, seed=0).write_h5ad(data)
    port = free_port("127.0.0.1", 18765)
    extra = getattr(request, "param", [])
    env = {
        **os.environ,
        "LUMENCELL_CACHE_DIR": str(tmp_path / "c"),
        "LUMENCELL_DATA_DIR": str(tmp_path / "s"),
    }
    proc = subprocess.Popen(
        [str(LUMENCELL), "launch", str(data), "--no-browser", "--port", str(port), *extra],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    url = f"http://127.0.0.1:{port}"
    try:
        deadline = time.time() + 60
        while time.time() < deadline:
            try:
                _get(f"{url}/api/v1/healthz")
                break
            except (urllib.error.URLError, ConnectionError):
                if proc.poll() is not None:
                    raise RuntimeError(
                        proc.stdout.read().decode() if proc.stdout else "launch exited"
                    ) from None
                time.sleep(0.25)
        else:
            raise TimeoutError("server did not start")
        yield url
    finally:
        proc.terminate()
        proc.wait(10)


def test_launch_serves_the_api_and_prepares_in_the_background(launched: str) -> None:
    info = _get(f"{launched}/api/v1/dataset")
    assert info["n_obs"] == 500
    deadline = time.time() + 120
    while _get(f"{launched}/api/v1/prepare")["state"] != "ready":
        assert time.time() < deadline, "prep did not finish"
        time.sleep(0.5)


@pytest.mark.slow
@pytest.mark.parametrize("launched", [["--root", "cluster:true_label=trunk_early"]], indirect=True)
def test_launch_with_a_root_has_a_run_ready(launched: str) -> None:
    deadline = time.time() + 180
    while True:
        runs = _get(f"{launched}/api/v1/trajectories")["runs"]
        if runs and runs[0]["state"] == "done":  # type: ignore[index]
            break
        assert time.time() < deadline, f"no finished run: {runs}"
        time.sleep(0.5)
    root = _get(f"{launched}/api/v1/root")["root"]
    assert root is not None and root["spec"]["value"] == "trunk_early"  # type: ignore[index]
