"""The assembled server: static web app, host checking, and `lumencell launch` end to end."""

import json
import os
import signal
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
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


def test_free_port_skips_a_busy_one_unless_strict() -> None:
    import socket

    with socket.socket() as busy:
        busy.bind(("127.0.0.1", 0))
        taken = busy.getsockname()[1]
        assert free_port("127.0.0.1", taken) != taken
        with pytest.raises(RuntimeError, match="no free port"):
            free_port("127.0.0.1", taken, attempts=1)


def _get(url: str) -> dict[str, object]:
    with urllib.request.urlopen(url, timeout=5) as res:
        data: dict[str, object] = json.loads(res.read())
        return data


def _processes() -> dict[int, tuple[int, str]]:
    """Every process on the machine: pid -> (parent pid, command line)."""
    out = subprocess.run(
        ["ps", "-Ao", "pid=,ppid=,command="], capture_output=True, text=True, check=True
    ).stdout
    procs = {}
    for line in out.splitlines():
        pid, ppid, *command = line.split(None, 2)
        procs[int(pid)] = (int(ppid), command[0] if command else "")
    return procs


def _children(pid: int) -> dict[int, str]:
    return {child: cmd for child, (ppid, cmd) in _processes().items() if ppid == pid}


def _survivors(children: dict[int, str], timeout: float) -> dict[int, str]:
    """The children still running after up to `timeout` seconds (same pid, same command)."""
    deadline = time.time() + timeout
    while True:
        running = _processes()
        alive = {pid: cmd for pid, cmd in children.items() if running.get(pid, (0, ""))[1] == cmd}
        if not alive or time.time() >= deadline:
            return alive
        time.sleep(0.2)


@contextmanager
def _launch(tmp_path: Path, extra: Sequence[str] = ()) -> Iterator[tuple["subprocess.Popen[bytes]", str]]:
    """`lumencell launch` on a small dataset, once it answers. On the way out the server is
    stopped and any child of it that outlives it (its job worker, the resource tracker) is
    killed, so a failing run never leaves processes behind."""
    data = tmp_path / "d.h5ad"
    generate("bifurcation", "realistic", n_cells=500, n_genes=300, seed=0).write_h5ad(data)
    port = free_port("127.0.0.1", 18765)
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
    children: dict[int, str] = {}
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
        children = _children(proc.pid)
        yield proc, url
    finally:
        proc.terminate()
        try:
            proc.wait(30)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        for pid in _survivors(children, timeout=10):
            os.kill(pid, signal.SIGKILL)


@pytest.fixture
def launched(tmp_path: Path, request: pytest.FixtureRequest) -> Iterator[str]:
    with _launch(tmp_path, getattr(request, "param", [])) as (_proc, url):
        yield url


def test_launch_serves_the_api_and_prepares_in_the_background(launched: str) -> None:
    info = _get(f"{launched}/api/v1/dataset")
    assert info["n_obs"] == 500
    deadline = time.time() + 120
    while _get(f"{launched}/api/v1/prepare")["state"] != "ready":
        assert time.time() < deadline, "prep did not finish"
        time.sleep(0.5)


def test_sigterm_stops_launch_and_its_job_worker(tmp_path: Path) -> None:
    with _launch(tmp_path) as (proc, _url):
        children = _children(proc.pid)
        assert any("spawn_main" in cmd for cmd in children.values()), f"no job worker among {children}"
        proc.terminate()
        code = proc.wait(30)
        survivors = _survivors(children, timeout=10)
        assert survivors == {}, f"still running after the server stopped: {survivors}"
        # 143 = 128 + SIGTERM, as a shell reports it: the server unwound through its cleanup
        # (which stops the job manager) instead of being ended by the signal on the spot (-15).
        assert code == 128 + signal.SIGTERM


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
