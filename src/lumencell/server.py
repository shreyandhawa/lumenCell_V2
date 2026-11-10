"""Assemble and run the local server: API, web app and background work."""

import socket
import threading
import time
import webbrowser
from importlib import resources
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from starlette.middleware.trustedhost import TrustedHostMiddleware

from lumencell.api.app import create_app
from lumencell.api.session import Session

LOCAL_HOSTS = ["127.0.0.1", "localhost"]
REPO_WEB_DIST = Path(__file__).resolve().parents[2] / "web" / "dist"

_NOT_BUILT = """<!doctype html><meta charset="utf-8"><title>lumenCell</title>
<body style="font-family:system-ui;margin:3rem;max-width:40rem">
<h1>lumenCell</h1><p>The API is running, but the web app has not been built.</p>
<p>From the repository, run <code>make web-dist</code> (or <code>make dev-web</code> for live
development), then reload. The API is documented at <a href="/api/docs">/api/docs</a>.</p></body>"""


def web_dist() -> Path | None:
    """The built web app: bundled in the package, or the repo's web/dist in development."""
    packaged = Path(str(resources.files("lumencell"))) / "web_dist"
    for candidate in (packaged, REPO_WEB_DIST):
        if (candidate / "index.html").exists():
            return candidate
    return None


def build_app(session: Session | None, trusted_hosts: list[str] | None = None) -> FastAPI:
    app = create_app(session)
    if trusted_hosts:
        # Only answer requests addressed to this machine by name: blocks DNS rebinding,
        # where a web page re-points its own domain at 127.0.0.1 to reach the API.
        app.add_middleware(TrustedHostMiddleware, allowed_hosts=trusted_hosts)
    dist = web_dist()
    if dist is not None:
        app.mount("/", StaticFiles(directory=dist, html=True), name="web")
    else:

        @app.get("/", include_in_schema=False)
        def not_built() -> HTMLResponse:
            return HTMLResponse(_NOT_BUILT)

    return app


def free_port(host: str, preferred: int, attempts: int = 20) -> int:
    for port in range(preferred, preferred + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((host, port))
            except OSError:
                continue
            return port
    raise RuntimeError(f"no free port in {preferred}-{preferred + attempts - 1}")


def open_when_ready(url: str, timeout: float = 30.0) -> None:
    """Open the browser once /healthz answers (from a background thread)."""
    import urllib.request

    def wait_and_open() -> None:
        deadline = time.time() + timeout
        while time.time() < deadline:
            try:
                with urllib.request.urlopen(f"{url}/api/v1/healthz", timeout=1):
                    webbrowser.open(url)
                    return
            except OSError:
                time.sleep(0.2)

    threading.Thread(target=wait_and_open, daemon=True).start()
