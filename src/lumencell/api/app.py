"""FastAPI application factory."""

import json
from pathlib import Path
from typing import Any

from fastapi import FastAPI

from lumencell import __version__
from lumencell.api.errors import install_error_handlers
from lumencell.api.routes import router


def create_app() -> FastAPI:
    app = FastAPI(
        title="lumenCell",
        version=__version__,
        summary="Local API for exploring single-cell pseudotime and lineage trajectories.",
        docs_url="/api/docs",
        redoc_url=None,
        openapi_url="/api/openapi.json",
    )
    install_error_handlers(app)
    app.include_router(router)
    return app


def openapi_document() -> dict[str, Any]:
    """The OpenAPI document, with the version pinned so the contract file only
    changes when the API does (not on every version bump)."""
    doc = create_app().openapi()
    doc["info"]["version"] = "v1"
    return doc


def write_openapi(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(openapi_document(), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    import sys

    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("contracts/openapi.json")
    write_openapi(target)
    print(f"wrote {target}")
