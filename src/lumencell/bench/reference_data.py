"""The reference benchmark rows packaged with lumenCell, read by the server.

Only json and importlib are imported here, so the installed wheel can read the rows without its
other dependencies. ``reference.py`` writes the file from the committed results.
"""

import json
from functools import cache
from importlib.resources import files
from typing import Any


@cache
def _methods() -> dict[str, list[dict[str, Any]]]:
    text = files("lumencell.bench").joinpath("reference_scores.json").read_text(encoding="utf-8")
    methods: dict[str, list[dict[str, Any]]] = json.loads(text)["methods"]
    return methods


def rows_for(method: str) -> list[dict[str, Any]]:
    """The reference rows of a method's defaults; empty for a method with no benchmark."""
    return [dict(row) for row in _methods().get(method, [])]
