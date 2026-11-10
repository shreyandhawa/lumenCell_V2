"""The user's saved choices for one dataset: the root and lineage renames.

They live in a small JSON file outside the cache, so they survive cache eviction.
"""

import contextlib
import json
import os
import tempfile
import threading
from collections.abc import Callable
from pathlib import Path
from typing import Any

import platformdirs


def default_state_dir() -> Path:
    override = os.environ.get("LUMENCELL_DATA_DIR")
    return Path(override) if override else Path(platformdirs.user_data_dir("lumencell")) / "state"


class StateStore:
    def __init__(self, path: Path) -> None:
        self.path = path
        self._lock = threading.Lock()

    def read(self) -> dict[str, Any]:
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            return {}
        except (OSError, ValueError):
            data = None
        if isinstance(data, dict):
            return data
        # Unreadable or not the expected shape. Set it aside (for inspection) and
        # start clean, rather than failing every route that reads the root.
        with contextlib.suppress(OSError):
            self.path.replace(self.path.with_suffix(".corrupt"))
        return {}

    def update(self, change: Callable[[dict[str, Any]], None]) -> None:
        """Read, change and write the state file as one step.

        Requests run on several threads; without the lock two of them could each
        read the old state and the second write would drop the first one's change.
        The temporary file has its own name so writers never share one either.
        """
        with self._lock:
            state = self.read()
            change(state)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(prefix=f".{self.path.stem}-", suffix=".tmp", dir=self.path.parent)
            tmp = Path(tmp_name)
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as f:
                    f.write(json.dumps(state, indent=2))
                tmp.replace(self.path)
            except BaseException:
                tmp.unlink(missing_ok=True)
                raise
