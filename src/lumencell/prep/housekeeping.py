"""What the on-disk cache holds, and removing what is no longer needed.

Layout (see ``prep/cache.py`` and ``engine/store.py``)::

    <cache root>/v1/<fingerprint>/    one dataset: manifest.json, stages/, runs/

Only those dataset folders, and leftovers of interrupted jobs inside them, are
ever removed. Everything else under the cache root (downloaded benchmark data in
``datasets/``, anything else put there) and the state directory holding the
saved roots are left alone, and every path is checked before anything goes.
"""

import json
import os
import re
import shutil
import time
from dataclasses import dataclass
from datetime import timedelta
from pathlib import Path

from lumencell.api.state import default_state_dir
from lumencell.prep.cache import CACHE_LAYOUT, LAST_USED
from lumencell.prep.stages import STAGES

FINGERPRINT = re.compile(r"[0-9a-f]{24}")
# Stages and runs sit in a temporary folder only while they are being saved, which
# takes seconds. One older than this belongs to a job that was cancelled or crashed.
TEMP_GRACE = timedelta(hours=1)
_AGE = re.compile(r"([1-9][0-9]*)([hdw])")
# The temporary names the cache and the run store use (tempfile adds 8 characters).
# Anything else, a .DS_Store say, is not ours to remove.
_RANDOM = r"[a-z0-9_]{8}"
_STAGE_NAMES = "|".join(s.name for s in STAGES)
_STAGE_TEMP = re.compile(rf"\.(?:{_STAGE_NAMES})-(?:{_RANDOM}|[0-9a-f]{{16}}\.damaged-[0-9a-f]{{8}})")
_RUN_TEMP = re.compile(rf"\.[A-Za-z0-9_-]{{1,64}}-{_RANDOM}")
_MANIFEST_TEMP = re.compile(rf"\.manifest-{_RANDOM}\.tmp")
_GRAPH_TEMP = "result.json.tmp"
_DATASET_PARTS = ("manifest.json", "stages", "runs")


class UnsafePath(ValueError):
    """A path that must not be removed: outside the cache, a symbolic link, or the saved state."""


@dataclass(frozen=True)
class CacheEntry:
    fingerprint: str
    path: Path
    name: str | None  # the dataset file's name, from the manifest (None in caches made before it was kept)
    source: str | None  # the dataset file's path
    size: int  # bytes
    last_used: float  # POSIX time
    is_link: bool = False  # a symbolic link, put there by hand: listed, never removed


@dataclass(frozen=True)
class Removal:
    path: Path
    size: int
    reason: str


def parse_age(text: str) -> timedelta:
    """``12h``, ``30d`` or ``2w``."""
    m = _AGE.fullmatch(text.strip())
    if m is None:
        raise ValueError(f"{text!r} is not an age like 12h, 30d or 2w")
    n, unit = int(m.group(1)), m.group(2)
    return {"h": timedelta(hours=n), "d": timedelta(days=n), "w": timedelta(weeks=n)}[unit]


def format_size(n: int) -> str:
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= scale:
            return f"{n / scale:.1f} {unit}"
    return f"{n} B"


def format_ago(then: float, now: float) -> str:
    days = int((now - then) // 86400)
    return "today" if days <= 0 else "1 day ago" if days == 1 else f"{days} days ago"


def _tree_size(path: Path) -> int:
    """Bytes in the files under ``path``, without following symbolic links."""
    if path.is_symlink() or not path.is_dir():
        return path.lstat().st_size
    total = 0
    for dirpath, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(dirpath, f)).st_size
            except FileNotFoundError:
                continue  # renamed or removed by a job since it was listed
    return total


def _last_used(entry: Path) -> float:
    try:
        return (entry / LAST_USED).lstat().st_mtime
    except FileNotFoundError:  # a cache from before the stamp existed
        return entry.lstat().st_mtime


def _dataset_of(entry: Path) -> tuple[str | None, str | None]:
    try:
        data = json.loads((entry / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None, None
    dataset = data.get("dataset") if isinstance(data, dict) else None
    if not isinstance(dataset, dict):
        return None, None
    name, source = dataset.get("name"), dataset.get("path")
    return (name if isinstance(name, str) else None), (source if isinstance(source, str) else None)


def list_cache(root: Path) -> list[CacheEntry]:
    """Every dataset in the cache, most recently used first."""
    base = root / CACHE_LAYOUT
    if not base.is_dir():
        return []
    out = []
    for p in base.iterdir():
        if not FINGERPRINT.fullmatch(p.name) or not p.is_dir():  # is_dir follows a link
            continue
        if not any((p / part).exists() for part in _DATASET_PARTS):
            continue  # a folder with a fingerprint's name, but nothing the cache wrote
        link = p.is_symlink()
        name, source = (None, None) if link else _dataset_of(p)
        try:
            out.append(CacheEntry(p.name, p, name, source, _tree_size(p), _last_used(p), link))
        except FileNotFoundError:
            continue  # removed since it was listed
    return sorted(out, key=lambda e: e.last_used, reverse=True)


def _leftovers(entry: Path, before: float) -> list[Path]:
    """Temporary files and folders of interrupted jobs, last touched before ``before``."""
    found = [p for p in entry.iterdir() if _MANIFEST_TEMP.fullmatch(p.name)]
    stages, runs = entry / "stages", entry / "runs"
    if stages.is_dir() and not stages.is_symlink():
        found += [p for p in stages.iterdir() if _STAGE_TEMP.fullmatch(p.name)]
    if runs.is_dir() and not runs.is_symlink():
        for p in runs.iterdir():
            if _RUN_TEMP.fullmatch(p.name):
                found.append(p)
            elif p.is_dir() and not p.is_symlink() and (p / _GRAPH_TEMP).exists():
                found.append(p / _GRAPH_TEMP)
    old = []
    for p in found:
        try:
            if p.lstat().st_mtime < before:
                old.append(p)
        except FileNotFoundError:
            continue  # the job finished with it meanwhile
    return old


def plan_prune(
    root: Path, *, older_than: timedelta | None, everything: bool, now: float | None = None
) -> list[Removal]:
    """What a prune would remove: whole datasets, then leftovers inside the datasets that stay."""
    now = time.time() if now is None else now
    out = []
    for e in list_cache(root):
        label = e.name or e.fingerprint
        if e.is_link:
            continue
        if everything:
            out.append(Removal(e.path, e.size, f"dataset {label}"))
        elif older_than is not None and e.last_used < now - older_than.total_seconds():
            out.append(Removal(e.path, e.size, f"dataset {label}, last used {format_ago(e.last_used, now)}"))
        else:
            for p in _leftovers(e.path, now - TEMP_GRACE.total_seconds()):
                out.append(Removal(p, _tree_size(p), f"leftover of an interrupted job in {label}"))
    return out


def check_prune(root: Path, removals: list[Removal]) -> None:
    """Raise UnsafePath unless every removal is a dataset folder, or something inside one."""
    base = (root / CACHE_LAYOUT).resolve()
    state = default_state_dir().resolve()
    for r in removals:
        if r.path.is_symlink():
            raise UnsafePath(
                f"refusing to remove {r.path}: it is a symbolic link, which may lead outside the cache"
            )
        resolved = r.path.resolve()
        if resolved == base or not resolved.is_relative_to(base):
            raise UnsafePath(f"refusing to remove {r.path}: it is outside the cache's datasets ({base})")
        entry = base / resolved.relative_to(base).parts[0]
        if not FINGERPRINT.fullmatch(entry.name) or not any((entry / p).exists() for p in _DATASET_PARTS):
            raise UnsafePath(f"refusing to remove {r.path}: {entry} is not a dataset the cache wrote")
        inside_state = state.is_relative_to(base) and resolved.is_relative_to(state)
        if inside_state or state.is_relative_to(resolved):
            raise UnsafePath(f"refusing to remove {r.path}: it holds the saved state ({state})")


def apply_prune(root: Path, removals: list[Removal]) -> list[tuple[Removal, str | None]]:
    """Check every removal first (nothing is removed if one is refused), then remove each.

    Returns each removal with what went wrong (None when it was removed). Removing a
    leftover does not count as using its dataset: the folder keeps its times.
    """
    check_prune(root, removals)
    base = (root / CACHE_LAYOUT).resolve()
    entries = {base / r.path.resolve().relative_to(base).parts[0] for r in removals}
    times = {e: (e.stat().st_atime, e.stat().st_mtime) for e in entries if e.is_dir()}
    outcomes: list[tuple[Removal, str | None]] = []
    for r in removals:
        try:
            if r.path.is_dir():
                shutil.rmtree(r.path)
            else:
                r.path.unlink()
            outcomes.append((r, None))
        except FileNotFoundError:
            outcomes.append((r, "already gone"))
        except OSError as exc:
            outcomes.append((r, f"could not remove: {exc}"))
    for e, (atime, mtime) in times.items():
        if e.is_dir():
            os.utime(e, (atime, mtime))
    return outcomes
