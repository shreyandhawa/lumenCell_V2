"""On-disk cache for preprocessing stages.

Layout, under ``<cache root>/v1/<dataset fingerprint>/``::

    manifest.json                 stage name -> {key, seconds, state}; used for ETAs.
                                  seconds is the time to compute the stage, kept across cache hits
    stages/<stage>-<key>/         one directory per computed stage
        <artifact>.npy            numpy arrays
        <artifact>.npz            scipy CSR matrices
        <artifact>.json           small JSON values

A stage's key hashes the pipeline version, the stage name, the parameters it
declares and the keys of the stages that produced its inputs. Changing
``n_top_genes`` therefore recomputes HVG selection and everything downstream,
and nothing upstream.

Each stage is written to a temporary directory and moved into place in one
rename, so a crash or a cancelled job can never leave a half-written stage that
later looks valid. A stage whose files go missing or cannot be read afterwards
(a disk problem, a file removed by hand) counts as not cached: it is removed
and computed again, never trusted.
"""

import hashlib
import json
import os
import shutil
import sys
import tempfile
import time
import uuid
import zipfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import platformdirs
from scipy import sparse

from lumencell.data.loader import Dataset
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import PIPELINE_VERSION, STAGES, PrepParams, Stage, StageRecord, run_stage

CACHE_LAYOUT = "v1"
ProgressFn = Callable[[str, int, int], None]
_SUFFIXES = (".npz", ".npy", ".json")
# What reading a damaged artifact raises: truncated arrays (ValueError, EOFError),
# broken or cut-short .npz archives (BadZipFile, KeyError), a missing file. Other
# OSErrors (permissions, I/O on an external disk) say nothing about the file being
# damaged, so they are raised as they are and the stage is kept.
_UNREADABLE = (
    FileNotFoundError,
    IsADirectoryError,
    NotADirectoryError,
    ValueError,
    EOFError,
    KeyError,
    zipfile.BadZipFile,
)


class CacheMiss(Exception):
    """A stage is not in the cache, or its files cannot be read."""

    def __init__(self, stage: str, reason: str) -> None:
        super().__init__(f"the {stage} stage is not in the cache ({reason})")
        self.stage = stage


def default_cache_root() -> Path:
    override = os.environ.get("LUMENCELL_CACHE_DIR")
    return Path(override) if override else Path(platformdirs.user_cache_dir("lumencell"))


def dataset_dir(fingerprint: str, root: Path | None = None) -> Path:
    return (root or default_cache_root()) / CACHE_LAYOUT / fingerprint


# ---- artifact files -------------------------------------------------------------------


def _save_artifact(directory: Path, name: str, value: Any) -> None:
    if sparse.issparse(value):
        sparse.save_npz(directory / f"{name}.npz", sparse.csr_matrix(value), compressed=False)
    elif isinstance(value, np.ndarray):
        np.save(directory / f"{name}.npy", value, allow_pickle=False)
    else:
        (directory / f"{name}.json").write_text(json.dumps(value), encoding="utf-8")


def _load_artifact(directory: Path, name: str) -> Any:
    if (p := directory / f"{name}.npz").exists():
        return sparse.csr_matrix(sparse.load_npz(p))
    if (p := directory / f"{name}.npy").exists():
        return np.load(p, allow_pickle=False)
    if (p := directory / f"{name}.json").exists():
        return json.loads(p.read_text(encoding="utf-8"))
    raise FileNotFoundError(f"artifact {name!r} missing in {directory}")


# ---- keys ------------------------------------------------------------------------------


def stage_key(stage: Stage, params: PrepParams, upstream: dict[str, str], dataset_fingerprint: str) -> str:
    """Hash of everything that determines this stage's output."""
    all_params = params.to_dict()
    payload = {
        "pipeline": PIPELINE_VERSION,
        "dataset": dataset_fingerprint,
        "stage": stage.name,
        "params": {k: all_params[k] for k in stage.params},
        "inputs": {name: upstream[name] for name in sorted(stage.inputs)},
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def stage_keys(ds: Dataset, params: PrepParams) -> dict[str, str]:
    """The cache key of every stage, without computing anything."""
    keys: dict[str, str] = {}
    producer: dict[str, str] = {}  # artifact -> key of the stage that makes it
    for stage in STAGES:
        key = stage_key(stage, params, producer, ds.fingerprint)
        keys[stage.name] = key
        for out in stage.outputs:
            producer[out] = key
    return keys


def prep_key(ds: Dataset, params: PrepParams) -> str:
    """One key for the whole prep: changes whenever any stage would."""
    keys = stage_keys(ds, params)
    return hashlib.sha256(json.dumps(keys, sort_keys=True).encode()).hexdigest()[:16]


# ---- the cache -------------------------------------------------------------------------


class PrepCache:
    def __init__(self, ds: Dataset, root: Path | None = None) -> None:
        self.ds = ds
        self.dir = dataset_dir(ds.fingerprint, root)
        self.stages_dir = self.dir / "stages"
        self._keys: dict[PrepParams, dict[str, str]] = {}

    def keys(self, params: PrepParams) -> dict[str, str]:
        """:func:`stage_keys` for this dataset, worked out once per set of parameters."""
        if params not in self._keys:
            self._keys[params] = stage_keys(self.ds, params)
        return self._keys[params]

    def _stage_dir(self, stage: Stage, key: str) -> Path:
        return self.stages_dir / f"{stage.name}-{key}"

    def has(self, stage: Stage, key: str) -> bool:
        """The stage finished writing and every artifact it declares is there."""
        try:
            names = set(os.listdir(self._stage_dir(stage, key)))
        except FileNotFoundError:
            return False
        return ".complete" in names and all(
            any(f"{out}{suffix}" in names for suffix in _SUFFIXES) for out in stage.outputs
        )

    def load(self, stage: Stage, key: str) -> dict[str, Any]:
        """The stage's artifacts. One that cannot be read invalidates the stage and raises CacheMiss."""
        d = self._stage_dir(stage, key)
        try:
            return {name: _load_artifact(d, name) for name in stage.outputs}
        except _UNREADABLE as exc:
            reason = f"{type(exc).__name__}: {exc}"
            self.invalidate(stage, key)
            print(
                f"lumencell: the cached {stage.name} stage at {d} cannot be read ({reason}); "
                "it was taken out of the cache and will be recomputed",
                file=sys.stderr,
            )
            raise CacheMiss(stage.name, reason) from exc

    def invalidate(self, stage: Stage, key: str) -> None:
        """Take a stage out of the cache, so the next run computes it again."""
        d = self._stage_dir(stage, key)
        # Out of the way in one rename (nobody can load it from then on), then deleted.
        trash = d.with_name(f".{d.name}.damaged-{uuid.uuid4().hex[:8]}")
        try:
            d.rename(trash)
        except FileNotFoundError:
            return  # another process got there first
        try:
            shutil.rmtree(trash)
        except OSError as exc:
            # The stage is already out of the cache; only the disk space is left behind
            # (`lumencell cache prune` removes it later).
            print(f"lumencell: could not delete {trash} ({exc}); it can be removed by hand", file=sys.stderr)

    def save(self, stage: Stage, key: str, outputs: dict[str, Any]) -> None:
        self.stages_dir.mkdir(parents=True, exist_ok=True)
        final = self._stage_dir(stage, key)
        tmp = Path(tempfile.mkdtemp(prefix=f".{stage.name}-", dir=self.stages_dir))
        try:
            for name, value in outputs.items():
                _save_artifact(tmp, name, value)
            (tmp / ".complete").write_text(str(time.time()), encoding="utf-8")
            if self.has(stage, key):  # a concurrent writer finished first; keep theirs
                shutil.rmtree(tmp)
                return
            if final.exists():  # an incomplete or damaged copy, which must not stay in the way
                print(f"lumencell: replacing the incomplete {stage.name} stage at {final}", file=sys.stderr)
                self.invalidate(stage, key)
            try:
                tmp.rename(final)
            except OSError:
                if not self.has(stage, key):
                    raise
                shutil.rmtree(tmp)  # a concurrent writer landed a complete copy just before us
        except BaseException:
            shutil.rmtree(tmp, ignore_errors=True)
            raise

    def _write_manifest(self, records: list[StageRecord], keys: dict[str, str], params: PrepParams) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)
        previous = self.manifest()
        known = previous.get("stages", {}) if previous else {}
        stages: dict[str, dict[str, Any]] = {}
        for r in records:
            before = known.get(r.name) if isinstance(known, dict) else None
            if r.state == "cached":
                if isinstance(before, dict) and before.get("key") == keys[r.name]:
                    # A cache hit: keep how long the stage took to compute. The time
                    # it took to load says nothing useful and would overwrite it.
                    stages[r.name] = before
                else:
                    stages[r.name] = {"key": keys[r.name], "state": "cached", "seconds": None}
            else:
                stages[r.name] = {"key": keys[r.name], "state": r.state, "seconds": round(r.seconds, 3)}
        manifest = {"pipeline": PIPELINE_VERSION, "params": params.to_dict(), "stages": stages}
        if manifest == previous:
            return  # nothing new to record; the server and the worker both load the cache often
        # Its own temporary file: the server and the worker can both be writing.
        fd, tmp_name = tempfile.mkstemp(prefix=".manifest-", suffix=".tmp", dir=self.dir)
        tmp = Path(tmp_name)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(json.dumps(manifest, indent=2))
            tmp.replace(self.dir / "manifest.json")
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise

    def manifest(self) -> dict[str, Any] | None:
        """The manifest, or None when there is none or it cannot be read (it only holds timings)."""
        try:
            data = json.loads((self.dir / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None
        return data if isinstance(data, dict) else None

    def status(self, params: PrepParams) -> dict[str, str]:
        """Per stage: 'done' (cached), 'skipped' (does not apply) or 'pending'."""
        keys = self.keys(params)
        out: dict[str, str] = {}
        for stage in STAGES:
            if not stage.applies(self.ds, params):
                out[stage.name] = "skipped"
            else:
                out[stage.name] = "done" if self.has(stage, keys[stage.name]) else "pending"
        return out

    def run(self, params: PrepParams | None = None, progress: ProgressFn | None = None) -> PrepResult:
        """Like :func:`run_pipeline`, but each stage is loaded from the cache when present.

        A stage that is missing, or whose files cannot be read, is computed.
        """
        return self._run(params or PrepParams(), progress, compute=True)

    def load_all(self, params: PrepParams) -> PrepResult:
        """Every stage from the cache, computing nothing: CacheMiss when one is not there.

        The server uses this; computing belongs to the worker.
        """
        return self._run(params, None, compute=False)

    def _run(self, params: PrepParams, progress: ProgressFn | None, compute: bool) -> PrepResult:
        keys = self.keys(params)
        artifacts: dict[str, Any] = {}
        records: list[StageRecord] = []
        for i, stage in enumerate(STAGES):
            if progress:
                progress(stage.name, i, len(STAGES))
            if not stage.applies(self.ds, params):
                records.append(StageRecord(stage.name, "skipped", 0.0))
                continue
            key = keys[stage.name]
            if self.has(stage, key):
                start = time.perf_counter()
                try:
                    artifacts.update(self.load(stage, key))
                    records.append(StageRecord(stage.name, "cached", time.perf_counter() - start))
                    continue
                except CacheMiss:
                    if not compute:
                        raise
            elif not compute:
                raise CacheMiss(stage.name, "missing or incomplete")
            out, seconds = run_stage(stage, self.ds, artifacts, params)
            self.save(stage, key, out)
            artifacts.update(out)
            records.append(StageRecord(stage.name, "done", seconds))
        if progress:
            progress("done", len(STAGES), len(STAGES))
        if compute:  # the server, which only loads, never writes to the cache
            self._write_manifest(records, keys, params)
        return PrepResult(artifacts=artifacts, stages=records)
