"""The worker process: runs one job at a time and reports progress as events.

It is started with ``spawn`` and kept alive between jobs, so scanpy and numba
stay imported and warm. Datasets it has loaded are kept by path. Nothing here
writes partial results: the stage cache and the run store both commit with an
atomic rename, so terminating the worker mid-job (to cancel it) is safe.

The worker ends itself as soon as the process that started it is gone, however
it went (killed, crashed). It would otherwise wait on the task queue for ever:
it holds a write end of that queue itself, so it never sees end-of-file.

Events sent back on the events queue (plain tuples, so they pickle cheaply)::

    ("started",  job_id)
    ("progress", job_id, stage, i, n)
    ("done",     job_id, result_id_or_None)
    ("failed",   job_id, message)
"""

import multiprocessing as mp
import os
import sys
import threading
import time
import traceback
from collections.abc import Callable
from dataclasses import dataclass, field
from multiprocessing.process import BaseProcess
from multiprocessing.queues import Queue
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Task:
    job_id: str
    kind: str
    payload: dict[str, Any] = field(default_factory=dict)


Emit = Callable[[str, int, int], None]
_datasets: dict[tuple[str, str | None], Any] = {}


def _dataset(path: str, x_is: str | None) -> Any:
    from lumencell.data.loader import load

    key = (path, x_is)
    if key not in _datasets:
        _datasets[key] = load(Path(path), x_is)  # type: ignore[arg-type]
    return _datasets[key]


def _cache_root(payload: dict[str, Any]) -> Path | None:
    root = payload.get("cache_root")
    return Path(root) if root else None


def handle_prepare(payload: dict[str, Any], emit: Emit) -> str | None:
    from lumencell.prep.cache import PrepCache

    ds = _dataset(payload["path"], payload.get("x_is"))
    PrepCache(ds, _cache_root(payload)).run(payload["prep_params"], progress=emit)
    return None


def handle_infer(payload: dict[str, Any], emit: Emit) -> str | None:
    from lumencell.engine.context import InferenceContext
    from lumencell.engine.registry import run_method
    from lumencell.engine.roots import parse_root
    from lumencell.engine.store import RunStore
    from lumencell.prep.cache import PrepCache, dataset_dir, prep_key

    # Required, with no default: a job that forgot it would let methods see clone columns.
    if "clonal_columns" not in payload:
        raise ValueError(
            "the infer job's payload has no 'clonal_columns'; every infer job must name the columns "
            "that are only for scoring (an empty list when there are none)"
        )
    clonal_columns = payload["clonal_columns"]
    if not isinstance(clonal_columns, list | tuple) or not all(isinstance(c, str) for c in clonal_columns):
        # A bare string would otherwise be read as one column per character.
        raise ValueError(
            f"the infer job's 'clonal_columns' must be a list of column names, not {clonal_columns!r}"
        )
    ds = _dataset(payload["path"], payload.get("x_is"))
    params = payload["prep_params"]
    emit("load", 0, 1)
    prep = PrepCache(ds, _cache_root(payload)).run(params)
    ctx = InferenceContext.from_prep(ds, prep, params, clonal_columns=clonal_columns)
    result = run_method(
        ctx,
        payload["method"],
        parse_root(payload["root"]),
        terminals=payload.get("terminals"),
        params=payload.get("params") or {},
        progress=emit,
        label_key=payload.get("label_key"),
    )
    rid: str = payload["run_id"]
    RunStore(dataset_dir(ds.fingerprint, _cache_root(payload))).save(
        rid, result, terminals=payload.get("terminals"), prep_key=prep_key(ds, params)
    )
    return rid


SCORE_CLONAL_KEYS = ("path", "prep_params", "run_id", "settings", "trust_name")


def handle_score_clonal(payload: dict[str, Any], emit: Emit) -> str | None:
    """Score a finished run against the dataset's clones and store it as the run's trust file.

    Only this job writes a run's clonal trust. A run that is gone or incomplete is never scored
    anywhere else: the job fails instead.
    """
    import numpy as np

    from lumencell.api.models import ClonalSettings
    from lumencell.engine.store import STORE_VERSION, RunStore
    from lumencell.prep.cache import PrepCache, dataset_dir, prep_key
    from lumencell.trust.clonal import agreement_for_run, not_scored, trust_name

    missing = [k for k in SCORE_CLONAL_KEYS if k not in payload]
    if missing:
        raise ValueError(f"the score_clonal job's payload has no {', '.join(map(repr, missing))}")
    ds = _dataset(payload["path"], payload.get("x_is"))
    params = payload["prep_params"]
    settings = ClonalSettings.model_validate(payload["settings"])
    rid: str = payload["run_id"]
    expected = trust_name(settings)
    if payload["trust_name"] != expected:
        raise ValueError(
            f"the score_clonal job's trust_name {payload['trust_name']!r} is not the one for its "
            f"settings ({expected!r}); nothing was written"
        )
    store = RunStore(dataset_dir(ds.fingerprint, _cache_root(payload)))
    gone = f"run {rid} was replaced or is incomplete; not scored"
    emit("load", 0, 2)
    try:
        if not (store.path(rid) / ".complete").exists():
            raise FileNotFoundError(rid)
        graph = store.graph(rid)
        current = graph.get("store_version") == STORE_VERSION
        same_prep = graph.get("prep_key") == prep_key(ds, params)
    except FileNotFoundError:
        raise RuntimeError(gone) from None
    if not current:
        agreement = not_scored("refused", settings, "made before inference isolation; run it again")
    elif not same_prep:
        # Not a saved verdict: the prep key includes settings the session can change back, and a
        # stored "failed" file would block rescoring once the matching prep is reopened.
        raise RuntimeError("the run was prepared with other settings; not scored")
    else:
        try:
            result = store.load(rid)
        except FileNotFoundError:
            raise RuntimeError(gone) from None
        prep = PrepCache(ds, _cache_root(payload)).run(params)
        emit("score", 1, 2)
        agreement = agreement_for_run(
            result,
            ds.obs,
            settings,
            np.asarray(prep.artifacts["pca"]),
            params.neighbors_pcs,
            graph["obs_keys"],
        )
    try:
        store.write_trust(rid, expected, agreement.model_dump(mode="json"))
    except FileNotFoundError:
        raise RuntimeError(gone) from None
    return rid


def handle_sleep(payload: dict[str, Any], emit: Emit) -> str | None:
    """Diagnostic job: sleeps in steps. Used to test progress and cancellation."""
    steps = int(payload.get("steps", 5))
    for i in range(steps):
        emit("sleep", i, steps)
        time.sleep(float(payload.get("seconds", 0.1)))
    if payload.get("fail"):
        raise RuntimeError("asked to fail")
    return None


HANDLERS: dict[str, Callable[[dict[str, Any], Emit], str | None]] = {
    "prepare": handle_prepare,
    "infer": handle_infer,
    "score_clonal": handle_score_clonal,
    "sleep": handle_sleep,
}


def _exit_with(parent: BaseProcess) -> None:
    parent.join()  # returns once the parent has exited
    os._exit(1)  # mid-job too: numba or igraph code cannot be interrupted any other way


def worker_main(tasks: "Queue[Task | None]", events: "Queue[tuple[Any, ...]]") -> None:
    os.environ.setdefault("PYTHONWARNINGS", "ignore")
    parent = mp.parent_process()
    if parent is not None:
        threading.Thread(target=_exit_with, args=(parent,), name="parent-watch", daemon=True).start()
    while True:
        task = tasks.get()
        if task is None:
            return
        events.put(("started", task.job_id))

        def emit(stage: str, i: int, n: int, _job: str = task.job_id) -> None:
            events.put(("progress", _job, stage, i, n))

        try:
            handler = HANDLERS[task.kind]
            events.put(("done", task.job_id, handler(task.payload, emit)))
        except Exception as exc:
            # The job keeps a one-line reason for the app; the full traceback goes to
            # the terminal running lumencell, which is the only place to debug from.
            print(f"lumencell: {task.kind} job {task.job_id} failed", file=sys.stderr)
            traceback.print_exc()
            detail = "".join(traceback.format_exception_only(type(exc), exc)).strip()
            events.put(("failed", task.job_id, detail))
