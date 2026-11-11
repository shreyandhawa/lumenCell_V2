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
    from lumencell.prep.cache import PrepCache, dataset_dir

    ds = _dataset(payload["path"], payload.get("x_is"))
    params = payload["prep_params"]
    emit("load", 0, 1)
    prep = PrepCache(ds, _cache_root(payload)).run(params)
    ctx = InferenceContext.from_prep(ds, prep, params)
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
        rid, result, terminals=payload.get("terminals")
    )
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
