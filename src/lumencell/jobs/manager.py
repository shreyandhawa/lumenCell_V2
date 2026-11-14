"""Job manager: queues jobs for the worker process and tracks their state.

* One worker process at a time; a dispatcher thread hands it the next queued
  job only when it is idle, so cancelling a queued job simply removes it.
* Cancelling the running job terminates the worker and starts a fresh one.
* A listener thread turns worker events into job state and fans them out to
  subscribers (the SSE endpoint reads from a subscriber queue).
* Submitting a job whose key matches a queued or running job returns that job
  instead of starting a second copy.
* If the worker dies mid-job (a crash in compiled code), the job fails with a
  clear message and a new worker is started. If it dies between jobs, a new
  one is started before the next job is handed over.
* A job whose payload cannot be pickled is refused at submit (PayloadError):
  the queue would drop it silently and the job would look running for ever.
* Two priorities. Prepare, infer and every other kind run ahead of the queued
  background kinds (BACKGROUND_KINDS: clonal scoring and the stability check).
  A running job is never pre-empted, so a background job that already started
  finishes first.
* add_done_callback registers a function that runs once for each job that ends
  "done" (not failed, cancelled, timed out or lost to a worker death). It runs
  on the listener thread after the lock is released. Callbacks must be quick and
  may only submit jobs: while one runs the listener reads no worker events, so a
  slow one delays the next dispatch, and calling wait() from one deadlocks the
  listener (it is the thread that ends jobs). wait() and subscribers are told a
  job is done before its callbacks run. Callbacks are never removed.
* Each kind of job has a time limit. A job that runs past it fails with a
  message saying so, and the worker is replaced, as for a cancel. The limits
  come from DEFAULT_TIME_LIMITS, then the LUMENCELL_JOB_TIME_LIMITS environment
  variable ("prepare=7200,infer=10800", in seconds), then the constructor.
"""

import math
import multiprocessing as mp
import os
import queue
import sys
import threading
import time
import traceback
import uuid
from collections import deque
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from multiprocessing.reduction import ForkingPickler
from typing import Any

from lumencell.jobs.worker import HANDLERS, Task, worker_main

ACTIVE = ("queued", "running")
# Queued behind every other kind: slow work that nobody is waiting at the screen for.
BACKGROUND_KINDS: frozenset[str] = frozenset({"score_clonal", "stability"})
# Generous: about twenty times the 100k-cell budgets (prep 8 min, slowest method 10 min).
DEFAULT_TIME_LIMITS: dict[str, float] = {
    "prepare": 3 * 3600,
    "infer": 3 * 3600,
    "score_clonal": 3600,
    "sleep": 600,
}
FALLBACK_TIME_LIMIT = 3 * 3600.0
TIME_LIMITS_ENV = "LUMENCELL_JOB_TIME_LIMITS"


class PayloadError(ValueError):
    """A job's payload cannot be sent to the worker process."""


def _checked_limits(limits: Mapping[str, object], where: str) -> dict[str, float]:
    """Every limit names a job kind and is a finite number of seconds above 0."""
    out: dict[str, float] = {}
    for kind, value in limits.items():
        number = isinstance(value, int | float) and not isinstance(value, bool)
        seconds = float(value) if number else math.nan  # type: ignore[arg-type]
        if kind not in HANDLERS or not (math.isfinite(seconds) and seconds > 0):
            raise ValueError(
                f"{where}: the time limit for {kind!r} must be a number of seconds above 0, not {value!r} "
                f"(job kinds: {', '.join(HANDLERS)})"
            )
        out[kind] = seconds
    return out


def _env_time_limits() -> dict[str, float]:
    raw = os.environ.get(TIME_LIMITS_ENV, "").strip()
    if not raw:
        return {}
    parsed: dict[str, object] = {}
    for part in raw.split(","):
        kind, eq, value = part.strip().partition("=")
        if not eq:
            raise ValueError(
                f"{TIME_LIMITS_ENV} must look like 'prepare=7200,infer=10800' (seconds), not {raw!r}"
            )
        try:
            parsed[kind.strip()] = float(value)
        except ValueError:
            parsed[kind.strip()] = value.strip()
    return _checked_limits(parsed, TIME_LIMITS_ENV)


def _duration(seconds: float) -> str:
    if seconds < 120:
        return f"{seconds:g} s"
    if seconds < 2 * 3600:
        return f"{seconds / 60:g} min"
    return f"{seconds / 3600:g} h"


@dataclass
class Job:
    job_id: str
    kind: str
    key: str
    payload: dict[str, Any]
    state: str = "queued"  # queued | running | done | failed | cancelled
    stage: str | None = None
    i: int = 0
    n: int = 0
    message: str | None = None
    error: str | None = None
    result_id: str | None = None
    created: float = field(default_factory=time.time)
    started: float | None = None
    finished: float | None = None

    @property
    def fraction(self) -> float:
        if self.state == "done":
            return 1.0
        return self.i / self.n if self.n else 0.0

    def event(self) -> dict[str, Any]:
        now = self.finished or time.time()
        return {
            "job_id": self.job_id,
            "state": self.state,
            "stage": self.stage,
            "i": self.i,
            "n": self.n,
            "fraction": round(self.fraction, 4),
            "message": self.error or self.message,
            "elapsed": round(now - (self.started or now), 3),
            "eta": None,
        }


class JobManager:
    def __init__(self, time_limits: Mapping[str, float] | None = None) -> None:
        self._time_limits = {
            **DEFAULT_TIME_LIMITS,
            **_env_time_limits(),
            **_checked_limits(time_limits or {}, "time_limits"),
        }
        self._ctx = mp.get_context("spawn")
        self._lock = threading.RLock()
        self._wake = threading.Condition(self._lock)
        self._jobs: dict[str, Job] = {}
        self._foreground: deque[str] = deque()
        self._background: deque[str] = deque()
        self._callbacks: list[Callable[[Job], None]] = []
        self._running: str | None = None
        self._subscribers: dict[str, list[queue.Queue[dict[str, Any]]]] = {}
        self._stopping = False
        self._worker: Any = None
        self._tasks: Any = None
        self._events: Any = None
        self._threads: list[threading.Thread] = []

    # ---- lifecycle ---------------------------------------------------------------------

    def start(self) -> "JobManager":
        self._spawn()
        for target in (self._dispatch_loop, self._listen_loop):
            t = threading.Thread(target=target, daemon=True, name=f"lumencell-{target.__name__}")
            t.start()
            self._threads.append(t)
        return self

    def _spawn(self) -> None:
        self._tasks = self._ctx.Queue()
        self._events = self._ctx.Queue()
        self._worker = self._ctx.Process(
            target=worker_main, args=(self._tasks, self._events), daemon=True, name="lumencell-worker"
        )
        self._worker.start()

    def shutdown(self, timeout: float = 5.0) -> None:
        with self._lock:
            self._stopping = True
            self._wake.notify_all()
        if self._worker is not None and self._worker.is_alive():
            self._tasks.put(None)
            self._worker.join(timeout)
            if self._worker.is_alive():
                self._stop_worker()
        for t in self._threads:
            t.join(timeout)

    # ---- public API -------------------------------------------------------------------

    def time_limit(self, kind: str) -> float:
        """Seconds a job of this kind may run before it is stopped."""
        return self._time_limits.get(kind, FALLBACK_TIME_LIMIT)

    def submit(self, kind: str, key: str, payload: dict[str, Any]) -> Job:
        try:
            # What the task queue does in a background thread, where a failure is
            # only printed and the job would never reach the worker.
            ForkingPickler.dumps(Task("", kind, payload))
        except Exception as exc:
            raise PayloadError(
                f"the {kind} job cannot be sent to the worker: its payload cannot be pickled ({exc})"
            ) from exc
        with self._lock:
            for job in self._jobs.values():
                if job.key == key and job.state in ACTIVE:
                    return job
            job = Job(job_id=uuid.uuid4().hex[:12], kind=kind, key=key, payload=payload)
            self._jobs[job.job_id] = job
            (self._background if kind in BACKGROUND_KINDS else self._foreground).append(job.job_id)
            self._wake.notify_all()
            return job

    def add_done_callback(self, callback: Callable[[Job], None]) -> None:
        """Call `callback(job)` once for every job that reaches "done" (not failed,
        cancelled, timed out or lost to a worker death). It runs on the listener thread
        with the lock released. Callbacks must be quick and may only submit jobs: never
        call wait() from one (it deadlocks the listener). They are never removed."""
        with self._lock:
            self._callbacks.append(callback)

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)

    def find(self, key: str) -> Job | None:
        """The most recent job with this key."""
        with self._lock:
            matches = [j for j in self._jobs.values() if j.key == key]
            return max(matches, key=lambda j: j.created) if matches else None

    def of_kind(self, kind: str) -> list[Job]:
        """Every job of this kind, in the order they were submitted."""
        with self._lock:
            return [j for j in self._jobs.values() if j.kind == kind]

    def cancel(self, job_id: str) -> Job | None:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None or job.state not in ACTIVE:
                return job
            if job.state == "queued":
                for waiting in (self._foreground, self._background):
                    if job_id in waiting:
                        waiting.remove(job_id)
                self._finish(job, "cancelled")
                return job
            self._replace_worker()
            self._finish(job, "cancelled")
            return job

    def subscribe(self, job_id: str) -> "queue.Queue[dict[str, Any]]":
        q: queue.Queue[dict[str, Any]] = queue.Queue()
        with self._lock:
            self._subscribers.setdefault(job_id, []).append(q)
            job = self._jobs.get(job_id)
            if job is not None:
                q.put(job.event())
        return q

    def unsubscribe(self, job_id: str, q: "queue.Queue[dict[str, Any]]") -> None:
        with self._lock:
            subs = self._subscribers.get(job_id, [])
            if q in subs:
                subs.remove(q)

    def wait(self, job_id: str, timeout: float = 60.0) -> Job:
        """Block until the job finishes (used by the CLI and tests)."""
        deadline = time.time() + timeout
        with self._lock:
            while True:
                job = self._jobs[job_id]
                if job.state not in ACTIVE:
                    return job
                remaining = deadline - time.time()
                if remaining <= 0:
                    raise TimeoutError(f"job {job_id} still {job.state} after {timeout}s")
                self._wake.wait(min(remaining, 0.5))

    # ---- internals ----------------------------------------------------------------------

    def _publish(self, job: Job) -> None:
        event = job.event()
        for q in self._subscribers.get(job.job_id, []):
            q.put(event)

    def _stop_worker(self) -> None:
        self._worker.terminate()
        self._worker.join(5)
        if self._worker.is_alive():  # it ignores SIGTERM (some libraries install their own handler)
            pid = self._worker.pid
            print(f"lumencell: the worker (pid {pid}) ignored SIGTERM; killing it", file=sys.stderr)
            self._worker.kill()
            self._worker.join(5)

    def _replace_worker(self) -> None:
        """Stop the running job by ending the worker (the only safe way to stop numba
        or igraph code) and start a fresh one. Call with the lock held."""
        self._stop_worker()
        self._running = None
        self._spawn()
        self._wake.notify_all()

    def _finish(self, job: Job, state: str, error: str | None = None) -> None:
        job.state = state
        job.error = error
        job.finished = time.time()
        self._publish(job)
        self._wake.notify_all()

    def _take_next(self) -> str | None:
        """The oldest queued foreground job, else the oldest background one. Call with
        the lock held."""
        if self._foreground:
            return self._foreground.popleft()
        if self._background:
            return self._background.popleft()
        return None

    def _dispatch_loop(self) -> None:
        with self._lock:
            while not self._stopping:
                if self._running is None and (self._foreground or self._background):
                    if self._worker is None or not self._worker.is_alive():
                        # It died between jobs (killed, or out of memory after the last
                        # one). Nothing was running, so nobody noticed until now.
                        self._spawn()
                    job_id = self._take_next()
                    assert job_id is not None
                    job = self._jobs[job_id]
                    self._running = job_id
                    job.state = "running"
                    job.started = time.time()
                    self._publish(job)
                    self._tasks.put(Task(job_id, job.kind, job.payload))
                    continue
                self._wake.wait(0.2)
                self._check_worker()
                self._check_time_limit()

    def _check_time_limit(self) -> None:
        if self._running is None:
            return
        job = self._jobs[self._running]
        limit = self.time_limit(job.kind)
        if job.started is not None and time.time() - job.started > limit:
            self._replace_worker()
            message = (
                f"timed out: the {job.kind} job ran longer than its time limit of {_duration(limit)} "
                "and was stopped"
            )
            print(f"lumencell: job {job.job_id} {message}", file=sys.stderr)
            self._finish(job, "failed", message)

    def _check_worker(self) -> None:
        if self._running is not None and self._worker is not None and not self._worker.is_alive():
            job = self._jobs[self._running]
            self._running = None
            self._finish(
                job, "failed", f"the worker process stopped unexpectedly (exit code {self._worker.exitcode})"
            )
            self._spawn()

    def _listen_loop(self) -> None:
        while True:
            with self._lock:
                if self._stopping:
                    return
                events = self._events
            try:
                msg = events.get(timeout=0.2)
            except (queue.Empty, OSError, ValueError, EOFError):
                continue
            finished: Job | None = None
            callbacks: list[Callable[[Job], None]] = []
            with self._lock:
                job = self._jobs.get(msg[1])
                if job is None or job.state not in ACTIVE:
                    continue  # e.g. a late event from a job that was cancelled
                kind = msg[0]
                if kind == "progress":
                    _, _, job.stage, job.i, job.n = msg
                    self._publish(job)
                elif kind == "done":
                    job.result_id = msg[2]
                    self._running = None
                    self._finish(job, "done")
                    finished = job
                    callbacks = list(self._callbacks)
                elif kind == "failed":
                    self._running = None
                    self._finish(job, "failed", msg[2])
            if finished is not None:
                for callback in callbacks:
                    try:
                        callback(finished)
                    except BaseException:  # SystemExit too: the listener must outlive any callback
                        print(f"lumencell: a callback after job {finished.job_id} failed", file=sys.stderr)
                        traceback.print_exc()
