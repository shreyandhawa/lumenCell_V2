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
"""

import multiprocessing as mp
import queue
import threading
import time
import uuid
from collections import deque
from dataclasses import dataclass, field
from typing import Any

from lumencell.jobs.worker import Task, worker_main

ACTIVE = ("queued", "running")


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
    def __init__(self) -> None:
        self._ctx = mp.get_context("spawn")
        self._lock = threading.RLock()
        self._wake = threading.Condition(self._lock)
        self._jobs: dict[str, Job] = {}
        self._pending: deque[str] = deque()
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
                self._worker.terminate()
        for t in self._threads:
            t.join(timeout)

    # ---- public API -------------------------------------------------------------------

    def submit(self, kind: str, key: str, payload: dict[str, Any]) -> Job:
        with self._lock:
            for job in self._jobs.values():
                if job.key == key and job.state in ACTIVE:
                    return job
            job = Job(job_id=uuid.uuid4().hex[:12], kind=kind, key=key, payload=payload)
            self._jobs[job.job_id] = job
            self._pending.append(job.job_id)
            self._wake.notify_all()
            return job

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
                self._pending.remove(job_id)
                self._finish(job, "cancelled")
                return job
            # Running: the only safe way to stop numba/igraph code is to end the process.
            self._worker.terminate()
            self._worker.join(5)
            self._running = None
            self._finish(job, "cancelled")
            self._spawn()
            self._wake.notify_all()
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

    def _finish(self, job: Job, state: str, error: str | None = None) -> None:
        job.state = state
        job.error = error
        job.finished = time.time()
        self._publish(job)
        self._wake.notify_all()

    def _dispatch_loop(self) -> None:
        with self._lock:
            while not self._stopping:
                if self._running is None and self._pending:
                    if self._worker is None or not self._worker.is_alive():
                        # It died between jobs (killed, or out of memory after the last
                        # one). Nothing was running, so nobody noticed until now.
                        self._spawn()
                    job_id = self._pending.popleft()
                    job = self._jobs[job_id]
                    self._running = job_id
                    job.state = "running"
                    job.started = time.time()
                    self._publish(job)
                    self._tasks.put(Task(job_id, job.kind, job.payload))
                    continue
                self._wake.wait(0.2)
                self._check_worker()

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
                elif kind == "failed":
                    self._running = None
                    self._finish(job, "failed", msg[2])
