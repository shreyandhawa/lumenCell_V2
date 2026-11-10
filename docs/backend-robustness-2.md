# Backend robustness, round 2

Notes for the `backend-robustness-2` branch: six backend issues, fixed test-first,
one commit each, in the order 1, 5, 3, 4, 2, 6. Each section says what changed, how
it was verified and what is still open.

## Issue 1: a damaged cache stage heals

**Before.** A stage counted as cached as soon as its `.complete` marker existed. A
stage folder with the marker but a missing or truncated artifact was reported as
`done` by `status()`, `run()` raised `FileNotFoundError`/`ValueError`/`BadZipFile`,
and every data route answered 500 until the cache was deleted by hand.

**Now** (`prep/cache.py`, `api/dataset_view.py`, `api/session.py`):

- `PrepCache.has()` also checks that every artifact the stage declares is on disk,
  so a missing artifact makes the stage `pending`.
- `PrepCache.load()` treats an artifact that cannot be read as damage. It prints
  one line to stderr naming the stage, its folder and the error, removes the stage
  (`invalidate()`: one rename out of the way, then delete), and raises `CacheMiss`.
- Only real damage counts. A missing file, or one that can't be parsed (cut short,
  not a valid archive), does. A `PermissionError` or another I/O error is raised as
  it is, and the stage stays on disk. If deleting a taken-out stage fails, that is
  reported, and the stage still counts as out of the cache.
- `PrepCache.run()` (the worker, and the CLI) computes a stage it could not load.
  `save()` replaces an incomplete copy instead of keeping it, and says so on
  stderr, so the healed stage stays healed on disk. If another writer's complete
  copy lands first, it is kept.
- The server never computes. It reads prep through `PrepCache.load_all()`, which
  raises `CacheMiss` rather than computing. `DatasetView.loaded_prep()` catches
  that, prints a line, clears the readiness latch and starts a prepare job.
  Routes that need prep answer the existing 409 `prep_not_ready`. `/dataset` and
  `/root` keep answering 200, with the prep state `running` and no root, so the
  web app keeps polling until the job has put the stage back. `load_all()` writes
  nothing to the cache, not even the manifest.
- A stage that is still unreadable after one heal is not sent back again in the
  same server session. Prep then shows as `missing`, and the app's Prepare button
  retries it. Without that limit, a disk that keeps damaging files would keep the
  worker busy for ever.

**Verified by** `tests/unit/test_prep_cache.py`:

- a missing artifact;
- a truncated `.npy`, and two kinds of truncated `.npz`. Only that stage is
  recomputed, the artifacts match the clean run, the damage is reported on
  stderr, and the next run is all cache hits;
- `load_all` never computes.

`tests/unit/test_prep_cache.py` also covers:

- a zero-byte artifact;
- `load_all` on a damaged stage;
- a permission error, which keeps the stage;
- a failed delete;
- losing the rename race to another writer.

And by `tests/api/test_cache_healing.py`, whose job manager starts only after the
damaged state has been checked, so the heal job can't race the assertions:

- a truncated or deleted artifact found by a data route answers 409, the server
  starts the prepare job, and the same route answers 200 after it;
- `/dataset` reports `running`;
- every data route answers 409 `prep_not_ready` while healing, and `/embeddings/dc12`
  answers 404 as on a dataset that was never prepared;
- four concurrent requests all get 409 and start exactly one job;
- a stage that stays unreadable is not prepared a second time automatically.

All of these failed before the fix: `status()` said `done`, `run()` raised, and
the routes answered 500.

**Still open.**

- Damage inside an artifact that keeps its length (flipped bytes in a `.npy`)
  loads without error, because nothing checksums artifacts.
- If the worker lands a fresh copy between the server's failed read and its
  rename, the server takes the fresh copy out. The cost is one more recompute,
  never wrong data.
- After a heal job fails, the 409 message still says to wait for the prepare job.
  `/prepare` shows `failed`, which is what the app displays.
- Once the server has loaded prep into memory, damage on disk only matters to the
  worker, which heals it on its next run.

## Issue 5: every root survives being written as text

**Before.** `format_root` wrote `cell:<name>` and `cluster:<key>=<value>`, and
`parse_root` stripped the text and split at the first `=`. So
`ClusterRoot("batch=1", "A")` came back as key `batch` and value `1=A`, and
`CellRoot(" AAAC-1 ")` lost its trailing space. Through the API, a cluster value
with a trailing space started a run from a different cluster without saying so.

**Now** (`engine/roots.py`). `format_root` keeps today's plain text whenever that
text already round-trips. That is: cell names and cluster values that are not
empty and don't end in whitespace, and cluster columns that are not empty and have
no `=`. Whitespace around a column name survives, because it sits mid-text. Everything else is written as `json:` followed by the API's
own shape, for example `json:{"kind": "cluster", "key": "batch=1", "value": "A"}`.
`parse_root` reads that form strictly: exactly those keys, string values, kind
`cell` or `cluster`. Today's parser rejects every text that starts with `json:`, so
no saved root, run id or CLI argument that worked before changes meaning. An
empty name now round-trips too, and is refused when it is resolved (still 422
`invalid_root`).

**Verified by** `tests/unit/test_roots.py`:

- hypothesis: `parse_root(format_root(spec)) == spec` for arbitrary text, with the
  reported cases pinned as examples (soaked once at 20,000 examples);
- today's ordinary forms are unchanged byte for byte;
- wherever today's form round-trips, the new text is identical;
- any text today's parser accepts parses to the same root;
- malformed `json:` texts are refused, absurdly deep JSON included, and any text
  gives a root or a `RootError`, nothing else.

A reviewer caught a first version that escaped cluster columns ending in a space,
which would have changed their run ids. The property now checks against today's
real parser instead of a copy of the new rule.

And by `tests/api/test_routes.py`: a cluster value with a trailing space is a 422
from resolving it, not a run from another cluster, and a root the plain text cannot
carry reaches the worker's payload unchanged. The round-trip property and the API test failed
before the fix.

**Fold into HANDOFF.** Root text forms: `cell:NAME`, `cluster:KEY=VALUE`, `file`,
and `json:{...}` for names the plain form can't carry.

## Issue 3: no job blocks the queue

**Before.** The task queue pickles each task in a background thread. A payload
that could not be pickled only printed a traceback there, so the job showed as
`running` for ever and every job behind it stayed `queued`. That was reproduced
before the fix: after 15 s, the job was still `running`, the next one was
`queued`, and the worker sat idle. No job had a time limit either.

**Now** (`jobs/manager.py`):

- `JobManager.submit` pickles the task exactly as the queue would, before
  queueing it. When that fails, it raises `PayloadError` ("the <kind> job cannot
  be sent to the worker: its payload cannot be pickled (...)") and queues nothing.
- Each job kind has a time limit: 3 h for `prepare` and `infer` (about twenty
  times the 100k-cell budgets), 10 min for `sleep`, and 3 h for any other kind.
  A running job past its limit fails with "timed out: the <kind> job ran longer
  than its time limit of <n> and was stopped". The worker is ended and replaced,
  as for a cancel, and the next job runs.
- Configuring the limits: `JobManager(time_limits={...})` beats the
  `LUMENCELL_JOB_TIME_LIMITS` environment variable, for example
  `prepare=7200,infer=10800` in seconds, which beats the defaults. Both sources
  are checked: the kind must be a real job kind, and the value a finite number of
  seconds above 0. A bad value raises at startup, and `lumencell launch` reports
  it as a one-line error instead of a traceback.
- A timeout is also printed to the terminal. A worker that ignores SIGTERM is
  killed after 5 s, on timeout, cancel or shutdown, so no second worker ever runs
  beside a stuck one.

**Verified by** `tests/unit/test_jobs.py`:

- a sleep job with a 2 s limit fails with the timeout message, its worker is
  gone, and a fresh worker runs the next job;
- a worker that ignores SIGTERM is still stopped, on timeout and on cancel;
- an unpicklable payload raises `PayloadError` at submit, nothing is queued, and
  the next job runs;
- the approved defaults: 3 h for prepare and infer, 10 min for sleep;
- the env var and constructor overrides work, and typos, unknown kinds,
  non-numbers, NaN and non-positive values are refused;
- `launch` reports a bad env value without a traceback.

All four failed before the fix (the constructor had no `time_limits`, and there
was no `PayloadError` or `time_limit`).

**Still open.**

- If a job's own `done` event is already waiting when its limit is hit, the job
  is reported as timed out although its run was saved. The run is still found
  on disk the next time it is asked for.
- Found while testing, older than this branch and not one of the six issues:
  stopping `lumencell launch` with SIGTERM leaves its worker running for ever.
  uvicorn re-raises the signal after shutting down, so `launch`'s
  `finally: jobs.shutdown()` never runs, and the worker keeps waiting on its
  queue. `tests/api/test_server.py` stops the server exactly this way, so every
  test run leaves two orphaned process pairs behind. Those from the main checkout
  and the frontend worktree were left alone, and the ones from this worktree
  were stopped.

**Fold into HANDOFF.** `LUMENCELL_JOB_TIME_LIMITS` sets the per-kind job time
limits (seconds). The defaults are 3 h for prepare and infer.

## Issue 4: `lumencell cache list` and `lumencell cache prune`

**Before.** Nothing ever removed anything from `<cache root>/v1/`. Cancelled jobs
left their temporary folders behind, and the cache could not tell which dataset a
folder belonged to.

**Now:**

- `prep/cache.py`: the manifest records the dataset's name and path. A
  `last-used` stamp file in the dataset folder is touched when a `run()` starts
  (so a long prepare counts as use) and when a `run()` or `load_all()` succeeds.
  A stamp file is used rather than the folder's mtime, so cleaning up inside the
  folder, or macOS dropping a `.DS_Store` there, never counts as use. Caches
  from before this branch fall back to the folder's mtime. The stamp is the only
  thing the server writes to the cache, and a failure there is reported, not
  raised.
- New `prep/housekeeping.py`:
  - `list_cache()`: entries, newest first, each with its fingerprint, name,
    source, size and last-used time.
  - `plan_prune()` works out what to remove: whole datasets with `--all` or past
    `--older-than`, and in the datasets that stay, leftovers of interrupted jobs
    that have been untouched for over an hour.
  - `check_prune()` and `apply_prune()`: every removal is checked before any is
    deleted. Each one is then removed on its own, and the outcome is reported
    (`removed`, `already gone`, `could not remove: <error>`, the last with exit
    code 1). Datasets keep their times when only leftovers are removed from them.
- `lumencell cache list` prints the root, then a row per dataset with its name,
  size, last-used date ("40 days ago"), fingerprint and source file, then a total.
  Caches written before this branch show their fingerprint, because their
  manifest has no name.
- `lumencell cache prune (--older-than 12h|30d|2w | --all) [--dry-run]` needs
  exactly one of the two selectors, prints what it removed or would remove and
  the space freed, and refuses unsafe paths with exit code 1 before anything is
  removed.

**Safety rules.**

- Only real folders directly under `<root>/v1/` are ever removed, and only when
  they have a 24-hex fingerprint name and something the cache wrote (a
  manifest, `stages/` or `runs/`). Things inside those folders can be removed
  too.
- A symbolic link put there by hand is listed but never removed, with a
  message. It doesn't block the rest of the prune.
- Anything whose resolved path is outside `<root>/v1/` is refused before anything
  is deleted. `datasets/` (benchmark downloads) and anything else under the root
  are never candidates.
- Nothing that contains the state directory, or is inside it, is ever removed. A
  state directory equal to the cache root is fine.
- Leftovers are matched by the exact temporary names the cache and the run store
  use, so unrelated dotfiles are never removed.
- Prune never ran against this machine's real cache.

**Verified by** `tests/unit/test_cache_admin.py`, on a fake cache under
`tmp_path` with three datasets: one 40 days old, one in use with old and fresh
temps, and one with an old-style manifest. The tests cover:

- what `list` reports;
- what `--older-than 30d` removes (the old dataset and the old temps, not the
  fresh temp);
- `--dry-run` changes nothing;
- `--all` leaves `datasets/` and the state file;
- a symlinked entry is refused and nothing is removed;
- `apply_prune` refuses anything outside the dataset folders;
- the flags and ages that are accepted;
- a missing cache root.

The review added these tests:

- removing a leftover keeps its dataset's last-used time;
- a failed removal is reported while the rest go ahead;
- the state directory is protected, whether it is the cache root or inside a
  dataset;
- `.DS_Store` files and a bare 24-hex folder survive;
- listing survives a file vanishing mid-walk;
- the dry-run snapshot includes mtimes.

`tests/unit/test_prep_cache.py` checks:

- the manifest names the dataset;
- `run()` and `load_all()` mark the cache as used;
- a dataset being prepared counts as in use before its first stage starts. All of these failed first:
`No such command 'cache'`, no `dataset` in the manifest, and an mtime that never
moved.

**Still open.**

- `--all` also removes a dataset that a running server is showing. The server
  keeps its prep in memory, but its runs are read from disk.
- The check and the delete are separate steps. `rmtree` refuses to follow links,
  but a parent folder swapped for a link between the two steps isn't guarded
  against.
- Sizes are file bytes, not disk blocks.

**Fold into HANDOFF.** `lumencell cache list` and `lumencell cache prune` exist.
Prune only touches `<cache root>/v1/<fingerprint>` folders. Don't prune with
`--all` while `lumencell launch` is running on one of those datasets.
