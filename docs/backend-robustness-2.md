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

**Fixed since, in PR #2** (`d660010`, `5c1410e`). Found while testing, older
than this branch and not one of the six issues: stopping `lumencell launch` with
SIGTERM left its worker running for ever. uvicorn re-raises the signal after
shutting down, so `launch`'s `finally: jobs.shutdown()` never ran, and the
worker kept waiting on its queue. `tests/api/test_server.py` stops the server
exactly this way, so every test run left two orphaned process pairs behind. Now
`launch` turns SIGTERM into an ordinary exit (code 143), so the shutdown runs,
and the worker also exits by itself when the process that started it dies, which
covers SIGKILL. `test_stopping_launch_stops_its_job_worker` checks both signals.

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

## Issue 2: memory at 100k cells

**Measured first.** I generated a 100,000-cell × 2,000-gene bifurcation into the
worktree with `lumencell synth .lumencell/perf/synth100k.h5ad --n-cells 100000`
(774 MB on disk, 864 MB in the cache once prepared).

Every run used `LUMENCELL_CACHE_DIR` and `LUMENCELL_DATA_DIR` inside the
worktree, on an Apple-silicon Mac, with the other agent's processes sharing the
machine:

- prep: `/usr/bin/time -l lumencell prepare`;
- `paga_dpt`: `lumencell run --method paga_dpt`, with the method time taken from
  the run's `runtime_seconds`;
- the server: a script that builds a `Session` the way `launch` does, then
  measures memory (RSS, and tracemalloc for the heatmap) and request timings in
  that process.

| | budget | before | after |
|---|---|---|---|
| prep, wall | < 8 min | 163 s | 135 s (see note) |
| prep, peak RSS | | 2.79 GB | 3.43 GB (see note) |
| `paga_dpt` method time | < 30 s | 0.09 s | 0.08 s |
| `lumencell run`, end to end | | 2.5 s | 2.4 s |
| server: expression copies held | one | 3 × 768 MB | **1 × 767 MB** |
| server: RSS after loading prep | | 1.93 GB | 2.22 GB |
| server: first gene (builds the by-gene copy) | | 305 ms | 1.1 ms (built at load) |
| server: one gene, warm (median of 20) | < 50 ms | 0.09 ms | 0.09 ms |
| server: uncached heatmap, top 50 | < 2 s | **2.97 s** | **0.79 s** |
| server: heatmap allocations (tracemalloc peak) | | **3,984 MB** | **21.7 MB** |
| server: trends, one gene | | 3.4 ms | 1.9 ms |
| server: RSS after the heatmap | | 3.77 GB | 2.09 GB |
| server: peak RSS of the process | | 4.45 GB | 2.23 GB |

Note on the prep and CLI rows: these paths load the full matrix, as they did
before, and this branch doesn't change their memory. Their peak RSS varies a lot
between identical runs. Three `lumencell run`s of the old code peaked at 1,754,
2,268 and 2,276 MB, and three of the new code at 2,126, 2,276 and 2,276 MB, so
the prep difference above is run-to-run spread, not this change. The `neighbors`
stage took 39 s in the first prep and 15 s in the second (compilation caches
warming up); `pca` took 71 s both times.

**Now:**

- `engine/trends.py` bins, averages and scores a block of genes at a time.
  Blocks are contiguous columns with at most `GENE_BLOCK_NNZ` = 500,000 nonzeros,
  and each block is converted to float64 on its own. No request copies or
  converts the whole matrix any more: no whole-matrix `astype`, no CSC→CSR
  conversion inside the sparse product, and no row subset. The results are the
  same operations in the same order. Input stored by cell is refused with a
  `TypeError` instead of silently losing the memory bound. A lineage whose cells
  have no pseudotime now has no trend (`null` in the API) and r2 0, where it used
  to raise `ZeroDivisionError` (a 500).
- `api/dataset_view.py` converts the cached lognorm to the by-gene (CSC) layout
  once, when prep is loaded, and drops the by-cell copy. `expression_matrix()`
  returns that one matrix.
- `data/loader.py`: `Dataset.matrix` is optional. `load(..., keep_matrix=False)`
  keeps the shape and a per-cell `genes_detected` and drops the matrix. The
  server's `Session` loads this way; the worker, the CLI and the benchmarks load
  in full, as before. Code that needs the matrix calls `require_matrix()`, which
  raises a clear `DatasetError` instead of failing somewhere downstream.
- `InferenceContext.counts` is replaced by `genes_detected`, which is all the
  gene-count root suggestion ever used.

**Verified by:**

- `tests/unit/test_trends.py`:
  - a frozen copy of the old implementation is reproduced on the shared
    synthetic fixture (with missing pseudotimes) and the toy data, at block sizes
    of 1 and 1,000 nonzeros and at the default: same genes in the same order,
    values within 1e-6;
  - the blocks cover every gene exactly once, within their budget;
  - a by-cell matrix is refused;
  - an empty lineage has no trend rather than an error;
  - a hypothesis property: results don't depend on the block size, over random
    matrices, missing pseudotimes, sparse lineages and named genes;
  - a slow test: an uncached heatmap on a 40,000 × 2,000 matrix allocates under
    0.2 copies. It allocated 883 MB against a 192 MB matrix before.
- `tests/api/test_routes.py`: walking everything reachable from the session
  finds one copy of expression, stored by gene in float32. Before there were 3.
  The server still offers the gene-count root suggestion without the matrix.
- `tests/unit/test_loader.py`: a dataset loaded without its matrix has the same
  shape, gene counts and fingerprint, and refuses to be prepared.

**Still open.**

- The server's peak RSS (2.2 GB) is now mostly the load itself: anndata reads
  `X` in full and then converts it to float32 CSR before the matrix is dropped.
  Reading the file in backed mode would avoid that.
- Making the by-gene copy at load takes a moment in the first request that loads
  prep. It was already in a request, only later.
- `_smooth` breaks for fewer than 5 bins. This is older than the branch, and
  every caller uses 50.

**Fold into HANDOFF.** The 100k-cell numbers above. The server's `Dataset` has no
matrix (`keep_matrix=False`), and the expression it serves is the cache's lognorm,
held once and by gene.

## Issue 6: CI on GitHub Actions

**Now.** `.github/workflows/ci.yml` runs `make ci` on every pull request and on
pushes to `main`, and nowhere else. It runs on `ubuntu-latest` and:

- installs Python 3.13 with uv (`astral-sh/setup-uv`, with the uv cache keyed
  on `uv.lock`) and Node 22 (`actions/setup-node`, with the npm cache keyed on
  `web/package-lock.json`);
- runs `make install`;
- restores Playwright's Chromium from a cache keyed on the lockfile. On a miss
  it downloads Chromium from Playwright's own server (no apt) and saves the cache
  straight away, so a later failure still leaves one. This is how `main` gets a
  cache: `actions/cache` saved only when the whole job passed, and a pull
  request's cache can't be used on `main` (PR #3);
- installs Chromium's system libraries with apt on every run: 3 tries, each
  capped at 3 minutes, then an error that says apt could not reach the mirrors
  (PR #3);
- runs `make ci`, and keeps `web/test-results` as an artifact when it fails.

- `UV_LOCKED=1`, so a `uv.lock` that is out of date fails instead of being
  resolved again.
- An explicit import of palantir and scFates. Their tests skip when the import
  fails, which would otherwise pass CI without them.
- Python comes from `.python-version` (3.13).
- A newer push to a pull request cancels its run in progress. Every push to
  `main` keeps its own result. The job has a 45-minute timeout.

**Suggested, in files this branch may not change.**

- `forbidOnly: !!process.env.CI` in `web/playwright.config.ts`, so that a
  committed `test.only` can't pass CI with fewer tests.
- The Makefile's `contract-check` prints "openapi.json is stale" even when the
  generator itself crashes.

**Verified by:**

- parsing the YAML and checking that its only triggers are `pull_request` and
  `push` to `main`;
- the full `make ci` in the worktree before the PR (exit 0);
- PR #1 on GitHub. The first run (43028851191) failed in one existing test,
  `test_unusable_file_is_reported_clearly`, and only there: with
  `GITHUB_ACTIONS` set, typer forces rich's terminal styling, which splits
  `--x-is` with colour codes, so the plain-text assertion missed it. Locally,
  `GITHUB_ACTIONS=true` reproduced it and `CI=true` did not. The `make ci`
  step now runs with the variable removed (`env -u GITHUB_ACTIONS make ci`).
  The test itself is unchanged. The second run (43030038492) passed in 6 min:
  lint, mypy, 351 Python tests, the contract check, 36 web tests, the build
  and the Playwright test.

**After the merge.** `main`'s first push run (43040969341) failed all 3
attempts in "Install Playwright Chromium", before `make ci` started. From about
13:48 UTC on 2026-11-11, GitHub's Ubuntu runners could not resolve Ubuntu's
package mirrors (actions/runner-images#15296, a Canonical outage that GitHub
called mitigated at 18:14 UTC). No test failed. PR #3 changed only `ci.yml`, to
the Playwright steps described above. Its run (43052546403) failed twice on the
outage, with the retries and the closing error working as intended, and passed
on its third attempt, after the outage. Then, on GitHub:

- 43130807421, `main` at `79e483e` (PR #3): the first green run on `main`. It
  downloaded Chromium and saved `main`'s first Playwright cache.
- 43130879682, PR #2 with `main` merged into it: green.
- 43131513810, `main` at `5c1410e` (PR #2): green in 5 min, with Chromium
  restored from `main`'s cache. 353 Python tests (12 deselected), 36 web tests
  and the Playwright test passed.
