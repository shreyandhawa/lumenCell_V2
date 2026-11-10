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
