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
