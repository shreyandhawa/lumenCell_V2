# Handoff: lumenCell V2, Milestone 1

Status as of 2026-11-13, commit `bfb2539` on `main` (67 commits, all pushed). Phase 1 Tasks 8 to 15 sit in open stacked PRs #14 to #21, not yet merged.
Read this first if you are picking the project up; the README covers usage.

## What exists

A local Python server plus a web app for exploring pseudotime and lineage
trajectories in `.h5ad` single-cell data. It is a ground-up rewrite of V1
(`/Users/shreyan/lumenCell`, a cellxgene fork, left untouched). V1 was dropped
because, measured against ground truth, it inferred on the 2-D UMAP and
auto-picked a branch tip as the root, so its default view had pseudotime
Spearman of about −0.2.

The full loop works and is covered end to end:
1. `lumencell launch data.h5ad` starts the server, and preprocessing runs in a
   background worker.
2. The user picks the root.
3. The user runs a method.
4. The user inspects lineages, gene trends and agreement between runs.

| Plan step | State |
|---|---|
| S0 scaffold, S1 API contract (pydantic → OpenAPI → TS types, drift check) | done |
| A1 synthetic generator + metrics, B1 loader + prep, A2 engine + roots + PAGA/DPT | done |
| B2 stage cache + run store + job worker, B3 API + `launch`/`prepare`/`roots`/`run`/`export` | done |
| C1–C3 web shell, root picker, runs, overlay, lineage tree, trust panel | done |
| I1 Playwright end-to-end test (in `make ci`) | done |
| A3 Palantir, A4 scFates adapters | done |
| A5 trends, heatmap, compare + C4/C5 web views | done (see deviations) |
| A6 real benchmark | done, acceptance **not** met (see results) |
| I2 packaging + 100k-cell performance budgets | budgets measured for prep, `paga_dpt` and the server, all within budget (`docs/backend-robustness-2.md`); Palantir, scFates and first paint are not in that table; packaging **not started** |

## How to work on it

```bash
make install        # uv sync --all-extras, npm ci (Python 3.13 via uv; Node 22+)
make ci             # ruff, mypy --strict, pytest (~65 s), contract check, web lint/types/tests/build, Playwright
make py-test-all    # also the slow tests (full synthetic matrix, calibration)
make bench-synth    # synthetic matrix vs benchmarks/baseline.json (fails on regressions)
make bench-baseline # re-measure and rewrite the baseline (commit it deliberately)
uv run lumencell bench fetch all   # real datasets, ~708 MB, into ~/Library/Caches/lumencell/datasets
make bench-real     # real-data release gate, writes benchmarks/results/real-<date>.json
uv run lumencell bench report --readme   # regenerate the README benchmark tables
make dev-server     # API on :8765 with a synthetic dataset; make dev-web for the app on :5173
make contract       # after changing api/models.py: regenerate OpenAPI and TS types
```

Test counts at handoff: 762 Python tests run in `make ci`, plus 24 slow or
real-data tests that are deselected there. The web side has 36 unit and
component tests and 1 end-to-end test.

## Working rules (agreed with the owner)

- **Commit as you go: plan, test, commit.** Each commit is one working piece,
  pushed to `main`, and nothing is committed with a failing check. Gate the
  commit on make's exit code, not on piped output: a `make … | grep` pipeline
  once let a red commit through.
- **No Claude attribution anywhere.** No `Co-Authored-By` lines and no
  "Generated with" notes. Commits are authored as `shreyandhawa`.
- **Benchmark thresholds are never lowered to fit results.** Every change to a
  threshold or a metric definition is logged, with its reason, in
  `benchmarks/thresholds.toml`. Defaults may be tuned on synthetic data and
  **pancreas** only; bone marrow and gastrulation are held out.
- **A method's `experimental` flag must match whether it misses any synthetic
  gate.** A slow golden test enforces this; update the method's
  `experimental_note` when its results change.

## Map of the code

- `src/lumencell/data/`: loader. Detects whether the file holds counts,
  log-normalised or scaled data; finds gene symbols, embeddings and a stored
  root; fingerprints the file. The server loads the dataset with
  `keep_matrix=False` and serves expression from the prep cache.
- `src/lumencell/prep/`: declared stages (normalize → HVG → PCA → kNN →
  Leiden → diffmap + multiscale → PAGA → UMAP only if the file has no
  embedding). `cache.py` caches each stage on disk, keyed by its parameters
  and its inputs.
- `src/lumencell/engine/`:
  - `schema.py`: the result type every method returns, with `validate()`
  - `roots.py`: explicit roots and evidence-backed suggestions
  - `methods/`: the three method adapters
  - `topology.py`: the oriented PAGA tree
  - `projection.py`: draws trajectories on any embedding; fitted curves if
    the method has them, else binned
  - `trends.py`, `compare.py`
  - `store.py`: runs on disk; `cells.lcb` is served as-is
- `src/lumencell/jobs/`: a persistent spawn worker. One job at a time, with
  progress events, deduplication of identical jobs, and cancel (which
  terminates and respawns the worker; safe because the cache and run store
  only ever write by atomic rename). Each job kind has a time limit
  (`DEFAULT_TIME_LIMITS` in `manager.py`, overridable with
  `LUMENCELL_JOB_TIME_LIMITS`); a job past its limit, or a worker that dies,
  fails and the worker is replaced. Background jobs (`score_clonal`,
  `stability`) queue behind every other kind, so `prepare` and `infer` run first.
- `src/lumencell/api/`:
  - `models.py`: the contract
  - `binary.py`: LCB1, the binary format for per-cell arrays
  - `session.py`: per-dataset server state
  - `routes.py`: thin route handlers
- `src/lumencell/bench/`: synthetic generator, truth, metrics, gates,
  real-dataset registry, runners, README report.
- `web/src/`:
  - `api/`: client, TanStack Query hooks, LCB1 decoder, server-sent events
  - `components/`: deck.gl scatter, sidebar, trajectory panel, gene panel
  - `state/ui.ts`: Zustand store for UI state
- `contracts/`: `openapi.json` and the LCB1 fixture. Both test suites decode the
  same fixture.

## Decisions worth knowing before changing things

- **Roots are always explicit.** There is no automatic root anywhere, and the
  API returns `409 root_required` when none is set.
- **Inference never uses the UMAP.** UMAP and the other 2-D embeddings are
  display only. The trust panel shows each curve's fidelity on the current
  embedding.
- **Palantir departs from its defaults in two places** (both measured, both
  documented in `palantir_adapter.py`):
  - `use_early_cell_as_start=True`: by default Palantir swaps the given start
    for a diffusion-component extreme, which can be a branch tip.
  - Terminal states come from the user, or else from the leaves of the PAGA
    tree. Palantir's automatic detection failed on most synthetic cases.
- **scFates settings** (`scfates_adapter.py`):
  - It is fitted in the multiscale diffusion space; on 30 PCs it grew many
    spurious tips.
  - Side branches shorter than 5 points are pruned.
  - The root follows `_root_point`: on a dead-end arm, the arm's tip; between
    forks, the root cell's own point. This rule was found on real bone marrow,
    where stem cells sit in the middle of the tree.
- **Branch metrics only count cells at least 0.15 past their fork** (in true
  pseudotime), where branches can be told apart at all. The root gate uses the
  root's excess pseudotime over the truth. Both are logged in
  `thresholds.toml`.

## Benchmark results (also in the README, generated)

**Synthetic** (12 cases). All three methods are flagged experimental:
- **PAGA + DPT** fails 7 cases: branch accuracy on three-way splits, and global
  timing when branches have different lengths.
- **Palantir** fails 2 cases: the root's excess pseudotime.
- **scFates** fails 2 cases: global timing on two-fork data.

**Real data**, one row per dataset:

| Dataset | Passes | Fails |
|---|---|---|
| pancreas | Palantir | PAGA + DPT (Alpha split across two lineages), scFates (branch accuracy 0.755 against 0.80) |
| bone marrow | PAGA + DPT, scFates | Palantir (its generous branch masks put monocytes on two lineages) |
| gastrulation | PAGA + DPT | Palantir and scFates (spurious branches split Erythroid3: terminal share 0.59 / 0.61) |

The plan's bar was two passing methods per dataset; only bone marrow meets it.
Pseudotime order against the stages is good everywhere (0.72–0.95). The
failures are all about branch structure.

## Deviations from the plan

- **Run comparison** is a panel (rank agreement, matched lineages, cells
  coloured by disagreement), not two synced scatter plots.
- **Terminal names and method parameters** can't be set in the web app yet.
  The API, the engine and `lumencell run` all support them; `/methods` already
  serves each method's JSON Schema for building a form.
- **Reference benchmark scores** are not yet shown in the trust panel
  (`reference_scores` is empty).
- **CI runs on GitHub Actions** (`.github/workflows/ci.yml`): the `make ci`
  steps on Ubuntu with Python 3.13 from `.python-version`. There is still no
  Python 3.14 job.
- **The browser mock layer (MSW)** is used only in component tests, not as a
  standalone mock server.

## Open work, in suggested order

1. **I2 packaging.** `src/lumencell/web_dist/` is gitignored, so the wheel needs
   it declared under `[tool.hatch.build.targets.wheel] artifacts`. Done when a
   clean venv `pip install ".[all]"` followed by `lumencell launch` serves the
   app.
2. **Performance at 100k cells.** Generate a large synthetic dataset with
   `lumencell synth --n-cells 100000` (the generator already works in chunks)
   and record against the plan's budgets: prep under 8 min, `paga_dpt` under
   30 s, Palantir and scFates under 10 min, first paint under 1.5 s, one gene
   under 50 ms, heatmap under 2 s.
3. **Branching, tuned on pancreas only.** Candidates:
   - Leiden resolution for PAGA, or refining cells at forks.
   - The `q` and `eps` settings of Palantir's `select_branch_cells`.
   - scFates `nodes` and `minbranchlength`.
   - Handling spurious tips on single-lineage data.

   Re-run `make bench-synth` (no regressions allowed) and `make bench-real`, and
   log every change.
4. **App inputs:** terminals and a parameter form, plus reference scores in the
   trust panel.
5. **Repo setup:** done. `.github/workflows/ci.yml` mirrors `make ci`.

## Gotchas met along the way

- **pandas 3 returns read-only arrays.** `to_numpy()` and `np.asarray` on pandas
  data can hand back read-only views; copy with `np.array(...)` before writing
  into them.
- **MSW 3 renamed an option.** Unhandled requests are configured with
  `onUnhandledFrame`, not `onUnhandledRequest`.
- **TypeScript is pinned to `~6.0`,** because typescript-eslint 8.72 needs a
  version below 6.1.
- **python.org Python on macOS has no CA certificates.** Downloads use
  `certifi`'s bundle; certificate verification is never turned off.
- **`lumencell launch` moves to the next port if one is busy.** Use
  `--strict-port` when something depends on the port, such as the Vite proxy or
  Playwright.
- **The 2019-era datasets trigger anndata `OldFormatWarning`.** It is harmless
  and filtered in the real-data tests.
- **Starlette suggests a different HTTP client package for TestClient.** It was
  not installed (the installer refused an unfamiliar package); the warning is
  filtered and `httpx` is kept.
