# Phase 2: T-cell Robustness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give every run a stability block (2 reseeds, 3 subsamples at 80%), a plain tree-support verdict whose thresholds come from synthetic data only, and a report-only cell-cycle check. The worker computes them, or the server reads them from cached prep, and they are served with the run.

**Architecture:** A report-only `cell_cycle` prep stage scores the Tirosh 2016 lists and stays out of `prep_key`, so no run id changes. `engine/stability.py` re-runs a finished run five times: reseeds rerun Leiden, PAGA and the method on all cells, and subsamples rerun kNN through PAGA, all in memory with `run_stage`. It compares each replicate with the reference on shared cells through an extended `engine/compare.py`. A background `stability` job stores the result through `RunStore.write_trust`. `engine/tree_support.py` turns coverage, PAGA tree-likeness, stability and the clonal score into `supported`, `weak`, `not_supported` or `provisional`. Its τs come from a fixed quantile rule on the synthetic matrix against four negative controls (`bench/synth_controls.py`), and they are verified on fresh seeds before any real-data use.

**Tech Stack:** Python 3.13 via uv; numpy, scipy (`scipy.sparse`, `scipy.optimize.linear_sum_assignment`, `scipy.stats.spearmanr`), pandas, scanpy (`sc.tl.score_genes_cell_cycle`); FastAPI + pydantic v2; pytest. `web/` gets only regenerated types and one fixture update.

**Spec:** `docs/superpowers/specs/2026-11-12-tcell-clonotypes-design.md`: §2 (constraints), §5.1–5.5 (this phase), and §7 (the units `bench/synth_controls.py`, `engine/stability.py`, `engine/tree_support.py`, `trust/constants.py` and the `cell_cycle` part of `prep/stages.py`, plus the `stability` branch of the data flow). Read §5 before any task. If the spec and this plan conflict, the spec wins. Report the conflict instead of choosing.

**Baseline:** `main` after every phase 1 task has merged (phase 1 plan `docs/superpowers/plans/2026-11-12-phase-1-clonotypes.md`, Task 17). Before starting, `git fetch origin && git rev-list --left-right --count main...origin/main` must print `0	0`, and phase 1's Done When commands must pass. Where phase 1's plan and the old code differ, phase 1's plan is the state you start from.

## Phase 1 Interfaces Consumed

Use these as phase 1 named them. Do not redefine them. The task number is phase 1's.
- `api/models.py`: `RunTrust(Model): clonal: ClonalAgreement | None` (T13); `ClonalAgreement.status`, `.primary: ClonalScore | None` with `.p` and `.delta_eta` (T13); `TrajectoryRun.trust: RunTrust | None`, `JobKind`, `JobInfo.kind: JobKind | str` (T14).
- `engine/store.py`: `STORE_VERSION = 2`, `graph["prep_key"]`, `graph["obs_keys"]`, `RunStore.save(rid, result, terminals=None, prep_key=None)` (T5); `TRUST_NAME`, `RunStore.trust_path`, `RunStore.write_trust(rid, name, payload) -> Path` (strict JSON, atomic, `ValueError` on NaN), `RunStore.read_trust(rid, name) -> dict | None`, and `save` keeps `trust/` (T11).
- `jobs/manager.py`: `BACKGROUND_KINDS = frozenset({"score_clonal", "stability"})`, `_take_next`, `add_done_callback`, `DEFAULT_TIME_LIMITS` (T12).
- `api/runs.py`, `api/session.py`: `TrustLookup(stored, job)`, `Runs._trust_lookup(rid, kind, name, h, payload) -> TrustLookup`, `Runs.submit_trust_jobs(rid) -> list[Job]`, `Runs.clonal_agreement(rid, graph) -> ClonalAgreement`, `Session._on_job_done` (T14).
- `engine/context.py`: `InferenceContext.from_prep(ds, prep, params, clonal_columns=())`, `ScoringOnlyColumnError`; `DatasetView.clonal_columns` (T5).
- `trust/constants.py` (`CLONAL_P_MAX = 0.01` is C5's bound) and `tests/unit/test_trust_constants.py::MIRRORED` with its equality test (T1).
- `bench/runner.py`: `ClonalRow`, `CaseResult.clonal`, and a `CaseResult.row()` that always emits its fields; `bench/gates.py` `check_gates` (`nan`, `status`), `GateFailure.reason` (T1).
- `clonal/scores.py`: `ClonalResult` with `.status` (`"ok" | "single_lineage" | "too_few_clones"`), `.p` (headline) and `.delta_eta` (T2). This is C5's input.
- `bench/synth.py` `generate(...)` output is pinned by phase 1 T3's digest test, and that test must keep passing.

## Global Constraints

Every task's requirements include these (spec §2).
- Clone labels are used only for scoring and never reach inference. Stability replicates get the same `clonal_columns` as inference and go through `run_method`, which refuses scoring-only keys.
- Every metric, threshold rule and control in this phase is written into `benchmarks/thresholds.toml` before the first run that measures it (Task 1). The τs come from the pre-declared rule, and they are never moved afterwards to fit a result.
- No output or text claims direction (which state came first) from clones. That includes C5's reasons.
- No method default and no threshold is tuned on T-cell data. The verdict τs come only from synthetic seeds 0–2 of the 12-case matrix and of the four controls.
- Thresholds are never lowered. A new gate or definition goes in with its reason and a date, in a task that lands before the first run that measures it.
- Roots are always explicit. Every stability replicate is rooted at `cell:<obs_name>` of the reference run's root cell. Inference never uses the UMAP. The input file is never modified.
- A method's `experimental` flag matches its synthetic gates (`tests/golden/test_synthetic.py`). The controls sit in `[synthetic.controls]`, outside the gated matrix, and never affect a flag.
- API changes start in `src/lumencell/api/models.py`, then `make contract`. In `web/`, only `make contract` output and the smallest update to the `dataset` fixture in `web/src/test/fixtures.ts`.
- The server never runs prep or methods. Stability runs in the worker. A GET that finds neither the trust file nor a job submits one.
- State-file writes are locked and atomic (`api/state.py`), and trust files are written only through `RunStore.write_trust`. Replicate prep is kept in memory and never enters the stage cache.
- `engine/` never imports `bench/`. No new runtime dependency. scikit-learn is never imported directly.
- Datasets and downloaded sources are never committed. A download asks the owner first and records its SHA-256. Real-data tests are marked `real_data` and slow tests `slow`. CI uses synthetic data only.
- Values read from a user's file (obs names, column names, lineage names) never reach a file path or a shell command. Trust files use hashed names.
- The runtime thresholds (cell-cycle ratios, stability parameters, C1's bound, the verdict τs) live in `src/lumencell/trust/constants.py`, and a test asserts they equal the TOML.
- Process: test-first; one branch and PR per task (`pivot/NN-<slug>`); squash-merge only on green CI after review by a different agent; nothing is pushed to `main` directly; no `Co-Authored-By` trailers or "Generated with" footers; each commit gated on the exit code of `make ci`.

## Review Focus

These are five inputs the spec implies but its own tests never exercise, most likely first. Each one has a test in the task named.
1. **Cell-cycle scoring that raises on real data** (duplicate gene symbols, list genes that are all zero, a scanpy error). A report-only stage must never fail prep. The stage completes with `not_scored`, and its reason names the error. Tests: Task 3, `test_a_scoring_error_is_reported_and_prep_still_completes` and `test_duplicate_symbols_are_scored_once`.
2. **Undefined stability values** (Spearman on fewer than 3 shared cells, a replicate that failed, no Hungarian pairs). Every such field is `null` in the API, and the trust file is strict JSON. `write_trust` refuses NaN, so a missed NaN would fail the job. Test: Task 10, `test_undefined_values_are_null_and_the_file_is_strict_json`.
3. **A replicate whose method raises** (for example a user-named terminal with no cells left in a subsample). The replicate is recorded with its error and counts as zero agreement. It is never dropped from a median, and the job still completes. Test: Task 6, `test_a_failing_replicate_counts_as_no_agreement`.
4. **Stability cancelled or failed.** The verdict never reads `supported`, and repeated GETs do not resubmit the job. Test: Task 11, `test_a_cancelled_stability_keeps_the_verdict_provisional_and_is_not_resubmitted`.
5. **A PAGA graph with no edges** (one Leiden cluster). C2 is undefined: `null` in the API, `not_supported` with a reason, and no division error. Test: Task 7, `test_a_graph_without_edges_fails_c2_with_a_reason`.

## Task Conventions

- **Branch:** each task header names `pivot/NN-<slug>`. NN is assigned at execution time.
- **Model:** each task names its implementer model and a reason.
- **TDD:** write the failing test, run it and see it fail for the stated reason, implement, then run it green.
- **Closing steps** ("C1–C4" below): C1, update `docs/superpowers/HANDOFF.md` with the current state, the decisions taken in this task, and the branch and PR row. C2, `make ci; echo $?` prints `0`; gate the commit on that exit code, never on piped output. C3, commit with the task's message, with no `Co-Authored-By` trailer and no "Generated with" footer. C4, `git push -u origin pivot/NN-<slug>`, then open the PR with `gh pr create`. Do not merge. The orchestrator squash-merges after review.
- **Background runs:** steps marked **(background)** take tens of minutes to hours. Start them with the Bash tool's `run_in_background: true` and keep the log under `$TMPDIR`. Task 10 does not depend on Tasks 8 and 9, so the orchestrator may run it while their long runs are going.
- **Order:** 1 → (2 → 3 → 4), 5, 6, 7 → 8 → 9 → 11, 12 → 13; 10 needs only 6. Tasks 11 and 12 use the verdict on data, so they wait for the τs (Task 8) and a passing verification (Task 9).
- **Heavy imports** inside CLI commands stay inside the function, as `cli.py` does today.

## Process Notes (not tasks)

- After each PR, run the pr-review-toolkit agents `code-reviewer`, `silent-failure-hunter` and `pr-test-analyzer`. Add `type-design-analyzer` on the PRs that add types: Tasks 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 and 12.
- At the end of the phase, run `/security-review` (Task 13).
- The calibration (Task 8) and the verification (Task 9) also go through the `data:statistical-analysis` and `data:validate-data` skills, and the real-data stability numbers (Task 12) through `data:validate-data`.
- DPT on a blob is likely stable, so C4 (pseudotime stability) may well prove non-discriminating. That is an acceptable outcome. It is recorded and reported, and C4 then drops out of the verdict. No rule, seed or control is changed to rescue it.

## Not in This Plan

- Spec §9: web views, the multi-user version, new trajectory methods, RNA velocity, receptor design, packaging, the 100k-cell performance budgets, cell-cycle regression, excluding immunoglobulin genes, AIRR or scirpy import in the app, similarity-based clonotype clusters, any direction inferred from clones, per-clone pseudotime coherence.
- Phase 1 work, which is consumed and not redone: clonal scores and models, `write_trust`, job priorities and done-callbacks, `_trust_lookup`, the constants test, and the bench row and gate machinery.
- Phase 3: Daniel 2022 and Sheih 2020. The case study's stability, verdict and cell-cycle report use what this phase builds.
- Stability in the synthetic `bench synth` matrix. Spec §5.2 asks for it only on the real datasets and in the calibration.
- C5 in the calibration and the verification. The controls have no clones and C5 is not calibrated, so both runs score C5 as not applicable.
- A web label for the `cell_cycle` stage (`PrepProgress.tsx` falls back to the stage name), a route to re-queue a cancelled stability job (a server restart retries it), CLI or API options for `cc_species`, per-cell cell-cycle arrays over the API, and stability in `lumencell run` (`launch` heals it on the first GET).

## File Map

| File | Responsibility | Task |
|---|---|---|
| `benchmarks/thresholds.toml` | cell-cycle, stability and verdict definitions, the calibration and verification rule, the controls; calibrated τs; changelog | 1, 8, 9, 12 |
| `src/lumencell/trust/constants.py` | runtime values mirrored from the TOML; τs | 1, 8 |
| `src/lumencell/prep/cell_cycle_data/tirosh2016.tsv` | the S and G2M lists, with source and verification header | 2 |
| `src/lumencell/prep/cell_cycle_data/human_mouse_orthologs.tsv` | the human-to-mouse table, with source, date and verification header | 2 |
| `src/lumencell/prep/cell_cycle_genes.py` | load the lists, resolve species, match symbols | 2 |
| `src/lumencell/prep/stages.py` | `Stage.report_only`, `PrepParams.cc_*`, the `cell_cycle` stage | 3 |
| `src/lumencell/prep/cache.py` | `prep_key` leaves report-only stages out | 3 |
| `src/lumencell/cli.py` | `prepare` cell-cycle line; `bench real --stability` | 3, 12 |
| `src/lumencell/trust/cell_cycle.py` | dataset cell-cycle info and the per-run lineage check | 4 |
| `src/lumencell/api/models.py` | cell-cycle models; stability models; verdict models, `RunTrust` and `DatasetInfo` fields, `JobKind` | 4, 10, 11 |
| `src/lumencell/bench/synth.py` | `_sample_counts` made public as `sample_counts` (rename only) | 5 |
| `src/lumencell/bench/synth_controls.py` | blob, two_groups, loop and sheet | 5 |
| `src/lumencell/engine/compare.py` | Jaccard matching with a shared-cell mask | 6 |
| `src/lumencell/engine/stability.py` | replicate plan, restricted prep, terminal mapping, padding, comparison, summaries | 6 |
| `src/lumencell/engine/tree_support.py` | C2 tree-likeness, the checks and the levels; `calibrated_thresholds` | 7, 8 |
| `src/lumencell/bench/verdict_calibration.py` | collect check values, derive τs, verify on fresh seeds | 8, 9 |
| `benchmarks/results/verdict-calibration-<date>.json` | recorded calibration | 8 |
| `benchmarks/results/verdict-verification-<date>.json` | recorded verification | 9 |
| `src/lumencell/trust/stability.py` | stability trust name and hash; API mapping | 10 |
| `src/lumencell/jobs/worker.py` | `stability` handler | 10 |
| `src/lumencell/jobs/manager.py` | `stability` time limit | 10 |
| `src/lumencell/engine/store.py` | `result.json` records the requested terminals | 10 |
| `src/lumencell/trust/verdict.py` | API trust blocks to verdict inputs, verdict to the API model | 11 |
| `src/lumencell/api/runs.py`, `api/session.py`, `api/dataset_view.py` | serve, submit and heal stability; verdict; cell-cycle check; `DatasetInfo` fields | 11 |
| `contracts/openapi.json`, `web/src/api/generated/*`, `web/src/test/fixtures.ts` | `make contract` output; `dataset` fixture | 11 |
| `src/lumencell/bench/runner.py`, `src/lumencell/bench/report.py` | `StabilityRow`, real-data stability and verdict, report column | 12 |
| `benchmarks/results/real-<date>.json`, `README.md` | real re-run with stability; regenerated tables | 12 |
| `docs/superpowers/HANDOFF.md` | state, decisions, PRs, verification records, security review | every task |
| tests: `tests/unit/test_gates.py`, `tests/unit/test_trust_constants.py`, `tests/unit/test_cell_cycle_genes.py`, `tests/unit/test_prep.py`, `tests/unit/test_prep_cache.py`, `tests/unit/test_commands.py`, `tests/api/test_cache_healing.py`, `tests/unit/test_cell_cycle_check.py`, `tests/unit/test_synth_controls.py`, `tests/unit/test_compare.py`, `tests/unit/test_stability.py`, `tests/unit/test_tree_support.py`, `tests/unit/test_verdict_calibration.py`, `tests/unit/test_jobs.py`, `tests/unit/test_store.py`, `tests/unit/test_stability_trust.py`, `tests/api/test_robustness_trust.py`, `tests/unit/test_runner_stability.py` | the tests named in each task | per task |

## Done When

Spec §5.5, as commands. Run them from the repo root on `main` after Task 13.
- `make ci; echo $?` → `0`.
- **(background)** `make bench-synth > "$TMPDIR/synth.log" 2>&1; echo $?` → `0`, and `grep -c '^REGRESSION' "$TMPDIR/synth.log"` → `0`.
- The cell-cycle stage is cached, tested and changes no run id: `uv run pytest tests/unit/test_prep.py tests/unit/test_prep_cache.py tests/unit/test_cell_cycle_genes.py tests/api/test_cache_healing.py -q; echo $?` → `0` (this includes the pinned `prep_key` test).
- Stability is queued at background priority, served, healed after a restart and cancellable: `uv run pytest tests/api/test_robustness_trust.py -q; echo $?` → `0`.
- The τs were derived and verified: `uv run python -c "import json,glob; c=sorted(glob.glob('benchmarks/results/verdict-calibration-*.json'))[-1]; v=sorted(glob.glob('benchmarks/results/verdict-verification-*.json'))[-1]; print(c, json.load(open(c))['complete'], v, json.load(open(v))['passed'])"` → both files, `True`, `True`. `grep -c '^\[tree_support.calibrated\]' benchmarks/thresholds.toml` → `1`, and `uv run pytest tests/unit/test_trust_constants.py tests/unit/test_verdict_calibration.py -q` passes.
- The definitions came before the measurements: `git merge-base --is-ancestor $(git log --format=%h -1 --grep='Declare the T-cell robustness') $(git log --format=%h --reverse -- 'benchmarks/results/verdict-*.json' | head -1); echo $?` → `0`.
- `docs/superpowers/HANDOFF.md` has a "Security review (phase 2)" section that records each finding as resolved or recorded.

---

### Task 1: Declare the robustness definitions and runtime constants

**Branch:** `pivot/NN-robustness-definitions` · **Model:** sonnet. Configuration plus constants, fully specified here.

**Files:**
- Modify: `benchmarks/thresholds.toml`, `src/lumencell/trust/constants.py`, `tests/unit/test_trust_constants.py`, `tests/unit/test_gates.py`

**Interfaces:**
- Consumes: phase 1 T1's `trust/constants.py` and `MIRRORED`.
- Produces, `trust/constants.py`: `CELL_CYCLE_MIN_GENES = 20`, `CELL_CYCLE_SEED = 0`, `CELL_CYCLE_RATIO_MIN = 2.0`, `CELL_CYCLE_SHARE_MIN = 0.5`, `STABILITY_METRIC_VERSION = 1`, `N_RESEEDS = 2`, `N_SUBSAMPLES = 3`, `SUBSAMPLE_FRACTION = 0.8`, `STABILITY_SEED = 0`, `MATCH_JACCARD = 0.6`, `COVERAGE_MIN = 0.95`. `MIRRORED` gains one entry each, e.g. `"CELL_CYCLE_RATIO_MIN": ("cell_cycle", "lineage_ratio_min")`, `"MATCH_JACCARD": ("stability", "match_jaccard")`, `"COVERAGE_MIN": ("tree_support", "coverage_min")`. Task 8 adds the τs.
- Produces, TOML tables `[cell_cycle]`, `[stability]`, `[tree_support]`, `[tree_support.calibration]`, `[tree_support.verification]` and `[synthetic.controls]`, read by Tasks 3–12. Task 8 appends `[tree_support.calibrated]`.

- [ ] **Step 1: Write the failing tests.** In `tests/unit/test_gates.py`:
  - `test_robustness_definitions_are_declared`: with `cfg = load_thresholds()`, assert `cfg["cell_cycle"] == {"gene_set": "tirosh2016", "min_genes": 20, "seed": 0, "lineage_ratio_min": 2.0, "lineage_share_min": 0.5}`; `cfg["stability"] == {"metric_version": 1, "n_reseeds": 2, "n_subsamples": 3, "subsample_fraction": 0.8, "seed": 0, "match_jaccard": 0.6}`; `cfg["tree_support"]["coverage_min"] == 0.95`; `cfg["tree_support"]["calibration"] == {"matrix_seeds": [0, 1, 2], "control_seeds": [0, 1, 2], "methods": ["paga_dpt", "palantir", "scfates_ppt"], "matrix_quantile": 0.10, "control_quantile": 0.90, "c2_resolution_factors": [0.5, 1.0, 2.0]}`; `cfg["tree_support"]["verification"] == {"matrix_seeds": [3, 4, 5], "control_seeds": [3, 4, 5], "max_controls_supported": 0}`; `"calibrated" not in cfg["tree_support"]`.
  - `test_controls_sit_outside_the_gated_matrix`: `cfg["synthetic"]["controls"]["kinds"] == ["blob", "two_groups", "loop", "sheet"]`; those kinds are disjoint from `cfg["synthetic"]["topologies"]`; `"controls" not in cfg["synthetic"]["gates"]`.
  - In `tests/unit/test_trust_constants.py`, add the eleven `MIRRORED` entries. The existing equality test now fails because the constants are missing.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_gates.py tests/unit/test_trust_constants.py -q`. Expected: FAIL (missing tables and constants).
- [ ] **Step 3: Add to `benchmarks/thresholds.toml`.** Put this right after `[synthetic.gates]`:

```toml
[synthetic.controls]
# Negative controls for the tree verdict (spec §5.3), built directly by bench/synth_controls.py,
# never by the tree generator. Outside the gated matrix: no gate, baseline, regression check or
# experimental flag reads this table. Counts use the [synthetic] preset's noise model.
kinds = ["blob", "two_groups", "loop", "sheet"]
n_cells = 3000
n_genes = 2000
programme_share = 0.10       # share of genes driven by each latent coordinate (two_groups, sheet)
loop_programme_share = 0.20  # loop: genes following the cycle
lfc_scale = 1.0              # median |log fold change| of a programme gene (as the preset's de_scale)
group_lfc_scale = 2.0        # two_groups: far enough apart to be unconnected in the kNN graph
```

  Put this block before `[regression]`. Its comments are the definitions the spec requires before any result:

```toml
# ---- T-cell robustness (spec docs/superpowers/specs/2026-11-12-tcell-clonotypes-design.md §5)
# Defined before any cell-cycle score, stability replicate or verdict value was computed.
# src/lumencell/trust/constants.py mirrors what the worker and server need at run time.

[cell_cycle]
# Report-only prep stage, left out of prep_key: scanpy score_genes_cell_cycle on the
# log-normalised matrix with the Tirosh et al. 2016 S and G2M lists (mouse through the shipped
# human-to-mouse ortholog table). Species "auto" goes by the case of the matching symbols.
gene_set = "tirosh2016"
min_genes = 20            # fewer list genes found (S and G2M together), or none in one list: not_scored
seed = 0                  # score_genes random_state
# Per-run check, a hint and never a verdict (proliferating effector cells are real biology): a
# lineage is flagged when the share of its single-lineage cells in S or G2M is at least
# lineage_ratio_min times the share of all scored cells, and at least lineage_share_min.
lineage_ratio_min = 2.0
lineage_share_min = 0.5

[stability]
metric_version = 1
n_reseeds = 2             # Leiden, PAGA and the method rerun on all cells with new seeds
n_subsamples = 3          # subsample_fraction of cells, uniform, the root cell always kept; kNN,
subsample_fraction = 0.8  #   Leiden, diffusion maps and PAGA recomputed; normalisation, HVG, PCA reused
seed = 0
# Reseed i (from 0) uses prep seed + 1 + i for Leiden and PAGA and, when the method has a seed
# parameter, the method's seed + 1 + i. Subsample j draws its cells with numpy
# default_rng(seed + 1 + j) and keeps the reference's prep and method seeds. Every replicate is
# rooted at cell:<obs name> of the reference run's root cell. Terminals named in the Leiden column
# map to the replicate's Leiden cluster with the largest overlap; other terminals pass unchanged.
# Per replicate, against the reference, on the cells both contain (subsamples padded to all cells):
#   spearman of pseudotime over cells finite in both; lineages paired by the Hungarian algorithm
#   on Jaccard over shared cells; share_matched = pairs with Jaccard >= match_jaccard / reference
#   lineages; mean_jaccard = mean Jaccard of the pairs; n_lineages of the replicate.
# A replicate that fails, or whose value is undefined, counts as 0 in every median and minimum.
match_jaccard = 0.6

[tree_support]
# C1 connectivity: share of cells with finite pseudotime >= coverage_min.
# C2 tree-likeness (dataset level): weight of PAGA's maximum spanning tree / total PAGA edge
#   weight, each edge once (upper triangle) >= tau_tree. No edges: undefined, which fails.
# C3 lineage stability: median over all R replicates of share_matched >= tau_lin.
# C4 pseudotime stability: median over all R replicates of spearman >= tau_pt.
# C5 clonal support (clones configured, L >= 2): status ok, headline p <= clonal.p_max, delta_eta > 0.
# Level: not_supported if C1 or C2 fails; provisional while C1 or C2 is unknown; weak if an
#   applicable check among C3-C5 fails (a failed clonal score fails C5); provisional while one of
#   them is pending or stability could not be measured; supported otherwise. A non-discriminating
#   check is left out and the verdict says so.
coverage_min = 0.95

[tree_support.calibration]
# tau = (quantile(matrix values, matrix_quantile) + quantile(control values, control_quantile)) / 2,
# numpy.quantile, method "linear". If the matrix quantile is not above the control quantile, the
# check is non-discriminating and left out of the verdict. C3 and C4 pool every (dataset, method)
# run; C2 takes one value per dataset (methods do not change it), and an undefined C2 counts as
# 0. C2 is also measured at Leiden resolution x0.5 and x2 and is kept only if it discriminates at
# all three; tau_tree is the one at the default resolution. A run whose reference fails is left
# out and counted. C1's bound and C5's rule are fixed, not calibrated. Matrix: the [synthetic]
# topologies, preset and size. Controls: every [synthetic.controls] kind.
matrix_seeds = [0, 1, 2]
control_seeds = [0, 1, 2]
methods = ["paga_dpt", "palantir", "scfates_ppt"]
matrix_quantile = 0.10
control_quantile = 0.90
c2_resolution_factors = [0.5, 1.0, 2.0]

[tree_support.verification]
# Fresh seeds, the calibrated taus, C5 not applicable (no clonal settings). No control run may be
# verdicted supported. The share of matrix runs verdicted supported is reported, not gated.
matrix_seeds = [3, 4, 5]
control_seeds = [3, 4, 5]
max_controls_supported = 0
```

  Add a changelog entry dated with this commit's date: "T-cell robustness definitions: the cell_cycle stage and its per-run hint, the stability replicates and comparison, the tree-support checks and levels, the τ calibration and verification rule, and four negative controls outside the gated matrix. All are set before any of them was measured. Nothing existing changed."
- [ ] **Step 4: Implement** the eleven constants in `trust/constants.py`.
- [ ] **Step 5: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 6: C1–C4.** Commit message: `Declare the T-cell robustness definitions, verdict calibration rule and controls before any measurement`.

---

### Task 2: Ship and verify the Tirosh 2016 lists and the mouse ortholog table

**Branch:** `pivot/NN-cell-cycle-genes` · **Model:** sonnet. Careful data work with a fixed procedure; no architecture.

**Files:**
- Create: `src/lumencell/prep/cell_cycle_data/tirosh2016.tsv`, `src/lumencell/prep/cell_cycle_data/human_mouse_orthologs.tsv`, `src/lumencell/prep/cell_cycle_genes.py`, `tests/unit/test_cell_cycle_genes.py`

**Interfaces:**
- Produces, `prep/cell_cycle_genes.py`:
  - `GENE_SETS = ("tirosh2016",)`; `Species = Literal["human", "mouse"]`.
  - `@dataclass(frozen=True) class CellCycleGenes: gene_set: str; species: Species; s: tuple[tuple[str, ...], ...]; g2m: tuple[tuple[str, ...], ...]; source: str; unmapped: tuple[str, ...]`. Each gene is a tuple of the symbols it may appear under. For human that is the published symbol, plus the current HGNC symbol if the gene was renamed (e.g. `("MLF1IP", "CENPU")`). For mouse it is the MGI symbol. `unmapped` lists the human genes without a 1:1 mouse ortholog, and it is empty for human.
  - `load_gene_set(gene_set: str, species: Species) -> CellCycleGenes` reads the TSVs with `importlib.resources.files("lumencell.prep") / "cell_cycle_data"`. An unknown `gene_set` raises `ValueError`.
  - `resolve_species(symbols: Sequence[str], requested: Literal["auto", "human", "mouse"], gene_set: str = "tirosh2016") -> Species | None`. An explicit species is returned as it is. For `auto`, take the symbols that match the human or the mouse list case-insensitively. Count those written all upper-case against those with any lower-case letter. The larger count gives `human` or `mouse`, and a tie, including no match at all, gives `None`.
  - `match_genes(symbols: Sequence[str], genes: tuple[tuple[str, ...], ...]) -> list[int]`. For each gene, the index of the first symbol in var order that matches any of its accepted symbols, case-insensitively. Genes without a match are left out, and each var index appears at most once.
- File format: `#` header lines, then a tab-separated table. `tirosh2016.tsv` has columns `phase` (`S` or `G2M`), `symbol` (as published) and `current_symbol` (HGNC today, empty if unchanged). `human_mouse_orthologs.tsv` has `human`, `mouse` and `mgi_class_key`. Each header states the source citation or URL, the retrieval date, the SHA-256 of every downloaded file, and the verification done, with its result.

- [ ] **Step 1: Stop and ask the owner** before any download. Give each file's URL and its size, taken from a `HEAD` request:
  - the Tirosh et al. 2016 supplementary tables (Science, doi:10.1126/science.aad0501) holding the G1/S and G2/M gene lists, from the paper's page or its PubMed Central copy;
  - the scanpy tutorial copy of the same list, `https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt`, as an independent cross-check;
  - MGI's `HOM_MouseHumanSequence.rpt` (`https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt`);
  - HGNC HCOP `human_mouse_hcop_fifteen_column.txt.gz` (`https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/`), for the cross-check and for current symbols.

  If a source needs a login, ask the owner to download it by hand. Wait for a clear yes.
- [ ] **Step 2: Write the failing tests** (`tests/unit/test_cell_cycle_genes.py`):
  - `test_tirosh_lists_have_43_s_and_54_g2m_genes`: human `len(s) == 43`, `len(g2m) == 54`, with no gene in both lists and no repeated symbol. These are the counts distributed by Seurat's `cc.genes` and the scanpy list. If the primary source gives other counts, stop and report to the owner; do not edit the numbers.
  - `test_known_members`: human S contains `("MCM5",)`, `("PCNA",)` and `("MLF1IP", "CENPU")`, and G2M contains `("TOP2A",)` and `("HN1", "JPT1")`. Mouse S contains `("Mcm5",)` and G2M contains `("Top2a",)`.
  - `test_every_human_gene_is_mapped_or_listed_as_unmapped`: `len(mouse.s) + len(mouse.g2m) + len(mouse.unmapped) == 97`.
  - `test_headers_name_their_sources`: the Tirosh header matches `10.1126/science.aad0501`, and the ortholog header matches `HOM_MouseHumanSequence.rpt`. Both contain an ISO date and at least one 64-hex SHA-256, plus a line starting `# verified:`.
  - `test_resolve_species_from_symbol_case`: upper-case Tirosh symbols give `human`; title-case symbols give `mouse`; `["G00000", "G00001"]` gives `None`; an equal mix gives `None`; `requested="mouse"` gives `mouse` for any input.
  - `test_matching_is_case_insensitive_and_counts_genes_once`: on `["pcna", "PCNA", "CENPU", "x"]`, matching the human S list returns `[0, 2]`, meaning one index per gene and the first occurrence wins.
  - `test_unknown_gene_set_raises`.
- [ ] **Step 3: Run** `uv run pytest tests/unit/test_cell_cycle_genes.py -q`. Expected: FAIL (module missing).
- [ ] **Step 4: Fetch and verify.** Download as the existing fetch pattern in `bench/datasets.py` does: `certifi`'s CA bundle, with certificate checks always on. Put the downloads in a new, empty directory `$TMPDIR/cc-sources/` and never commit them. Record `shasum -a 256` of each file. Keep the extraction script outside that directory and run it with `uv run python -I script.py <file>`. Verify both sources and record how:
  - **Lists:** compare the paper's table with the scanpy copy, symbol for symbol and phase for phase. Write the result into the `# verified:` line, either "identical" or each difference with how it was resolved; the paper wins. Fill `current_symbol` from HCOP's current HGNC symbol for each published symbol that is no longer approved.
  - **Orthologs:** for each human gene, take the MGI homology class (DB Class Key) holding its current symbol. The mapping is kept only if that class has exactly one human gene and exactly one mouse gene; anything else is unmapped and listed in the header. Cross-check every kept mapping against HCOP's mouse symbol. Record the agreement count and each disagreement, resolved in favour of MGI's 1:1 class, in `# verified:`.
- [ ] **Step 5: Implement** `cell_cycle_genes.py` as in Interfaces.
- [ ] **Step 6: Run** the Step 3 command. Expected: PASS.
- [ ] **Step 7: C1–C4.** In HANDOFF, record each source with its SHA-256, the date and the two verification results. Commit message: `Ship the Tirosh 2016 cell-cycle lists and a verified human-to-mouse ortholog table`.

---

### Task 3: The report-only `cell_cycle` stage

**Branch:** `pivot/NN-cell-cycle-stage` · **Model:** sonnet. One new stage and one key rule in the existing pattern.

**Files:**
- Modify: `src/lumencell/prep/stages.py`, `src/lumencell/prep/cache.py`, `src/lumencell/cli.py` (`prepare`)
- Test: `tests/unit/test_prep.py`, `tests/unit/test_prep_cache.py`, `tests/unit/test_commands.py`, `tests/api/test_cache_healing.py`

**Interfaces:**
- Consumes: Task 2 (`load_gene_set`, `resolve_species`, `match_genes`); `CELL_CYCLE_MIN_GENES`, `CELL_CYCLE_SEED`.
- Produces, `prep/stages.py`:
  - `Stage.report_only: bool = False`, as the new last field.
  - `PrepParams.cc_species: Literal["auto", "human", "mouse"] = "auto"` and `PrepParams.cc_gene_set: Literal["tirosh2016"] = "tirosh2016"`. A `__post_init__` raises `ValueError` for any other value.
  - `STAGES` gains, last: `Stage("cell_cycle", ("lognorm",), ("cc_s", "cc_g2m", "cc_phase", "cc_info"), ("cc_species", "cc_gene_set"), _cell_cycle, report_only=True)`.
  - The outputs are always all written. `cc_s` and `cc_g2m` are float32 `(n,)`, all NaN when not scored. `cc_phase` is a `<U8` array of `"G1"`, `"S"` and `"G2M"`, all `"unscored"` when not scored. `cc_info = {"status": "scored" | "not_scored", "reason": str | None, "requested_species": str, "species": "human" | "mouse" | None, "gene_set": str, "genes_found": {"s": int, "g2m": int}, "genes_listed": {"s": int, "g2m": int}, "source": str, "seed": CELL_CYCLE_SEED}`.
  - The checks run in this order, and the first that applies gives the `not_scored` reason.
    1. Count the list genes that the symbols match case-insensitively: by their human or mouse symbols for `auto`, and by the requested species' symbols otherwise. Fewer than `CELL_CYCLE_MIN_GENES` gives "only {k} of {total} Tirosh 2016 genes were found; at least 20 are needed". Synthetic `G00000` symbols always stop here.
    2. Only then resolve the species. If it is unresolved, the reason is "the species could not be told from the gene symbols".
    3. Match the resolved species' lists. Fewer than 20 across both lists, or none in one list, gives the reason from step 1 with the new count.
    4. If scoring raises any `Exception`, the reason is `f"scoring failed: {type(exc).__name__}: {exc}"`, with one stderr line as well.
  - Scoring runs `sc.tl.score_genes_cell_cycle(view, s_genes=…, g2m_genes=…, random_state=CELL_CYCLE_SEED)`. The view is an AnnData with `X = art["lognorm"]` and var names `anndata.utils.make_index_unique(pd.Index(ds.gene_symbols))`. The matched genes are passed by those unique names.
- Produces, `prep/cache.py`: `prep_key(ds, params)` hashes only the keys of stages with `report_only=False`. `stage_keys`, `status`, `run` and `load_all` still cover every stage, so the report-only stage counts towards readiness.
- Produces, `cli.prepare`: one line after the stage table. It reads `cell cycle: scored as {species} from {k} of {total} Tirosh 2016 genes` or `cell cycle: not scored ({reason})`.

- [ ] **Step 0: Pin the current key.** On the branch's starting commit, before any change, print `prep_key(cast(Dataset, SimpleNamespace(fingerprint="f" * 24)), PrepParams())`. Paste it into Step 1's test as a literal. It is the key every existing run id was built from.
- [ ] **Step 1: Write the failing tests.**
  - `tests/unit/test_prep_cache.py`:
    - `test_prep_key_is_what_it_was_before_the_cell_cycle_stage`: the key above equals the Step 0 literal.
    - `test_report_only_stages_are_left_out_of_prep_key`: `prep_key` is equal for `PrepParams()` and `PrepParams(cc_species="mouse")`, while `stage_keys(...)["cell_cycle"]` differs.
    - `test_changing_cell_cycle_params_recomputes_only_that_stage`: after `run()`, `run(PrepParams(cc_species="human"))` has every stage `cached` except `cell_cycle`, which is `done`.
    - `test_report_only_stages_count_towards_readiness`: invalidate the `cell_cycle` stage, then `status(...)["cell_cycle"] == "pending"` and `load_all` raises `CacheMiss` with `.stage == "cell_cycle"`.
  - `tests/unit/test_prep.py`:
    - `test_no_inference_stage_reads_a_report_only_output`: the outputs of report-only stages are disjoint from the inputs of every other stage.
    - `test_cell_cycle_is_not_scored_on_synthetic_genes`: `cc_info["status"] == "not_scored"`, `cc_info["species"] is None`, and the reason contains `"0 of 97"`, not the species reason.
    - `test_an_unresolved_species_is_reported`: 30 list symbols, 15 upper-case and 15 title-case, give the species reason. `cc_s` and `cc_g2m` are float32 of length n and all NaN, and `cc_phase` is all `"unscored"`.
    - `test_cell_cycle_scores_human_symbols`: start from a 600-cell, 400-gene bifurcation (seed 0). Rename var 0–96 to the human S symbols then the G2M symbols, and multiply the S genes' counts by 8 in cells 0–149 and the G2M genes' counts by 8 in cells 150–299. Then `species == "human"` and `genes_found == {"s": 43, "g2m": 54}`. At least 80% of cells 0–149 are `"S"` and at least 80% of cells 150–299 are `"G2M"`. The share of `"G1"` in cells 300–599 is above that in cells 0–299.
    - `test_cell_cycle_resolves_mouse_from_symbol_case`: the same with the mouse symbols gives `species == "mouse"`.
    - `test_duplicate_symbols_are_scored_once` (Review Focus 1): when `MCM5` appears twice, the stage scores, and `genes_found["s"] == 43`.
    - `test_a_scoring_error_is_reported_and_prep_still_completes` (Review Focus 1): monkeypatch `scanpy.tl.score_genes_cell_cycle` to raise `RuntimeError("boom")`. `run_pipeline` completes, `status == "not_scored"`, and the reason contains `"RuntimeError: boom"`.
    - `test_unknown_cell_cycle_params_are_refused`: `PrepParams(cc_species="dog")` raises `ValueError`.
  - `tests/unit/test_commands.py`, `test_prepare_reports_the_cell_cycle_status`: on the synthetic file, the output contains `cell cycle: not scored (only 0 of 97`.
  - `tests/api/test_cache_healing.py`, `test_a_cache_from_before_the_cell_cycle_stage_is_completed_by_one_prepare`: remove the `cell_cycle-*` stage folder. A new `Session` then reports `prep_ready() is False`. Start the manager and run `start_prep()`, then wait. Afterwards `prep_ready()` is true, the manifest shows every other stage with its earlier key, and `session.prep_key` is unchanged.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_prep.py tests/unit/test_prep_cache.py tests/unit/test_commands.py tests/api/test_cache_healing.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** as in Interfaces.
- [ ] **Step 4: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add a report-only cell_cycle prep stage that is left out of prep_key, so no run id changes`.

---

### Task 4: The per-run cell-cycle check

**Branch:** `pivot/NN-cell-cycle-check` · **Model:** sonnet. A small pure function and its models.

**Files:**
- Modify: `src/lumencell/api/models.py`
- Create: `src/lumencell/trust/cell_cycle.py`, `tests/unit/test_cell_cycle_check.py`

**Interfaces:**
- Consumes: Task 3's `cc_phase` and `cc_info`; `CELL_CYCLE_RATIO_MIN`, `CELL_CYCLE_SHARE_MIN`.
- Produces, `api/models.py`. These are defined here, and Task 11 references them from a route, so `make contract` shows no change yet.
  - `CellCycleInfo(Model): status: Literal["scored", "not_scored"]; reason: str | None; requested_species: Literal["auto", "human", "mouse"]; species: Literal["human", "mouse"] | None; gene_set: str; genes_found: int; genes_listed: int; source: str; share_cycling: float | None`
  - `LineageCycling(Model): lineage: int; name: str; n_cells: int; share_cycling: float | None; flagged: bool`
  - `CellCycleCheck(Model): status: Literal["scored", "not_scored"]; dataset_share: float | None; ratio_min: float; share_min: float; lineages: list[LineageCycling]; warnings: list[Notice]; note: str | None`
- Produces, `trust/cell_cycle.py`:
  - `CYCLING = frozenset({"S", "G2M"})`.
  - `cell_cycle_info(cc_info: dict[str, Any], cc_phase: npt.NDArray[np.str_]) -> models.CellCycleInfo`. `genes_found` and `genes_listed` are summed over both lists. `share_cycling` is the share of all cells in `CYCLING`, and `None` when not scored.
  - `cell_cycle_check(primary_lineage: npt.NDArray[np.int16], names: Sequence[str], cc_phase: npt.NDArray[np.str_], cc_info: dict[str, Any]) -> models.CellCycleCheck`. `dataset_share` is the share of all cells in `CYCLING`. A lineage's share is taken over its single-lineage cells (`primary_lineage == l`), and it is `None` when it has none. `flagged` means `share >= CELL_CYCLE_RATIO_MIN * dataset_share and share >= CELL_CYCLE_SHARE_MIN`. Each flagged lineage adds `Notice("lineage_mostly_cycling", "warning", f"{label} is mostly cycling cells; it may be a cell-cycle state rather than a fate.")`, where `label` is the name with "Lineage " in front unless the name already starts with "Lineage ". When not scored, `lineages` and `warnings` are empty and `note = cc_info["reason"]`.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_cell_cycle_check.py`):
  - `test_a_mostly_cycling_lineage_is_flagged`: 10 cells with `primary_lineage = [0,0,0,0,1,1,1,1,-1,-2]`, and S or G2M at cells 0, 1, 2. That gives `dataset_share == 0.3` and lineage 0 a share of 0.75, so lineage 0 is flagged and lineage 1 (share 0.0) is not. There is exactly one warning, reading exactly `"Lineage Alpha is mostly cycling cells; it may be a cell-cycle state rather than a fate."` for names `["Alpha", "Beta"]`.
  - `test_the_bounds_are_inclusive`: a share of exactly 0.5 with `dataset_share` exactly 0.25 is flagged.
  - `test_a_high_share_below_twice_the_dataset_is_not_flagged`: `dataset_share` 0.4 and lineage share 0.6 is not flagged.
  - `test_shared_and_unassigned_cells_count_only_for_the_dataset`: cells at −1 and −2 change `dataset_share` but no lineage's `n_cells`.
  - `test_a_lineage_without_single_lineage_cells_has_no_share`: `share_cycling is None`, `flagged is False`, and no exception.
  - `test_names_starting_with_lineage_are_not_doubled`: name `"Lineage 2"` gives a message starting `"Lineage 2 is mostly"`.
  - `test_not_scored_passes_the_reason_on`: `status == "not_scored"`, `note` equals the reason, and there are no lineages or warnings.
  - `test_dataset_info_sums_both_lists`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_cell_cycle_check.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the models first, then `make contract` (expected: no diff), then `trust/cell_cycle.py`.
- [ ] **Step 4: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add the per-run cell-cycle check that flags lineages made mostly of cycling cells`.

---

### Task 5: Negative controls

**Branch:** `pivot/NN-synth-controls` · **Model:** sonnet. The generative models are fixed here.

**Files:**
- Create: `src/lumencell/bench/synth_controls.py`, `tests/unit/test_synth_controls.py`
- Modify: `src/lumencell/bench/synth.py` (rename `_sample_counts` to `sample_counts` and nothing else)

**Interfaces:**
- Consumes: `[synthetic.controls]` and `[synthetic].preset` (Task 1); `PRESETS`, `SynthParams` and `sample_counts` from `bench/synth.py`.
- Produces, `bench/synth_controls.py` (read by Tasks 8 and 9):
  - `ControlKind = Literal["blob", "two_groups", "loop", "sheet"]`; `KINDS: tuple[ControlKind, ...]`; module constants equal to `[synthetic.controls]` (`PROGRAMME_SHARE = 0.10`, `LOOP_PROGRAMME_SHARE = 0.20`, `LFC_SCALE = 1.0`, `GROUP_LFC_SCALE = 2.0`); `CONTROL_STREAM = 7100502`.
  - `generate_control(kind: ControlKind, n_cells: int = 3000, n_genes: int = 2000, seed: int = 0, preset: str = "realistic") -> ad.AnnData`. It returns raw counts with var names `G00000…` and obs names `cell000000…`, and no embedding and no `uns["lumencell_truth"]`. `uns["lumencell_control"] = {"kind", "seed", "preset", "root_cell", "programme_genes"}`.
  - `control_root_cell(adata: ad.AnnData) -> str`, the obs name the run is rooted at with `cell:`.
- The model. Use one seeded generator per (kind, seed). The base log-means (`log_base`) and programme-gene fold changes follow the tree generator's conventions, with the scales and shares from the TOML. Counts come from `sample_counts(log_mu, PRESETS[preset], rng)`. What each control is:
  - `blob`: `log_mu = log_base` for every cell, so there are no programme genes. The root is `cell000000`.
  - `two_groups`: `group ~ Bernoulli(0.5)`; `log_mu = log_base + group · lfc` over a `PROGRAMME_SHARE` gene set with scale `GROUP_LFC_SCALE`. obs gains `control_group`. The root is the first cell of group 0.
  - `loop`: `θ ~ U(0, 2π)`; over a `LOOP_PROGRAMME_SHARE` gene set with magnitude `m_g` (scale `LFC_SCALE`) and phase `φ_g ~ U(0, 2π)`, `log_mu = log_base + m_g·cos(θ − φ_g)`. obs gains `control_angle`. The root is the cell with the smallest `min(θ, 2π − θ)`.
  - `sheet`: `x, y ~ U(0, 1)` independently; two disjoint `PROGRAMME_SHARE` gene sets with `log_mu = log_base + (2x − 1)·lfc_x + (2y − 1)·lfc_y`. obs gains `control_x` and `control_y`. The root is the cell with the smallest `x² + y²`.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_synth_controls.py`; 1000 cells × 500 genes unless stated; prep with `run_pipeline` where needed):
  - `test_controls_are_deterministic_and_differ_by_seed` (parametrized over the kinds).
  - `test_controls_are_counts_with_synthetic_gene_names`: `X` is integer-valued; the var names start `G00000`; `"lumencell_truth" not in uns`; `generate("blob")` from `bench/synth.py` raises `ValueError`.
  - `test_blob_has_no_programme`: `uns["lumencell_control"]["programme_genes"] == 0`, and no `control_*` obs column exists.
  - `test_two_groups_are_unconnected`: less than 0.1% of kNN connectivity edges join the two groups.
  - `test_the_loop_closes`: the median circular angle gap between kNN neighbours is under 0.3 rad, and at least one kNN edge joins a cell with θ < 0.2 to one with θ > 2π − 0.2.
  - `test_the_sheet_resolves_both_gradients`: the median |Δx| and the median |Δy| between kNN neighbours are both under 0.15.
  - `test_root_cells_sit_where_declared` (parametrized).
  - `test_constants_equal_the_thresholds_file`.
  - `test_the_generator_is_unchanged_by_the_rename`: phase 1 T3's `test_generator_output_is_unchanged_byte_for_byte` still passes. Run it by node id in Step 4.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_synth_controls.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the rename and `synth_controls.py`.
- [ ] **Step 4: Run** `uv run pytest tests/unit/test_synth_controls.py tests/unit/test_synth.py tests/unit/test_synth_clones.py -q`. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add four negative controls for the tree verdict, outside the gated matrix`.

---

### Task 6: Stability replicates and their comparison

**Branch:** `pivot/NN-stability-engine` · **Model:** opus. It restricts datasets, recomputes prep and compares on shared cells; the seed and root invariants are easy to break.

**Files:**
- Modify: `src/lumencell/engine/compare.py`
- Create: `src/lumencell/engine/stability.py`, `tests/unit/test_compare.py`, `tests/unit/test_stability.py`

**Interfaces:**
- Consumes: `run_stage`, `STAGES`, `PrepParams`, `PrepResult`; `InferenceContext.from_prep(..., clonal_columns=…)`; `run_method`; `CellRoot`; `N_RESEEDS`, `N_SUBSAMPLES`, `SUBSAMPLE_FRACTION`, `STABILITY_SEED`, `MATCH_JACCARD`.
- Produces, `engine/compare.py`. With `shared=None` everything behaves exactly as today.
  - `jaccard_matrix(ma: BoolArray, mb: BoolArray, shared: BoolArray | None = None) -> npt.NDArray[np.float64]`, with Jaccard computed over `shared` columns only.
  - `match_lineages(ma: BoolArray, mb: BoolArray, shared: BoolArray | None = None) -> list[LineagePair]` (Hungarian on that matrix).
  - `compare(a: TrajectoryResult, b: TrajectoryResult, shared: BoolArray | None = None) -> RunComparison`. With a mask, `both = finite(a) & finite(b) & shared`, and the lineages are paired with `match_lineages(..., shared)`.
- Produces, `engine/stability.py`. `run_method` is imported into the module namespace, and tests patch `lumencell.engine.stability.run_method`.
  - `class StabilityError(ValueError)`; `ReplicateKind = Literal["reseed", "subsample"]`.
  - `@dataclass(frozen=True) class StabilityConfig: n_reseeds: int = N_RESEEDS; n_subsamples: int = N_SUBSAMPLES; fraction: float = SUBSAMPLE_FRACTION; seed: int = STABILITY_SEED; match_jaccard: float = MATCH_JACCARD`
  - `@dataclass(frozen=True) class RunSpec: method: str; params: dict[str, Any]; terminals: list[str] | None; label_key: str | None; clonal_columns: tuple[str, ...] = ()`
  - `@dataclass(frozen=True) class ReplicatePlan: kind: ReplicateKind; index: int; seed: int; method_seed: int | None`; `replicate_plan(config: StabilityConfig, prep_seed: int, method_params: Mapping[str, Any]) -> tuple[ReplicatePlan, ...]`. The reseeds come first, then the subsamples, following the `[stability]` seed scheme. `method_seed` is `None` when the params have no `"seed"`. For a subsample it is the reference's own method seed.
  - `subsample_cells(n_obs: int, root: int, fraction: float, seed: int) -> npt.NDArray[np.int64]`. It returns `round(fraction·n_obs)` sorted, unique indices that always include `root`, drawn with `default_rng(seed)` from the other cells.
  - `restrict(ds: Dataset, cells: npt.NDArray[np.int64]) -> Dataset`: a `Dataset` of the subset, with no matrix and a fingerprint distinct from the full dataset's.
  - `reseed_prep(ds: Dataset, prep: PrepResult, params: PrepParams, seed: int) -> tuple[PrepResult, PrepParams]`. It reruns the `leiden` and `paga` stages with `replace(params, seed=seed)` and reuses every other artifact.
  - `subsample_prep(ds: Dataset, prep: PrepResult, params: PrepParams, cells: npt.NDArray[np.int64]) -> tuple[Dataset, PrepResult]`. It reuses `hvg`, `pca[cells]` and `pca_variance_ratio`, and reruns `neighbors`, `leiden`, `diffmap` and `paga` through `run_stage` on `restrict(ds, cells)`. Nothing is written to disk.
  - `map_leiden_terminals(terminals: Sequence[str], reference_leiden: IntArray, replicate_leiden: IntArray, cells: IntArray) -> list[str]`. Each terminal maps to the replicate label holding the most of that terminal's reference cells among `cells`; a tie goes to the lowest label. A terminal with no cells left is dropped, and duplicates collapse in order.
  - `pad_result(result: TrajectoryResult, cells: IntArray, n_obs: int) -> TrajectoryResult`: a full-length result in which cells outside `cells` are unassigned. The root index is mapped back, and the padded result passes `validate()`.
  - `@dataclass(frozen=True) class ReplicateScores: n_compared: int; n_lineages: int; spearman: float; share_matched: float; mean_jaccard: float`; `score_replicate(reference: TrajectoryResult, replicate: TrajectoryResult, shared: BoolArray, match_jaccard: float) -> ReplicateScores`. A value that cannot be computed is NaN.
  - `@dataclass(frozen=True) class ReplicateRow: kind: ReplicateKind; index: int; seed: int; method_seed: int | None; n_cells: int; n_compared: int; n_lineages: int; spearman: float; share_matched: float; mean_jaccard: float; terminals: tuple[str, ...] | None; error: str | None = None`
  - `@dataclass(frozen=True) class KindSummary: kind: ReplicateKind; r: int; fraction: float; seeds: tuple[int, ...]; spearman_median: float; spearman_min: float; share_matched_median: float; share_matched_min: float; mean_jaccard_median: float; mean_jaccard_min: float; n_lineages_median: float; n_lineages_min: int; rows: tuple[ReplicateRow, ...]`; `summarise(kind: ReplicateKind, fraction: float, rows: Sequence[ReplicateRow]) -> KindSummary`. NaN counts as 0 in every median and minimum, and a failed row counts as `n_lineages = 0`.
  - `@dataclass(frozen=True) class StabilityResult: n_obs: int; reference_lineages: int; match_jaccard: float; reseeds: KindSummary; subsamples: KindSummary; lineage_stability: float; pseudotime_stability: float`. The last two are the medians over all R rows, NaN counted as 0; they are the C3 and C4 values.
  - `run_stability(ds: Dataset, prep: PrepResult, params: PrepParams, reference: TrajectoryResult, spec: RunSpec, config: StabilityConfig = StabilityConfig(), progress: ProgressFn = no_progress) -> StabilityResult`. It first raises `StabilityError("the root cell's name {name!r} occurs more than once in the file, so replicates cannot be rooted at it")` when the name is not unique. For each `plan` in `replicate_plan(config, params.seed, spec.params)`, a reseed uses `reseed_prep(ds, prep, params, plan.seed)`, and a subsample uses `subsample_prep(ds, prep, params, subsample_cells(n, root, config.fraction, plan.seed))`. The method params are `{**spec.params, "seed": plan.method_seed}` when `plan.method_seed` is not `None`, and `spec.params` otherwise. The replicate then builds `InferenceContext.from_prep(..., clonal_columns=spec.clonal_columns)` and runs `run_method(ctx, spec.method, CellRoot(reference.root.obs_name), terminals, method_params, label_key=spec.label_key)`. When `params.get("cluster_key", LEIDEN_KEY) == LEIDEN_KEY`, the terminals are mapped with `map_leiden_terminals`; otherwise they pass unchanged. Subsample results are padded, and every replicate is scored against `reference` with `shared` as the replicate's cells (all of them for reseeds). An `Exception` in a replicate becomes `error=f"{type(exc).__name__}: {exc}"`, with its traceback on stderr; NaN values; and `n_lineages=0`. `progress(f"{kind} {index + 1}", i, R)` is called before each replicate.

- [ ] **Step 1: Write the failing tests.** `tests/unit/test_compare.py`:
  - `test_a_shared_mask_ignores_cells_outside_it`: if lineage masks agree on the shared cells and one run has extra `True` cells outside the mask, Jaccard is 1.0 with the mask and below 1.0 without it.
  - `test_without_a_mask_compare_is_unchanged`: equal to today's result field by field on the `bifurcation` fixture (Spearman, pairs, `rank_difference`).

  `tests/unit/test_stability.py` (600-cell, 400-gene bifurcation, seed 0, `paga_dpt` from `cluster:true_label=trunk_early`):
  - `test_replicate_plan_follows_the_declared_seeds`: `replicate_plan(StabilityConfig(), 0, {"seed": 20})` gives the reseeds `(1, 21), (2, 22)` and the subsamples `(1, 20), (2, 20), (3, 20)` as (seed, method_seed); with `{}`, every `method_seed` is `None`.
  - `test_subsample_keeps_the_root_and_draws_the_fraction`: `subsample_cells(1000, 17, 0.8, 1)` has 800 sorted unique values including 17, is deterministic, and differs for seed 2.
  - `test_reseed_reruns_only_leiden_and_paga`: `pca`, `connectivities` and `diffmap` are the same arrays, and `leiden` equals `run_stage(leiden, ds, art, replace(params, seed=1))`.
  - `test_subsample_prep_reuses_pca_and_recomputes_the_graph`: `pca == prep["pca"][cells]`, `hvg` is equal, `connectivities.shape == (k, k)`, `"lognorm" not in` the result, and the restricted obs index equals `ds.obs.index[cells]`.
  - `test_leiden_terminals_map_to_the_largest_overlap`: for `reference_leiden=[0,0,0,1,1,2]`, `replicate_leiden=[4,4,3,3,3]` and `cells=[0..4]`, terminals `["0","1"]` give `["4","3"]` and `["2"]` gives `[]`. When two terminals map to one cluster they collapse.
  - `test_pad_result_validates_and_matches_on_the_subset`.
  - `test_a_replicate_identical_to_the_reference_scores_one`: `score_replicate(r, r, all_true, 0.6) == ReplicateScores(n_finite, L, 1.0, 1.0, 1.0)`.
  - `test_every_replicate_is_rooted_at_the_reference_root_cell`: a spy around `run_method` sees `CellRoot(reference.root.obs_name)` five times.
  - `test_user_column_terminals_pass_through`: with `params={"cluster_key": "true_label"}` and terminals `["A_late", "B_late"]`, every call receives those terminals.
  - `test_a_failing_replicate_counts_as_no_agreement` (Review Focus 3): patch `run_method` to raise `RuntimeError` on its 2nd call only. Then `reseeds.rows[1].error` starts with `"RuntimeError"`, `math.isnan(reseeds.rows[1].share_matched)`, `reseeds.share_matched_min == 0.0`, and `lineage_stability` equals the median of the five `share_matched` values with that one as 0.
  - `test_summaries_count_undefined_values_as_zero`: `summarise` over shares `[1.0, nan, 1.0]` gives a median of 1.0 and a minimum of 0.0.
  - `test_a_repeated_root_name_is_refused`: rename a second cell to the root's obs name, and `run_stability` raises `StabilityError`.
  - `test_stability_runs_two_reseeds_and_three_subsamples`: `reseeds.r == 2`, `subsamples.r == 3`, `subsamples.fraction == 0.8`, `reseeds.fraction == 1.0`; every subsample row has `n_cells == 480`; `0 <= lineage_stability <= 1` and `-1 <= pseudotime_stability <= 1`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_compare.py tests/unit/test_stability.py tests/unit/test_trends.py -q`. Expected: FAIL (`lumencell.engine.stability` missing).
- [ ] **Step 3: Implement** the `compare.py` extension (keep `compare(a, b)` call sites as they are), then `stability.py`.
- [ ] **Step 4: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add stability replicates (two reseeds, three subsamples) compared with the reference on shared cells`.

---

### Task 7: The tree-support verdict

**Branch:** `pivot/NN-tree-support` · **Model:** sonnet. Fixed rules on given values.

**Files:**
- Create: `src/lumencell/engine/tree_support.py`, `tests/unit/test_tree_support.py`

**Interfaces:**
- Consumes: `COVERAGE_MIN`, `CLONAL_P_MAX`; `ClonalResult` (phase 1 T2).
- Produces, `engine/tree_support.py`:
  - `CheckId = Literal["C1", "C2", "C3", "C4", "C5"]`; `CheckState = Literal["pass", "fail", "not_applicable", "pending", "unavailable", "excluded"]`; `Level = Literal["supported", "weak", "not_supported", "provisional"]`; `CHECK_NAMES = {"C1": "connectivity", "C2": "tree-likeness", "C3": "lineage stability", "C4": "pseudotime stability", "C5": "clonal support"}`.
  - `tree_likeness(paga_connectivities: sparse.spmatrix, paga_tree: sparse.spmatrix) -> float`. Take `sym = conn.maximum(conn.T)`. The total is the sum of `triu(sym, k=1)`. The tree edges are the nonzero upper-triangle positions of `paga_tree + paga_tree.T`, weighted from `sym`. The result is the tree weight divided by the total, or NaN when the total is 0.
  - `coverage(pseudotime: FloatArray) -> float`, the share of finite values.
  - `@dataclass(frozen=True) class Thresholds: tau_tree: float | None; tau_lin: float | None; tau_pt: float | None; coverage_min: float = COVERAGE_MIN; clonal_p_max: float = CLONAL_P_MAX`. A `None` τ means the check was non-discriminating, so its state is `excluded`.
  - `@dataclass(frozen=True) class StabilityEvidence: state: Literal["done", "pending", "unavailable"]; lineage: float = math.nan; pseudotime: float = math.nan; note: str = ""`
  - `@dataclass(frozen=True) class ClonalEvidence: state: Literal["not_configured", "pending", "scored"]; status: str = ""; p: float = math.nan; delta_eta: float = math.nan; note: str = ""` with `@classmethod from_result(cls, r: ClonalResult) -> ClonalEvidence` (`state="scored"`).
  - `@dataclass(frozen=True) class Check: id: CheckId; name: str; state: CheckState; value: float; threshold: float | None; reason: str | None`. `reason` is set for `fail`, `pending`, `unavailable` and `excluded`, and is `None` otherwise. C5's `value` is the headline p.
  - `@dataclass(frozen=True) class Verdict: level: Level; checks: tuple[Check, ...]; reasons: tuple[str, ...]` (the non-`None` reasons in check order).
  - `verdict(coverage: float, tree_likeness: float | None, stability: StabilityEvidence, clonal: ClonalEvidence, thresholds: Thresholds) -> Verdict`. A `tree_likeness` of `None` means prep has not loaded, so C2 is `pending`; NaN means no edges, which is a `fail`. C3 and C4 are `pending` or `unavailable` (with the note) unless stability is `done`. C5 is `not_applicable` for `not_configured` or a scored `single_lineage`, and `pending` only while the clonal score is pending. A scored `ok` passes only if `p <= clonal_p_max and delta_eta > 0` (NaN fails). Any other scored status fails, with the status in the reason, `failed` included (spec: C5 applies and is not `ok`). The levels follow `[tree_support]` in that order. C5's reasons say whether clones support the lineages beyond chance and beyond expression clustering, never which state came first.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_tree_support.py`):
  - `test_tree_likeness_by_hand`: `conn = [[0,.5,.2],[.5,0,.4],[.2,.4,0]]` with a tree holding the edges (0,1) and (1,2) stored once gives `0.9 / 1.1`, to 1e-12. The same with a symmetric tree gives the same value.
  - `test_a_path_graph_scores_one`.
  - `test_a_graph_without_edges_fails_c2_with_a_reason` (Review Focus 5): a 1×1 zero matrix gives NaN tree-likeness. The verdict is `not_supported`, C2 is `fail`, and its reason is non-empty.
  - `test_tree_likeness_of_the_prep_is_in_the_unit_interval`, on the `bifurcation` fixture.
  - `test_levels` (parametrized; τs 0.7/0.8/0.85; values given as a pass or a fail):
    - C1 fails: `not_supported`.
    - C2 fails while C3–C5 pass: `not_supported`.
    - C2 is `None`: `provisional`.
    - C3 fails: `weak`.
    - C5 fails while stability is pending: `weak`.
    - Stability is pending and C5 does not apply: `provisional`.
    - Stability is unavailable: `provisional`.
    - Everything passes: `supported`.
    - `tau_pt=None` with the rest passing: `supported`, with C4 `excluded` and a reason saying it was left out.
  - `test_c1_bound_is_inclusive`: coverage exactly 0.95 passes.
  - `test_c5_rules` (parametrized): `ok` with p 0.01 and Δη 0.001 passes; p 0.011 fails; Δη 0 fails; Δη NaN fails; `too_few_clones` fails; `refused` fails; `failed` fails, so the level is `weak` even while stability is done and passes; `not_configured` and `single_lineage` give `not_applicable`; `pending` gives `pending`, so the level is `provisional`.
  - `test_from_result_reads_the_clonal_result`.
  - `test_reasons_one_per_check_that_did_not_pass_or_apply`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_tree_support.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** as in Interfaces.
- [ ] **Step 4: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add the tree-support verdict: five checks, three levels and provisional`.

---

### Task 8: Calibrate the verdict thresholds (slow, run once and recorded)

**Branch:** `pivot/NN-verdict-calibration` · **Model:** opus. Statistical design and interpretation; the rule must be applied exactly as declared.

**Files:**
- Create: `src/lumencell/bench/verdict_calibration.py`, `tests/unit/test_verdict_calibration.py`, `benchmarks/results/verdict-calibration-<date>.json`
- Modify: `benchmarks/thresholds.toml` (`[tree_support.calibrated]`, changelog), `src/lumencell/trust/constants.py`, `src/lumencell/engine/tree_support.py`, `tests/unit/test_trust_constants.py`

**Interfaces:**
- Consumes: `[tree_support.calibration]`, `[synthetic]`, `[synthetic.controls]` (Task 1); `generate_control`, `control_root_cell` (Task 5); `run_stability`, `RunSpec` (Task 6); `tree_likeness`, `coverage` (Task 7); `generate`, `run_pipeline`, `run_stage`, `run_method`.
- Produces, `bench/verdict_calibration.py` (Task 9 adds `verify`):
  - `Group = Literal["matrix", "control"]`.
  - `@dataclass(frozen=True) class RunRecord: group: Group; dataset: str; seed: int; method: str; coverage: float; tree_likeness: tuple[tuple[float, float], ...]; lineage_stability: float; pseudotime_stability: float; n_lineages: int; error: str | None = None`. `dataset` is `"<topology or kind>/seed<s>"`, and `tree_likeness` holds (resolution factor, value) pairs.
  - `collect_records(cfg: dict[str, Any], seeds: Mapping[Group, Sequence[int]], methods: Sequence[str], workdir: Path) -> list[RunRecord]`. Per dataset: write it, `load` it, and run `run_pipeline(ds, PrepParams())`. C2 is taken at each factor by rerunning the `leiden` and `paga` stages with `leiden_resolution = 1.0 × factor` on the same kNN; factor 1.0 uses the prep as it is. Per installed method: run the reference from `cluster:true_label=<root_label>` (matrix) or `cell:<control_root_cell>` (control), then `coverage`, then `run_stability(..., RunSpec(method, result.params, None, None))`. A reference that raises gives a record with `error` set and NaN values, plus one stderr line.
  - `@dataclass(frozen=True) class CheckCalibration: check: str; factor: float | None; matrix_quantile: float; control_quantile: float; n_matrix: int; n_control: int; discriminating: bool; tau: float | None`.
  - `@dataclass(frozen=True) class Calibration: checks: tuple[CheckCalibration, ...]; tau_tree: float | None; tau_lin: float | None; tau_pt: float | None; non_discriminating: tuple[str, ...]; n_failed: int`.
  - `derive(records: Sequence[RunRecord], cal: dict[str, Any]) -> Calibration` applies the `[tree_support.calibration]` rule exactly. That means `np.quantile(..., method="linear")`; C2 at one value per dataset per factor, with NaN taken as 0.0; and errored records left out and counted. C2 is kept only if every factor discriminates.
  - `calibration_payload(records, calibration, cfg, commit: str) -> dict[str, Any]` has top-level `"generated"`, `"commit"`, `"config"`, `"records"`, `"checks"`, `"taus"`, `"non_discriminating"`, `"n_failed"`, and `"complete"` (every expected (dataset, method) pair is present). `write_json(path: Path, payload: dict[str, Any]) -> None` writes NaN as `null`, with `allow_nan=False`.
- Produces, constants: `TAU_TREE: float | None`, `TAU_LIN: float | None`, `TAU_PT: float | None`. `engine/tree_support.py` gains `calibrated_thresholds() -> Thresholds` built from them.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_verdict_calibration.py`):
  - `test_tau_is_the_midpoint_of_the_declared_quantiles`: matrix C3 values `[0.6, 0.7, 0.8, 0.9, 1.0]` and control values `[0.1, 0.2, 0.3, 0.4, 0.5]` give `matrix_quantile = 0.64`, `control_quantile = 0.46` and `tau = 0.55`.
  - `test_ranges_that_overlap_or_touch_are_non_discriminating`: when the matrix quantile is ≤ the control quantile, `discriminating is False`, `tau is None`, and the check is listed in `non_discriminating`.
  - `test_c2_is_kept_only_if_it_discriminates_at_all_three_resolutions`: if factor 2.0 overlaps, then `tau_tree is None` even when 0.5 and 1.0 discriminate. Otherwise `tau_tree` is the factor-1.0 τ.
  - `test_c2_counts_each_dataset_once`: with three methods on one dataset, `n_matrix` equals the number of datasets.
  - `test_failed_runs_are_left_out_and_counted`, and `test_an_undefined_tree_likeness_counts_as_zero`.
  - `test_the_payload_is_strict_json`.
  - `test_calibration_on_the_declared_seeds` (`@pytest.mark.slow`). A module fixture runs `collect_records` on seeds 0–2 for every matrix topology and every control kind, with the three methods. When `LUMENCELL_VERDICT_CALIBRATION_OUT` is set, it writes the payload there before any assertion, so an honest result is always recorded. Assert `complete` and that there are 72 records (24 datasets × 3 methods), and print every `CheckCalibration`.
  - `test_the_committed_taus_follow_from_the_recorded_calibration`: rebuild the records from the newest `verdict-calibration-*.json`, then `derive` equals `[tree_support.calibrated]`.
  - In `tests/unit/test_trust_constants.py`, `CALIBRATED = {"TAU_TREE": ("C2", "tau_tree"), "TAU_LIN": ("C3", "tau_lin"), "TAU_PT": ("C4", "tau_pt")}`. The set-equality check becomes `set(MIRRORED) | set(CALIBRATED)`, and `test_calibrated_taus_equal_the_thresholds_file` asserts that each constant is `None` exactly when its check is in `non_discriminating`, and otherwise equals the TOML value.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_verdict_calibration.py tests/unit/test_trust_constants.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** `verdict_calibration.py`.
- [ ] **Step 3a: Smoke-run** without recording. Use a short `uv run python -` script on a copy of the config narrowed to `topologies = ["bifurcation"]` and `controls.kinds = ["blob"]`, with seeds `{"matrix": [0], "control": [0]}` and `paga_dpt` only, through `collect_records` and `derive`. Expected: two records, no `error`, and no exception.
- [ ] **Step 4 (background, hours):** `LUMENCELL_VERDICT_CALIBRATION_OUT=benchmarks/results/verdict-calibration-$(date +%F).json uv run pytest -m slow tests/unit/test_verdict_calibration.py -q > "$TMPDIR/verdict-cal.log" 2>&1; echo $?`.
- [ ] **Step 5:** Run the `data:statistical-analysis` and `data:validate-data` skills on the recorded file. Check that the per-check distributions of matrix against controls are plausible and that no record is NaN where it should not be. Check also that the quantile arithmetic in `checks` reproduces from `records`, and how close each τ sits to the two quantiles. Record their conclusions in HANDOFF.
- [ ] **Step 6: Commit the τs** exactly as derived. Append `[tree_support.calibrated]` with `tau_tree`, `tau_lin` and `tau_pt` written with `repr(float)` (a non-discriminating check's key is omitted), plus `non_discriminating = [...]`. Add the dated changelog line: "Verdict calibration, results in benchmarks/results/verdict-calibration-<date>.json: <one line per check: the two quantiles, τ, or non-discriminating>. Rule, seeds and controls unchanged." Then set the three constants and add `calibrated_thresholds()`. If C4 (or any other check) is non-discriminating, record it as such. That is acceptable. Never change the rule, the seeds or the controls to rescue a check.
- [ ] **Step 7: Run** `uv run pytest tests/unit/test_verdict_calibration.py tests/unit/test_trust_constants.py -m "not slow" -q`. Expected: PASS.
- [ ] **Step 8: C1–C4.** Commit message: `Calibrate the tree-verdict thresholds on synthetic seeds 0-2 against four negative controls`.

---

### Task 9: Verify the verdict on fresh seeds (slow, run once and recorded)

**Branch:** `pivot/NN-verdict-verification` · **Model:** sonnet. It applies an existing verdict to recorded runs and checks one condition.

**Files:**
- Modify: `src/lumencell/bench/verdict_calibration.py`, `tests/unit/test_verdict_calibration.py`, `benchmarks/thresholds.toml` (changelog only)
- Create: `benchmarks/results/verdict-verification-<date>.json`

**Interfaces:**
- Consumes: `collect_records`, `RunRecord`, `write_json` (Task 8); `verdict`, `StabilityEvidence`, `ClonalEvidence`, `calibrated_thresholds` (Tasks 7–8); `[tree_support.verification]`.
- Produces: `@dataclass(frozen=True) class Verification: verdicts: tuple[tuple[str, str, str], ...]; controls_supported: tuple[str, ...]; matrix_supported_share: float; n_failed: int; passed: bool`. `verdicts` holds (dataset, method, level). `verify(records: Sequence[RunRecord], thresholds: Thresholds, max_controls_supported: int) -> Verification` calls `verdict(r.coverage, <factor-1.0 tree_likeness>, StabilityEvidence("done", r.lineage_stability, r.pseudotime_stability), ClonalEvidence("not_configured"), thresholds)` for each record without an error. `passed` means `len(controls_supported) <= max_controls_supported`. `verification_payload(...)` has a top-level `"passed"`.

- [ ] **Step 1: Write the failing tests:**
  - `test_a_supported_control_fails_verification`.
  - `test_the_matrix_share_is_reported_not_gated`: every matrix record is `weak` and no control is supported, so `passed is True` and `matrix_supported_share == 0.0`.
  - `test_failed_runs_get_no_verdict_and_are_counted`.
  - `test_verification_on_fresh_seeds` (`@pytest.mark.slow`). Run seeds 3–5. Write `LUMENCELL_VERDICT_VERIFICATION_OUT` before any assertion, then assert `passed` and print the verdict counts by group and level.
  - `test_the_recorded_verification_follows_from_its_records`: rebuild the records from the newest `verdict-verification-*.json`, then `verify` gives the recorded `passed` and `controls_supported`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_verdict_calibration.py -m "not slow" -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `verify` and the payload.
- [ ] **Step 4 (background, hours):** `LUMENCELL_VERDICT_VERIFICATION_OUT=benchmarks/results/verdict-verification-$(date +%F).json uv run pytest -m slow tests/unit/test_verdict_calibration.py -k fresh_seeds -q > "$TMPDIR/verdict-ver.log" 2>&1; echo $?`.
- [ ] **Step 5:** Run `data:statistical-analysis` and `data:validate-data` on the file, and record the conclusions in HANDOFF. Add a dated changelog line: "Verdict verification on seeds 3–5, benchmarks/results/verdict-verification-<date>.json: <controls supported> of <n> control runs supported; <share> of matrix runs supported (reported, not gated)."
- [ ] **Step 6:** If verification fails, do not change the τs, the rule or the seeds. Commit the honest record, write the failure under "Open questions for Shreyan" in HANDOFF, open the PR, and stop for the owner. Tasks 11 and 12 wait for the owner's decision.
- [ ] **Step 7: C1–C4.** Commit message: `Verify the tree verdict on fresh seeds 3-5 and record the result`.

---

### Task 10: The `stability` job

**Branch:** `pivot/NN-stability-job` · **Model:** sonnet. Worker glue in phase 1's `score_clonal` pattern.

**Files:**
- Modify: `src/lumencell/api/models.py`, `src/lumencell/jobs/worker.py`, `src/lumencell/jobs/manager.py`, `src/lumencell/engine/store.py`
- Create: `src/lumencell/trust/stability.py`, `tests/unit/test_stability_trust.py`
- Test: `tests/unit/test_jobs.py`, `tests/unit/test_store.py`

**Interfaces:**
- Consumes: `run_stability`, `RunSpec`, `StabilityError`, `StabilityResult` (Task 6); `RunStore.write_trust`, `STORE_VERSION`, `graph["prep_key"]` (phase 1); stability constants (Task 1).
- Produces, `api/models.py`. These are defined here, and Task 11 references them from a route, so `make contract` shows no change yet.
  - `ReplicateKind = Literal["reseed", "subsample"]`
  - `StabilityReplicate(Model): kind: ReplicateKind; index: int; seed: int; method_seed: int | None; n_cells: int; n_compared: int; n_lineages: int; spearman: float | None; share_matched: float | None; mean_jaccard: float | None; terminals: list[str] | None; error: str | None`
  - `StabilitySummary(Model): kind: ReplicateKind; r: int; fraction: float; seeds: list[int]; spearman_median: float; spearman_min: float; share_matched_median: float; share_matched_min: float; mean_jaccard_median: float; mean_jaccard_min: float; n_lineages_median: float; n_lineages_min: int; replicates: list[StabilityReplicate]`
  - `StabilityStatus = Literal["ok", "pending", "failed", "not_available"]`
  - `Stability(Model): status: StabilityStatus; match_jaccard: float; reference_lineages: int | None; reseeds: StabilitySummary | None; subsamples: StabilitySummary | None; lineage_stability: float | None; pseudotime_stability: float | None; notes: list[str]`
- Produces, `trust/stability.py`:
  - `STABILITY_NOTE = "Reseeds rerun Leiden, PAGA and the method with new seeds on all cells. Subsamples rerun the neighbour graph, Leiden, diffusion maps, PAGA and the method on 80% of cells, reusing normalisation, HVG and PCA, which slightly understates subsampling variance."`
  - `stability_hash() -> str`: the first 16 hex characters of SHA-256 over `json.dumps({"metric_version": STABILITY_METRIC_VERSION, "params": {"n_reseeds", "n_subsamples", "subsample_fraction", "seed", "match_jaccard"}}, sort_keys=True)`. `stability_name() -> str` is `f"stability-{stability_hash()}"`.
  - `to_api_stability(r: StabilityResult) -> models.Stability`. The status is `ok`, every NaN becomes `None`, and `notes` holds `STABILITY_NOTE` plus `"<kind> <index + 1> failed: <error>"` for each failed row.
  - `not_measured(status: Literal["failed", "not_available"], note: str) -> models.Stability`.
- Produces, `engine/store.py`: `result.json` gains `"terminals"`, the requested terminals (`list[str] | None`), with no store-version bump.
- Produces, `jobs/manager.py`: `DEFAULT_TIME_LIMITS["stability"] = 6 * 3600`.
- Produces, `jobs/worker.py`: `HANDLERS["stability"] = handle_stability`. Its payload is the session's `path`, `x_is`, `prep_params` and `cache_root`, plus `run_id`, `trust_name`, `label_key` and `clonal_columns`. It writes `not_measured` in these cases:
  - `store_version != STORE_VERSION`: `not_available`, "made before inference isolation; run it again".
  - `graph["prep_key"] != prep_key(ds, params)`: `failed`, "the run was prepared with other settings".
  - `terminals_source == "user"` with no `"terminals"` key: `not_available`, "the terminals this run was asked for were not recorded; run it again".
  - `StabilityError`: `failed` with its message.

  Otherwise it loads the prep (`PrepCache(ds, cache_root).run(params)`) and the run (`RunStore.load(rid)`). It runs `run_stability` with `RunSpec(graph["method"], graph["params"], graph.get("terminals"), payload["label_key"], tuple(payload["clonal_columns"]))` and `progress=emit`, then writes `to_api_stability(...)`. It always writes with `write_trust(rid, payload["trust_name"], …)` and returns `rid`. Any other exception fails the job.

- [ ] **Step 1: Write the failing tests.**
  - `tests/unit/test_stability_trust.py`:
    - `test_the_name_is_a_valid_trust_name_and_follows_the_parameters`: it fully matches `TRUST_NAME`, and it changes when `MATCH_JACCARD` is monkeypatched.
    - `test_undefined_values_are_null_and_the_file_is_strict_json` (Review Focus 2): a `StabilityResult` with a failed row and NaN Spearman maps to `None` fields. `RunStore.write_trust` succeeds, and the file parses with `json.loads(text, parse_constant=_raise)`.
    - `test_failed_rows_are_noted`.
  - `tests/unit/test_store.py`, `test_result_json_records_the_requested_terminals`: `["3", "5"]` round-trips through `graph`, and `None` is stored for auto.
  - `tests/unit/test_jobs.py`:
    - `test_stability_job_writes_its_trust_file`: as in `test_prepare_then_infer_through_the_worker`, prepare and infer a 400-cell bifurcation, then submit `stability`. The job finishes `done`. `read_trust(rid, stability_name())["status"] == "ok"`, `reseeds.r == 2` and `subsamples.r == 3`. The `stages/` listing of the prep cache is identical before and after.
    - `test_a_cancelled_stability_job_leaves_nothing`: on a 1500-cell bifurcation, so the job runs for several seconds, `cancel` it once it is `running` and has sent its first progress event. Its state is `cancelled`, and `trust/` holds no `stability-*.json` and no `.stability-*.tmp`.
    - `test_stability_has_its_own_time_limit`: `JobManager().time_limit("stability") == 21600`.
    - `test_a_stale_run_is_reported_not_measured`: a run whose recorded `prep_key` is edited is written as `failed` with "prepared with other settings".
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_stability_trust.py tests/unit/test_store.py tests/unit/test_jobs.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the models (then `make contract`; expected: no diff), then `trust/stability.py`, the store field, the time limit and the handler.
- [ ] **Step 4: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add the background stability job, stored per run as a trust file`.

---

### Task 11: Serve stability, the verdict and the cell-cycle check

**Branch:** `pivot/NN-robustness-api` · **Model:** opus. The GET-healing, cancel and verdict paths are concurrency-sensitive and join four blocks.

**Files:**
- Modify: `src/lumencell/api/models.py`, `src/lumencell/api/runs.py`, `src/lumencell/api/session.py`, `src/lumencell/api/dataset_view.py`, `contracts/openapi.json`, `web/src/api/generated/*`, `web/src/test/fixtures.ts`
- Create: `src/lumencell/trust/verdict.py`, `tests/api/test_robustness_trust.py`

**Interfaces:**
- Consumes: Tasks 4, 7, 8 (`calibrated_thresholds`), 10 (models, `stability_name`, `stability_hash`, the handler); phase 1's `_trust_lookup`, `submit_trust_jobs`, `clonal_agreement` and `_on_job_done`.
- Produces, `api/models.py`:
  - `SupportCheck(Model): id: Literal["C1", "C2", "C3", "C4", "C5"]; name: str; state: Literal["pass", "fail", "not_applicable", "pending", "unavailable", "excluded"]; value: float | None; threshold: float | None; reason: str | None`
  - `TreeSupport(Model): level: Literal["supported", "weak", "not_supported", "provisional"]; checks: list[SupportCheck]; reasons: list[str]`
  - `RunTrust` gains `stability: Stability | None`, `tree_support: TreeSupport | None` and `cell_cycle: CellCycleCheck | None`, with no defaults. All three are filled for a done run, except `cell_cycle`, which is `None` while prep is not loaded.
  - `DatasetInfo` gains `tree_likeness: float | None` (C2 at the default resolution; `None` while prep is not loaded or when there are no edges) and `cell_cycle: CellCycleInfo | None`, with no defaults.
  - `JobKind` gains `"stability"`.
- Produces, `trust/verdict.py`: `stability_evidence(s: models.Stability) -> StabilityEvidence` (`ok` maps to `done`, `pending` to `pending`, and `failed` or `not_available` to `unavailable` with the first note); `clonal_evidence(c: models.ClonalAgreement | None) -> ClonalEvidence` (`None` or `not_configured` maps to `not_configured`, `pending` to `pending`, and every other status, `failed` included, to `scored` with that status and `primary.p` and `primary.delta_eta`, where `None` becomes NaN); `to_api_tree_support(v: Verdict) -> models.TreeSupport` (NaN becomes `None`); `tree_support_for_run(graph: dict[str, Any], tree_likeness: float | None, stability: models.Stability, clonal: models.ClonalAgreement | None, thresholds: Thresholds | None = None) -> models.TreeSupport`, where coverage is `1 − n_unassigned / n_obs` from the graph and the default is `calibrated_thresholds()`.
- Produces, `api/dataset_view.py`: `DatasetView.tree_likeness() -> float | None` returns `None` only while prep is not loaded, and NaN when PAGA has no edges, so `verdict` fails C2 rather than reading it as pending. `DatasetView.cell_cycle_info() -> models.CellCycleInfo | None` comes from `loaded_prep()` too. Only `Session.dataset_info()` maps NaN to `None` for the API.
- Produces, `api/runs.py`:
  - `Runs.stability(self, rid: str, graph: dict[str, Any]) -> models.Stability`. A legacy `store_version` gives `not_available` with no job. Otherwise it calls `_trust_lookup(rid, "stability", stability_name(), stability_hash(), payload)`, and a stored file is read with `Stability.model_validate` (a file that does not validate counts as missing, with one stderr line). An active job gives `pending`, and a failed or cancelled job gives `failed` with the note `job.error or "cancelled"`, and is not resubmitted. The dedup key is `f"stability:{rid}:{stability_hash()}"`.
  - `Runs.submit_trust_jobs(rid)` also submits stability, unless its file exists. The payload is `self._job_payload()` plus `run_id`, `trust_name=stability_name()`, `label_key=self.view.label_key` and `clonal_columns=list(self.view.clonal_columns)`, and the same builder serves the GET-healing path.
  - `Runs.cell_cycle(self, rid: str, names: list[str]) -> models.CellCycleCheck | None` runs `cell_cycle_check(self._load(rid).primary_lineage, names, …)` on the loaded prep.
  - `Runs.get_run` fills `trust=RunTrust(clonal=…, stability=…, tree_support=tree_support_for_run(g, self.view.tree_likeness(), stability, clonal), cell_cycle=…)`, using the lineage names after renames. It appends the cell-cycle warnings to `TrajectoryRun.warnings`.
- Produces, `api/session.py`: `dataset_info()` fills the two new fields.

- [ ] **Step 1: Write the failing tests** (`tests/api/test_robustness_trust.py`). A 400-cell bifurcation is prepared into a temporary cache, and each test builds its own `JobManager` and `Session`, as `tests/api/test_cache_healing.py` does.
  - `test_a_finished_infer_job_queues_stability_at_background_priority`: with the manager started, POST a `paga_dpt` run and wait for the infer job. Exactly one `stability` job for the run exists (its key ends with `stability_hash()`), and `"stability" in BACKGROUND_KINDS`. Wait for it, then GET: `stability.status == "ok"`, `stability.reseeds.r == 2`, `stability.subsamples.r == 3`, `stability.subsamples.fraction == 0.8`, and `tree_support.level` is one of `supported`, `weak` or `not_supported`. C5 is `not_applicable`, since there are no settings.
  - `test_foreground_jobs_overtake_a_queued_stability_job`: on an unstarted manager, submit `stability` and then `prepare`; `_take_next()` yields the prepare first.
  - `test_a_restart_heals_a_missing_stability_job`: save a run in-process, as `lumencell run` does. A new `Session` on an unstarted manager answers the GET with `stability.status == "pending"`, a `tree_support.level` other than `supported`, and exactly one queued `stability` job. Start, wait and GET: `ok`. A third GET submits nothing.
  - `test_a_cancelled_stability_keeps_the_verdict_provisional_and_is_not_resubmitted` (Review Focus 4): on an unstarted manager, a GET queues the job; `DELETE /api/v1/jobs/{id}` then returns `kind == "stability"` and state `cancelled` (Task 10 covers cancelling a running job). Two later GETs both show `stability.status == "failed"` with the note `"cancelled"`, `tree_support.level == "provisional"` with a reason naming stability, and `len(jobs.of_kind("stability")) == 1`.
  - `test_legacy_runs_get_no_stability_job`: a run folder with `store_version = 1` gives `not_available`, and no job is submitted.
  - `test_dataset_info_carries_tree_likeness_and_cell_cycle`: after prep, `0 < tree_likeness <= 1` and `cell_cycle.status == "not_scored"`. Before prep, both are `None`.
  - `test_a_cycling_lineage_is_warned_on_the_run`: set the loaded prep's `cc_phase` to `"S"` on lineage 0's single-lineage cells and `"G1"` elsewhere, and its `cc_info` to scored. The GET then has a `lineage_mostly_cycling` notice in both `trust.cell_cycle.warnings` and `warnings`.
  - `test_clonal_support_enters_the_verdict`: with clonal settings and both trust files done, C5 is not `not_applicable`.
  - `test_a_dataset_without_paga_edges_is_not_supported`: set the loaded prep's `paga_connectivities` and `paga_tree` to 1×1 zero matrices. The GET then gives `tree_support.level == "not_supported"` and C2 `fail` with `value is None` and a reason, and `GET /dataset` gives `tree_likeness is None`.
- [ ] **Step 2: Run** `uv run pytest tests/api/test_robustness_trust.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the models first, then `make contract`, then `trust/verdict.py`, the view, runs and session. Add `tree_likeness: null, cell_cycle: null` to the `dataset()` fixture in `web/src/test/fixtures.ts` and nothing else. (The `run` fixture already has `trust: null`.)
- [ ] **Step 4: Run** `uv run pytest tests/api -q`. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Serve stability, the tree verdict and the cell-cycle check with every run, healed on GET`.

---

### Task 12: Stability and the verdict on the real benchmarks

**Branch:** `pivot/NN-bench-real-stability` · **Model:** sonnet. Glue over specified functions, then a run and a logged report.

**Files:**
- Modify: `src/lumencell/bench/runner.py`, `src/lumencell/bench/report.py`, `src/lumencell/cli.py` (`bench real`), `benchmarks/thresholds.toml` (changelog), `README.md`
- Create: `tests/unit/test_runner_stability.py`, `benchmarks/results/real-<date>.json`

**Interfaces:**
- Consumes: `run_stability`, `RunSpec` (Task 6); `tree_likeness`, `coverage`, `verdict`, `StabilityEvidence`, `ClonalEvidence`, `calibrated_thresholds` (Tasks 7–8); phase 1's `CaseResult`/`row()` pattern.
- Produces, `bench/runner.py`:
  - `StabilityRowStatus = Literal["ok", "failed", "not_measured"]`.
  - `@dataclass(frozen=True) class StabilityRow: status: StabilityRowStatus; lineage: float = nan; pseudotime: float = nan; tree_likeness: float = nan; verdict: str = "not_measured"; reason: str | None = None`.
  - `CaseResult.stability: StabilityRow = StabilityRow("not_measured")`, as the new last field. `row()` always adds `stability_status`, `stability_lineage`, `stability_pt`, `tree_likeness` (rounded to 4 places) and `tree_support` (the level).
  - `stability_row(data: Dataset, prep: PrepResult, params: PrepParams, result: TrajectoryResult, label_key: str | None, clonal_columns: tuple[str, ...] = ()) -> StabilityRow`. It forwards `clonal_columns` to `RunSpec`, so replicates hide clone columns from methods as inference does; phase 3 calls it by this name. C5 is `ClonalEvidence("not_configured")`, because the developmental datasets carry no clones. Any exception gives `StabilityRow("failed", reason=f"{type(exc).__name__}: {exc}")`, also printed to stderr.
  - `run_real_case(dataset_id, method="paga_dpt", params=None, stability: bool = False)`.
- Produces, CLI: `bench real` gains `--stability/--no-stability` (default on), and each printed line gains `verdict=<level>`. Produces, report: each real-table cell gains `, verdict <level>` when the row has a `tree_support` field other than `not_measured`.

- [ ] **Step 0: Check the datasets are cached** with `uv run python -c "from lumencell.bench.datasets import DATASETS, is_cached; print({k: is_cached(v) for k, v in DATASETS.items()})"`. If any is `False`, stop and ask the owner before `uv run lumencell bench fetch all` (about 708 MB into `~/Library/Caches/lumencell/datasets`).
- [ ] **Step 1: Write the failing tests** (`tests/unit/test_runner_stability.py`):
  - `test_rows_always_carry_the_stability_fields`: a synthetic `CaseResult(...).row()` has `stability_status == "not_measured"` and a NaN `stability_lineage`.
  - `test_stability_row_on_a_synthetic_run`: for a 600-cell bifurcation with `paga_dpt`, the status is `ok`, the verdict is in the four levels, and `0 <= lineage <= 1`.
  - `test_clonal_columns_reach_every_replicate`: a spy on `run_stability` sees `spec.clonal_columns == ("clone_id", "donor")`.
  - `test_a_stability_error_is_recorded_as_failed`: monkeypatch `run_stability` to raise, and the status is `failed` with a reason naming the error.
  - `test_bench_real_has_a_stability_switch`: `--help` lists `--no-stability`.
  - `test_the_report_shows_the_verdict`: `render()`, given a temporary real file whose rows carry `tree_support = "weak"`, contains `verdict weak`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_runner_stability.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** as in Interfaces.
- [ ] **Step 4: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 5 (background):** `make bench-real > "$TMPDIR/real.log" 2>&1; echo $?` → `0`. Gate failures do not fail the command.
- [ ] **Step 6:** Compare with phase 1's `real-<date>.json`. Every gated metric and `n_clusters` must be identical, because stability changes no reference run. Any difference is an error: stop and report it. Check the stability and verdict columns with `data:validate-data`. Add a dated changelog line: "Stability and the tree verdict reported on the real benchmarks (not gated): <per dataset and method, C3, C4 and the level>. No gate changed." Then run `uv run lumencell bench report --readme`.
- [ ] **Step 7: C1–C4.** Commit message: `Report stability and the tree verdict for every method on the real benchmarks`.

---

### Task 13: Security review of phase 2

**Branch:** `pivot/NN-phase2-security` · **Model:** sonnet. It runs the review and fixes small findings test-first.

**Files:**
- Modify: whatever the findings require, and `docs/superpowers/HANDOFF.md`.

- [ ] **Step 1:** Invoke the `security-review` skill on the diff from phase 1's last merge to `main`. Point it at the untrusted inputs: obs names in replicate roots (`cell:<obs_name>` text, the repeated-name refusal), lineage names in notices, gene symbols in cell-cycle matching, the package-data TSV parsing, the trust file names and job keys, the `LUMENCELL_VERDICT_*_OUT` paths, the cancel path, and the handling of the downloaded sources (never committed, parsed with `-I`).
- [ ] **Step 2:** For each finding, write a failing test that shows it, fix it, and see the test pass. A finding too large for this PR gets its own task in HANDOFF instead.
- [ ] **Step 3:** Add a "Security review (phase 2)" section to `docs/superpowers/HANDOFF.md` listing every finding and whether it was fixed (with its test) or recorded, and why.
- [ ] **Step 4 (background):** Run the Done When commands at the top of this plan and record each result in HANDOFF.
- [ ] **Step 5: C1–C4.** Commit message: `Resolve or record the phase 2 security review findings`.
