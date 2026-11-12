# Phase 1: Clonotypes End to End — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Score every trajectory run against TCR clonotypes (clonal agreement against three nulls and a lineage-free baseline), keep clone identity out of inference, gate the synthetic benchmark on it, and serve the result through the API and the CLI.

**Architecture:** A pure package `src/lumencell/clonal/` (labels to clone ids; scores, nulls, baseline, coupling) is called by the benchmark in oracle-cells mode and by a background `score_clonal` worker job in own-cells mode. Clone and receptor columns become scoring-only (removed from the obs methods see, refused at every key), and receptor genes are dropped from HVG selection. Results are stored per run under `<run>/trust/clonal-<h>.json`, served as `TrajectoryRun.trust`, and healed by a GET that finds neither the file nor a job.

**Tech Stack:** Python 3.13 via uv; numpy, scipy (`scipy.cluster.vq.kmeans2`, `scipy.stats.false_discovery_control`), pandas, scanpy; FastAPI + pydantic v2; typer; pytest + Hypothesis. `web/` gets only regenerated types and fixture updates.

**Spec:** `docs/superpowers/specs/2026-11-12-tcell-clonotypes-design.md` — §2 (constraints), §3 (terms), §4.1–4.7 (this phase), §7 (interfaces). Read §3 and §4 before any task. The spec wins over this plan on any conflict; report the conflict instead of choosing.

**Baseline:** `main` at 142a788. Before starting, `git fetch origin && git rev-list --left-right --count main...origin/main` must print `0	0`.

## Global Constraints

Every task's requirements include these (spec §2).
- Clone labels are used only for scoring: never through obs (a prep stage, a method, a root, a terminal, `cluster_key`, `label_key`) and never through expression (genes matching `^TR[ABDG]([VDJ]\d|C\d*$)`, case-insensitive, are excluded from HVG selection).
- Each metric, null and the baseline is defined in `benchmarks/thresholds.toml` before any result is seen; every score is reported against its nulls.
- Clone identity is the pair (donor, label). The same label string in two donors is two clones. A file without a donor column must say `single_donor = true`.
- The minimum clone size is stated (m = 2 primary, m = 3 sensitivity); every score reports cells, clones and the effective clone count; single-cell clones never contribute.
- No output or text claims direction (which state came first) from clones.
- No method default and no threshold is tuned on T-cell data.
- Thresholds are never lowered. New gates and definitions go into `thresholds.toml`, with their reasons and a date, before the first run that measures them.
- Roots are always explicit. Inference never uses the UMAP. The input file is never modified.
- A method's `experimental` flag matches its synthetic gates (`tests/golden/test_synthetic.py`).
- API changes start in `src/lumencell/api/models.py`, then `make contract`. In `web/`, only `make contract` output and the smallest updates to the `dataset` and `run` fixtures in `web/src/test/fixtures.ts`.
- The server process never runs prep or methods; heavy work runs in the worker. A GET that finds a missing artefact may submit a job.
- State-file writes are locked and atomic (`api/state.py`); run, cache and trust files are written by atomic rename.
- `engine/` never imports `bench/`. `clonal/` imports nothing from `engine/`, `api/`, `bench/` or anndata (only `lumencell.trust.constants`).
- No new runtime dependency. scikit-learn is never imported directly.
- Datasets are never committed. A download asks the owner first. Real-data tests are marked `real_data`; slow tests `slow`. CI uses synthetic data only.
- Values read from a user's file (column names, labels) never reach a file path or a shell command; derived files use hashed names.
- Runtime thresholds live in `src/lumencell/trust/constants.py`; a test asserts they equal the TOML.
- Process: test-first; one branch and PR per task (`pivot/NN-<slug>`); squash-merge only on green CI after review by a different agent; nothing is pushed to `main` directly; no `Co-Authored-By` trailers or "Generated with" footers; each commit gated on the exit code of `make ci`.

## Review Focus

Five inputs the spec implies but no spec test exercises, most likely first. Each has a test in the task named.
1. **Undefined statistics** (zero scored clones, one clone, a null with zero spread, a lineage with no scored cells): every affected field is `null` in the API and the trust file is strict JSON, never `NaN`. Test: Task 13, `test_undefined_values_are_null_and_the_file_is_strict_json`.
2. **Messy label columns** (categorical with unused categories, `None`/`NaN`/`pd.NA`/`""`/`"nan"`/`"None"`, a labelled cell with a missing donor): missing labels are unlabelled, unused categories are not clones, a missing donor drops the cell from scoring and is counted in a note, never pooled into a `"nan"` donor. Tests: Task 2 `test_missing_donor_drops_the_cell_and_counts_it`, Task 13 `test_missing_donors_are_noted`.
3. **Clonal settings that name the saved root's column** (for example a donor column the user had used as `cluster:` root): `PUT /dataset/clonal` answers 422 `scoring_only_column` naming the root, rather than the root silently disappearing. Test: Task 10, `test_settings_that_capture_the_saved_root_are_refused`.
4. **A scoring job that fails**: the run shows `status = "failed"` with the job's reason, and repeated GETs do not resubmit it. Test: Task 14, `test_a_failed_score_is_reported_and_not_resubmitted_on_every_get`.
5. **Settings changed while a score is queued or running**: a result computed for old settings is never served for new ones. Test: Task 14, `test_a_result_for_old_settings_is_never_served`.

## Task Conventions

- **Branch:** each task header names `pivot/NN-<slug>`; NN is assigned at execution time.
- **Model:** each task names its implementer model with a reason.
- **TDD:** write the failing test, run it and see it fail for the stated reason, implement, run it green.
- **Closing steps** (every task ends with these; "C1–C4" below means this list):
  - C1. Update `docs/superpowers/HANDOFF.md`: current state, decisions taken in this task, the branch and PR row.
  - C2. `make ci; echo $?` prints `0`. Gate the commit on that exit code, never on piped output.
  - C3. Commit with the task's message. No `Co-Authored-By` trailer, no "Generated with" footer.
  - C4. `git push -u origin pivot/NN-<slug>` and open the PR with `gh pr create`. Do not merge; the orchestrator squash-merges after review.
- **Background runs:** steps marked **(background)** take tens of minutes to hours. Start them with the Bash tool's `run_in_background: true` and keep the output log under `$TMPDIR`; the orchestrator may start the next independent task meanwhile.
- **Imports:** heavy imports inside CLI commands stay inside the function, as `cli.py` does today.

## Process Notes (not tasks)

- After each PR: run the pr-review-toolkit agents `code-reviewer`, `silent-failure-hunter` and `pr-test-analyzer`. Add `type-design-analyzer` on PRs that add types: Tasks 1, 2, 3, 8, 10, 13, 14.
- At the end of the phase: `/security-review` (Task 17).
- Statistical work goes through the `data:statistical-analysis` and `data:validate-data` skills: the calibration run (Task 8) and the bench re-measures' comparisons (Tasks 7 and 9).
- README benchmark tables are regenerated by Tasks 7 and 9; whichever merges second rebases and regenerates again.

## Not in This Plan

- Spec §9: web views, the multi-user version, new trajectory methods, RNA velocity, receptor design, packaging, the 100k-cell performance budgets, cell-cycle regression, excluding immunoglobulin genes, AIRR or scirpy import in the app, similarity-based clonotype clusters, any direction inferred from clones, per-clone pseudotime coherence.
- Phase 2: the `cell_cycle` stage and `Stage.report_only`, stability replicates (`engine/stability.py`, the `stability` job handler and its time limit), the tree-support verdict (`engine/tree_support.py`), negative controls (`bench/synth_controls.py`), `RunTrust.stability` / `tree_support` / `cell_cycle`, verdict τs and cell-cycle ratios in `trust/constants.py`. Phase 1 only reserves `"stability"` in `BACKGROUND_KINDS` and leaves the `_trust_lookup` and `submit_trust_jobs` hooks for it.
- Phase 3: Daniel 2022 and Sheih 2020 (fetch, builders, `RealDataset.sources`/`licence`, truth, namespace, clonal gates on real data, the case study).
- A web form for clonal settings, and scoring runs with a clone column per run.

## File Map

| File | Responsibility | Task |
|---|---|---|
| `benchmarks/thresholds.toml` | clonal definitions, gates, calibration plan, changelog | 1, 7, 8, 9 |
| `src/lumencell/trust/__init__.py` | package marker (empty) | 1 |
| `src/lumencell/trust/constants.py` | runtime thresholds mirrored from the TOML | 1 |
| `src/lumencell/bench/gates.py` | per-gate NaN and status policy; `clonal_excess_ratio` in `HIGHER_IS_BETTER` | 1 |
| `src/lumencell/bench/runner.py` | `ClonalRow`, always-emitted row fields; clonal layer and oracle scoring | 1, 6 |
| `src/lumencell/clonal/__init__.py` | package marker (empty) | 2 |
| `src/lumencell/clonal/labels.py` | raw labels and donors to clone, donor and stratum codes | 2 |
| `src/lumencell/clonal/scores.py` | entropy, η, nulls, p, baseline, coupling | 2 |
| `src/lumencell/clonal/isolation.py` | receptor gene and column patterns, scoring-only set, isolation version | 4, 5 |
| `src/lumencell/bench/truth.py` | `Truth.clone_ids`, `Truth.donors` | 3 |
| `src/lumencell/bench/synth_clones.py` | synthetic clone layer and its two controls | 3 |
| `src/lumencell/prep/stages.py` | `exclude_receptor_genes`, `hvg_info` | 4 |
| `src/lumencell/cli.py` | `prepare` receptor line; `run --clone-key …`; `synth --fate-bias` | 4, 15 |
| `src/lumencell/engine/context.py` | filtered obs, `scoring_only`, `ScoringOnlyColumnError` | 5 |
| `src/lumencell/engine/roots.py` | refuse a scoring-only `cluster:` root | 5 |
| `src/lumencell/engine/registry.py` | refuse scoring-only `cluster_key`/`label_key`; record `obs_keys` | 5 |
| `src/lumencell/engine/schema.py` | `TrajectoryResult.obs_keys` | 5 |
| `src/lumencell/engine/store.py` | store v2 fields; `write_trust`/`read_trust`; `save` keeps `trust/` | 5, 11 |
| `src/lumencell/jobs/worker.py` | infer passes clonal columns and prep key; `score_clonal` handler | 5, 13 |
| `src/lumencell/api/runs.py` | create-run refusal; trust and GET healing | 5, 14 |
| `src/lumencell/bench/report.py` | clonal ratio column and the "not a branch gate" note | 6 |
| `Makefile` | `BENCH_METHODS` for `bench-synth` and `bench-baseline` | 7 |
| `benchmarks/baseline.json` | re-measured baseline with clonal fields | 7 |
| `src/lumencell/bench/calibration.py` | calibration runs and their report | 8 |
| `benchmarks/results/clonal-calibration-<date>.json` | recorded calibration | 8 |
| `benchmarks/results/real-<date>.json` | real re-run after the TR-gene exclusion | 9 |
| `README.md` | regenerated benchmark tables | 7, 9 |
| `src/lumencell/api/models.py` | `ClonalSettings`, `ClonalSummary`, `DatasetInfo` fields; trust models; `JobKind` | 10, 13, 14 |
| `src/lumencell/trust/clonal.py` | settings checks and summary; `agreement_for_run`; trust hash | 10, 13 |
| `src/lumencell/api/dataset_view.py` | one column-kind rule; clonal columns and context rebuild | 5, 10 |
| `src/lumencell/api/session.py` | clonal settings in state; `job_info`; done-callback | 10, 14 |
| `src/lumencell/api/routes.py` | `PUT`/`DELETE /dataset/clonal`; `/obs` description | 10 |
| `src/lumencell/prep/housekeeping.py` | trust temp files count as leftovers | 11 |
| `src/lumencell/jobs/manager.py` | two priorities, done-callbacks, `score_clonal` time limit | 12 |
| `contracts/openapi.json`, `web/src/api/generated/*` | `make contract` output | 10, 14 |
| `web/src/test/fixtures.ts` | `dataset` and `run` fixtures | 10, 14 |
| `docs/HANDOFF.md` | stale lines corrected | 16 |
| `docs/superpowers/HANDOFF.md` | state, decisions, PRs | every task |
| tests: `tests/unit/test_trust_constants.py`, `tests/unit/test_gates.py`, `tests/clonal/{__init__,test_labels,test_scores,test_properties,test_calibration}.py`, `tests/unit/test_synth_clones.py`, `tests/unit/test_prep.py`, `tests/unit/test_prep_cache.py`, `tests/unit/test_isolation.py`, `tests/unit/test_roots.py`, `tests/unit/test_store.py`, `tests/unit/test_cache_admin.py`, `tests/unit/test_jobs.py`, `tests/unit/test_clonal_agreement.py`, `tests/unit/test_commands.py`, `tests/unit/test_cli.py`, `tests/unit/test_makefile.py`, `tests/unit/test_runner_clonal.py`, `tests/golden/test_synthetic.py`, `tests/api/test_clonal_settings.py`, `tests/api/test_clonal_trust.py`, `tests/api/test_contract.py` | the tests named in each task | per task |

## Done When

Spec §4.7, as commands. Run from the repo root on `main` after Task 17.
- `make ci; echo $?` → `0`.
- **(background)** `make bench-synth > "$TMPDIR/synth.log" 2>&1; echo $?` → `0`, and `grep -c '^REGRESSION' "$TMPDIR/synth.log"` → `0`. (`BENCH_METHODS` covers all three methods since Task 7.)
- `ls benchmarks/results/real-*.json | tail -1` names the Task 9 file, and `grep -n 'TR-gene exclusion' benchmarks/thresholds.toml` shows its dated changelog entry.
- `uv run pytest tests/unit/test_isolation.py tests/clonal/test_labels.py tests/clonal/test_scores.py tests/clonal/test_properties.py -q; echo $?` → `0` (includes the slow Palantir and scFates isolation cases).
- `uv run python -c "import json,glob; f=sorted(glob.glob('benchmarks/results/clonal-calibration-*.json'))[-1]; print(f, json.load(open(f))['passed'])"` → the file and `True`.
- `uv run pytest tests/api/test_clonal_trust.py -q; echo $?` → `0` (trust through the API; restart heals).
- `d=$TMPDIR/p1; LUMENCELL_CACHE_DIR=$d/cache LUMENCELL_DATA_DIR=$d/state uv run lumencell synth $d/c.h5ad --fate-bias 0.7 && LUMENCELL_CACHE_DIR=$d/cache LUMENCELL_DATA_DIR=$d/state uv run lumencell run $d/c.h5ad --root cluster:true_label=trunk_early --clone-key clone_id --donor-key donor` → output contains `clonal agreement: ok`.
- `uv run pytest tests/unit/test_trust_constants.py -q` → passes, and `git log --reverse --format='%h %ad %s' --date=iso 142a788.. -- benchmarks/thresholds.toml benchmarks/baseline.json benchmarks/results | head -1` is the Task 1 commit.
- `git diff 142a788 -- docs/HANDOFF.md` corrects lines 3, 32, 49, 71–73, 89–90, 159–160 and 185.
- `docs/superpowers/HANDOFF.md` has a "Security review (phase 1)" section with each finding resolved or recorded.

---

### Task 1: Declare the clonal definitions, gates and runtime constants

**Branch:** `pivot/NN-clonal-definitions` · **Model:** sonnet — configuration plus small, fully specified code.

**Files:**
- Modify: `benchmarks/thresholds.toml`, `src/lumencell/bench/gates.py`, `src/lumencell/bench/runner.py`, `tests/unit/test_gates.py`
- Create: `src/lumencell/trust/__init__.py` (empty), `src/lumencell/trust/constants.py`, `tests/unit/test_trust_constants.py`

**Interfaces:**
- Produces, `trust/constants.py` (phase 2 appends its τs and cell-cycle ratios here):
  `CLONAL_METRIC_VERSION = 1`, `MIN_CLONE_SIZE = 2`, `SENSITIVITY_CLONE_SIZE = 3`, `MIN_CLONES = 20`, `N_PERMUTATIONS = 1000`, `CLONAL_SEED = 0`, `PSEUDOTIME_BINS = 10`, `CLONAL_P_MAX = 0.01`, `ORACLE_ETA_MIN = 0.05`, `ORACLE_P_MAX = 0.01`, and `MAX_DONORS = 1000` (a validation limit from spec §4.1, not mirrored in the TOML).
- Produces, `tests/unit/test_trust_constants.py`: `MIRRORED: dict[str, tuple[str, ...]]`, constant name → TOML key path, e.g. `"MIN_CLONES": ("clonal", "min_clones")`, `"ORACLE_ETA_MIN": ("clonal", "synthetic", "oracle_eta_min")`. Phase 2 appends entries.
- Produces, `bench/gates.py`:
  - `HIGHER_IS_BETTER` gains `"clonal_excess_ratio"`.
  - `GateFailure(case: str, metric: str, value: float, bound: str, reason: str | None = None)`; `str()` appends ` [<reason>]` when set.
  - `check_gates(row: dict[str, Any], gates: dict[str, dict[str, Any]]) -> list[GateFailure]`. A gate entry may carry `nan` (`"skip"`, the default, or `"fail"`; anything else raises `ValueError`, checked for every gate whether or not its value is NaN) and `status` (the name of a row field). With `status` set: status `"not_applicable"` skips the gate; any status other than `"ok"` fails it with `reason=<status>`; a missing status field counts as `"not_measured"`; a missing metric counts as NaN. With status `"ok"`, a NaN value fails under `nan = "fail"` with `reason="nan"`. Gates without `status` behave exactly as today (missing key raises `KeyError`, NaN skips).
- Produces, `bench/runner.py`:
  - `ClonalRowStatus = Literal["ok", "single_lineage", "too_few_clones", "not_assessable", "failed", "not_applicable", "not_measured"]`
  - `@dataclass(frozen=True) class ClonalRow: status: ClonalRowStatus; excess_ratio: float = nan; p: float = nan; delta_base: float = nan; eta_oracle: float = nan; eta_method: float = nan; n_clones: int = 0; n_cells: int = 0; effective_clones: float = nan; eta_m3: float = nan; p_m3: float = nan; reason: str | None = None`
  - `CaseResult.clonal: ClonalRow = ClonalRow("not_measured")` (new last field).
  - `CaseResult.row()` always adds `clonal_excess_ratio`, `clonal_p`, `clonal_delta_base` (each `round(·, 4)`, NaN stays NaN), `clonal_status` (str), `clonal_eta_oracle`, `clonal_eta_method`, `clonal_effective_clones`, `clonal_eta_m3`, `clonal_p_m3` (each rounded 4), `clonal_n_clones` and `clonal_n_cells` (int). The first four are the spec's row fields; the other seven are the evidence counts of §2 and §4.3, reported only and never gated (NaN or 0 when not applicable). Phases 2 and 3 read these exact names.

- [ ] **Step 1: Write the failing tests.** In `tests/unit/test_gates.py`, extend `_row()` with `clonal_excess_ratio=0.8, clonal_p=0.001, clonal_status="ok"`, update the expected gate set in `test_thresholds_file_has_the_documented_gates` to add `clonal_excess_ratio` and `clonal_p`, and add:

```python
GATES = load_thresholds()["synthetic"]["gates"]

def test_clonal_gate_entries_are_declared() -> None:
    assert GATES["clonal_excess_ratio"] == {"min": 0.5, "nan": "fail", "status": "clonal_status"}
    assert GATES["clonal_p"] == {"max": 0.01, "nan": "fail", "status": "clonal_status"}

def test_clonal_nan_fails_with_its_reason() -> None:
    (f,) = check_gates(_row(clonal_excess_ratio=math.nan), GATES)
    assert (f.metric, f.reason) == ("clonal_excess_ratio", "nan")

@pytest.mark.parametrize("status", ["too_few_clones", "failed", "not_assessable", "single_lineage", "not_measured"])
def test_any_status_but_ok_fails_both_clonal_gates(status: str) -> None:
    failed = check_gates(_row(clonal_status=status), GATES)
    assert {(f.metric, f.reason) for f in failed} == {("clonal_excess_ratio", status), ("clonal_p", status)}

def test_not_applicable_skips_the_clonal_gates() -> None:
    assert check_gates(_row(clonal_status="not_applicable", clonal_excess_ratio=math.nan, clonal_p=math.nan), GATES) == []

def test_a_row_from_before_the_clonal_fields_fails_as_not_measured() -> None:
    old = {k: v for k, v in _row().items() if not k.startswith("clonal_")}
    assert {f.reason for f in check_gates(old, GATES)} == {"not_measured"}

def test_existing_gates_still_skip_nan_and_unknown_policies_raise() -> None:
    assert check_gates(_row(branch_accuracy=math.nan), GATES) == []
    with pytest.raises(ValueError, match="nan"):
        check_gates(_row(), {"branch_accuracy": {"min": 0.9, "nan": "maybe"}})

def test_regression_guard_covers_the_clonal_ratio() -> None:
    assert [f.metric for f in check_regression(_row(clonal_excess_ratio=0.6), {"clonal_excess_ratio": 0.7}, 0.03)] == ["clonal_excess_ratio"]

def test_case_rows_always_carry_the_clonal_fields() -> None:
    row = CaseResult("c", "m", scores=Scores(<every float 0.0, ints 0>, pt_spearman_by_lineage={}, matches=[]), n_clusters=1, runtime_seconds=0.0).row()
    assert row["clonal_status"] == "not_measured" and math.isnan(row["clonal_excess_ratio"])
    assert {"clonal_p", "clonal_delta_base", "clonal_eta_oracle", "clonal_eta_method", "clonal_n_clones",
            "clonal_n_cells", "clonal_effective_clones", "clonal_eta_m3", "clonal_p_m3"} <= row.keys()
    assert row["clonal_n_cells"] == 0 and math.isnan(row["clonal_effective_clones"]) and math.isnan(row["clonal_p_m3"])
    assert "clonal_n_cells" not in load_thresholds()["synthetic"]["gates"]  # reported only
```

  In `tests/unit/test_trust_constants.py`, `test_runtime_constants_equal_the_thresholds_file`: for each `MIRRORED` entry, walk the TOML path and assert `getattr(constants, name) == value`; also assert `set(MIRRORED)` equals every upper-case name in `constants` except `MAX_DONORS`.

- [ ] **Step 2: Run** `uv run pytest tests/unit/test_gates.py tests/unit/test_trust_constants.py -q`. Expected: FAIL (no `clonal_*` gates, no `lumencell.trust`).

- [ ] **Step 2a: Verify the method references** that spec §10 lists as unverified and phase 1 cites: Phipson & Smyth 2010 and Vinh, Epps & Bailey 2010 (the definitions), and Sturm et al. 2020 (scirpy) and Schattgen et al. 2022 (CoNGA) (the reason for the baseline). For each, read its Crossref record (WebFetch on `https://api.crossref.org/works?query.bibliographic=<title and year>&rows=3`, then `https://api.crossref.org/works/<doi>`) or the publisher's page; nothing is downloaded. Record title, authors, journal, year and DOI (JMLR has no DOI: record volume and pages from jmlr.org). Put them in the Sources lines below and in HANDOFF, with the date checked. A reference that cannot be verified stays marked "unverified" in both. The spec is not edited.
- [ ] **Step 3: Add to `benchmarks/thresholds.toml`.** Inside `[synthetic.gates]`, after `root_pt_excess`:

```toml
clonal_excess_ratio = { min = 0.5, nan = "fail", status = "clonal_status" }  # see [clonal]
clonal_p = { max = 0.01, nan = "fail", status = "clonal_status" }
```

  Before `[regression]`, this block (comments are the definitions the spec requires before any result):

```toml
# ---- clonal agreement (spec docs/superpowers/specs/2026-11-12-tcell-clonotypes-design.md §4.3-4.5)
# Defined before any clonal score was computed. src/lumencell/trust/constants.py mirrors what the
# worker needs at run time; tests/unit/test_trust_constants.py keeps the two equal.
# Clone: the pair (donor, label). A missing, empty, "nan" or "None" label means no receptor.
# Scored cells: single-lineage cells with a label whose clone has >= m such cells.
#   own-cells mode (API): the run's single-lineage cells, counted in their lineage.
#   oracle-cells mode (gates): the true lineages' single-lineage cells, weighted by the method's
#   masks: one lineage 1, k lineages 1/k each, none 1/L each (L = the method's lineages).
# H_k = -sum_l p_kl log p_kl / log L; H = sum_k n_k H_k / sum_k n_k; purity P = 1 - H.
# eta = 1 - H_obs / E_null[H], reported against each null; headline eta is against N1.
# Nulls shuffle clone identities among the scored cells, B times, seeded, within:
#   N1 donor; N2 donor x decile of pseudotime (the run's own; the truth's order for the oracle;
#   cells without pseudotime form their own stratum); N3 donor x strata_keys, when declared.
# p = (1 + #{H_b <= H_obs}) / (B + 1), one-sided; headline p = the largest over computed nulls.
# Lineage-free baseline: scipy.cluster.vq.kmeans2 with k = L, seeded, on the PCs the kNN graph
#   uses, scored on the same cells: eta_base, delta_eta = eta - eta_base, and a paired
#   permutation p for delta_eta over the same N1 shuffles.
# Coupling (reported, never gated): cosine similarity of the raw scored-cell count columns of
#   each lineage pair, against N1: z, p and Benjamini-Hochberg p_adj. Not CoSpar's normalised
#   default. It says two lineages share clones, never which came first.
# Every score also reports cells and clones scored, the effective clone count (inverse Simpson
#   index of scored mass, 1 / sum_k (n_k / N)^2) and the same block at m = 3.
# Sources: STARTRAC, Zhang et al., Nature 2018, doi:10.1038/s41586-018-0694-x (purity);
#   LRT, Xie et al., PLoS Comput Biol 2023, doi:10.1371/journal.pcbi.1011300 (lineage bias);
#   CoSpar, Wang et al., Nat Biotechnol 2022, doi:10.1038/s41587-022-01209-1 (coupling);
#   Vinh, Epps & Bailey, JMLR 2010, <reference verified in Step 2a> (chance-corrected eta);
#   Phipson & Smyth, Stat Appl Genet Mol Biol 2010, <DOI verified in Step 2a> (p never zero);
#   Sturm et al. 2020 and Schattgen et al. 2022, <DOIs verified in Step 2a> (clones are
#   transcriptionally coherent, hence the lineage-free baseline).
[clonal]
metric_version = 1
min_clone_size = 2          # primary m; single-cell clones never contribute
sensitivity_clone_size = 3  # always reported alongside
min_clones = 20             # fewer scored clones: too_few_clones, values reported, never a pass
n_permutations = 1000       # B per null
seed = 0
pseudotime_bins = 10
p_max = 0.01                # the clonal_p gate, and phase 2's C5

[clonal.synthetic]
fate_bias = 0.7             # every branching matrix case gets this clonal layer
n_donors = 3
label_rate = 0.85           # share of cells with a receptor, near 10x capture rates
singleton_share = 0.5       # clone size: 1 with this probability, else 1 + Zipf(zipf_a)
zipf_a = 2.2
donor_confound_share = 0.8  # donor-confound control: share of lineage j's cells in donor j mod D
oracle_eta_min = 0.05       # below either oracle floor the case is not_assessable (and fails)
oracle_p_max = 0.01

[clonal.calibration]
# Slow, run once and recorded (tests/clonal/test_calibration.py), at the [synthetic] preset and
# size. Seeds are disjoint from the gated matrix. Seed s uses topologies[(s - seed_start) % 3].
seed_start = 100
n_seeds = 200
topologies = ["bifurcation", "trifurcation", "binary2"]
depth_only_topology = "binary2"   # the topology whose lineages occupy different depths
alpha = 0.05
binomial_quantile = 0.99          # binom.ppf(0.99, 200, 0.05) = 18 rejections at most
monotonic_fate_bias = [0.0, 0.35, 0.7, 1.0]
monotonic_seeds = [0, 1, 2]
monotonic_topologies = ["bifurcation", "trifurcation", "binary2"]
```

  Changelog entry at the end, dated with this commit's date:

```toml
# <YYYY-MM-DD>  Clonal agreement: definitions above and two synthetic gates, set before any
#   clonal score was computed. clonal_excess_ratio = eta_method / eta_oracle >= 0.5: a method
#   whose lineages capture less than half of the chance-corrected clonal signal of the true
#   lineages is not separating the fates the clones show. clonal_p <= 0.01: with B = 1000 the
#   smallest p is about 0.001 and its Monte Carlo error at the boundary about 0.003. A case is
#   assessable only if the oracle has eta >= 0.05 and headline p <= 0.01. These two gates fail
#   on NaN and on every status but ok (too_few_clones, failed, not_assessable, single_lineage),
#   except not_applicable (fewer than 2 true lineages); every other gate still skips NaN. The
#   clonal gate is not a branch gate: merging two true lineages can push the ratio above 1, and
#   lineage_recall catches that.
```

- [ ] **Step 4: Implement** `trust/constants.py` (module docstring: the TOML is not in the wheel, so the worker reads these; a test keeps them equal), `check_gates`/`GateFailure`/`HIGHER_IS_BETTER` in `gates.py`, and `ClonalRowStatus`/`ClonalRow`/`CaseResult.clonal`/`row()` in `runner.py` as in Interfaces.
- [ ] **Step 5: Run** `uv run pytest tests/unit/test_gates.py tests/unit/test_trust_constants.py -q`. Expected: PASS.
- [ ] **Step 6: C1–C4.** Commit message: `Declare clonal agreement definitions, gates and runtime constants before any measurement`.

---

### Task 2: Clone labels and clonal scores

**Branch:** `pivot/NN-clonal-scores` · **Model:** opus — statistics-heavy; the shuffle scheme decides which invariances hold exactly.

**Files:**
- Create: `src/lumencell/clonal/__init__.py` (docstring only), `src/lumencell/clonal/labels.py`, `src/lumencell/clonal/scores.py`
- Test: `tests/clonal/__init__.py`, `tests/clonal/test_labels.py`, `tests/clonal/test_scores.py`, `tests/clonal/test_properties.py`

**Interfaces:**
- Consumes: Task 1 constants (`MIN_CLONE_SIZE`, `MIN_CLONES`, `N_PERMUTATIONS`, `CLONAL_SEED`, `PSEUDOTIME_BINS`), and the verified Sources lines in `thresholds.toml`.
- Produces, `clonal/labels.py`:
  - `is_missing(values: npt.NDArray[np.object_]) -> npt.NDArray[np.bool_]` — `None`, float NaN, `pd.NA`, or `str(v)` in `{"", "nan", "None"}`.
  - `@dataclass(frozen=True) class CloneCodes: clone: IntArray; donor: IntArray; stratum: IntArray | None; n_missing_donor: int` — all `(n,)` int64. `clone` is −1 for an unlabelled cell or a labelled cell with a missing donor, otherwise a dense id over distinct (donor code, label string) pairs. `donor` is −1 where the donor is missing. `stratum` codes the tuple (donor, each strata value) when strata are given, a strata value that `is_missing` flags being its own level distinct from every string; `None` otherwise. `is_missing` defines missing for labels, donors and strata alike. `n_missing_donor` counts labelled cells whose donor is missing.
  - `encode_clones(labels: npt.NDArray[np.object_], donors: npt.NDArray[np.object_] | None, strata: Sequence[npt.NDArray[np.object_]] = ()) -> CloneCodes` — `donors=None` means one donor (all 0). Callers pass `obs[col].to_numpy(dtype=object)`.
- Produces, `clonal/scores.py`:
  - `Mode = Literal["own_cells", "oracle_cells"]`, `NullName = Literal["N1", "N2", "N3"]`, `ScoreStatus = Literal["ok", "single_lineage", "too_few_clones"]`
  - `@dataclass(frozen=True) class NullComparison: null: NullName; mean_entropy_mean: float; mean_entropy_sd: float; eta: float; p: float`
  - `@dataclass(frozen=True) class Coupling: a: int; b: int; observed: float; null_mean: float; z: float; p: float; p_adj: float`
  - `@dataclass(frozen=True) class ClonalResult:` `status: ScoreStatus; mode: Mode; n_lineages: int; n_cells: int; n_clones: int; effective_clones: float; n_donors: int; n_strata: int | None; min_clone_size: int; n_permutations: int; seed: int; purity: float; mean_entropy: float; nulls: tuple[NullComparison, ...]; eta: float; p: float; eta_base: float; delta_eta: float; p_delta: float; clone_weighted_purity: float; purity_without_largest: float; coupling: tuple[Coupling, ...]` — `nulls` are in order N1, N2, then N3 when strata exist; `n_donors` and `n_strata` count values among the scored cells.
  - `baseline_space(pca: npt.NDArray[np.floating], n_pcs: int) -> npt.NDArray[np.float64]` — `pca[:, :min(n_pcs, pca.shape[1])]` as float64; callers pass `PrepParams.neighbors_pcs`.
  - `kmeans_partition(space: npt.NDArray[np.float64], k: int, seed: int) -> IntArray` — `kmeans2(space, k, minit="++", seed=seed)` labels as int64, over all cells.
  - `score_clonal(membership: BoolArray, codes: CloneCodes, pseudotime: FloatArray, baseline: IntArray, *, oracle_cells: BoolArray | None = None, min_clone_size: int = MIN_CLONE_SIZE, min_clones: int = MIN_CLONES, n_permutations: int = N_PERMUTATIONS, seed: int = CLONAL_SEED, n_bins: int = PSEUDOTIME_BINS) -> ClonalResult` — `membership` is `(L, n)`; `pseudotime` and `baseline` are `(n,)`. `oracle_cells=None` is own-cells mode (candidates: cells in exactly one lineage of `membership`); otherwise oracle-cells mode (candidates: `oracle_cells`, weighted from `membership` as spec §3 says). Scored cells: candidates with `clone >= 0` whose clone has at least `min_clone_size` candidates.
  - Status: L < 2 → `single_lineage` (counts computed, every float NaN, `nulls` and `coupling` empty). `n_clones < min_clones` → `too_few_clones`, values still computed. `n_clones == 0` → floats NaN, `nulls` and `coupling` empty.

- [ ] **Step 1: Write the failing label tests** (`tests/clonal/test_labels.py`):
  - `test_missing_values_are_unlabelled`: labels `[None, np.nan, "", "nan", "None", pd.NA, "c1"]` → `clone[:6] == -1`, `clone[6] >= 0`.
  - `test_clone_identity_is_donor_and_label`: labels `["a","a","a","b"]`, donors `["d1","d2","d1","d1"]` → `clone[0] == clone[2] != clone[1]`, three distinct ids.
  - `test_missing_donor_drops_the_cell_and_counts_it` (Review Focus 2): labels `["a","a","a"]`, donors `[None, "d1", "d1"]` → `clone[0] == -1`, `donor[0] == -1`, `n_missing_donor == 1`, `clone[1] == clone[2]`.
  - `test_unused_categories_are_not_clones`: `pd.Series(pd.Categorical(["a","a"], categories=["a","b","zz"])).to_numpy(dtype=object)` → one distinct id.
  - `test_strata_combine_with_donor_and_missing_is_its_own_level`: donors `["d1","d1","d1","d2"]`, strata `["t1","t1",None,"t1"]` → `stratum[0] == stratum[1]`; `stratum[2]` differs from both; `stratum[3] != stratum[0]` (another donor).
  - `test_no_donor_column_means_one_donor`: `donors=None` → `donor` all 0.

- [ ] **Step 2: Write the failing score tests** (`tests/clonal/test_scores.py`), hand-computed:

```python
def test_entropy_purity_and_diagnostics_by_hand() -> None:
    # L=2, one donor. Clone A: two cells in lineage 0. Clone B: one cell in each lineage.
    r = score_clonal(
        own([[1, 1, 1, 0], [0, 0, 0, 1]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        min_clones=1,
        n_permutations=99,
    )
    assert r.mean_entropy == pytest.approx(0.5) and r.purity == pytest.approx(0.5)
    assert r.clone_weighted_purity == pytest.approx(0.5)
    assert r.purity_without_largest == pytest.approx(0.0)  # tie on size: the purer clone (A) is left out
    assert r.effective_clones == pytest.approx(2.0) and (r.n_cells, r.n_clones) == (4, 2)


def test_oracle_cells_use_fractional_counts() -> None:
    # cells: c0 in both method lineages, c1 in lineage 0, c2 in none, c3 in lineage 1; clones X={c0,c1}, Y={c2,c3}
    r = score_clonal(
        np.array([[1, 1, 0, 0], [1, 0, 0, 1]], bool),
        codes(["X", "X", "Y", "Y"]),
        PT0,
        BASE0,
        oracle_cells=np.ones(4, bool),
        min_clones=1,
        n_permutations=99,
    )
    assert r.mode == "oracle_cells" and r.mean_entropy == pytest.approx(0.811278, abs=1e-6)


def test_null_mean_and_p_match_enumeration() -> None:
    # 4 cells (L0,L0,L1,L1), clones {0,1},{2,3}: 2 of the 6 equally likely pairings give H=0, 4 give H=1.
    r = score_clonal(
        own([[1, 1, 0, 0], [0, 0, 1, 1]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        min_clones=1,
        n_permutations=9999,
    )
    n1 = r.nulls[0]
    assert n1.null == "N1" and n1.mean_entropy_mean == pytest.approx(2 / 3, abs=0.01)
    assert r.eta == 1.0 and n1.p == pytest.approx(1 / 3, abs=0.02)


def test_p_is_never_below_one_over_b_plus_one(): ...  # every null p >= 1/(B+1)
def test_headline_p_is_the_largest_null_p_and_eta_is_against_n1(): ...
def test_n3_exists_only_with_strata(): ...  # [n.null for n in r.nulls] == ["N1","N2"] / ["N1","N2","N3"]
def test_baseline_equal_to_the_lineages_gives_zero_delta(): ...  # baseline == lineage index: delta_eta == 0 (abs 1e-12), eta_base == eta
def test_one_lineage_is_single_lineage_with_nan_values(): ...  # status, math.isnan(r.purity), r.nulls == ()
def test_too_few_clones_still_reports_values(): ...  # 5 clones, min_clones=20: status too_few_clones, math.isfinite(r.purity)
def test_coupling_by_hand() -> None:
    # same data as the enumeration test: X = [[2,0],[0,2]]; null cosine is 1 in 4/6 pairings, 0 in 2/6
    (c,) = r.coupling
    assert (
        (c.a, c.b, c.observed) == (0, 1, 0.0) and c.null_mean == pytest.approx(2 / 3, abs=0.01) and c.p == 1.0
    )


def test_cells_without_pseudotime_form_their_own_n2_stratum(): ...
def test_kmeans_partition_is_seeded(): ...  # same seed, same labels; values in range(k)
def test_sources_are_cited(): ...  # clonal.__doc__ holds the three spec DOIs and the two Step 2a references
```

  `own(rows)` builds a bool membership; `codes(labels)` calls `encode_clones` with `donors=None`; `PT0` is all-zero pseudotime; `BASE0` an all-zero baseline.

- [ ] **Step 3: Write the failing properties** (`tests/clonal/test_properties.py`, Hypothesis, `@settings(max_examples=30, deadline=None)`, `n_permutations=499`). A strategy draws L in 2..4, n in 20..120, 1..3 donors, labels from a 12-string alphabet with some missing, membership, pseudotime with some NaN, and a baseline in `range(L)`. A helper `same(a, b)` compares two `ClonalResult`s field by field: ints and strings exactly, floats with `math.isclose(abs_tol=1e-12)`, NaN equal to NaN.
  - `test_relabelling_clones_or_lineages_changes_nothing`: a bijection on label strings, and a permutation of membership rows (coupling pairs mapped through it) → `same(...)`.
  - `test_reordering_cells_keeps_deterministic_fields_and_null_fields_within_monte_carlo_error`: n_cells, n_clones, effective_clones, purity, mean_entropy, clone_weighted_purity, purity_without_largest equal within 1e-12; for each null, `|Δp| <= 5*sqrt(0.5/B) + 2/(B+1)` and `|Δeta| <= 5*sqrt(2/B)*mean_entropy*sd/mean**2 + 1e-9` (from that null's reported mean and sd).
  - `test_singleton_clones_never_change_a_score`: append cells each with a new unique label (any donor, membership, pseudotime, baseline) → `same(...)` exactly.
  - `test_purities_lie_in_the_unit_interval`: P, clone-weighted purity and purity without the largest clone are in [0, 1] when finite.
  - `test_splitting_clones_across_donors_never_merges_them`: number of distinct `clone >= 0` ids equals the number of distinct (donor, label) pairs among labelled cells with a donor.
  - `test_ignoring_donors_changes_the_score_when_labels_collide`: donor A's copies of a shared label all in lineage 0, donor B's in lineage 1 → purity with donors > purity with `donors=None`.

- [ ] **Step 4: Run** `uv run pytest tests/clonal -q`. Expected: FAIL (`lumencell.clonal` missing).

- [ ] **Step 5: Implement `labels.py` and `scores.py`.** The shuffle is fixed, because it is what makes relabelling and added singletons exact: null i (N1 = 1, N2 = 2, N3 = 3) has its own `np.random.default_rng([seed, i])`; each of the B shuffles draws `u = rng.random(n_scored)` over the scored cells in their original order and permutes clone ids within strata `g` by `lexsort((u, g))`; the method weights, the baseline and (N1 only) the coupling all use those same shuffled ids. Strata: N1 donor; N2 donor × decile, with decile `min(n_bins − 1, ⌊n_bins·(average rank − 1)/n_finite⌋)` over the finite scored pseudotime and NaN in a bin of its own; N3 `codes.stratum`. The docstrings of `clonal/__init__.py` (all five) and `scores.py` (each where used: STARTRAC and LRT for purity, Vinh et al. for η, Phipson & Smyth for p, CoSpar for coupling) cite the sources with the references recorded in `thresholds.toml`. Comparisons use a 1e-12 slack: `H_b <= H_obs + 1e-12`, `cos_b >= cos_obs - 1e-12`, `Δ_b >= Δη - 1e-12`. `eta = 1 - H_obs / mean(H_b)` (NaN when the mean is 0); sd uses `ddof=1`. η_base and Δη are against N1; `p_delta = (1 + #{(η_b − η_base,b) >= Δη − 1e-12}) / (B + 1)` with `η_b = 1 − H_b / mean(H_b)` and `η_base,b = 1 − Hbase_b / mean(Hbase_b)` over the N1 shuffles. Coupling: z uses the N1 null sd (`z` NaN when sd is 0); a pair with an all-zero column gets NaN fields and is left out of the Benjamini–Hochberg adjustment (`scipy.stats.false_discovery_control(p, method="bh")`). Purity without the largest clone: NaN with one clone. Effective clones: `1 / sum((n_k / N) ** 2)`.

- [ ] **Step 6: Run** `uv run pytest tests/clonal -q`. Expected: PASS.
- [ ] **Step 7: C1–C4.** Commit message: `Add the clonal package: clone codes, entropy-based agreement, three nulls, k-means baseline and coupling`.

---

### Task 3: Synthetic clones

**Branch:** `pivot/NN-synth-clones` · **Model:** sonnet — the model and draw order are fixed here.

**Files:**
- Create: `src/lumencell/bench/synth_clones.py`, `tests/unit/test_synth_clones.py`
- Modify: `src/lumencell/bench/truth.py`

**Interfaces:**
- Consumes: `encode_clones`, `score_clonal` (Task 2); `generate`, `truth_from_adata` (public, unchanged).
- Produces, `bench/truth.py`: `Truth.clone_ids: npt.NDArray[np.object_] | None = None` and `Truth.donors: npt.NDArray[np.object_] | None = None`, appended after `terminals`. `__post_init__` checks both are `(n,)` and that both or neither are set.
- Produces, `bench/synth_clones.py` (phase 2's controls may reuse it on any AnnData carrying the same truth):
  - Input contract: obs `true_segment` and `true_pseudotime`, and `uns["lumencell_truth"]["segments"]` (name → `{"parent": str}`); a leaf is a segment that is no one's parent. Nothing private from `synth.py` is imported, and `synth.py` is not modified.
  - Module constants equal to `[clonal.synthetic]`: `FATE_BIAS = 0.7`, `N_DONORS = 3`, `LABEL_RATE = 0.85`, `SINGLETON_SHARE = 0.5`, `ZIPF_A = 2.2`, `DONOR_CONFOUND_SHARE = 0.8`; plus `CLONE_STREAM = 7100501`.
  - `Control = Literal["none", "depth_only", "donor_confound"]`
  - `@dataclass(frozen=True) class ClonalLayer: fate_bias: float = FATE_BIAS; n_donors: int = N_DONORS; label_rate: float = LABEL_RATE; control: Control = "none"`; `DEFAULT_LAYER = ClonalLayer()`.
  - `add_clones(adata: ad.AnnData, layer: ClonalLayer, seed: int) -> ad.AnnData` — returns a copy; the input is untouched. obs gains `clone_id` (categorical; NaN for unlabelled) and `donor` (categorical `d0`…, every cell); `uns["lumencell_truth"]["clones"]` records the layer, `seed` and `CLONE_STREAM`.
  - `clone_truth(truth: Truth, adata: ad.AnnData) -> Truth` — `dataclasses.replace` with `clone_ids` (`obs["clone_id"].to_numpy(dtype=object)`) and `donors`.
  - `draw_clone_size(rng: np.random.Generator) -> int` — `1` with probability `SINGLETON_SHARE`, else `1 + int(rng.zipf(ZIPF_A))`.
- Model and draw order (from `np.random.default_rng([seed, CLONE_STREAM])`): (1) `labelled = rng.random(n) < label_rate`; (2) donors: uniform `rng.integers(n_donors, size=n)`, or for `donor_confound` a cell exclusive to true lineage j (sorted leaf names) goes to donor `j % D` with probability `DONOR_CONFOUND_SHARE`, else uniformly to another donor, and all other cells uniformly; (3) `join_leaf = rng.random(n) < fate_bias` (control `none`); (4) per donor in ascending order, pools in the order trunk pool then leaves by sorted name: a labelled leaf cell with `join_leaf` joins its leaf's pool, every other labelled cell the trunk pool; each pool is permuted with `rng.permutation`, then cut into consecutive clones of `draw_clone_size` cells, only the last clone truncated. Label strings are `c0, c1, …`, counting across a donor's pools and restarting in each donor. `depth_only`: per donor, one pool of all labelled cells sorted by true pseudotime (stable), cut consecutively. `donor_confound`: per donor, one permuted pool of all labelled cells.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_synth_clones.py`; data from `generate(..., n_cells=3000, n_genes=50, seed=0)` unless stated):

```python
def _digest(adata: ad.AnnData) -> str:
    h = hashlib.sha256()
    for a in (adata.X.data, adata.X.indices, adata.X.indptr):
        h.update(np.ascontiguousarray(a).tobytes())
    h.update(adata.obs.to_csv().encode())
    h.update(json.dumps(adata.uns["lumencell_truth"], sort_keys=True, default=str).encode())
    return h.hexdigest()


def test_generator_output_is_unchanged_byte_for_byte() -> None:
    # Computed at 142a788. If only Linux CI disagrees, stop and report; never edit this value.
    a = generate("binary2", "realistic", n_cells=400, n_genes=120, seed=7)
    assert _digest(a) == "87d55724becfbc9790e0b8697e42f39077d4c14feed972da02a1f41d45e78e3e"
    before = _digest(a)
    add_clones(a, DEFAULT_LAYER, seed=7)
    assert _digest(a) == before
```

  Plus:
  - `test_the_layer_adds_two_columns_and_its_spec_only`: X, var and the original obs columns of the output equal the input's; obs gains exactly `clone_id`, `donor`; `uns["lumencell_truth"]["clones"]["fate_bias"] == 0.7`.
  - `test_same_seed_same_clones` and a different seed differs.
  - `test_labelled_share_and_sizes_do_not_depend_on_fate_bias`: for fate_bias in `[0, 0.35, 0.7, 1]` on bifurcation, labelled share in `[0.83, 0.87]`; singleton shares within 0.06 of each other; the 75th and 90th percentile clone sizes differ by at most 1 across levels.
  - `test_full_fate_bias_confines_every_leaf_cells_clone_to_its_leaf` (parametrized over bifurcation, trifurcation, binary2): at fate_bias 1, every clone holding a leaf-segment cell has all its cells on that one leaf segment.
  - `test_label_strings_collide_across_donors`: `c0` occurs in every donor.
  - `test_depth_only_clones_sit_in_a_pseudotime_window`: median over clones of (max − min true pseudotime) `< 0.05`.
  - `test_donor_confound_skews_donor_composition`: among labelled cells exclusive to true lineage j, the share in donor `j % 3` is `>= 0.6`.
  - `test_clone_truth_fills_and_checks_the_truth`: `clone_truth` sets both arrays aligned with obs; `Truth(..., clone_ids=<wrong length>)` raises `ValueError`; setting only one raises `ValueError`.
  - `test_layer_defaults_equal_the_thresholds_file`: module constants equal `[clonal.synthetic]`.
  - `test_oracle_eta_rises_from_no_bias_to_full_bias` (CI-sized monotonicity): bifurcation 1500 cells seed 0, oracle-cells mode on the truth (`oracle_cells = truth.exclusive() >= 0`), true pseudotime, a constant baseline, `n_permutations=199` → `eta(fate_bias=1) > eta(fate_bias=0)`.

- [ ] **Step 2: Run** `uv run pytest tests/unit/test_synth_clones.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** the `Truth` fields and `bench/synth_clones.py` as in Interfaces.
- [ ] **Step 4: Run** `uv run pytest tests/unit/test_synth_clones.py tests/unit/test_synth.py -q`. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add a synthetic clone layer on its own random stream, with depth-only and donor-confound controls`.

---

### Task 4: Keep receptor genes out of HVG selection

**Branch:** `pivot/NN-receptor-genes` · **Model:** sonnet — one stage change with a clear rule.

**Files:**
- Create: `src/lumencell/clonal/isolation.py`
- Modify: `src/lumencell/prep/stages.py`, `src/lumencell/cli.py` (`prepare`)
- Test: `tests/unit/test_prep.py`, `tests/unit/test_prep_cache.py`, `tests/unit/test_commands.py`

**Interfaces:**
- Produces, `clonal/isolation.py`: `RECEPTOR_GENE = re.compile(r"^TR[ABDG]([VDJ]\d|C\d*$)", re.IGNORECASE)` (spec pattern); `receptor_gene_mask(symbols: npt.NDArray[np.str_]) -> npt.NDArray[np.bool_]`.
- Produces, `prep/stages.py`: `PrepParams.exclude_receptor_genes: bool = True`; the `hvg` stage's params become `("n_top_genes", "exclude_receptor_genes")` and its outputs `("hvg", "hvg_info")`. Matching is on `ds.gene_symbols`. `hvg_info = {"receptor_genes": int, "excluded": bool, "receptor_hvg_without_exclusion": int}`. With exclusion on and at least one receptor gene: HVGs are selected on the non-receptor genes with `n_top = min(n_top_genes, n_non_receptor)` (receptor genes are `False` in `hvg`), and HVG selection is also run once on all genes only to count `receptor_hvg_without_exclusion`. With no receptor genes, or exclusion off, the stage runs exactly as today (one call, all genes), so synthetic output is unchanged.

- [ ] **Step 1: Write the failing tests.** In `tests/unit/test_prep.py` (a fixture appends eight receptor genes `TRAV1-2, TRAV12-1, TRAJ33, TRAC, TRBV5-1, TRBV20-1, TRBJ2-7, TRBC2` to a 1500-cell, 1000-gene synthetic bifurcation; gene i's counts are Poisson(30) on the true lineage it marks (i mod the `true_lineage` categories) and Poisson(0.5) elsewhere, `rng = np.random.default_rng(0)`, so they rank among the top HVGs). Every test here uses `PrepParams(n_top_genes=200)`: with the default 2000 every one of the 1,008 genes would be an HVG and the tests would pass vacuously.
  - `test_receptor_genes_are_never_hvgs`: `art["hvg"].sum() == 200` exactly (a select-then-mask implementation leaves fewer), `art["hvg"][tr].sum() == 0`, `hvg_info == {"receptor_genes": 8, "excluded": True, "receptor_hvg_without_exclusion": <k>}` with `k >= 1`.
  - `test_receptor_genes_can_be_kept_on_request`: `PrepParams(n_top_genes=200, exclude_receptor_genes=False)` → `art["hvg"].sum() == 200`, `art["hvg"][tr].sum() == hvg_info["receptor_hvg_without_exclusion"] >= 1` and `excluded is False`.
  - `test_synthetic_hvgs_are_unchanged_by_the_exclusion`: on the plain synthetic dataset with `n_top_genes=200`, `hvg` with exclusion on equals `hvg` with it off, `hvg.sum() == 200`, and `hvg_info["receptor_genes"] == 0`.
  - `test_the_pattern_covers_human_and_mouse_segments`: `receptor_gene_mask` is True for `TRAV1-2, TRAV14DV4, Trbc1, trgv9, TRDC, TRAC, TRBD1, TRBJ1-1`; False for `TRADD, TRAT1, TRAPPC1, TRDN, TRDMT1`.
  In `tests/unit/test_prep_cache.py`: `test_changing_the_receptor_exclusion_recomputes_hvg_and_downstream`: the `normalize` key is unchanged; `hvg` and every later stage key changes.
  In `tests/unit/test_commands.py`: `test_prepare_reports_receptor_genes` on a file with the eight genes: output contains `receptor genes: 8`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_prep.py tests/unit/test_prep_cache.py tests/unit/test_commands.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `isolation.py`, the HVG change, and in `cli.prepare` one line after the stage table when `hvg_info["receptor_genes"] > 0`: `receptor genes: {n} matched ^TR[ABDG]([VDJ]\d|C\d*$); {excluded from|kept in} HVG selection ({m} would be HVGs without the exclusion)`.
- [ ] **Step 4: Run** the same tests. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Exclude TCR segment genes from HVG selection so clone identity cannot reach the PCA`.

---

### Task 5: Keep clone columns out of inference

**Branch:** `pivot/NN-inference-isolation` · **Model:** opus — touches the engine, store and worker invariants at once.

**Files:**
- Modify: `src/lumencell/clonal/isolation.py`, `src/lumencell/engine/context.py`, `src/lumencell/engine/roots.py`, `src/lumencell/engine/registry.py`, `src/lumencell/engine/schema.py`, `src/lumencell/engine/store.py`, `src/lumencell/jobs/worker.py`, `src/lumencell/api/dataset_view.py`, `src/lumencell/api/runs.py`, `src/lumencell/cli.py` (`run` saves with `prep_key`)
- Test: `tests/unit/test_isolation.py` (new), `tests/unit/test_roots.py`, `tests/unit/test_store.py`, `tests/api/test_routes.py`

**Interfaces:**
- Consumes: `RECEPTOR_GENE` and `PrepParams.exclude_receptor_genes` (Task 4).
- Produces, `clonal/isolation.py`: `ISOLATION_VERSION = 1`; `RECEPTOR_COLUMN = re.compile(r"^(?:clone_id|clonotype|cdr3|ir_|airr:|receptor_)|^(?:raw_clonotype_id|chain_pairing)$", re.IGNORECASE)`; `scoring_only_columns(columns: Iterable[object], named: Iterable[str] = ()) -> frozenset[str]` — every column whose `str` name matches `RECEPTOR_COLUMN`, plus every name in `named`.
- Produces, `engine/context.py`:
  - `class ScoringOnlyColumnError(ValueError)`.
  - `inference_obs(obs: pd.DataFrame, clonal_columns: Iterable[str] = ()) -> tuple[pd.DataFrame, frozenset[str]]` — a copy without the scoring-only columns, and that set.
  - `InferenceContext.scoring_only: frozenset[str] = frozenset()` (new last field); `InferenceContext.from_prep(ds, prep, params, clonal_columns: Iterable[str] = ())` builds `obs` and `label_key` from `inference_obs`.
  - `check_inference_column(ctx: InferenceContext, column: str, role: str) -> None` raises `ScoringOnlyColumnError(f"{column!r} holds clone or receptor data, which lumenCell uses only to score runs; it cannot be the {role}")`.
- Produces, refusals: `resolve_root` raises `RootError` with that message for a `ClusterRoot` whose key is in `ctx.scoring_only`, checked before the existence check. `run_method` calls `check_inference_column` for the method's `cluster_key` (role "cluster column", which terminals refer to) and for the `label_key` it uses (role "label column").
- Produces, `engine/schema.py`: `TrajectoryResult.obs_keys: list[str] = field(default_factory=list)` — set by `run_method` to the sorted, non-empty set of the resolved root's `cluster_key`, the method's `cluster_key` and the label key used.
- Produces, `engine/store.py`: `STORE_VERSION = 2`; `RunStore.save(rid, result, terminals=None, prep_key: str | None = None)`; `result.json` gains `"prep_key"`, `"inference_isolation": ISOLATION_VERSION` and `"obs_keys"`; `load` restores `obs_keys` (`[]` when absent).
- Produces, worker and CLI: the infer payload carries `"clonal_columns": list[str]`; `handle_infer` passes it to `from_prep` and saves with `prep_key=prep_key(ds, params)`; `cli.run` saves with `prep_key`.
- Produces, `api/dataset_view.py`: `DatasetView.clonal_columns: tuple[str, ...]` (empty until Task 10); `label_key` is detected on `inference_obs(ds.obs, clonal_columns)`; `ctx()` passes `clonal_columns`.
- Produces, `api/runs.py`: `create_run` checks, before resolving the root, the `cluster:` root key, `params["cluster_key"]` and `view.label_key` against `scoring_only_columns(view.ds.obs.columns, view.clonal_columns)`, answering `ApiError(422, "scoring_only_column", <the message above>, {"columns": [...]})`; the infer payload gains `clonal_columns`.

- [ ] **Step 1: Write the failing isolation test** (`tests/unit/test_isolation.py`). Build the plain `Dataset` from a 1500-cell, 1000-gene synthetic bifurcation (seed 0) with the eight receptor genes of Task 4 appended, each cell's eight counts drawn `rng.poisson(4)` independently of lineage (`rng = np.random.default_rng(0)`), written and `load`ed once. Build the clonal variant in memory with `dataclasses.replace(plain, obs=…, matrix=…)`: obs gains `clone_id` = a copy of `true_lineage`; each cell's eight receptor counts are the same multiset reordered so its largest count sits at index `(true_lineage code * 3) % 8` (the receptor genes encode the clone; library sizes are unchanged).

```python
@pytest.mark.parametrize(
    "method",
    [
        "paga_dpt",
        pytest.param("palantir", marks=pytest.mark.slow),
        pytest.param("scfates_ppt", marks=pytest.mark.slow),
    ],
)
def test_clone_columns_and_receptor_genes_change_no_result(
    method: str, plain: Dataset, clonal: Dataset
) -> None:
    a, b = (_run(ds, method, clonal_columns=("clone_id",)) for ds in (plain, clonal))
    _assert_identical(
        a, b
    )  # every TrajectoryResult field except runtime_seconds; arrays NaN-aware, value for value


def test_the_test_has_teeth(plain: Dataset, clonal: Dataset) -> None:
    off = PrepParams(exclude_receptor_genes=False)
    assert not np.array_equal(
        run_pipeline(plain, off).artifacts["pca"], run_pipeline(clonal, off).artifacts["pca"]
    )
```

  `_run` is `run_pipeline` then `InferenceContext.from_prep(..., clonal_columns=…)` then `run_method(ctx, method, ClusterRoot("true_label", "trunk_early"))`.
  Also in this file:
  - `test_receptor_columns_never_reach_methods`: parametrize over `clone_id, clone_id_size, clonotype, raw_clonotype_id, IR_VJ_1_junction_aa, airr:clone_id, receptor_type, chain_pairing, cdr3_nt` → absent from `ctx.obs` and present in `ctx.scoring_only`; `cell_type, donor_age, tissue` stay.
  - `test_named_columns_are_scoring_only`: `from_prep(..., clonal_columns=("donor",))` drops `donor`.
  - `test_scoring_only_cluster_key_and_label_key_are_refused`: `run_method(..., params={"cluster_key": "clone_id"})` and `run_method(..., label_key="clone_id")` raise `ScoringOnlyColumnError`.
  - `test_runs_record_what_steered_them`: `result.obs_keys == sorted({"true_label", "lumencell_leiden"})` for the `cluster:true_label=…` root with default params and detected label key.
- [ ] **Step 2: Write the other failing tests.** `tests/unit/test_roots.py`: `test_a_scoring_only_column_is_refused_as_a_root` — `resolve_root(ctx, ClusterRoot("clone_id", "A"))` raises `RootError` matching `"only to score"`. `tests/unit/test_store.py`: `test_result_json_records_isolation` — `store_version == 2`, `prep_key` as passed, `inference_isolation == 1`, `obs_keys` round-trips through `load`; a run with `store_version = 1` is not listed. `tests/api/test_routes.py`: `test_create_run_refuses_a_scoring_only_root_or_cluster_key` — the dataset gains a `clone_id` column; `POST /trajectories` with root `cluster:clone_id=…`, and with `params={"cluster_key": "clone_id"}`, answer 422 with `error.code == "scoring_only_column"` and `detail.columns == ["clone_id"]`.
- [ ] **Step 3: Run** `uv run pytest tests/unit/test_isolation.py tests/unit/test_roots.py tests/unit/test_store.py tests/api/test_routes.py -q`. Expected: FAIL.
- [ ] **Step 4: Implement** as in Interfaces.
- [ ] **Step 5: Run** the same command, then `uv run pytest tests/unit/test_isolation.py -q` (includes the slow cases). Expected: PASS.
- [ ] **Step 6: C1–C4.** Commit message: `Make clone and receptor columns scoring-only and record what steered each run (store version 2)`.

---

### Task 6: Clonal scoring in the synthetic benchmark

**Branch:** `pivot/NN-bench-clonal` · **Model:** sonnet — glue over specified functions.

**Files:**
- Modify: `src/lumencell/bench/runner.py`, `src/lumencell/bench/report.py`, `tests/golden/test_synthetic.py`
- Test: `tests/unit/test_runner_clonal.py` (new)

**Interfaces:**
- Consumes: `ClonalRow`/`CaseResult.clonal` (Task 1); `encode_clones`, `score_clonal`, `baseline_space`, `kmeans_partition` (Task 2); `DEFAULT_LAYER`, `add_clones`, `clone_truth` (Task 3); `from_prep(..., clonal_columns=…)` (Task 5); `ORACLE_ETA_MIN`, `ORACLE_P_MAX`, `CLONAL_SEED`, `SENSITIVITY_CLONE_SIZE`.
- Produces, `bench/runner.py` (phase 3 reuses `clonal_row` for real data):
  - `run_synthetic_case(..., clonal_layer: ClonalLayer | None = DEFAULT_LAYER) -> CaseResult` — the layer is applied after `generate` and before writing, only when the truth has at least 2 lineages; inference gets `clonal_columns=("clone_id", "donor")` when it was applied; `CaseResult.clonal` is `ClonalRow("not_applicable")` otherwise.
  - `clonal_row(truth: Truth, result: TrajectoryResult, pca: npt.NDArray[np.floating], neighbors_pcs: int) -> ClonalRow`, in order: truth with fewer than 2 lineages or no `clone_ids` → `not_applicable`; the oracle = `score_clonal(truth.membership(), codes, truth.order, kmeans k = truth L, oracle_cells=truth.exclusive() >= 0)`; method L < 2 → `single_lineage` (with `eta_oracle`, `n_clones`); the method = the same call with `result.membership`, `result.pseudotime` and k = method L; method status `too_few_clones` → that status; oracle `eta < ORACLE_ETA_MIN` or oracle `p > ORACLE_P_MAX` → `not_assessable` (with `eta_oracle`, `eta_method`); else `ok` with `excess_ratio = method.eta / oracle.eta`, `p = method.p`, `delta_base = method.delta_eta`. Evidence counts: `n_cells`, `n_clones` and `effective_clones` come from the oracle call (the same cells the method is scored on) for every status after `not_applicable`; `eta_m3` and `p_m3` are the method's η and headline p from the same call at `min_clone_size=SENSITIVITY_CLONE_SIZE`, filled whenever the method is scored (`too_few_clones`, `not_assessable`, `ok`), NaN otherwise. Any exception inside → `ClonalRow("failed", reason=f"{type(exc).__name__}: {exc}")`, also printed to stderr.
  - `run_real_case` sets `ClonalRow("not_applicable")` (no clone truth on the developmental datasets).
- Produces, `bench/report.py`: the synthetic table gains a "Clonal excess ratio" column (`_range(rows, "clonal_excess_ratio")`) and, under it, the sentence: "The clonal gate is not a branch gate: a method that merges two true lineages can score a clonal ratio above 1, and lineage recall catches that."

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_runner_clonal.py`, 1500-cell cases):
  - `test_branching_cases_carry_a_clonal_score`: `run_synthetic_case(tmp, "bifurcation", n_cells=1500, n_genes=1000).row()` has `clonal_status` in `ClonalRowStatus`, `clonal_eta_oracle >= 0.05`, `clonal_n_clones >= 20`, `clonal_n_cells >= 2 * clonal_n_clones`, `clonal_effective_clones > 0`, `math.isfinite(clonal_eta_m3)` and `0 < clonal_p_m3 <= 1`.
  - `test_linear_cases_are_not_applicable`: status `not_applicable`, ratio NaN, `clonal_n_cells == 0`, `clonal_effective_clones` and `clonal_eta_m3` NaN, and `check_gates` raises no clonal failure.
  - `test_existing_metrics_do_not_move_with_the_clonal_layer`: for bifurcation seed 0, the row with `clonal_layer=DEFAULT_LAYER` equals the row with `clonal_layer=None` on every key except `runtime_seconds` and `clonal_*`.
  - `test_a_single_lineage_method_fails_the_gate`: `clonal_row` with a one-lineage result → status `single_lineage`.
  - `test_an_oracle_below_its_floor_is_not_assessable`: a truth built with `fate_bias=0` → status `not_assessable`.
  - `test_an_error_is_recorded_as_failed`: monkeypatch `score_clonal` to raise → status `failed`, `reason` names the error.
  - `test_the_report_says_the_clonal_gate_is_not_a_branch_gate`: `render()` contains "not a branch gate".
  In `tests/golden/test_synthetic.py`, `test_v1_trap_bifurcation_holds_its_baseline` also asserts `row["clonal_status"] != "not_measured"`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_runner_clonal.py tests/golden/test_synthetic.py -m "not slow" -q`. Expected: FAIL.
- [ ] **Step 3: Implement** as in Interfaces. The PCs are `prep.artifacts["pca"]` and `pp.neighbors_pcs`.
- [ ] **Step 4: Run** the same command. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Score clonal agreement in oracle-cells mode on every branching synthetic case`.

---

### Task 7: Re-measure the synthetic baseline, on purpose

**Branch:** `pivot/NN-bench-baseline` · **Model:** sonnet — runs plus a careful comparison.

**Files:**
- Modify: `Makefile`, `benchmarks/baseline.json`, `README.md`, possibly `experimental_note` in `src/lumencell/engine/methods/*.py`
- Test: `tests/unit/test_makefile.py` (new)

**Interfaces:**
- Consumes: the complete bench path (Tasks 1–6 merged).
- Produces: `Makefile` `BENCH_METHODS ?= paga_dpt,palantir,scfates_ppt`, passed as `--methods $(BENCH_METHODS)` by `bench-synth` and `bench-baseline` (today both default to `paga_dpt` only, and a re-measure would drop the other two methods from `baseline.json`). A `baseline.json` whose rows carry the `clonal_*` fields.

- [ ] **Step 1: Write the failing test** `tests/unit/test_makefile.py::test_bench_targets_cover_every_method`: `subprocess.run(["make", "-n", target], capture_output=True, text=True).stdout` contains `--methods paga_dpt,palantir,scfates_ppt` for `bench-synth` and `bench-baseline`.
- [ ] **Step 2: Run** it. Expected: FAIL. **Step 3:** add `BENCH_METHODS`. **Step 4:** run it. Expected: PASS.
- [ ] **Step 5 (background):** `make bench-baseline > "$TMPDIR/baseline.log" 2>&1; echo $?` → `0`.
- [ ] **Step 6: Check that no existing metric moved** (it is an error if one did: the clone layer or the exclusion would be leaking into synthetic inference):

```bash
git show 142a788:benchmarks/baseline.json > "$TMPDIR/old-baseline.json"
uv run python - "$TMPDIR/old-baseline.json" <<'EOF'
import json, math, sys
old = json.load(open(sys.argv[1]))["methods"]; new = json.load(open("benchmarks/baseline.json"))["methods"]
bad = [(m, c, k, v, new[m][c].get(k)) for m, cases in old.items() for c, row in cases.items()
       for k, v in row.items() if k != "runtime_seconds"
       and not (v == new[m][c].get(k) or (isinstance(v, float) and math.isnan(v) and math.isnan(new[m][c][k])))]
print(len(bad), bad[:10])
EOF
```

  Expected: `0 []`. If not, stop, do not commit, and report the differing rows.
- [ ] **Step 7 (background):** `uv run pytest -m slow tests/golden/test_synthetic.py -q; echo $?` → `0` (the `experimental` flags are still honest; all three methods already fail a gate). If a method now also fails a clonal gate, update its `experimental_note` to name it.
- [ ] **Step 8:** `uv run lumencell bench report --readme`; review the clonal columns with `data:validate-data` (ranges plausible, no NaN where a status is `ok`).
- [ ] **Step 9: C1–C4.** Commit message: `Re-measure the synthetic baseline with the clonal layer (existing metrics unchanged)`.

---

### Task 8: Calibrate the nulls (slow, run once and recorded)

**Branch:** `pivot/NN-clonal-calibration` · **Model:** opus — statistical design and interpretation.

**Files:**
- Create: `src/lumencell/bench/calibration.py`, `tests/clonal/test_calibration.py`, `benchmarks/results/clonal-calibration-<date>.json`
- Modify: `benchmarks/thresholds.toml` (changelog line only)

**Interfaces:**
- Consumes: `[clonal.calibration]` and `[synthetic]` (Task 1); `score_clonal`, `encode_clones` (Task 2); `add_clones`, `ClonalLayer`, `clone_truth` (Task 3); `from_prep(..., clonal_columns=…)` (Task 5).
- Produces, `bench/calibration.py`:
  - `@dataclass(frozen=True) class RejectionCount: name: str; null: str; n: int; rejections: int; bound: int; expect_calibrated: bool; p_values: tuple[float, ...]` with property `passed` (`rejections <= bound` when calibrated is expected, `> bound` otherwise).
  - `@dataclass(frozen=True) class MonotonicCase: topology: str; seed: int; fate_bias: tuple[float, ...]; eta: tuple[float, ...]; p_at_zero: float`
  - `@dataclass(frozen=True) class CalibrationReport: checks: tuple[RejectionCount, ...]; monotonic: tuple[MonotonicCase, ...]; config: dict[str, Any]; commit: str` with property `passed`.
  - `run_calibration(cfg: dict[str, Any], methods: Sequence[str], workdir: Path) -> CalibrationReport`; `write_report(report: CalibrationReport, path: Path) -> None` (JSON, `allow_nan=False`, `"passed"` at the top level).
- The checks (bound = `int(binom.ppf(0.99, 200, 0.05)) == 18`; a rejection is p ≤ 0.05):
  - `null/<method>`: fate_bias 0, own-cells mode on each method's own masks and pseudotime, headline p; expect calibrated. Prep runs once per seed, shared by the methods.
  - `donor_confound/N1`: the `donor_confound` control, oracle (truth lineages, true order, oracle cells), the N1 p; expect calibrated.
  - `donor_confound/N1_donors_ignored`: the same data with `codes.donor` replaced by zeros (clone ids still keyed by donor), the N1 p; expect *not* calibrated.
  - `depth_only/N2` on `depth_only_topology`, oracle, the N2 p; expect calibrated. `depth_only/N1`, same data, the N1 p; expect *not* calibrated (the control has teeth).
  - Monotonicity: for each monotonic topology and seed, oracle η strictly increases over `monotonic_fate_bias`, and the oracle's headline p at fate_bias 0 is above 0.01.
  - Controls and monotonicity use a constant baseline (their baseline fields are not examined) and need no prep.

- [ ] **Step 1: Write the tests** (`tests/clonal/test_calibration.py`, every test `@pytest.mark.slow`). A module fixture runs `run_calibration` once and, when `LUMENCELL_CALIBRATION_OUT` is set, writes the report there *before* any assertion, so an honest failure is still recorded. Tests: `test_binomial_bound_is_18`; `test_the_null_is_calibrated_on_each_methods_own_masks[method]`; `test_donor_confound_is_calibrated_under_n1`; `test_donor_confound_with_donors_ignored_is_not`; `test_depth_only_is_calibrated_under_n2`; `test_depth_only_is_not_calibrated_under_n1`; `test_oracle_eta_rises_with_fate_bias`; `test_oracle_p_at_no_bias_is_above_0_01`. Each asserts the matching `RejectionCount.passed` or monotonic condition and prints the counts.
- [ ] **Step 2: Run** `uv run pytest -m slow tests/clonal/test_calibration.py -q -k bound`. Expected: FAIL (`lumencell.bench.calibration` missing).
- [ ] **Step 3: Implement** `bench/calibration.py`.
- [ ] **Step 3a: Smoke-run** on 3 seeds and one monotonic case, without recording: a short `uv run python -` script that loads the thresholds, sets `cfg["clonal"]["calibration"]["n_seeds"] = 3`, `monotonic_seeds = [0]` and `monotonic_topologies = ["bifurcation"]`, calls `run_calibration(cfg, ["paga_dpt", "palantir", "scfates_ppt"], workdir)` and prints each check's `n` and `rejections`. Expected: every check has `n == 3`; no exception. The full 200-seed counts come only from Step 4.
- [ ] **Step 4 (background, hours):** `LUMENCELL_CALIBRATION_OUT=benchmarks/results/clonal-calibration-$(date +%F).json uv run pytest -m slow tests/clonal/test_calibration.py -q > "$TMPDIR/calibration.log" 2>&1; echo $?`.
- [ ] **Step 5:** Run the `data:statistical-analysis` and `data:validate-data` skills on the recorded JSON: is the rejection rate plausible for discrete, conservative permutation p-values; do the per-seed p distributions look uniform or conservative under the null; are the controls' teeth real. Record their conclusions in HANDOFF.
- [ ] **Step 6:** Add a dated changelog line to `thresholds.toml`: "Calibration run, results in benchmarks/results/clonal-calibration-<date>.json: <one line per check>. No threshold or definition changed." If any check failed, do not change seeds, definitions or thresholds; record the failure under "Open questions for Shreyan" in HANDOFF and stop for the owner after opening the PR.
- [ ] **Step 7: C1–C4.** Commit message: `Calibrate the clonal nulls on 200 pre-declared seeds and record the result`.

---

### Task 9: Re-run the real benchmarks after the TR-gene exclusion

**Branch:** `pivot/NN-bench-real-rerun` · **Model:** sonnet — a run and a logged comparison.

**Files:**
- Create: `benchmarks/results/real-<date>.json`
- Modify: `benchmarks/thresholds.toml` (changelog), `README.md`

**Interfaces:**
- Consumes: Tasks 4, 5 and 6 merged (receptor exclusion; real rows carry `clonal_status = "not_applicable"`).
- Produces: a results file and a changelog entry that phase 3 reads as the pre-T-cell real-data reference.

- [ ] **Step 1: Check the datasets are cached:** `uv run python -c "from lumencell.bench.datasets import DATASETS, is_cached; print({k: is_cached(v) for k, v in DATASETS.items()})"`. If any is `False`, **stop and ask the owner** before `uv run lumencell bench fetch all` (about 708 MB into `~/Library/Caches/lumencell/datasets`).
- [ ] **Step 2 (background):** `make bench-real > "$TMPDIR/real.log" 2>&1; echo $?` → `0` (gate failures do not fail the command; read them in the log).
- [ ] **Step 3: Compare** the new file with `benchmarks/results/real-2026-11-10.json`: for each method and dataset, the five gated metrics and `n_clusters`, old and new, side by side. Record each dataset's `hvg_info` from its prep cache (`PrepCache(load(local_path(ds))).run(PrepParams()).artifacts["hvg_info"]`). Check the comparison with `data:validate-data`.
- [ ] **Step 4:** Add a dated changelog entry to `thresholds.toml` titled "TR-gene exclusion re-run": per dataset, the receptor genes matched and how many would have been HVGs, and every metric that moved (or "no metric moved"). Thresholds unchanged.
- [ ] **Step 5:** `uv run lumencell bench report --readme`.
- [ ] **Step 6: C1–C4.** Commit message: `Re-run the real benchmarks after excluding TCR genes from HVG selection, and log the change`.

---

### Task 10: Clonal settings for a dataset

**Branch:** `pivot/NN-clonal-settings` · **Model:** sonnet — API work in the existing pattern.

**Files:**
- Modify: `src/lumencell/api/models.py`, `src/lumencell/api/dataset_view.py`, `src/lumencell/api/session.py`, `src/lumencell/api/routes.py`, `contracts/openapi.json`, `web/src/api/generated/*`, `web/src/test/fixtures.ts`, `tests/api/test_contract.py`
- Create: `src/lumencell/trust/clonal.py`, `tests/api/test_clonal_settings.py`

**Interfaces:**
- Consumes: `encode_clones` (Task 2); `scoring_only_columns`, `inference_obs`, `DatasetView.clonal_columns` (Task 5); `MAX_DONORS`.
- Produces, `api/models.py`:
  - `ClonalSettings(Model): clone_key: str = Field(min_length=1); donor_key: str | None = None; single_donor: bool = False; strata_keys: list[str] = []` with a `model_validator(mode="after")`: exactly one of `donor_key` or `single_donor = True`.
  - `ClonalSummary(Model): settings: ClonalSettings; n_cells_labelled: int; n_clones: int; n_clones_2plus: int; n_donors: int; size_median: float; size_p90: float; size_p99: float; size_max: int` — over every labelled cell with a donor, singletons included.
  - `DatasetInfo` gains `clonal: ClonalSettings | None` and `clonal_candidates: list[str]` (no defaults).
- Produces, `trust/clonal.py`:
  - `CANDIDATE = re.compile(r"^(?:airr:)?(?:clone_id|raw_clonotype_id|clonotype)$", re.IGNORECASE)`; `clonal_candidates(columns: Iterable[object]) -> list[str]` in obs order. Suggestions only, never applied.
  - `settings_columns(settings: ClonalSettings) -> tuple[str, ...]` — the clone and donor (if any) columns, which become scoring-only. Never the strata keys: strata columns such as a time point stay visible to methods and may define a root.
  - `check_settings(obs: pd.DataFrame, settings: ClonalSettings) -> list[str]` — problems, empty when valid: each named column exists; the clone and donor columns differ; strata keys are distinct and differ from both; the clone column is categorical, pandas string, or object with `infer_dtype` in `{"string", "empty"}`; the donor column has 1 to `MAX_DONORS` distinct non-missing values; at least one cell has a clone label.
  - `settings_summary(obs: pd.DataFrame, settings: ClonalSettings) -> ClonalSummary`.
  - `settings_from_state(state: dict[str, Any]) -> ClonalSettings | None` — `state["clonal"]` validated; an invalid value gives `None` and one stderr line.
- Produces, `api/dataset_view.py`: `column_kind(s: pd.Series, categories: list[str]) -> ObsKind` (numeric → `continuous`; boolean → `boolean`; at most `MAX_CATEGORIES_LISTED` categories → `categorical`; otherwise `string`), used by both `_file_columns` and `obs_payload` (whose meta `kind` is now `column_kind`; codes and categories are still sent for `string`). `DatasetView.set_clonal(columns: Sequence[str]) -> None` — under the lock: stores `clonal_columns`, drops the cached context, re-detects `label_key`.
- Produces, `api/session.py`: `Session.clonal_settings() -> ClonalSettings | None`; `Session.set_clonal(settings: ClonalSettings) -> ClonalSummary` (422 `invalid_clonal_settings` with `detail = {"problems": [...]}`; 422 `scoring_only_column` when the saved root is a `cluster:` root on the clone or donor column; else save to state through `StateStore.update`, then `view.set_clonal`); `Session.clear_clonal() -> None` (removes the setting from state and calls `view.set_clonal(())`). `__init__` applies saved settings to the view. `dataset_info()` fills the two new fields.
- Produces, routes: `PUT /api/v1/dataset/clonal` (body `ClonalSettings`, 200 `ClonalSummary`) and `DELETE /api/v1/dataset/clonal` (204). The `/obs/{column}` description says `kind` is `string` above 1,000 categories.

- [ ] **Step 1: Write the failing tests** (`tests/api/test_clonal_settings.py`; a module fixture writes a 1000-cell bifurcation with `add_clones(DEFAULT_LAYER)`, a 1,001-value string column `barcode2`, and a 1,001-value column `many_donors`, then prepares it into a temporary cache as `test_cache_healing.py` does; the job manager is not started):
  - `test_put_returns_the_summary`: `PUT {"clone_key": "clone_id", "donor_key": "donor"}` → 200; `n_donors == 3`; `n_cells_labelled` equals the count of non-missing `clone_id`; `n_clones` equals the distinct (donor, label) pairs; `size_max >= size_p99 >= size_p90 >= size_median >= 1`.
  - `test_donor_or_single_donor_exactly_one`: both and neither → 422 `invalid_request`.
  - `test_invalid_settings_answer_422_with_the_error_body` (parametrized): unknown column, clone = donor, numeric clone column (`true_pseudotime`), `many_donors` as donor, unknown strata key, strata key equal to the clone key → 422 `invalid_clonal_settings`, `detail.problems` non-empty.
  - `test_dataset_info_lists_candidates_and_never_selects_one`: before any PUT, `clonal is None` and `clonal_candidates == ["clone_id"]`; after PUT, `clonal` equals the settings; after DELETE (204), `clonal is None`.
  - `test_settings_survive_a_new_session`.
  - `test_put_rebuilds_the_inference_context`: after PUT, `session.view.ctx().scoring_only >= {"clone_id", "donor"}` and `"donor" not in ctx.obs`.
  - `test_settings_that_capture_the_saved_root_are_refused` (Review Focus 3): `PUT /root` `cluster:true_segment=trunk`, then PUT clonal with `{"clone_key": "clone_id", "donor_key": "true_segment"}` → 422 `scoring_only_column`; the settings are not saved.
  - `test_strata_columns_stay_visible_to_methods`: PUT clonal with `{"clone_key": "clone_id", "donor_key": "donor", "strata_keys": ["true_segment"]}` → 200; `"true_segment" in session.view.ctx().obs` and `"true_segment" not in ctx.scoring_only`; `PUT /root` `cluster:true_segment=trunk` → 200, and `POST /trajectories` from that root is not refused.
  - `test_create_run_refuses_a_settings_column`: with settings saved, `POST /trajectories` rooted on `cluster:donor=d0` → 422 `scoring_only_column`.
  - `test_listing_and_serving_agree_on_every_column_kind`: for every obs column, the listed `kind` equals the `kind` in `GET /obs/{column}` meta; `barcode2` is `string` in both.
  In `tests/api/test_contract.py`, add the two routes to the required set.
- [ ] **Step 2: Run** `uv run pytest tests/api/test_clonal_settings.py tests/api/test_contract.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** models first, then `make contract`; then `trust/clonal.py`, the view, the session and the routes. Update the `dataset()` fixture in `web/src/test/fixtures.ts` with `clonal: null, clonal_candidates: []` and nothing else.
- [ ] **Step 4: Run** the same tests. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add dataset-level clonal settings, and list and serve obs columns by one rule`.

---

### Task 11: Trust files in the run store

**Branch:** `pivot/NN-run-trust-store` · **Model:** sonnet — small storage change with clear atomicity rules.

**Files:**
- Modify: `src/lumencell/engine/store.py`, `src/lumencell/prep/housekeeping.py`
- Test: `tests/unit/test_store.py`, `tests/unit/test_cache_admin.py`

**Interfaces:**
- Produces, `engine/store.py` (phase 2 stores `stability-<h>.json` the same way):
  - `TRUST_NAME = re.compile(r"[a-z][a-z_]{0,31}-[0-9a-f]{16}")`
  - `RunStore.trust_path(rid: str, name: str) -> Path` — `<run>/trust/<name>.json`; `ValueError` unless `name` fully matches `TRUST_NAME`.
  - `RunStore.write_trust(rid: str, name: str, payload: dict[str, Any]) -> Path` — `FileNotFoundError` unless the run is complete; writes `json.dumps(payload, indent=2, allow_nan=False)` to `tempfile.mkstemp(prefix=f".{name}-", suffix=".tmp", dir=<run>/trust)` and `os.replace`s it; the temporary file is removed on any error (a NaN raises `ValueError` and leaves nothing).
  - `RunStore.read_trust(rid: str, name: str) -> dict[str, Any] | None` — `None` when missing; an unreadable or non-object file also gives `None`, with one stderr line naming the file.
  - `RunStore.save` of an existing complete run moves its `trust/` folder into the new run folder before the swap.
- Produces, `prep/housekeeping.py`: `_TRUST_TEMP = re.compile(rf"\.[a-z][a-z_]{{0,31}}-[0-9a-f]{{16}}-{_RANDOM}\.tmp")`; `_leftovers` also lists matching files in `runs/<rid>/trust/` (not a symlink), older than `TEMP_GRACE`.

- [ ] **Step 1: Write the failing tests.** `tests/unit/test_store.py`: `test_trust_round_trips_atomically` (no `.tmp` left in `trust/`); `test_trust_names_are_checked` (`"../x"`, `"clonal-xyz"`, `"Clonal-0123456789abcdef"` raise `ValueError`); `test_trust_refuses_nan_and_leaves_nothing`; `test_trust_needs_a_complete_run`; `test_saving_a_run_again_keeps_its_trust`; `test_a_damaged_trust_file_reads_as_missing` (stderr names it). `tests/unit/test_cache_admin.py`: `test_prune_removes_old_trust_temp_files_only` — an old `.clonal-<16 hex>-abcd1234.tmp` is removed, a fresh one and `clonal-<16 hex>.json` stay.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_store.py tests/unit/test_cache_admin.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement.** **Step 4: Run** again. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Store per-run trust results atomically, keep them when a run is saved again, and prune their leftovers`.

---

### Task 12: Two job priorities and done-callbacks

**Branch:** `pivot/NN-job-priorities` · **Model:** opus — concurrency in the job manager.

**Files:**
- Modify: `src/lumencell/jobs/manager.py`
- Test: `tests/unit/test_jobs.py`

**Interfaces:**
- Produces, `jobs/manager.py` (phase 2 relies on all of it):
  - `BACKGROUND_KINDS: frozenset[str] = frozenset({"score_clonal", "stability"})`. Every other kind is foreground.
  - `DEFAULT_TIME_LIMITS` gains `"score_clonal": 3600`.
  - The single `_pending` deque becomes `_foreground` and `_background`. `submit` appends by kind. `_take_next(self) -> str | None` (call with the lock held) pops the oldest foreground job, else the oldest background job. Running jobs are never pre-empted. `cancel` removes a queued job from whichever deque holds it.
  - `add_done_callback(self, callback: Callable[[Job], None]) -> None`. Each callback runs once per job that reaches `done` (not `failed` or `cancelled`), on the listener thread, **after** the manager lock is released. An exception in a callback is printed to stderr with its traceback (`lumencell: a callback after job <id> failed`) and the listener carries on.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_jobs.py`):
  - `test_foreground_jobs_overtake_queued_background_jobs`: started manager; submit `sleep` `"block"` (5 steps × 0.1 s); while it runs submit `score_clonal` `"bg"` with payload `{}` and then `sleep` `"fg"`; wait for both (the `score_clonal` job fails, having no handler yet); assert `fg.started < bg.started`.
  - `test_take_next_order` (unstarted manager): submit bg1, fg1, bg2, fg2 → `_take_next()` yields fg1, fg2, bg1, bg2, then `None`.
  - `test_a_queued_background_job_can_be_cancelled`.
  - `test_done_callbacks_run_outside_the_lock`: the callback starts a thread that calls `manager.get(job_id)` and joins it with a 2 s timeout; assert the thread finished and the callback saw the job `done`.
  - `test_callbacks_see_only_done_jobs`: a `sleep` with `fail=True` triggers no callback.
  - `test_a_failing_callback_is_reported_and_the_listener_keeps_going`: stderr contains `callback after job`; the next job completes.
  - `test_scoring_has_its_own_time_limit`: `JobManager().time_limit("score_clonal") == 3600`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_jobs.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement.** **Step 4: Run** again. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Run prepare and infer ahead of background jobs, and let callers react when a job finishes`.

---

### Task 13: The `score_clonal` job

**Branch:** `pivot/NN-score-clonal-job` · **Model:** opus — joins scores, settings, refusal and storage; defines the trust models.

**Files:**
- Modify: `src/lumencell/api/models.py`, `src/lumencell/trust/clonal.py`, `src/lumencell/jobs/worker.py`
- Test: `tests/unit/test_clonal_agreement.py` (new), `tests/unit/test_jobs.py`

**Interfaces:**
- Consumes: `score_clonal`, `ClonalResult`, `encode_clones`, `baseline_space`, `kmeans_partition` (Task 2); `scoring_only_columns` (Task 5); `ClonalSettings`, `settings_columns` (Task 10); `RunStore.write_trust`, `TRUST_NAME` (Task 11); `STORE_VERSION`, `graph["obs_keys"]`, `graph["prep_key"]` (Task 5).
- Produces, `api/models.py` (defined here, referenced by a route only in Task 14, so `make contract` shows no change yet):
  - `NullComparison(Model): null: Literal["N1", "N2", "N3"]; mean_entropy_mean: float | None; mean_entropy_sd: float | None; eta: float | None; p: float | None`
  - `ClonalScore(Model): mode: Literal["own_cells", "oracle_cells"]; n_lineages: int; min_clone_size: int; n_cells: int; n_clones: int; effective_clones: float | None; n_donors: int; n_strata: int | None; n_permutations: int; seed: int; purity: float | None; mean_entropy: float | None; nulls: list[NullComparison]; eta: float | None; p: float | None; eta_base: float | None; delta_eta: float | None; p_delta: float | None; clone_weighted_purity: float | None; purity_without_largest: float | None` (`eta` is against N1; `p` is the headline, the largest over the nulls).
  - `LineageCoupling(Model): a: int; b: int; observed: float | None; null_mean: float | None; z: float | None; p: float | None; p_adj: float | None`
  - `ClonalStatus = Literal["ok", "not_configured", "single_lineage", "too_few_clones", "refused", "pending", "failed"]`
  - `ClonalAgreement(Model): status: ClonalStatus; settings: ClonalSettings | None; primary: ClonalScore | None; sensitivity_m3: ClonalScore | None; coupling: list[LineageCoupling]; notes: list[str]`
  - `RunTrust(Model): clonal: ClonalAgreement | None` (phase 2 adds `stability`, `tree_support`, `cell_cycle`).
- Produces, `trust/clonal.py`:
  - `COUPLING_NOTE = "Coupling is the cosine similarity of raw scored-cell counts per clone, not CoSpar's normalised default. A high z says two lineages share more clones than chance; it does not say which came first."`
  - `trust_hash(settings: ClonalSettings) -> str` — the first 16 hex characters of SHA-256 over `json.dumps({"settings": <model_dump(mode="json") with strata_keys sorted>, "metric_version": CLONAL_METRIC_VERSION, "params": {"min_clone_size", "sensitivity_clone_size", "min_clones", "n_permutations", "seed", "pseudotime_bins"}}, sort_keys=True)`. `trust_name(settings) -> str` is `f"clonal-{trust_hash(settings)}"`. Column names never appear in a path.
  - `to_api_score(r: ClonalResult) -> ClonalScore` and `to_api_coupling(c: Coupling) -> LineageCoupling` map every NaN to `None`.
  - `not_scored(status: ClonalStatus, settings: ClonalSettings | None, note: str) -> ClonalAgreement`.
  - `agreement_for_run(result: TrajectoryResult, obs: pd.DataFrame, settings: ClonalSettings, pca: npt.NDArray[np.floating], neighbors_pcs: int, obs_keys: Sequence[str]) -> ClonalAgreement` — own-cells mode. `refused` (note names the column) when `obs_keys` meets `scoring_only_columns(obs.columns, settings_columns(settings))`; `single_lineage` when L < 2; else `primary` at m = 2 and `sensitivity_m3` at m = 3, baseline k-means with k = L on `baseline_space(pca, neighbors_pcs)` and seed `CLONAL_SEED`, `status = primary.status`, coupling from `primary`. Notes always hold `COUPLING_NOTE`; add `"<n> cells with a receptor have no donor and were not scored."` when `n_missing_donor > 0`, and `"Only <n> clones were scored; at least 20 are needed before the score counts."` for `too_few_clones`.
- Produces, `jobs/worker.py`: `HANDLERS["score_clonal"] = handle_score_clonal`. Payload: the session's `path`, `x_is`, `prep_params`, `cache_root`, plus `run_id`, `settings` (`ClonalSettings.model_dump(mode="json")`) and `trust_name`. The handler: loads the dataset and the run graph; `store_version != STORE_VERSION` → `not_scored("refused", …, "made before inference isolation; run it again")`; `graph["prep_key"] != prep_key(ds, params)` → `not_scored("failed", …, "the run was prepared with other settings")`; otherwise loads the run and the prep (`PrepCache(...).run(params)`) and calls `agreement_for_run` with `graph["obs_keys"]`. It writes `agreement.model_dump(mode="json")` with `write_trust(rid, payload["trust_name"], …)` and returns `rid`.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_clonal_agreement.py`; a module fixture: 1500-cell bifurcation with `add_clones(DEFAULT_LAYER, seed=0)`, `run_pipeline`, `paga_dpt` from `cluster:true_label=trunk_early` with `clonal_columns=("clone_id", "donor")`):
  - `test_a_run_with_clones_scores_ok`: `status == "ok"`; `primary.mode == "own_cells"`, `primary.min_clone_size == 2`, `sensitivity_m3.min_clone_size == 3`; `primary.n_clones >= 20`; `[n.null for n in primary.nulls] == ["N1", "N2"]`; `len(coupling) == L * (L - 1) // 2`; `COUPLING_NOTE in notes`.
  - `test_strata_add_the_third_null`: settings with `strata_keys=["true_segment"]` give nulls `["N1", "N2", "N3"]` and `primary.n_strata >= 3`.
  - `test_a_run_steered_by_a_clone_column_is_refused`: `obs_keys=["clone_id"]` → `status == "refused"`, `primary is None`, a note names `clone_id`.
  - `test_one_lineage_is_single_lineage`.
  - `test_missing_donors_are_noted` (Review Focus 2): blank the donor of 10 labelled cells → a note starts with `10 cells with a receptor have no donor`.
  - `test_undefined_values_are_null_and_the_file_is_strict_json` (Review Focus 1): a run whose clone column yields 0 scored clones, and one with 1 → `status == "too_few_clones"`; `RunStore.write_trust` succeeds and `json.loads(text, parse_constant=_raise)` parses the file; `primary.purity_without_largest is None` for the one-clone case.
  - `test_the_trust_name_hides_column_names`: settings with `clone_key="../../etc"` (an obs column of that name) → `trust_name` fully matches `TRUST_NAME`; the hash is the same for strata `["a","b"]` and `["b","a"]` and differs when `donor_key` changes.
  In `tests/unit/test_jobs.py`, `test_score_clonal_job_writes_the_trust_file`: prepare and infer through the worker on a file with `add_clones(DEFAULT_LAYER)` applied (as `test_prepare_then_infer_through_the_worker` does), then submit `score_clonal` with the payload above → `done`, and `RunStore(...).read_trust(rid, name)["status"]` is `ok` or `too_few_clones`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_clonal_agreement.py tests/unit/test_jobs.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** models (then `make contract`; expected: no diff), `trust/clonal.py`, the handler.
- [ ] **Step 4: Run** the same tests. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add the score_clonal worker job and the clonal agreement models`.

---

### Task 14: Clonal agreement on every run, with healing

**Branch:** `pivot/NN-run-trust-api` · **Model:** opus — the GET-healing and settings-change paths are concurrency-sensitive.

**Files:**
- Modify: `src/lumencell/api/models.py`, `src/lumencell/api/runs.py`, `src/lumencell/api/session.py`, `contracts/openapi.json`, `web/src/api/generated/*`, `web/src/test/fixtures.ts`
- Create: `tests/api/test_clonal_trust.py`

**Interfaces:**
- Consumes: `add_done_callback`, `ACTIVE` (Task 12); trust models, `trust_hash`, `trust_name`, `not_scored` (Task 13); `RunStore.read_trust` (Task 11); `settings_from_state` (Task 10).
- Produces, `api/models.py`: `TrajectoryRun.trust: RunTrust | None` (no default; always filled for a done run); `JobKind = Literal["prepare", "infer", "heatmap", "score_clonal"]`; `JobInfo.kind: JobKind | str` (a kind the contract does not list is reported by its own name).
- Produces, `api/session.py`: `job_info` passes `job.kind` through unchanged; `__init__` calls `jobs.add_done_callback(self._on_job_done)`; `_on_job_done(job)` acts only when `job.kind == "infer"`, `job.payload.get("path") == str(self.path)` and `job.result_id`, calling `self.runs.submit_trust_jobs(job.result_id)`; `set_clonal` ends with `self.runs.submit_clonal_all()`.
- Produces, `api/runs.py` (phase 2 reuses `_trust_lookup` and extends `submit_trust_jobs`):
  - `@dataclass(frozen=True) class TrustLookup: stored: dict[str, Any] | None; job: Job | None`
  - `Runs._trust_lookup(self, rid: str, kind: str, name: str, h: str, payload: Callable[[], dict[str, Any]]) -> TrustLookup` — the stored file if `read_trust(rid, name)` finds it; else the job `jobs.find(f"{kind}:{rid}:{h}")`; when there is none, or it is `done` but its file is gone, a newly submitted job with that key. A `failed` or `cancelled` job is returned as it is and not resubmitted (a settings PUT or a server restart retries it).
  - `Runs.clonal_agreement(self, rid: str, graph: dict[str, Any]) -> ClonalAgreement`: no settings → `not_configured`; `graph["store_version"] != STORE_VERSION` → `refused` with the note "made before inference isolation; run it again"; a stored file → `ClonalAgreement.model_validate` of it (a file that does not validate counts as missing, with one stderr line); an active job → `pending`; a failed or cancelled job → `failed` with `job.error or "cancelled"` as the note.
  - `Runs.get_run` fills `trust=RunTrust(clonal=self.clonal_agreement(rid, g))`.
  - `Runs.submit_clonal(self, rid: str, settings: ClonalSettings) -> Job | None` (`None` when the file already exists); `Runs.submit_clonal_all(self) -> int` (every run in `store.list()`, under the current settings); `Runs.submit_trust_jobs(self, rid: str) -> list[Job]` (clonal when settings exist).
  - Dedup key: `f"score_clonal:{rid}:{trust_hash(settings)}"`.

- [ ] **Step 1: Write the failing tests** (`tests/api/test_clonal_trust.py`; a 1000-cell bifurcation with clones, prepared into a temporary cache; each test builds its own `JobManager` and `Session`):
  - `test_a_finished_infer_job_queues_scoring`: started manager; PUT settings; POST a `paga_dpt` run; wait for infer; a `score_clonal` job for that run appears in `jobs.of_kind("score_clonal")`; wait for it; `GET /trajectories/{rid}` → `trust.clonal.status in ("ok", "too_few_clones")` and `trust.clonal.primary.n_permutations == 1000`.
  - `test_without_settings_the_status_is_not_configured`.
  - `test_changing_settings_rescores_every_complete_run`: two done runs; PUT new settings → two queued `score_clonal` jobs with the new hash.
  - `test_a_restart_while_scoring_is_pending_heals_itself`: save a run to the store in-process (as `lumencell run` does) with settings in the state file and no job; a new `Session` on an unstarted manager → GET gives `pending` and exactly one queued job; start the manager, wait, GET again → `ok` or `too_few_clones`; a third GET submits nothing.
  - `test_a_failed_score_is_reported_and_not_resubmitted_on_every_get` (Review Focus 4): delete the run's `cells.lcb` after saving; GET, start, wait; two more GETs → `status == "failed"`, a note carries the error, `len(jobs.of_kind("score_clonal")) == 1`.
  - `test_a_result_for_old_settings_is_never_served` (Review Focus 5): unstarted manager; PUT settings A then B; GET → `pending` with `settings == B`; start and wait for both jobs; GET → `settings == B` and the stored file is `clonal-<hash of B>.json`.
  - `test_legacy_runs_are_never_scored`: a run folder with `store_version = 1` → `refused`, no job.
  - `test_job_info_reports_kinds_by_name`: `GET /jobs/{id}` of a `score_clonal` job has `kind == "score_clonal"`, and of a `sleep` job `kind == "sleep"`.
- [ ] **Step 2: Run** `uv run pytest tests/api/test_clonal_trust.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** models first, `make contract`, then runs and session. Add `trust: null` to the `run` fixture in `web/src/test/fixtures.ts` and nothing else.
- [ ] **Step 4: Run** `uv run pytest tests/api -q`. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Serve clonal agreement with every run, scored in the background and healed on GET`.

---

### Task 15: Clonal agreement from the CLI

**Branch:** `pivot/NN-cli-clonal` · **Model:** sonnet — CLI options over existing functions.

**Files:**
- Modify: `src/lumencell/cli.py`, `tests/unit/test_commands.py`, `tests/unit/test_cli.py`

**Interfaces:**
- Consumes: `check_settings`, `settings_columns`, `agreement_for_run`, `trust_name` (Tasks 10, 13); `RunStore.write_trust` (Task 11); `add_clones`, `ClonalLayer` (Task 3).
- Produces: `lumencell run` gains `--clone-key TEXT`, `--donor-key TEXT`, `--single-donor`, `--strata TEXT` (comma-separated). With `--clone-key`, the settings are validated (problems → `typer.BadParameter`), the clone and donor columns are scoring-only for this run (`from_prep(..., clonal_columns=settings_columns(settings))`; strata columns stay visible), and after saving, the agreement is computed in-process, written with `write_trust(rid, trust_name(settings), …)` (so `launch` with the same settings serves it without rescoring), and printed. The block starts with `clonal agreement: <status> (own cells, m = 2)` and then gives cells, clones (effective), donors, purity, η against N1, the headline p with the nulls it is the largest of, η_base, Δη and its p, an `m = 3:` line with η and headline p, and each note. `--donor-key` and `--single-donor` without `--clone-key` → `BadParameter`.
- Produces: `lumencell synth` gains `--fate-bias FLOAT` (0 to 1); when given, `add_clones(adata, ClonalLayer(fate_bias=…), seed=seed)` is applied before writing.

- [ ] **Step 1: Write the failing tests.** `tests/unit/test_commands.py` (the `data` fixture's file with `add_clones` applied): `test_run_prints_clonal_agreement` (exit 0; output contains `clonal agreement: ok (own cells, m = 2)`, `m = 3:` and `delta`; the trust file exists in the run folder); `test_run_refuses_invalid_clonal_options` (both `--donor-key` and `--single-donor`; an unknown column; `--donor-key` without `--clone-key`) → exit 2 with the problem in the output; `test_run_refuses_a_clone_column_as_root` (`--root cluster:clone_id=c0 --clone-key clone_id --donor-key donor`) → exit 2 mentioning "only to score". `tests/unit/test_cli.py`: `test_synth_can_add_clones` (`--fate-bias 0.7` writes `clone_id` and `donor`; without it, neither).
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_commands.py tests/unit/test_cli.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement.** **Step 4: Run** again. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Add clonal options to lumencell run and lumencell synth`.

---

### Task 16: Correct the stale lines in docs/HANDOFF.md

**Branch:** `pivot/NN-handoff-corrections` · **Model:** haiku — a mechanical documentation edit from listed sources.

**Files:**
- Modify: `docs/HANDOFF.md` (only the lines below)

**Interfaces:**
- Consumes: `main` after Tasks 1–15; `docs/backend-robustness-2.md`; `.github/workflows/ci.yml`; `make ci` output.

- [ ] **Step 1:** Run `make ci > "$TMPDIR/ci.log" 2>&1; echo $?`; it must print `0`. Read the test counts from the log.
- [ ] **Step 2:** Edit, each once, keeping the file's style:
  - line 3: the date and `main`'s current short commit.
  - line 32: the 100k-cell budgets were measured for prep, `paga_dpt` and the server in `docs/backend-robustness-2.md` (all within budget); Palantir, scFates and first paint were not; packaging is still not started.
  - line 49: the Python, deselected, web and end-to-end test counts from Step 1.
  - lines 71–73: the server loads the dataset with `keep_matrix=False` and serves expression from the prep cache.
  - lines 89–90: each job kind has a time limit (`DEFAULT_TIME_LIMITS`, `LUMENCELL_JOB_TIME_LIMITS`); a job past its limit, or a worker that dies, fails the job and the worker is replaced; `prepare` and `infer` run ahead of background jobs.
  - lines 159–160 and 185: GitHub Actions CI exists (`.github/workflows/ci.yml` runs the `make ci` steps on Ubuntu with Python 3.13); there is still no Python 3.14 job; item 5 of the open work is done.
- [ ] **Step 3:** `git diff docs/HANDOFF.md` touches only those lines.
- [ ] **Step 4: C1–C4.** Commit message: `Correct the stale lines in the milestone 1 handoff`.

---

### Task 17: Security review of phase 1

**Branch:** `pivot/NN-phase1-security` · **Model:** sonnet — runs the review and fixes small findings test-first.

**Files:**
- Modify: whatever the findings require; `docs/superpowers/HANDOFF.md`

- [ ] **Step 1:** Invoke the `security-review` skill on the diff `142a788..main`. Point it at the untrusted-input paths: obs column names and labels (settings, `clonal_candidates`, notes, trust file names, job keys), the `PUT /dataset/clonal` body, `LUMENCELL_CALIBRATION_OUT`, the trust store's file writes and the housekeeping patterns.
- [ ] **Step 2:** For each finding, write a failing test that shows it, fix it, and see the test pass. A finding too large for this PR gets its own task in HANDOFF instead.
- [ ] **Step 3:** Add a "Security review (phase 1)" section to `docs/superpowers/HANDOFF.md`: every finding, and whether it was fixed (with its test) or recorded, and why.
- [ ] **Step 4 (background):** the Done When commands at the top of this plan; record each result in HANDOFF.
- [ ] **Step 5: C1–C4.** Commit message: `Resolve or record the phase 1 security review findings`.
