# Handoff: T-cell clonotype work

A new session starts with: "Read docs/superpowers/HANDOFF.md and continue."
This file is updated inside every task's PR. `docs/HANDOFF.md` is the milestone 1 note; its
stale lines are corrected once, by phase 1's handoff-corrections task.

## Current state (2026-11-13)

- **Phase:** phase 1 has started. Both approval gates are passed: Shreyan approved the spec
  (PR #5) and the plans (PR #6, merged as 153a7f0) on 2026-11-12. Task 1 (clonal definitions,
  gates and runtime constants) merged as 6ff4676 (PR #7). Task 2 (the `clonal/` package: clone
  codes, scores, nulls, baseline, coupling) merged as f0715b6 (PR #8). Task 3 (synthetic clone
  layer and its two controls) merged as 744a0a8 (PR #9). Task 4 (receptor genes out of HVG
  selection) merged as 7b6254b (PR #10). Task 5 (clone and receptor columns scoring-only,
  store version 2) merged as 3f7e93d (PR #11). Task 6 (clonal scoring in the synthetic
  benchmark) merged as 857d03d (PR #12). Task 7 (re-measured baseline) merged as bfb2539
  (PR #13). Task 8 (calibration of the nulls) is on branch `pivot/10-clonal-calibration`; one
  check failed (see Open questions), so it stops for Shreyan after the PR is opened. Task 9 (real benchmarks re-run after the TR-gene
  exclusion) is on branch `pivot/11-bench-real-rerun`, stacked on `pivot/10-clonal-calibration`.
- **Spec:** approved and merged in PR #5 (142a788):
  [specs/2026-11-12-tcell-clonotypes-design.md](specs/2026-11-12-tcell-clonotypes-design.md).
  Branch `pivot/02-phase-plans` amends it in two places (see Decisions).
- **Plans**, merged in PR #6:
  - [phase 1, clonotypes end to end](plans/2026-11-12-phase-1-clonotypes.md): 17 tasks
  - [phase 2, T-cell robustness](plans/2026-11-12-phase-2-robustness.md): 13 tasks
  - [phase 3, real T-cell data](plans/2026-11-12-phase-3-tcell-data.md): 15 tasks
- **Next:** Shreyan's call on the `depth_only/N2` failure (Open questions), then phase 1 Task 9.
  - Execute phase 1 with `superpowers:subagent-driven-development`, which Shreyan chose:
    - one task per subagent, on the model the task names
    - TDD throughout
    - a fresh reviewer per task, plus the pr-review-toolkit agents on each PR
    - a squash-merge once CI is green
    - branches `pivot/NN-<slug>`, numbered from 03 in execution order
- **Handoffs:** at the end of every phase, and whenever the session nears 35% of its context.

## Controller handoff (2026-11-13, after Task 4)

The orchestrating session reached 33% of its context after Task 4 and handed off here.
- **Read the ledger first:** `.superpowers/sdd/2026-11-12-phase-1-clonotypes/progress.md`. It
  is git-ignored and exists only on this machine. It holds the pre-flight rulings for Tasks
  5–17 (F1, F2, F4, F9, F12–F15, F18, F20, F21, F23, F24, F27), the things carried into Task 13,
  and every deferred minor finding for the final review. `preflight.md` in the same folder has
  the full conflict table. Each task's brief, report and review diff are there too. The Task 5
  brief is already extracted.
- **Rulings that change later tasks** (the ledger has the reasons and costs):
  - **F4, Task 5:** `clonal_columns` is a required key in the infer payload, and Task 5 updates
    `tests/unit/test_jobs.py` to send it. The scan proposed a `()` default, which was rejected:
    it would silently let methods see clone columns.
  - **F24, Task 5:** the refusal message names the column's real role (clone, receptor, donor
    or strata).
  - **F2, Task 6:** the report's clonal column reads rows with `.get(key, nan)`.
  - **F11, Task 6:** a synthetic branching truth with no clone data is `not_measured`.
  - **F23, Task 6:** add one helper in `clonal/scores.py` that runs baseline_space, then
    kmeans_partition, then score_clonal. Tasks 8 and 13 reuse it.
  - **F18, Task 7:** first compare runs with and without the clone layer on the same code,
    then compare against baseline.json.
  - **F12, Task 8:** a NaN p is stored as null and counted apart, and this rule goes into
    thresholds.toml before the calibration run.
  - **F1 and F9, Task 10:** fixtures need at least 1,100 cells. Add a case with unused
    categories.
  - **F20 and F27, Task 13:** `trust_hash` params are a dict, not a set. `ClonalScore` gains a
    `status` field.
  - **Task 2 carry-over, Task 13:** a one-clone test case must spread its clone over two
    lineages, or it reads `single_lineage`.
  - **Task 13:** surface `CloneCodes.n_missing_donor` in the trust note.
  - **F13, Task 14:** a trust file that fails validation is treated as missing.
  - **F14, Task 15:** BadParameter is raised for conflicting donor flags.
  - **F15, Task 16:** edits may span whole bullets.
- **Process as practised.** Each task goes through these steps:
  - one implementer;
  - in parallel, a task reviewer plus the pr-review-toolkit code-reviewer,
    silent-failure-hunter and pr-test-analyzer, plus type-design-analyzer on Tasks 5, 8, 10, 13
    and 14 (Sonnet for all of them, Opus for the statistics-heavy task review);
  - one combined fix round, then a scoped re-review;
  - squash-merge once CI is green.

  The test analyzer mutation-tests a scratch copy, and every task so far had surviving mutants
  that the fix round closed.
- **Process change (2026-11-13), from Task 6 on.** Shreyan asked: "Use token efficient agents
  for the tasks so that it can finish quickly on less usage". So:
  - reviews are two seats: a Sonnet task reviewer that also covers silent failures and
    NaN/fallback paths, and pr-test-analyzer on Sonnet, mutation-testing a scratch copy;
  - no separate code-reviewer or type-design-analyzer; an Opus reviewer only for Task 8;
  - implementers run on Sonnet by default, on Opus only for Tasks 13 and 14, and on Haiku for
    mechanical or docs-only work;
  - scoped re-reviews run on Haiku when the fix is mostly tests;
  - subagents stop after replying, never watch CI, and keep replies under 15 lines.
- **Merges need Shreyan's word.** The auto-mode permission check blocks `gh pr merge`
  ("Merge Without Review"), so Shreyan approved each of #7, #8 and #9 in chat. Ask before every
  merge, or have Shreyan add an allow rule for `Bash(gh pr merge:*)`.
- **Check the 5-hour usage window** (`get_usage`) before starting a five-agent review set. On
  2026-11-12 it reached 91%, and the Task 2 review waited for the reset.
- **Next:** review Task 5 (PR #11, `pivot/07-inference-isolation`; type-design-analyzer is in its
  set), merge on Shreyan's word, then Task 6 on `pivot/08-<slug>`.

## Done

- Reconnaissance (loader, bench, prep, API, robustness-2, web research). The notes are in
  `research/`.
- Baseline at 865ce22: `make ci` exit 0, with 353 Python tests passed and 12 deselected, 36
  web tests and 1 Playwright test.
- Spec revisions 1 and 2. Revision 2 went through an independent method review and an
  architecture review. Shreyan approved it, and it merged as PR #5.
- Three phase plans, each written by its own Opus subagent with the phase 1 interfaces fixed
  first. One independent reviewer then cross-checked them (one blocker, eight should-fixes),
  and every finding was fixed.

## Task 4 decisions

- `PrepParams.exclude_receptor_genes` (default True) is a param of the `hvg` stage, so toggling it
  changes the `hvg` key and every later key; `normalize` is untouched. The `hvg` stage now also
  writes a small `hvg_info` artifact (`receptor_genes`, `excluded`, `receptor_hvg_without_exclusion`).
- With receptor genes present and exclusion on, HVGs are chosen among the non-receptor genes with
  `n_top = min(n_top_genes, n_non_receptor)`; a second selection over all genes only counts the
  receptor genes that would have been HVGs. Without receptor genes the output is identical to before.
- Existing cache keys for every dataset change (the `hvg` params gained a field), so old prep caches
  recompute from `hvg` on. No test pinned a literal key, so no assertion changed.
- `hvg_info` also records `symbol_source` (the var column the loader took symbols from, else
  `var_names`) and `symbols_look_like_ensembl`. `cli prepare` always prints a `receptor genes: ...`
  line naming the source, including `0 matched`. If 0 matched and the symbols look like Ensembl IDs
  (the loader's own test, `looks_like_ensembl`), the hvg stage logs a warning and prepare prints one:
  clone identity may reach the PCA.
- A dataset whose every gene matches the pattern raises `DatasetError` from the hvg stage.
- The regex is unchanged; genome-prefixed or versioned symbols (e.g. `GRCh38_TRAC`) stay unmatched.
  Deferred.

## Task 5 decisions

- **Scoring-only columns** are those matching `RECEPTOR_COLUMN` (clone_id, clonotype, cdr3,
  ir_, airr:, receptor_ prefixes; raw_clonotype_id, chain_pairing), plus any the caller names
  in `clonal_columns`. `InferenceContext.obs` leaves them out; `ctx.scoring_only` lists them.
  The label key is detected without them.
- **One refusal helper**, `scoring_only_message(column, role)` in `clonal/isolation.py`
  (ruling F24): a pattern match "holds clone or receptor data"; a named column "is set as a
  clone or donor column". Roles: root, cluster, label column. Engine and API both use it.
- **Refused before any work.** `resolve_root` checks before the existence check. `run_method`
  checks the method's `cluster_key` and the label key it uses before resolving the root.
  `create_run` checks root, `cluster_key`, label key (422 `scoring_only_column`); `detail.columns`
  names each offender once, in that order.
- **The label key used** is `label_key`, else the `cluster_key` when terminals are named, else
  `ctx.label_key`, which `name_lineages` already fell back to. It is checked and recorded.
- **`obs_keys`** is the sorted set of the resolved root's `cluster_key` (lumencell_leiden for a
  cell or file root), the method's `cluster_key` and the label key used.
- **`clonal_columns` is required in the infer payload** (ruling F4). Without it, or as anything
  but a list of names, the job fails with a message naming the key. `create_run` sends
  `list(view.clonal_columns)`, empty until Task 10.
- **Store version 2.** `result.json` gains `prep_key`, `inference_isolation` (1) and
  `obs_keys`; version 1 runs drop out of the list. Run ids change (`STORE_VERSION` is hashed).
  `clonal_columns` is not part of the run id.
- **Fix round 1:**
  - `InferenceContext.__post_init__` raises `ScoringOnlyColumnError` when `obs` holds a column in
    `scoring_only` or one matching `RECEPTOR_COLUMN`, so `dataclasses.replace` cannot leak one.
  - A bare string as `clonal_columns` (in `scoring_only_columns`, `inference_obs`, `from_prep`,
    `DatasetView`) is a `TypeError`, via `column_names`.
  - `inference_obs` raises `ValueError` for a named column the obs lacks (a typo such as
    "Donor" would leave the real column visible). The isolation test now names `clone_id` only
    for the clonal file.
  - A run from another store version is absent: `exists()` is False and `load()` raises
    `FileNotFoundError`, so it is recomputed. `list()` logs one warning with the count it skips.

## Task 6 decisions

- **`clonal_row`** (`bench/runner.py`) never raises: an error is a `failed` row with
  `reason`, also printed to stderr. `run_real_case` sets `not_applicable`.
- **Rulings F11 and F23.**
  - Fewer than 2 true lineages: `not_applicable`. Branching truth with no `clone_ids`:
    `not_measured`, so both gates fail. `clonal_layer=None` on a branching case gives this.
  - `score_with_baseline` in `clonal/scores.py` runs `baseline_space`, `kmeans_partition`,
    `score_clonal`. The oracle, the method and the m = 3 call all use it.
- **Method status `single_lineage`** (L >= 2 but all scored mass in one lineage) is the row
  status too, like L < 2. The m = 3 fields stay NaN there.
- **An oracle with NaN eta or p is `not_assessable`.** NaN is not above the floor.
- **Evidence counts** (`n_cells`, `n_clones`, `effective_clones`) come from the oracle call.
- **The report** reads rows with `.get(key, nan)`: committed `baseline.json` has no `clonal_*`
  keys until Task 7, so the new column shows the missing mark. It also carries the "not a
  branch gate" sentence.
- **Inference** gets `clonal_columns=("clone_id", "donor")` only when the layer was applied.
- **Gates now measure.** `test_full_matrix_holds_baseline_and_experimental_flag_is_honest`
  (slow, not run here) now sees real clonal rows. Its regression check skips the clonal ratio
  (no baseline value yet); the flag check may move until Task 7 re-baselines.

## Task 7 decisions

- **Makefile.** `BENCH_METHODS ?= paga_dpt,palantir,scfates_ppt` is passed as `--methods` by
  `bench-synth` and `bench-baseline` (both defaulted to `paga_dpt` only). Tested by
  `tests/unit/test_makefile.py` (`make -n`, cwd = repo root).
- **Existing metrics did not move.** Compared with main's `baseline.json` (unchanged since
  142a788), every non-runtime field is identical on all 36 rows (`0 []`). So neither the clone
  layer nor the receptor exclusion leaks into synthetic inference.
- **Clonal ranges, 27 branching rows** (linear: 9 rows, `not_applicable`):
  - excess ratio 0.30-0.62 (gate 0.5): paga_dpt 0.30-0.54, palantir 0.37-0.56, scfates_ppt 0.44-0.62;
  - below the gate: paga_dpt 6 of 9, palantir 7 of 9, scfates_ppt 5 of 9;
  - eta method 0.22-0.47, eta oracle 0.73-0.83, eta m=3 0.21-0.49;
  - effective clones 68-207; p = 0.001 on every row (the B = 1000 floor);
  - no NaN where status is `ok`; cells >= 2 x clones everywhere.
- **Experimental flags unchanged**; each method already failed a branch gate. The three
  `experimental_note`s now also name the clonal gate. No method code, default or threshold changed.
- **README bench table** regenerated: new "Clonal excess ratio" column (printed to 3 decimals,
  rounded down so a failing value never shows as 0.500) and a per-gate failure count. Cases
  failing any gate are now 9, 7, 7 of 12 for paga_dpt, palantir, scfates_ppt (were 7, 2, 2); the
  rise is the clonal gate alone.
- **`terminal_share: NaN`** is on every synthetic baseline row. Expected: the key was added to
  `runner.py` in 0621edc (real datasets); main's baseline predates it.

## Task 9 result

`make bench-real` exit 0 at 4b0b1d5 (3 real_data tests passed); record
`benchmarks/results/real-2026-11-13.json`. Real rows are `clonal_status = not_applicable`.
- `hvg_info` receptor_genes / receptor_hvg_without_exclusion: pancreas 0 / 0, bonemarrow 11 / 6,
  gastrulation_erythroid 258 / 3.
- Moved: pancreas nothing. Bone marrow: paga_dpt branch_accuracy 0.917 -> 0.617 and lineage_recall
  1.0 -> 0.5 (gates now fail), scfates_ppt pt_spearman_lineage 0.781 -> 0.668 (fails), palantir
  branch_accuracy 0.639 -> 0.943 (improves). Gastrulation: paga_dpt lineage_recall 0.0 -> 1.0.
  Eight gates fail in total, as the log lists. No threshold changed.
- The 11-10 file predates other prep/engine commits, so not every move is the exclusion alone.
  Full table: `.superpowers/sdd/2026-11-12-phase-1-clonotypes/task-9-report.md`.

## Task 8 decisions and result

Run 2026-11-13 at cb5c9c3 (the pre-declaration commit), 200 seeds (100-299), 39 min. Nothing was
tuned. Record: `benchmarks/results/clonal-calibration-2026-11-13.json`. Bound 18 (n scored 200).
- **NaN rule**, declared in `thresholds.toml` before the run: a NaN p is `null`, counted in
  `n_nan`, left out of the bound's n. There were none, and no method raised.
- **Result per check (rejections of 200, p <= 0.05):**

| Check | Rej. | Expected | Result |
|---|---|---|---|
| null/paga_dpt | 6 | calibrated | pass |
| null/palantir | 7 | calibrated | pass |
| null/scfates_ppt | 6 | calibrated | pass |
| donor_confound/N1 | 7 | calibrated | pass |
| donor_confound/N1_donors_ignored | 200 | not calibrated | pass (has teeth) |
| depth_only/N1 | 200 | not calibrated | pass (has teeth) |
| depth_only/N2 | 57 | calibrated | **FAIL** |
| monotonicity, 9 cases | all rise; p at 0 is 0.21-0.998 | | pass |

- **Reading (Step 5; script and output in the task report).**
  - The three method nulls reject 3.0-3.5% (exact 95% CI 1.1-7.1%): conservative, as discrete
    permutation p-values should be, and all inside the bound. KS against U(0,1) is 0.03 / 0.007
    / 0.11 for the three methods: a mild departure for paga_dpt and palantir, from having too
    few small p, not too many.
  - donor_confound/N1 rejects 3.5% with a KS p of 0.61: well calibrated.
  - Both controls without their stratum reject 200 of 200 at the p floor 0.001, so the controls
    have teeth: shuffling without donors (or without depth) is wrong on these data.
  - **depth_only/N2 failed its calibration check:** 57 of 200 reject (28.5%, CI 22.4-35.3%;
    bound 18). Cause, reproduced by a reviewer on seeds 100-139 (9 of the first 40 reject, as
    recorded): binary2's lineage mix steps at t = 0.65. Oracle cells are single-lineage cells, so
    below 0.65 every scored cell is lineage B, and above it B, A1 and A2 in about equal shares.
    That step lies inside one decile, so the within-decile shuffle mixes pure-B clones with
    mixed-region cells and the null entropy comes out too high. About 96% of clones lie in one
    decile, which does not by itself tell the two readings apart. Exploratory probe on 40 seeds,
    not a result: splitting at t >= 0.65 gave 1 of 40 rejections, 40 quantile bins 3 of 40.
  - Monotonicity is a weak check: eta at fate_bias 1.0 is 1.0 by construction and the steps are
    about 0.35, so it fails only if something is grossly broken.
- **Why no gate moves:** under depth confounding N1 sits at the floor (0.001 on all 200 seeds), so
  the headline p (max over nulls) equals N2's p and inherits its rejection rate; the max protects
  only when a null is valid for the confound present. The gate is `clonal_excess_ratio`, and p
  enters only as the oracle's assessability test `oracle.p <= 0.01` on layers with control `none`.
  Nothing was re-run.
- **Provenance:** the recorded run's commit, cb5c9c3, is HEAD at the time. The calibration code
  and tests were uncommitted then and first committed in 52f739d; reviewers reproduced seeds
  100-139 exactly on that code. The recorded JSON carries this as a `provenance` note, numbers
  unchanged; the recorder now appends `-dirty` to a dirty commit.
- **Fix round 1:** a check with no scored seed never passes, and a report with any method error
  does not pass.

## Verified references

Checked on 2026-11-12 against Crossref (`api.crossref.org`) and jmlr.org; nothing downloaded.
The spec's §10 still lists them as unverified (the spec is not edited). They are cited in the
Sources lines of `[clonal]` in `benchmarks/thresholds.toml`.

| Reference | Journal, year | DOI or volume and pages |
|---|---|---|
| Phipson B, Smyth GK. "Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn" | Stat Appl Genet Mol Biol 9(1), 2010 | 10.2202/1544-6115.1585 |
| Vinh NX, Epps J, Bailey J. "Information Theoretic Measures for Clusterings Comparison: Variants, Properties, Normalization and Correction for Chance" | J Mach Learn Res 11(95), 2010 | no DOI; pp. 2837-2854 (jmlr.org/papers/v11/vinh10a.html) |
| Sturm G, Szabo T, Fotakis G, Haider M, Rieder D, Trajanoski Z, Finotello F. "Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data" | Bioinformatics 36(18), 2020 | 10.1093/bioinformatics/btaa611; pp. 4817-4818 |
| Schattgen SA, Guion K, Crawford JC, Souquette A, et al. "Integrating T cell receptor sequences and transcriptional profiles by clonotype neighbor graph analysis (CoNGA)" | Nat Biotechnol 40, 2022 (online 2021-08-23) | 10.1038/s41587-021-00989-2; pp. 54-63 |

## Decisions

The spec's decisions are in its §8. Decisions made at plan time:
- **Spec amendment 1:** the receptor-gene regex is `^TR[ABDG]([VDJ]\d|C\d*$)`. The original
  also caught TRADD.
- **Spec amendment 2:** `strata_keys` columns are not scoring-only, so Sheih's time point may
  define a root.
- **Sheih runs per patient.** Pooled patients would make batch structure look like "no tree".
  A pooled run is a labelled diagnostic only. Cost: 12 runs instead of 3, with fewer cells
  each.
- **Daniel's "not assessable" is shown explicitly** in the results JSON and the README, never
  omitted.
- **Each plan lists its own resolved ambiguities**, each with its cost. Read them before
  executing a task.

Decisions taken in phase 1 Task 2 (`clonal/` package):
- **Undefined statistics are NaN.** η against a null is NaN when that null's mean H is 0, and
  so is η_base. `p_delta` is NaN whenever Δη is NaN (a constant baseline, or every scored cell
  in one k-means cluster), never the most significant value. A coupling pair with an all-zero
  column has NaN fields and is left out of the Benjamini-Hochberg adjustment. z is NaN when the
  N1 null sd is 0, to within the 1e-12 comparison slack, so rounding noise never yields a huge z.
- **`single_lineage` also covers one lineage holding every scored cell** (controller ruling,
  fix round 1): the status is `single_lineage` when L < 2, or when the scored weights put mass
  in fewer than two lineages. Before, such a run read as perfect agreement (purity 1, η NaN,
  p 1). The counts and `effective_clones` are reported, because the effective count depends only
  on clone sizes; every other float is NaN, and `nulls` and `coupling` are empty. With no scored
  clone the status stays `too_few_clones`, as the plan says, and the floats are NaN.
- **`min_clone_size` below 2 is refused** (ValueError), because single-cell clones never
  contribute (spec §2). `min_clones`, `n_permutations` and `n_bins` must be at least 1.
- **Normalised entropies are clipped to [0, 1]** against rounding, so P and both purity
  diagnostics stay in the unit interval exactly.
- **The baseline's entropy is divided by log max(L, number of baseline groups)**; η_base is a
  ratio, so the constant cancels, and the clip never bites.
- **Shuffles are scored in blocks of rows.** A block draws `rng.random((rows, n_scored))`, the
  same stream as `rows` draws of `rng.random(n_scored)`, so the plan's shuffle is unchanged. The
  block keeps rows × cells × columns near 2^18 elements (one shuffle per block when a single
  shuffle is larger).

Decisions taken in phase 1 Task 3 (`bench/synth_clones.py`):
- **"Lineage j" is the index in sorted leaf-name order** (controller ruling F16), in `add_clones`
  and in the tests. `Truth.names` follows topology order, which differs for binary2 (B, A1, A2
  against sorted A1, A2, B); the donor-confound test has a binary2 case so a mix-up fails.
- **Donor-confound draws.** After the uniform donor draw, the control draws `rng.random(n)`
  (keep the home donor, probability 0.8) and `rng.integers(D - 1, size=n)` (which other donor),
  used only for cells on exactly one lineage; the brief fixed the model but not these draws.
- **Exclusive cells** are found from the segment tree (a segment on exactly one leaf's path), so
  `add_clones` needs only `true_segment`, `true_pseudotime` and `segments`.
- `Truth.clone_ids` and `Truth.donors` are set together or not at all, each of shape `(n,)`.

## Found along the way

- **Local `main` was 13 commits behind `origin/main`** on 2026-11-12. Check
  `git rev-list --left-right --count main...origin/main` before reading code.
- **`main` has no branch protection.** Keeping work on PRs is discipline, not enforcement.
- **The `need_from_me` section of the brief was missing.** It is taken to mean dataset
  downloads and new runtime dependencies (spec §8, decision 20).
- **`make bench-baseline` measures only paga_dpt today.** Phase 1 Task 7 adds
  `BENCH_METHODS`.
- **ruff 0.17 formats Python code blocks inside Markdown files.** `make ci` fails on an
  unformatted ```` ```python ```` block in a doc, so run `uv run ruff format <file>.md` on new
  plans and docs.
- **Clonal gates failed by design until Task 6.** Before it every synthetic row was `not_measured`.
  Task 6 measures them; the slow golden test and `make bench-synth` have not been run on it.
- **Gate failures don't fail the benchmark commands.** `lumencell bench real` and
  `bench synth` exit 0 on gate failures; only regressions fail `bench-synth`.

## Open questions for Shreyan

- **depth_only/N2 is a failed calibration check (Task 8, ruled by the controller).** 57 of 200
  seeds reject at 0.05 (bound 18); N2 does not protect against depth confounding across a branch
  point that falls inside a decile. It is recorded as failed. No N2 definition, seed or threshold
  was changed. The slow test `test_depth_only_is_calibrated_under_n2` stays red.
  Phase-2 candidate (2026-11-13), for your sign-off: a pre-declared N2 amendment, either finer
  bins or a branch-point split, defined in thresholds.toml before a run and calibrated on a fresh
  seed range, never on 100-299. The reviewer's exploratory 40-seed numbers (3 of 40 for 40 bins, 1
  of 40 for the branch-point split) are exploratory, not a result.
- **Reordering cells: exact or within Monte Carlo error?** Spec §4.5 says η, P and the p-values
  are "unchanged" when cells are reordered. Task 2 keeps the plan's shuffle: each null draws one
  random key per scored cell, in the cells' original order. That makes relabelling clones or
  lineages and adding single-cell clones leave every field exactly unchanged, but reordering
  cells moves the keys between cells. Reordering then leaves the deterministic fields (cells,
  clones, effective clones, P, H, both purity diagnostics) exactly unchanged and the null-based
  fields (null mean and sd, η, p, η_base, Δη, p_delta, coupling) within Monte Carlo error;
  `test_reordering_cells_keeps_deterministic_fields_and_null_fields_within_monte_carlo_error`
  bounds each null's p and η. Exact invariance to both relabelling and reordering can't hold together with this
  shuffle. Needs your sign-off. If you want exact reordering invariance, the cost is reworking
  the shuffle to draw keys in a canonical cell order (for example by cell barcode, which callers
  would then pass in), then re-checking the relabelling and singleton properties.
- Commit 33c281d and PR #5 were created at 07:57 under your account, not by the orchestrating
  session. Was that you, or a tool?

## Branches and PRs

| Branch | PR | State |
|---|---|---|
| `pivot/01-tcell-spec` | #5 | merged (squash, 142a788) |
| `pivot/02-phase-plans` | #6 | merged (squash, 153a7f0) |
| `pivot/03-clonal-definitions` | #7 | merged (squash, 6ff4676) |
| `pivot/04-clonal-scores` | #8 | merged (squash, f0715b6) |
| `pivot/05-synth-clones` | #9 | merged (squash, 744a0a8) |
| `pivot/06-receptor-genes` | #10 | merged (squash, 7b6254b) |
| `pivot/07-inference-isolation` | #11 | merged (squash, 3f7e93d) |
| `pivot/08-bench-clonal` | #12 | merged (squash, 857d03d) |
| `pivot/09-bench-remeasure` | #13 | merged (squash, bfb2539) |
| `pivot/10-clonal-calibration` | #14 | phase 1 Task 8: one check failed, awaiting Shreyan |
| `pivot/11-bench-real-rerun` | (this PR) | phase 1 Task 9: real re-run, stacked on #14 |

## Verify the current state

```bash
git fetch origin && git rev-list --left-right --count main...origin/main
make install
make ci
make bench-synth
```
