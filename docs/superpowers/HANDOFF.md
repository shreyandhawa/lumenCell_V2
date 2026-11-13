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
  store version 2) is on branch `pivot/07-inference-isolation` (PR #11), awaiting review.
- **Spec:** approved and merged in PR #5 (142a788):
  [specs/2026-11-12-tcell-clonotypes-design.md](specs/2026-11-12-tcell-clonotypes-design.md).
  Branch `pivot/02-phase-plans` amends it in two places (see Decisions).
- **Plans**, merged in PR #6:
  - [phase 1, clonotypes end to end](plans/2026-11-12-phase-1-clonotypes.md): 17 tasks
  - [phase 2, T-cell robustness](plans/2026-11-12-phase-2-robustness.md): 13 tasks
  - [phase 3, real T-cell data](plans/2026-11-12-phase-3-tcell-data.md): 15 tasks
- **Next:** phase 1 Task 6, after Task 5's PR is reviewed and squash-merged.
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
- **Clonal gates fail by design until Task 6 merges.** From the Task 1 commit until then every
  synthetic row is `clonal_status` not_measured, so both clonal gates fail on every method.
  The slow golden test `test_full_matrix_holds_baseline_and_experimental_flag_is_honest`
  (deselected in `make ci`) and the `make bench-synth` gate report are red by design.
- **Gate failures don't fail the benchmark commands.** `lumencell bench real` and
  `bench synth` exit 0 on gate failures; only regressions fail `bench-synth`.

## Open questions for Shreyan

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
| `pivot/07-inference-isolation` | #11 | phase 1 Task 5: awaiting review |

## Verify the current state

```bash
git fetch origin && git rev-list --left-right --count main...origin/main
make install
make ci
make bench-synth
```
