# Handoff: T-cell clonotype work

A new session starts with: "Read docs/superpowers/HANDOFF.md and continue."
This file is updated inside every task's PR. `docs/HANDOFF.md` is the milestone 1 note; its
stale lines are corrected once, by phase 1's handoff-corrections task.

## Current state (2026-11-20)

- **Phase 2 is complete on `main` (570fb5e).** Merged: Tasks 1 (#26), 5 (#27), 6 (#28), 7 (#29), 8 (#32), 9 (#33),
  10 (#31), 11 (#34), 12 (#35), 13 (#36), N2a (#30), and the late Tasks 2 (#39), 3 (#42) and 4 (#46); the start (#25),
  the final review (#37) and the follow-ups H-1 (#38), the disconnected-kNN error (#40) and the X-1 phase 3 plan
  amendment (#41). Phase 2 PRs are #25 to #40, #42 and #46 (#41 is the phase 3 plan amendment; #43 to #45, #47 and #48 are phase 3). This PR (`pivot/45-phase2-close`)
  records the Done When on 570fb5e (see "Done When (phase 2)") and this state.
- **What failed and was accepted:** the verdict verification (Task 9) failed on fresh seeds 3-5 (28 of 31 control runs
  supported, max 0). Shreyan ruled (a): accept a provisional-only verdict; a supported verdict is capped at
  `provisional`, and a failed check is served `weak` (H-1, 2026-11-19). Two Done When commands therefore fail
  (the verification's `passed` is `False` and the slow test that asserts it); a pre-declared C1/C2 redesign is the only
  way to pass them and is not planned. C2, C3 and C4 were non-discriminating in Task 8, so no tau is set for them. The
  calibrated N2a split is adopted (clonal `metric_version` 2).
- **Phase 3 is under way** (base 35ac459). Merged: Task 1 (#43, multi-file sources, licences, build step), Task 2 (#44,
  the nine Daniel 2022 GEO files, sizes and SHA-256), Task 3 (#48, Daniel sample records and metadata; 23 captures),
  Task 4 (#47, clone calling from 10x contigs) and Task 9 (#45, Sheih 2020 joins at 0.99, so the Cao branch is not
  taken). Task 10 (#51, merged) declares `[case_study.sheih2020]` before any analysis. Task 5a (#50, merged): Daniel
  file formats and the capture map, all 23 TCR members joined to their cells. Task 11 (#52, merged) fetches the rest
  of GSE125881 and builds the case-study dataset (62,167 cells, no analysis yet). Task 5b (#53, merged): the Daniel
  builder, tested on a synthetic layout only; the real derived file is not built. Task 6 (#54, merged): the Daniel
  declaration (8 chronic spleen D8/D21 captures, namespace = capture, root `TEx_prog`, no gated lineage, clonal
  declared not assessable), before any clone was counted. Task 7 (#55, merged): clonal scores and reported fields for real
  datasets, tested on a synthetic `fake` fixture only; a dataset with no gated lineage is shown as not measured, never
  PASS. Task 8 (this PR, open): the first real run with Daniel rows, `benchmarks/results/real-2026-11-20.json`; no Daniel gate measured. Next: Task 12 after 7 and 11. See each task's
  decisions section.
- **Spec:** approved and merged in PR #5 (142a788):
  [specs/2026-11-12-tcell-clonotypes-design.md](specs/2026-11-12-tcell-clonotypes-design.md),
  with the two plan-time amendments under Decisions.
- **Plans**, merged in PR #6:
  - [phase 1, clonotypes end to end](plans/2026-11-12-phase-1-clonotypes.md): 17 tasks, done
  - [phase 2, T-cell robustness](plans/2026-11-12-phase-2-robustness.md): 13 tasks plus N2a, done
  - [phase 3, real T-cell data](plans/2026-11-12-phase-3-tcell-data.md): 15 tasks, amended by X-1 (#41); in progress
- **Phase 2 process:** executed task by task with subagents: one implementer per task, TDD, two review seats, one fix round
  and a scoped re-review, then a squash-merge on green CI. Branches `pivot/NN-<slug>`.
  Ledger (git-ignored, this machine only): `.superpowers/sdd/2026-11-17-phase-2-robustness/progress.md`.
- **Handoffs:** at the end of every phase, and whenever the session nears 35% of its context.

## Phase 2 rulings (Shreyan, 2026-11-17: "approve both")

- **N2 amendment approved for phase 2.** The phase 1 `depth_only/N2` failure (57/200, bound 18)
  stays recorded as failed, and `clonal-calibration-2026-11-13.json` is never edited. Phase 2 adds
  a pre-declared amendment to N2: its definition, reasons and date go into
  `benchmarks/thresholds.toml` in a commit that lands before any run that measures it, and it is
  calibrated on a fresh seed range, never 100-299. Where and how it fits is ruled in "Phase 2 pre-flight rulings" (task N2a).
- **F10 accepted: reordering cells is within Monte Carlo error, not exact.** Reordering cells
  leaves the deterministic fields exactly unchanged and the null-based fields (null mean and sd,
  eta, p, eta_base, delta_eta, p_delta, coupling) within Monte Carlo error, as
  `test_reordering_cells_keeps_deterministic_fields_and_null_fields_within_monte_carlo_error`
  bounds. Spec §4.5's "unchanged" is read that way for those fields. The shuffle is not reworked.
  Closed.

## Phase 2 pre-flight rulings (controller, 2026-11-17)

One Sonnet scan checked the plan against `main`; its full table is in the ledger folder
(`preflight.md`, findings F1-F22). Every consumed interface exists. Rulings that change a task:
- **N2a, the N2 amendment (its own task and PR, Opus, after Task 5).** A branch-point split, chosen
  over finer bins because it targets the diagnosed mechanism (lineage composition steps inside a
  decile), where finer bins only shrink the step and cost power on every run.
  - Definition: for each lineage, its onset is the smallest pseudotime among the scored cells whose
    whole weight is on that lineage. A lineage with no such cell adds no cut. N2's strata become donor
    x decile x segment, where segment counts the onsets at or below a cell's pseudotime. Cells
    without pseudotime keep their own stratum. Near-equal onsets are not merged. The strata use
    only lineage weights and pseudotime, never clone labels, so the shuffle stays a valid null.
  - It is the same rule in own-cells and oracle-cells modes, built from the weights and pseudotime
    each call already passes to `score_clonal`. N2 is redefined in place, so N3's RNG stream and the
    API's null names do not move.
  - Order inside the PR: (1) a commit that declares the definition, its reasons and the date in
    `thresholds.toml`, plus a new `[clonal.calibration_n2_amended]` table with seeds 300-499. That
    table keeps the phase 1 checks, alpha, bound rule and NaN rule. `[clonal.calibration]` and
    `clonal-calibration-2026-11-13.json` stay untouched. (2) The implementation, behind the pinned
    definition. (3) The 200-seed run in the background, about 40 min, recorded as a new file.
  - Adoption is decided now, before the run. If every check passes, the amended N2 becomes the
    runtime N2 and `CLONAL_METRIC_VERSION` goes from 1 to 2. That changes `trust_hash`, so every
    stored clonal agreement is re-scored on its next GET. If any check fails, it is recorded as
    failed, and the runtime stays on metric version 1.
  - Cost: an Opus implementer, a 40-minute run, and one `make bench-synth` to confirm nothing
    regresses.
- **F14, verification failure is planned for now.** If Task 9 fails, the verdict machinery still
  ships (Tasks 10-12 go ahead), but every served verdict is capped at `provisional` (amended 2026-11-19 by H-1: only `supported` is capped), with a note
  that the thresholds failed verification on fresh seeds. Task 1 declares this as
  `[tree_support.verification] on_failure = "provisional"`; Task 9 records the outcome as a constant
  (`VERDICT_VERIFIED`), mirrored from the TOML. No rule, seed or control is changed to rescue a
  failed verification.
- **F2, Tasks 1 and 8.** Task 1 writes its test as planned. Task 8 updates that test: it replaces
  `"calibrated" not in ...` with checks on the calibrated table.
- **F12, C5 (Task 7), the spec wins.** C5 is `not_applicable` only when clones are not configured or
  L < 2. With L >= 2 and every scored clone in one lineage (`single_lineage`), C5 fails.
- **F3, F4 and F5, Tasks 10 and 11.** `result.json` records the requested `terminals` (as sent) and
  the run's `clonal_columns`. The stability payload takes both from the run, never from the current
  settings, and the server reads the settings through `view.clonal_snapshot()`.
- **F1, Task 11.** Stability is queued for every finished run whether or not clones are
  configured. Task 11 updates
  `tests/api/test_clonal_trust.py::test_submit_trust_jobs_queues_clonal_only_with_settings` and any
  API test that now sees a stability job.
- **F6 and F7, Task 6.** `restrict()` slices `genes_detected`, `embeddings` and `obs`, and maps
  `uns_root` into the subset (dropped when the subset lacks it). `STAGE_BY_NAME` goes in `stages.py`.
  The Jaccard matrix and the Hungarian pairing are pulled out of `compare.py` into named functions.
- **F13.** Replicate prep stays in memory (plan), not in the job's temporary directory (spec §5.2).
  Neither touches the stage cache.
- **F16, Task 10.** The cancel test uses 3000 cells, so it does not hinge on timing.
- **F17, Task 12.** `runtime_seconds` stops before stability runs.
- **F8, Task 3, accepted.** The `cell_cycle` key covers only its two params, so a later edit to the
  list or seed leaves a stale, report-only hint until `PIPELINE_VERSION` moves.
- **F22.** After Task 3 every existing prep cache lacks `cell_cycle`, so the first server start
  re-runs that one cheap stage. That is not a regression.
- **Execution order:** 1, 5, N2a, 6, 7, then 2-4 (once the downloads are approved), 8 with 10
  alongside, 9, 11, 12, 13.

## Phase 2 N2a decisions and result

Three commits, in the declared order: b0e235a declares the amendment, its reasons, the adoption rule
and `[clonal.calibration_n2_amended]` (seeds 300-499) before any measurement; 2b15c77 implements it
(`n2_split` on `score_clonal`, `_onset_segments`, `run_calibration(table=, n2_split=)`) and the slow
entry point `tests/clonal/test_calibration_n2_amended.py`; 7adf711 records the run and adopts it, then
`origin/main` (#26 and #27, squash-merged) is merged in.
- **Run:** 2026-11-17 at 2b15c77 (clean tree), 200 seeds, 30.4 min, run once. Record:
  `benchmarks/results/clonal-calibration-n2-amended-2026-11-17.json`. Bound 18, no NaN p, no
  method error. Every null p, so every headline p, used the amended N2.

| Check | Rej. | Expected | Result |
|---|---|---|---|
| null/paga_dpt | 4 | calibrated | pass |
| null/palantir | 3 | calibrated | pass |
| null/scfates_ppt | 3 | calibrated | pass |
| donor_confound/N1 | 15 | calibrated | pass |
| donor_confound/N1_donors_ignored | 200 | not calibrated | pass (has teeth) |
| depth_only/N1 | 200 | not calibrated | pass (has teeth) |
| depth_only/N2 (amended) | 8 | calibrated | pass |
| monotonicity, 9 cases | all rise; p at 0 is 0.13-0.996 | | pass |

- **Reading.** The amended N2 brings depth_only from 57/200 (phase 1, seeds 100-299) to 8/200
  (4%), so the split removes the diagnosed excess. The method nulls stay conservative (1.5-2%).
  donor_confound/N1 is 15/200 (7.5%) against 7/200 in phase 1: inside the bound, but above
  nominal. N1 does not use N2 at all, so this is seed-range variation, not the amendment; it is
  worth watching if a later run lands near 18. Monotonicity remains a weak check (see Task 8).
- **Adopted, by the pre-declared rule:** `RUNTIME_N2_SPLIT = "lineage_onset"` (scores.py),
  `CLONAL_METRIC_VERSION` and `[clonal] metric_version` 1 -> 2. `trust_hash` changes, so every stored
  clonal score is re-scored on its next GET. No test pinned a v1 hash. Phase 1's calibration entry
  point keeps `n2_split="none"`, so it still describes its recorded file.
- **Implementation notes.** Onsets come from scored cells whose weight on the lineage is exactly 1
  (own-cells: every scored cell; oracle-cells: cells in exactly one of the method's lineages).
  NaN-pseudotime cells are segment 0, already apart through decile bin `n_bins`. The stratum id is
  `(donor * (n_bins + 1) + decile) * (L + 1) + segment`; with "none" it is the v1 id exactly.
  Cells after the step but before the later lineages' onsets form a thin segment (as declared:
  near-equal onsets are not merged).
- **bench-synth** at metric_version 2: exit 0, 0 `REGRESSION` lines. 29 gate lines are reported
  (as before: clonal_excess_ratio and pt_spearman rows); no clonal_p gate row fails.
- Caveat: each onset is the earliest exclusive scored cell, so one stray early cell moves a cut (as declared; the calibration passed).

## Phase 2 Task 8 decisions and result

Two commits: 0a584d5 adds `bench/verdict_calibration.py` and its tests (make ci 0); the run was made on
that clean tree; the second commit records the file and applies the rule. Run 2026-11-17, 6.3 min,
run once: `benchmarks/results/verdict-calibration-2026-11-17.json` (commit 0a584d5, `complete`, 72 of 72
records). The quantile arithmetic was recomputed independently from the records and matches exactly.

| Check | Matrix q10 (n) | Control q90 (n) | Result |
|---|---|---|---|
| C2 x0.5 | 0.908 (12) | 0.852 (12) | discriminates, tau 0.880 |
| C2 x1.0 | 0.833 (12) | 1.000 (12) | non-discriminating |
| C2 x2.0 | 0.714 (12) | 0.919 (12) | non-discriminating |
| C3 | 1.000 (36) | 1.000 (31) | non-discriminating (touch) |
| C4 | 0.998 (36) | 0.999 (31) | non-discriminating |

- **By the declared rule:** C2 must discriminate at all three resolutions and fails at x1.0 and x2.0, so
  C2, C3 and C4 are all non-discriminating. `[tree_support.calibrated]` holds only
  `non_discriminating = ["C2", "C3", "C4"]`; `TAU_TREE`, `TAU_LIN` and `TAU_PT` are `None`;
  `calibrated_thresholds()` returns them. Nothing was re-run; rule, seeds and controls unchanged (F14).
- **Runs left out (5, reference failed, all on `two_groups`):** palantir seeds 0 and 2
  (`InvariantError: fate probabilities must sum to 1`) and scfates_ppt seeds 0-2
  (`IndexError: index 0 is out of bounds for axis 0 with size 0`). Both methods fail on a kNN graph
  in two unconnected parts. That is an engine finding (an unhandled disconnected graph gives an
  internal error, not a clear message), not fixed here. Their datasets still give C2 through the
  paga_dpt record. No stability call raised; no prep failed.
- **NaN only where expected:** no NaN in coverage, C3 or C4 of any kept record. C2 is NaN (no PAGA edges,
  counted as 0) for blob seeds 1-2 and two_groups seeds 0-2 at x0.5 only.
- **Reading, plainly.** On these synthetic trees every method is almost perfectly stable (C3 exactly 1.0 on
  all 36 matrix runs, C4 0.993-0.9996), so the 10% matrix quantile sits at the top of the scale and any
  control that is also stable reaches it. Controls that are stable are common: 11 of 31 control runs
  find one lineage, which scores C3 = 1.0 trivially, and paga_dpt and palantir give C4 0.996-1.0 on
  loop and sheet. C4 was expected to drop out (Process Notes). C2 at x1.0 is beaten by two_groups,
  whose two unconnected components form a forest, so all their weight is on the tree (1.0), and by
  loop (0.91-0.93): PAGA keeps the loop's closing edge weak. two_groups would be caught by C1 instead
  (paga_dpt coverage 0.50-0.51). blob is clearly not tree-like at x1.0 and x2.0 (0.05-0.23) but is a
  single cluster or a pair (1.0 or NaN) at x0.5. So the checks as defined do not separate trees from
  these controls on this matrix; the controls show where they fail (single lineages, forests, loops).
- **Consequence for Task 9.** With every tau `None`, a run's level depends on C1 (and C5, not applicable
  here). Every control run with coverage >= 0.95 will be verdicted `supported`, so the verification is
  expected to fail and the declared `on_failure = "provisional"` cap to apply. Not rescued here.
- **Checks applied (data:validate-data, statistics by hand):** quantiles reproduced; denominators
  (12 datasets per group for C2, 36 and 31 runs for C3/C4) match the records; the left-out runs are all one
  control kind, and leaving them out cannot change any result (without two_groups the control q90 is
  still 1.0 for C3 and >= 0.998 for C4, and loop alone beats the matrix at C2 x1.0 and x2.0).
  Caveat: each control quantile rests on 12 datasets or 31 runs, so it is set by its top two to four values.

## Phase 2 Task 9 decisions and result

Two commits: 37197c8 adds `verify`, `Verification` and `verification_payload` (make ci 0); the run was made
on that clean tree, 2026-11-17, 6.2 min, once: `benchmarks/results/verdict-verification-2026-11-17.json`
(commit 37197c8, `complete`, 72 of 72 records). The second commit records the file and applies F14.

- **Result: verification FAILED, as expected.** Controls verdicted supported: 28 of 31 control runs (max 0):
  blob 9, loop 9, sheet 9, two_groups 1 (palantir seed 5). The 3 not supported are two_groups paga_dpt
  seeds 3-5 (coverage 0.48-0.51, C1 fails). `passed` is false.
- **Matrix:** 36 of 36 runs verdicted supported (share 1.0, reported, not gated).
- **Runs left out (5, reference failed, all two_groups):** palantir seeds 3 and 4, scfates_ppt seeds 3-5,
  the same two engine findings as Task 8 (a kNN graph in two parts gives an internal error). One palantir
  subsample replicate raised (non-finite input) on a dataset that still ran; it counts as zero agreement.
- **Reading, plainly.** With C2, C3 and C4 all non-discriminating (Task 8) the verdict is C1 alone, and C1
  only asks that 95% of cells have a pseudotime. Every single-lineage, loop or sheet control reaches
  that, so the verdict cannot tell a tree from these controls: supported means "connected", not "a tree".
  Nothing was tuned or re-run; rule, taus (all None), seeds and controls unchanged.
- **F14 applied.** `[tree_support.verified]` (passed = false, date, source) mirrors into
  `VERDICT_VERIFIED = False` and `VERDICT_VERIFIED_ON`. In `_level`, while it is False, supported
  becomes provisional with a reason naming the failed verification and its date; not_supported is kept.
  (Amended 2026-11-19 by H-1: weak is kept too, see "Phase 2 H-1 decision".)
  `verdict(..., verified=None)` reads the constant; `verify` passes `verified=True` so the rule is measured,
  not the cap. The share of matrix runs supported is reported, not gated. The existing verdict tests run
  with the constant patched True; the capped state has its own tests.
- **Tasks 11 and 12 go ahead with the cap (ruling F14, which supersedes plan Step 6).** Shreyan ruled on
  2026-11-17 ("yes downloads and a"): option (a), accept a provisional-only verdict. The F14 cap stays
  (`VERDICT_VERIFIED = False`); wherever the verdict is served or documented, the failed verification is
  stated plainly. Nothing is re-run or tuned to rescue it.
- **Follow-up, a pre-declared redesign of C1/C2 (not in phase 2).** Task 8 and this run show C1 (coverage)
  passes on every control except two_groups paga_dpt (0.48-0.51), and C2 (PAGA tree-likeness) does not
  separate loops and forests from trees; C3 reads 1.0 on many controls, mostly because a run with one
  lineage matches itself trivially (11 of 31 control runs in Task 8, 13 in Task 9), not because of the
  Jaccard >= 0.6 cut (multi-lineage blob runs score 0-0.2), and stable non-trees with two or more
  lineages (sheet, two_groups palantir) also reach 1.0 (see "Phase 2 final review"). A redesign must be
  declared in `thresholds.toml` with reasons and a date before any run that measures it, calibrated on a
  fresh seed range (not 0-5), and verified on another; until then the served verdict stays provisional.
- Test fix round (after review): survivor tests for a weak control, the factor-1.0 tree-likeness, n_failed,
  the matrix share, a failed payload by kind, `complete`, and no cap note on a pending C2. The recorded
  run and the measured commit 37197c8 are unchanged.

## Phase 2 Task 11 decisions

- **R7, no cell cycle here.** Tasks 2-4 are blocked on a manual download, so `CellCycleCheck`,
  `CellCycleInfo` and `cell_cycle_check` do not exist yet. `RunTrust.cell_cycle`,
  `DatasetInfo.cell_cycle`, `Runs.cell_cycle`, `DatasetView.cell_cycle_info`, the run's cell-cycle
  warnings and `test_a_cycling_lineage_is_warned_on_the_run` are left to Task 4. The dataset test is
  `test_dataset_info_carries_tree_likeness`, without its cell-cycle half. Task 4 also adds
  `cell_cycle: null` to the web `dataset()` fixture.
- **The capped verdict is served as is (F14, ruling (a)).** Nothing in `src/` passes `verified=` to
  `verdict`; it reads `VERDICT_VERIFIED` (False). Tests that need the uncapped level patch
  `lumencell.engine.tree_support.VERDICT_VERIFIED` to True and say so in their docstring;
  `test_the_served_verdict_is_capped_while_the_verification_has_failed` checks, unpatched, that a
  finished run is provisional (or not_supported) with the failed-verification reason.
- **All taus are None, so tests that exercise C2-C4 patch the thresholds.** Under the calibrated
  thresholds C2, C3 and C4 are `excluded`, so the verdict is C1 and C5. The Review Focus 4 test
  (cancelled stability) gives C3/C4 a tau of 0.5 and lifts the cap, so stability really decides the
  level; a second test checks the served defaults never reach supported. The Review Focus 5 test
  (no PAGA edges) first checks C2 `excluded` with a null value under the served defaults, then gives
  C2 a tau and checks `fail`, null value, a reason and `not_supported`. Patched at
  `lumencell.trust.verdict.calibrated_thresholds`.
- **F1.** Stability is queued for every finished run, with or without clonal settings.
  `submit_trust_jobs` queues clonal scoring first, then stability; `get_run` looks up clonal first
  too, so a GET that queues both runs the cheap one first.
  `test_submit_trust_jobs_queues_clonal_only_with_settings` is now
  `test_submit_trust_jobs_queues_clonal_only_with_settings_and_stability_always`.
- **The payload** is `_job_payload()` plus `run_id` and `trust_name` only (Task 10: the handler
  reads terminals, clone columns and the label column from `result.json`). One builder,
  `Runs._stability_payload`, serves submit and GET healing; dedup key `stability:<rid>:<hash>`.
- **`clonal_evidence(c, n_lineages)`** takes the run's lineage count, not only `c` as the plan lists:
  a failed or refused agreement has no `primary` to read L from, and defaulting it would make C5
  `not_applicable` instead of the plan's `scored`/fail. The run's count equals `primary.n_lineages`
  whenever a primary exists. A non-ok scored status carries the agreement's first note into C5's reason.
- **A pending stability** is `not_measured("pending", "being scored in the background")`;
  `not_measured` now accepts `pending`. Legacy runs (`store_version` other than 2) are
  `not_available` with no job; the store does not serve them on GET at all, so that test calls
  `Runs.stability` directly, as the clonal legacy test does.
- **`DatasetView.tree_likeness()`** is computed from the loaded prep on each call (cheap; no cache,
  so it follows the loaded artifacts). `Session.dataset_info()` maps NaN to null.
- Extra tests beyond the plan: a failed stability job is served with its error and not resubmitted;
  a stability file that does not validate counts as missing with one stderr line.

Fix round 1 (PR #34 review):
- **R10.** A check left out (tau None) stays `excluded` and the level is unchanged, but its reason now
  also says when its evidence is missing: C3/C4 say they wait for the stability run (pending) or carry
  the stability note, e.g. "cancelled" (unavailable); C2 says it waits for the dataset's prep.
- **R11.** `DatasetView.tree_likeness()` is cached per loaded prep object, under the view's lock, so a
  run GET no longer recomputes it; a different prep object recomputes.
- **Tests.** Pinned the served mappings that mutation found unpinned: `stability_evidence` (done, counts,
  values, NaN, failed note), `clonal_evidence` (failed and ok), a NaN check value served as null, the
  no-edges GET's null C2 value and dataset tree-likeness, and submit skipping an existing stability file.

## Phase 2 Task 12 decisions and result

`bench real` (default `--stability`) now runs `stability_row` after each method; the run is
`benchmarks/results/real-2026-11-17.json` (measured on commit 99c3786, a clean tree, re-run for provenance; `make bench-real` exit 0,
2 min, 9 of 9 rows, stability `ok` on all, none failed).

- **Identity check: PASS.** Every gated metric and `n_clusters` is identical to `real-2026-11-13.json` for
  all 9 dataset/method rows (NaN equals NaN). `runtime_seconds` stops before stability (F17) and differs from 2026-11-13 because of a warm cache and
  machine load (not gated). C3, C4, tree-likeness and levels are unchanged from the first run.
- **Verdict cap (F14, ruling (a)).** `verdict(...)` is called without `verified=` and with
  `calibrated_thresholds()` (all taus None), so a supported verdict is capped at provisional
  because the verification failed on fresh seeds (2026-11-17), and C2-C4 are non-discriminating. A
  not_supported row is a real failure (bonemarrow paga_dpt: coverage 0.914 < 0.95) and is not capped
  (nor is weak since H-1; no recorded row is weak).
  The changelog and README say so.
  C3/C4 are still numbers. C2 is `tree_likeness` of the prep's PAGA at the default resolution; C5 is
  `ClonalEvidence("not_configured", n_lineages)`.
- **`stability_row`** turns any exception into `StabilityRow("failed", reason=...)` plus one stderr line.
  `row()` also carries `stability_reason`, so a NaN always has its reason beside it.
- **Report** cells gain `, verdict <level>`, and a footnote (only when some verdict was measured) says which verdicts are capped.

| dataset | method | C3 | C4 | tree-likeness | level |
|---|---|---|---|---|---|
| bonemarrow | paga_dpt | 0.200 | 0.998 | 0.703 | not_supported (coverage 0.914) |
| bonemarrow | palantir | 0.200 | 0.996 | 0.703 | provisional |
| bonemarrow | scfates_ppt | 0.857 | 0.952 | 0.703 | provisional |
| gastrulation_erythroid | paga_dpt | 0.333 | 0.999 | 0.625 | provisional |
| gastrulation_erythroid | palantir | 0.667 | 0.996 | 0.625 | provisional |
| gastrulation_erythroid | scfates_ppt | 1.000 | 0.998 | 0.625 | provisional |
| pancreas | paga_dpt | 0.200 | 0.997 | 0.666 | provisional |
| pancreas | palantir | 0.400 | 0.999 | 0.666 | provisional |
| pancreas | scfates_ppt | 0.800 | 0.988 | 0.666 | provisional |

C4 sits at 0.95-1.00 everywhere, as Task 8 found. C3 is low (0.2-0.4) for PAGA+DPT and Palantir on
bonemarrow and pancreas: their lineages rarely match at Jaccard 0.6 across replicates. Real gate failures
are unchanged from phase 1 and no gate changed.

## Phase 2 Task 2 decisions and result

Sources (all untrusted data kept outside the repo, never committed; scripts ran with `python -I`, the xlsx was read with
stdlib `zipfile` and `ElementTree`; SHA-256 recomputed and matched on 2026-11-19):
- Tirosh et al. 2016, doi:10.1126/science.aad0501, Table S5, PMC author manuscript NIHMS791837 supplement 5
  (`NIHMS791837-supplement-5.xlsx`), downloaded by Shreyan by hand on 2026-11-19, the primary source.
  SHA-256 `3b51070702338214bf8bb824a397dd1e04e8a57e2abefc442e15c1182922bb93`.
- scanpy/Seurat copy `regev_lab_cell_cycle_genes.txt` (cross-check), extracted 2026-11-17 from
  `https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip?dl=1` (scanpy's docs fetch it there; the GitHub path in
  the plan no longer exists). Zip `6557fe2a2761d1550d41edf61b26c6d5ec9b42b4933d56c3161164a595f145ab`, file
  `eb666ba4f7b7a8b2845e1586d42687feaf813e9fa1fdadba89dacd4093411ac0`.
- MGI `HOM_MouseHumanSequence.rpt`, `https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt`, 2026-11-17T13:51Z,
  `0b8eb746eb41eb9505b3bf3c5c8b986aba8a51e6a75df77642dcac34322aeb06`.
- HCOP `human_mouse_hcop_fifteen_column.txt.gz`, `https://storage.googleapis.com/public-download-files/hcop/human_mouse_hcop_fifteen_column.txt.gz`,
  2026-11-17T13:51Z, `0cfb78e4eb273751e0557e8230445a60e844bfaf958e32740d4a266fb41cdb35` (current HGNC symbols only).
- HGNC `hgnc_complete_set.txt`, `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`,
  2026-11-19T19:23Z, `6df77ac3cbaeb86a6dce5e17d851ebb27b6b3cdcd2471a211a3d8464f00c54c4` (used for `prev_symbol`).

Lists verification: the paper's G1/S column has 43 cells (43 unique); its G2/M column has 55 cells and 54 unique genes, because
HJURP is entered twice (cells B35 and B36). Ruling R18: the plan's counts are gene counts, so the test stays at 43 and 54 and the
duplicate cell is recorded, not counted. The S list matches the scanpy copy in content and order; the G2M list matches in content,
and in order once the duplicate is dropped. No gene is in both lists. A test asserts no symbol is repeated in the shipped TSVs.

Renames since 2016, found by Approved-entry `prev_symbol` (alias_symbol was not needed): MLF1IP -> CENPU (S), FAM64A -> PIMREG
(G2M), HN1 -> JPT1 (G2M). All three are in `current_symbol`. No symbol maps to more than one gene or to none.

Orthologs verification: all 97 genes (by current HGNC symbol) sit in an MGI class with exactly one human and one mouse gene, so 97
are mapped and 0 unmapped. HCOP cross-check: 97 agree, 0 disagree, 0 absent.

Package data: hatchling includes the TSVs in the wheel (checked with `uv build --wheel`); `importlib.resources` reads them.

- **Fix round 1 (after review).** `_rows` now checks the header row and the column count of every
  row and strips fields, so a missing header or a short row is refused naming the line instead of
  silently dropping or misreading a row. Tests pin the renamed genes' mouse forms (Cenpu, Pimreg, Jpt1),
  `unmapped == ()`, an unknown species raising, and the first var index of any accepted symbol winning.
  An independent re-derivation from the five sources matched both TSVs row for row.

## Phase 2 Task 4 decisions

- **Scope.** The plan's Task 4 plus, by ruling R7, the cell-cycle parts Task 11 deferred: `RunTrust.cell_cycle` and `DatasetInfo.cell_cycle`
  (no defaults), `DatasetView.cell_cycle_info()`, `Runs.cell_cycle(rid, names)`, the run's cell-cycle warnings appended to
  `TrajectoryRun.warnings`, `Session.dataset_info()`, `make contract`, and `cell_cycle: null` in the web `dataset()` fixture.
- **`None` means prep is not loaded.** A prep whose stage was not scored gives a `CellCycleCheck` with `status="not_scored"` and the reason as `note`.
- **Names.** The check takes the lineage names after renames (list index = lineage id). They reach only the served notice text, never a path or a log line.
- **`CellCycleInfo.genes_found` and `genes_listed`** sum `cc_info`'s per-list dicts; `source` is `""` if the stage recorded none.
- **Doc fixes from the final review.** The README/report sentence now says C2 separates blob and sheet from the matrix but not loop or two_groups at x1.0,
  citing the recorded calibration (README regenerated; no recorded result changed). The phase 2 security review was done by hand.

- **Fix round 1 (ruling R20).** An unreadable run cell file no longer gives `cell_cycle: null` (which the API
  documents as "prep is not loaded"): `Runs.cell_cycle` serves `not_scored` with the note "the run's cell file
  could not be read (<ExcType>)" and one stderr line; any non-OSError still propagates. Test:
  `test_an_unreadable_cell_file_is_reported_not_hidden`.

## Phase 2 Task 3 decisions

- **Stage.** `cell_cycle` is last in `STAGES` with `report_only=True`; `prep_key` hashes only stages with `report_only=False`, and
  `prep_key(fingerprint "f"*24, PrepParams())` is still `ee3f123928ac0d1c` (checked on 3f5ab9a before and after). `stage_keys`, `status`, `run`
  and `load_all` cover it, so a cache without it is not ready and one prepare job adds it.
- **F9.** scanpy returns float64; `cc_s` and `cc_g2m` are cast to float32. When not scored they are all NaN (float32) and `cc_phase` is all `"unscored"`.
- **Check order** is the binding one in the plan: gene count first (by either species' symbols for `auto`, a gene counted once), then species,
  then the resolved lists. `total` in "only {k} of {total}" is the number of genes listed (97 = 43 + 54), so synthetic `G00000` symbols read "0 of 97".
- **Never fails prep.** Any `Exception` while checking or scoring becomes `not_scored` with `scoring failed: {Type}: {msg}` and one stderr line.
  Duplicate symbols are scored once: `match_genes` takes the first, and the view's var names go through `make_index_unique`.
- **Renamed genes.** In the human test the paper's symbols (for example MLF1IP) are used; Task 2's two accepted symbols make them match.
- **F8 (accepted).** The stage key covers only `cc_species` and `cc_gene_set`, so a later edit to the list or seed leaves a stale, report-only
  result until `PIPELINE_VERSION` moves.
- **F22.** Every existing prep cache lacks `cell_cycle`, so the first server start re-runs that one cheap stage (`test_a_cache_from_before_the_cell_cycle_stage_is_completed_by_one_prepare`).
- Two existing tests were adjusted to the new last stage: `test_existing_embedding_skips_umap` looks `umap` up by name, and the receptor-exclusion
  key test exempts `cell_cycle`, which reads only `lognorm`.
- **Fix round 1.** The stage setup (`gene_symbols`, both `load_gene_set` calls, the default info) now sits inside the `try`, so a broken
  shipped TSV gives `not_scored` with the error named and one stderr line instead of failing prep. Pinned by tests: the 10 S + 10 G2M versus
  9 S + 10 G2M boundary ("only 19 of 97"), an empty phase list (25 S, 0 G2M: the gene-count reason, no stderr), `random_state == CELL_CYCLE_SEED`,
  the never-raise test over four exception types, and the CLI "scored as human from 97 of 97" line.

## Phase 2 Task 1 decisions

- **F14 declared.** `[tree_support.verification] on_failure = "provisional"` sits after
  `max_controls_supported`, dated 2026-11-17. It has no runtime constant; Task 9 adds `VERDICT_VERIFIED`.
- **`test_every_scalar_clonal_key_is_mirrored` narrowed.** It took every two-part `MIRRORED` path as a
  clonal key; it now counts only paths starting with `clonal`, since phase 2 adds `cell_cycle`,
  `stability` and `tree_support` paths.

## Phase 2 Task 5 decisions

- **Rename only.** `bench/synth.py`'s `_sample_counts` is now `sample_counts`; the phase 1 digest test passes unchanged.
- **RNG.** One generator per (kind, seed): `default_rng([CONTROL_STREAM, KINDS.index(kind), seed])`.
- **Gene sets.** Programme genes are Bernoulli draws per gene (as the tree generator's `de_prob`), so the
  count varies with the seed. `sheet` draws one uniform per gene: below the share is an x gene, the next
  band of the same width is a y gene, so the sets are disjoint. Fold changes are signed log-normal with
  sigma 0.5 around the scale. `loop` magnitudes are `|lfc|`, the sign carried by the phase.
- **`two_groups` root** is the first cell (lowest index) of group 0.
- **Tests stop prep at the kNN graph.** At 1000 x 500, `two_groups` gives exactly two Leiden clusters with
  no edge between them, and the PAGA stage raises (igraph builds an empty adjacency). At 3000 x 2000 all
  four kinds pass full `run_pipeline` for seeds 0 and 1. Tasks 8 and 9 run at that size; if a control
  ever fails prep there, that is an engine finding, not a control to change.

## Phase 2 Task 10 decisions

- **Both inputs come from the run.** `result.json` now records `terminals` (as sent, `null` for
  auto) and `clonal_columns` (sorted; `[]` when none), written by `RunStore.save(..., clonal_columns=)`
  and passed by `handle_infer`. `handle_stability` reads them from the run, never from the payload:
  a `clonal_columns` key in the payload is ignored, so Task 11 need not send it. The payload needs
  `path`, `prep_params`, `run_id`, `trust_name` and `label_key` (may be `None`); a missing key fails
  the job naming it. A `trust_name` that is not `stability_name()` fails it and writes nothing.
- **Old runs** (no `clonal_columns` key, made before this PR) are written as `not_available`
  ("the clone columns this run used were not recorded; run it again"), not replayed with `()`:
  an empty set could let clone labels into a replicate. Plan cases kept: store version other than 2
  is `not_available`; another prep key is `failed` (stored, unlike `score_clonal`'s raise, as the
  plan says); `terminals_source == "user"` without `terminals` is `not_available`; `StabilityError` is `failed`.
- **Models.** `Stability` gains `n_failed: int` (controller ruling). `StabilitySummary` float fields
  are `float | None`, not `float` as the plan lists them: an empty kind has NaN medians, and Review
  Focus 2 needs them null. `make contract` shows no diff.
- **Trust hash** covers `STABILITY_METRIC_VERSION` and the five `[stability]` params, imported from
  `trust/constants.py`. A run's stored `failed` for other prep settings is permanent per name, as the plan wrote.
- A run that is gone or incomplete (at read or at `write_trust`) fails the job, as in `score_clonal`.
- The cancel test uses 3000 cells (F16).
- **Fix round 1.** `label_key` is recorded in `result.json` (the `used_label_key` from `run_method`,
  carried on `TrajectoryResult.label_key`, may be `null`; read with `.get`, so run ids and `STORE_VERSION`
  are unchanged). `handle_stability` replays it from the run and no longer takes it in the payload; a run
  with no `label_key` key is `not_available` ("the label column this run used was not recorded").
  `api/runs.py` is not changed: its `label_key` is the *infer* job's payload (on main), which the infer worker still reads.
- **Fix round 1.** `terminals` is stored as `terminals or None`, so `[]` reads as auto.
- **Fix round 1.** A stale prep key now RAISES and writes nothing, like `score_clonal` (a stored `failed`
  would block measuring once the session reopens with the matching prep). This replaces the stored `failed` above.
- **Fix round 1.** Tests close survivors: every `[stability]` parameter changes the name, a wrong
  `trust_name`, an incomplete run, user terminals replayed in order, and `run_stability` never called on a `not_available` run.

## Phase 2 Task 6 decisions

- **A method's default seed is reseeded.** `replicate_plan` reads `"seed"` from the params it is
  given (its test: `{}` gives `None`). `run_stability` hands it the run's params with the method's
  defaults filled in (`Params.model_validate(...).model_dump()`), so a Palantir run that left its
  seed at the default 20 still gets 21 and 22, as `[stability]` asks ("when the method has a seed
  parameter"). `paga_dpt` has none, so its reseeds vary only Leiden and PAGA (F11). The method is
  called with the run's own params plus `seed`, never with the filled-in defaults.
- **Every Leiden terminal dropped fails the replicate.** When the run named Leiden terminals and
  none has cells left, passing `[]` would let the method pick its own leaves, a different question.
  The replicate is recorded with a `StabilityError` and counts as 0, like any failure. Terminals
  in a user column pass unchanged; a label with no cells left makes the method raise, which is
  the same recorded failure (Review Focus 3).
- **A root name absent from the file** is refused with `StabilityError` too, next to the
  repeated-name refusal the plan names.
- **`restrict()` fingerprint** is `subset-<sha256 of the full fingerprint and the cell indices>`.
  `uns_root` maps to its position in the subset, or `None` when the subset lacks it (F6).
- **Subsample prep reuses** `hvg`, `pca_variance_ratio` and `pca[cells]` only; `hvg_info`, `lognorm`
  and `umap` are left out (nothing downstream of PCA reads them). Reseed prep keeps every artifact
  except `leiden` and the three PAGA ones, which are rerun.
- **Padding.** Cells outside a subsample get NaN pseudotime, progress, fate and extras, no
  membership, segment -1 and empty milestone-weight and curve-node rows; `pad_result` calls
  `validate()` before returning.
- **`jaccard_matrix` keeps the old double loop**, so `compare(a, b)` gives bit-identical values;
  `test_without_a_mask_compare_is_unchanged` checks it against a frozen copy of the old code.
- **Identical-run Spearman is 1 to rounding.** `spearmanr` of a vector with itself gives
  0.9999999999999999 here, so `test_a_replicate_identical_to_the_reference_scores_one` compares
  Spearman with `approx(1, abs=1e-12)` and the other four fields exactly.
- **Empty kinds** (`n_reseeds=0`, say) summarise to NaN medians and minimums and `n_lineages_min` 0.
- **Fix round 1:**
  - `StabilityResult.n_failed` counts replicates that raised (each still counts as 0); the recorded
    error is `"<ExceptionType>: <message>"`. Tested for one failure and for all five.
  - The Leiden-terminal test reads the method's default `cluster_key` from its Params model
    (`_terminal_key_default`), not a hard-coded `LEIDEN_KEY`.
  - `_method_params` carries a comment: a bad param aborts the whole job, as the reference run would have.
  - New tests close surviving mutants: subsample draw seed `seed + 1 + j`, `share_matched` over the
    reference lineage count (2 vs 3 lineages gives 0.5, `n_lineages` 3), Jaccard equal to the
    threshold is matched, `pad_result` leaves progress/fate/extras NaN, membership False, segment -1
    outside, and `compare()` counts cells finite in both (it already did: no bug found).
- On the 600-cell bifurcation, `paga_dpt` scores lineage stability 1.0 and pseudotime stability
  0.99 (one subsample finds 2 of 3 lineages). Not a measurement for any τ.

## Phase 2 Task 7 decisions

- **C5 and `single_lineage` (F12, spec wins).** `not_applicable` only when clones are not configured.
  A scored `single_lineage` fails C5 with the status in its reason; the plan's step 1 text that makes it
  not applicable is superseded. Test: `test_a_single_lineage_clonal_score_fails_c5`. L < 2 reaches the
  verdict as `not_configured` (the caller decides), since a scored result always has L >= 2 or `single_lineage`.
- **R2, PAGA with no edges.** `prep/stages.py::_paga` returns an all-zero (k, k) connectivity and tree and
  circle positions when no kNN edge joins two Leiden clusters; scanpy raised there. A single cluster did
  the same before and now takes this path too. Graphs with edges are untouched. `engine/topology.py`
  `_paga` (non-Leiden grouping) still calls scanpy and is not changed.
- **One place decides the level:** `_level` in `engine/tree_support.py`, so Task 9's cap on a failed
  verification goes there. It is not implemented.
- **All replicates failed.** `StabilityEvidence` takes `n_replicates` and `n_failed` (default 0); with
  `n_failed >= n_replicates > 0`, C3 and C4 are `unavailable` and the level is `provisional`. Partial
  failures stay in the medians as zero. To be fed from `StabilityResult.n_failed` once the Task 6 fix lands.
- **Level order.** C1/C2 fail beats everything; C1/C2 unknown (C2 `None`) gives `provisional`; then a C3-C5
  fail gives `weak` even while others are pending; then pending/unavailable gives `provisional`.
  An `excluded` or `not_applicable` check never changes the level. A NaN C3/C4/C5 value fails.
  C2 `excluded` (tau_tree `None`) wins over pending.
- **Fix round 1.** (1) F12: `ClonalEvidence.n_lineages` is required (no default), set by `from_result`;
  `_c5` is `not_applicable` with a reason when scored with `n_lineages < 2` (pending is checked first, so a
  pending caller passes any value); L >= 2 `single_lineage` still fails.
- **Fix round 1.** (2) Survivor tests: exact-threshold passes (C2-C4), max-symmetrisation of PAGA
  connectivity (upper-only, lower-only tree), excluded C2 with `None`, and not_supported over other fails.
- **Fix round 1.** (3) Removed the dead `unavailable` branch from the C1/C2 level handling.
- **Fix round 1.** (4) Comment in `engine/topology.py::_paga`: non-Leiden branch still raises on an edgeless graph.

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
  `clonal_columns` is part of the run id when any are named (see Task 10 decisions).
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

## Task 14 decisions

- `TrajectoryRun.trust: RunTrust | None` (no default; always filled for a done run).
  `JobKind` gains `score_clonal`, and `JobInfo.kind` is `JobKind | str`: a kind the contract does
  not list (`sleep`) is reported by its own name, no longer as `prepare`.
- `Runs.clonal_agreement`: no settings gives `not_configured`; a store version other than 2 gives
  `refused` and queues nothing; a stored file gives its `ClonalAgreement`; else `pending` while the
  job is queued or running, `failed` with the job's error (or `cancelled`) after it ends badly.
  It is looked up by the current settings' `trust_name` only, so a superseded job's file, which
  still gets written (F21, parked), is never served.
- `Runs._trust_lookup` (phase 2 reuses it): the file, else `jobs.find(f"{kind}:{rid}:{h}")`, else a
  new job. A done job whose file is missing is resubmitted (after one re-read, as the job may have
  written it since). A failed or cancelled job is not resubmitted by a GET; a settings PUT or a
  restart retries it. It takes an optional `validate`: a file that parses but does not validate
  counts as missing, with one stderr line naming it (F13).
- Background submission: `Session._on_job_done` (a done-callback, listener thread) calls
  `submit_trust_jobs(rid)` for an `infer` job on this session's file; `set_clonal` ends with
  `submit_clonal_all()`. `clear_clonal` submits nothing (controller ruling). Dedup is the job key
  `score_clonal:<rid>:<trust_hash>`, decided under the manager's lock.
- No deadlock: none of the trust code takes a lock of its own. The manager's lock is taken per call
  (`find`, `submit`) and never held while another lock is taken; the state-file lock is released
  before `set_clonal` submits; `Runs._lock` only guards the in-memory caches. The callback never
  waits.
- `Runs` takes a seventh argument, `clonal` (a getter of the current settings).

## Task 13 decisions

- `trust/clonal.py`: `agreement_for_run` (own-cells mode, m = 2 primary and m = 3 sensitivity, via
  `score_with_baseline`), `to_api_score` and `to_api_coupling` (NaN to `None`), `not_scored`, and
  `trust_hash`/`trust_name` (`clonal-<16 hex>`). The hash covers the settings (strata sorted), the
  metric version and a params dict built from `trust/constants.py`, never hand-typed (F20).
- `api/models.py`: `NullComparison`, `ClonalScore` (with its own `status`, F27), `LineageCoupling`,
  `ClonalStatus`, `ClonalAgreement`, `RunTrust`. No route uses them yet; `make contract` has no diff.
- Notes: `COUPLING_NOTE` only on an agreement that carries a coupling (controller ruling, fix
  round 1); a reason note for `refused`, `failed` and both `single_lineage` cases (one lineage, or
  scored clones in one lineage); the missing-donor count; the too-few-clones count from the primary.
  None claims direction. Non-finite values (NaN, infinity) map to null; `p_delta` is null whenever
  `delta_eta` is.
- Settings that `check_settings` rejects against the file give a `failed` agreement naming the
  problems, never a bare `KeyError`.
- `handle_score_clonal`: a payload missing a key fails naming it (`{}` names `'path'`). A missing or
  incomplete run, at read or at `write_trust`, fails the job with
  `run <rid> was replaced or is incomplete; not scored`; nothing is written elsewhere. Store version
  other than 2 writes `refused`; another prep key fails the job (`RuntimeError`) and writes nothing.
  This departs from plan l.808, which wrote a permanent `failed` file: `recompute_umap` is part of
  the prep key (`stages.py:325`), so the session can reopen with the matching prep, and a saved
  `failed` file would then block rescoring. A payload `trust_name` that is
  not `trust_name(settings)` fails the job and writes nothing.
- `JobKind` in the API still lacks `score_clonal`; Task 14 adds it with the route.

## Task 12 decisions

- `BACKGROUND_KINDS` (`score_clonal`, `stability`) queue in `_background`; every other kind goes to
  `_foreground`. `_take_next` (lock held) takes the oldest foreground job, else the oldest background
  one. A running job is never pre-empted.
- `add_done_callback`: runs once per job that ends `done`, on the listener thread, after the lock is
  released. Callbacks must be quick and may only submit jobs: while one runs the listener reads no
  worker events, so it delays the next dispatch. An exception is printed to stderr with its traceback.
- Callbacks are never removed, so a superseded scoring job still runs (ruling F21, parked).
- `score_clonal` had no worker handler until Task 13.

## Task 11 decisions

- `<run>/trust/<name>.json`, `name` fully matching `TRUST_NAME`; `write_trust` needs a complete run
  (`.complete`), writes via `mkstemp` and `os.replace`, and refuses NaN and infinity.
- `read_trust` gives `None` for a missing run, missing file, unreadable file or non-object file
  (the last two with one stderr line). Task 14 still validates content.
- `save` of an existing run moves its `trust/` folder into the new folder before the swap, and
  moves it back if the swap fails. `cache prune` removes `trust/` temp files older than the grace period.
- Task 13: only the worker that owns a run writes its trust, and it treats `FileNotFoundError` from
  `write_trust` as "run replaced": the job fails with the run-replaced message and nothing is written.

## Task 10 decisions

- Settings live in the state file under `clonal`; `Session` validates them once at start and keeps
  them (an invalid saved value is ignored with one stderr line). The clone and donor columns become
  scoring-only; strata columns never do.
- `run_id` hashes the sorted clonal columns only when there are any, so every earlier run id is
  unchanged; naming columns changes the id because it can change the detected label key.
- `api/runs.py` takes the scoring-only set from `view.ctx().scoring_only` (needs prep ready).
- `column_kind` is the one rule for listing and serving obs columns: above 1,000 categories a
  column is `string` in both places. Counts of clones and donors use values present, never unused
  categories.
- Not done here: the web UI for these settings (a later task).

## Task 9 result

`make bench-real` exit 0 at 4b0b1d5 (3 real_data tests passed); record
`benchmarks/results/real-2026-11-13.json`. Real rows are `clonal_status = not_applicable`.
- `hvg_info` receptor_genes / receptor_hvg_without_exclusion: pancreas 0 / 0, bonemarrow 11 / 6,
  gastrulation_erythroid 258 / 3.
- Moved: pancreas nothing. Bone marrow: paga_dpt branch_accuracy 0.917 -> 0.617 and lineage_recall
  1.0 -> 0.5 (gates now fail), scfates_ppt pt_spearman_lineage 0.781 -> 0.668 (fails), palantir
  branch_accuracy 0.639 -> 0.943 (improves). Gastrulation: paga_dpt lineage_recall 0.0 -> 1.0, a
  threshold flip, not a gain.
  Eight gates fail in total, as the log lists. No threshold changed.
- **Attribution: the exclusion alone causes every move.** Exclusion off at HEAD reproduces the
  11-10 file to 3 dp, on reproduces the new file (reviewer's check); cold-cache bone marrow
  paga_dpt at 744a0a8 vs 7b6254b is 0.9171 -> 0.6171, 12 -> 13 clusters, identical repeats. Only
  7b6254b changed the expression path since abe3f73. The would-be HVGs in bone marrow are TRAC,
  TRDC, TRGC1, TRGC2, TRBC1, TRBC2.
- **Fragility:** swapping 6 (bone marrow) or 3 (gastrulation) of 2000 HVGs flips real-data gates;
  in bone marrow Leiden goes 12 -> 13 and one lineage fails to match. No bug found; failing gates
  stay failing; thresholds unchanged.
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

- **Some existing gates are pinned only by the benchmark runs, not by unit tests** (found by PR #51's mutation seat, 2026-11-19): changing `[regression] max_drop`, `pt_spearman_global` or the pancreas `branch_accuracy` in `thresholds.toml` passes the unit tests. Phase 3 changes none of them; noted for the phase 3 security review.
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
  Task 6 measures them. Task 7 ran the slow golden tests (5 passed) and the Done When section
  records `make bench-synth` (exit 0, 0 regressions).
- **Gate failures don't fail the benchmark commands.** `lumencell bench real` and
  `bench synth` exit 0 on gate failures; only regressions fail `bench-synth`.

## Phase 2 final review (2026-11-17)

One Opus reviewer read the spec, plan, every phase 2 section and recorded file and the code diff
`2217f76..1167b2d`. Full report: the controller's notes folder, `phase2-final-review-2026-11-17.md`
(outside the repo). No blocker. Declare-before-measure holds for every recorded file (N2a 2b15c77 after
b0e235a; calibration 0a584d5 after 59d8213; verification 37197c8 after the taus in 0d0608a; real run
99c3786 after both). The fast suite on 1167b2d passes (1009).
- **H-1, closed 2026-11-19 (see "Phase 2 H-1 decision").** `engine/tree_support.py` `_level` capped `weak` at
  `provisional` as well as `supported`. A failed C5 (clones that do not agree, p 0.5) or a scored `single_lineage` is then
  served as provisional. C5's rule passed its own calibration (N2a); the verification failed because the
  rule is too lenient, which says nothing against a weak verdict. The reviewer proposes keeping `weak`
  when a check fails, as `not_supported` is kept. Amended F14 and ruling (a), so it went to Shreyan, who ruled to keep `weak`.
- **X-1, phase 3 plan amendment needed.** It calls `clonal_evidence(agreement)` (now `(c, n_lineages)`),
  and its slow test and case-study text expect supported/weak levels the cap makes impossible.
- **X-2, X-4.** `verify` passes no `n_replicates`/`n_failed` and `RunRecord` does not record replicate
  failures, so the "every replicate failed" path the API serves is never measured; real rows carry no
  `n_failed` or reference lineage count, so "none failed" in Task 12 cannot be checked from the file.
- **Minor:** the README sentence (from `bench/report.py`) says C2 does not separate trees from sheets or blobs,
  but at x1.0 it does separate blob (0.19-0.23) and sheet (<= 0.66) from the matrix (>= 0.816); loop and
  two_groups beat it. The cap note is missing while C5 is pending. Coverage has two equivalent formulas.
  "80%" and `1.0 * factor` are hard-coded. The served verdict carries no rule version. The security
  section does not say the review was by hand (the skill was not invoked).
- **Lessons (records, seeds 0-5).** C3: one-lineage runs score 1.0 trivially; even with L >= 2, sheet and
  two_groups palantir reach 1.0, so stability measures reproducibility, not tree shape. C4 separates only
  blob (max 0.94 vs matrix min 0.993); paga_dpt reseeds reuse the diffusion map, so 2 of 5 replicates
  are copies. C1: palantir gives finite pseudotime across two_groups' disconnected parts.
- **Candidate redesign (to pre-declare, not adopted):** C1' share of cells in the root's kNN component;
  C2a the current ratio with a 3-cluster floor (blob, sheet); C2b the longest cycle an off-tree PAGA edge
  closes (loop); C3' symmetric branch-only Jaccard with L >= 2 plus a "random terminals" control; C4'
  calibrated against blob only. Seeds never used: 6-99 and 500 up. Proposed: 500-509 calibrate,
  510-519 verify, 900-904 a declared, burned pilot range.
- **Follow-up tasks (not phase 2):** (1) the C1-C4 redesign, folding in X-2, X-4 and a
  `STABILITY_METRIC_VERSION` bump; (2) a clear engine error for a disconnected kNN graph (palantir
  "fate probabilities must sum to 1", scFates IndexError on `two_groups`) before phase 3's T-cell data;
  (3) the phase 3 plan amendment (X-1, done on `pivot/38-phase3-plan-x1`); (4) Shreyan's ruling on H-1; (5) the minor doc fixes, in Task 4's PR.

## Phase 2 H-1 decision (Shreyan, 2026-11-19)

- **Ruling.** Keep `weak`. While `VERDICT_VERIFIED = False`, only a supported verdict is capped at
  provisional (with the dated failed-verification note). A run whose check fails is served `weak` with its
  failing checks' reasons and no cap note, as `not_supported` already was. Unmeasured or pending runs stay
  `provisional`. It amends ruling (a) of 2026-11-17 and F14.
- **Changed.** `_level` in `engine/tree_support.py`; the unit test that pinned weak -> provisional now
  checks weak is kept, and a new test serves a failed C5 (p 0.5) as weak under the unverified constant;
  `TreeSupport` and constant comments; the bench report footnote and README; a dated line in
  `benchmarks/thresholds.toml`. PR `pivot/34-weak-under-cap`.
- **Recorded results unchanged.** `real-2026-11-17.json` has no weak row (one not_supported, eight
  provisional): with every tau None and no clones, nothing in C3-C5 can fail. Nothing was re-run. No rule,
  tau or seed changed.

## Disconnected kNN graph error (phase 2 follow-up)

- **Reproduction** (`two_groups`, seed 0, root = `control_root_cell`). The kNN graph splits into two parts at
  1000 x 500 (470 + 530 cells), 1500 x 1000 (743 + 757) and 3000 x 2000 (1514 + 1486); full prep passes at
  all three (Task 5's note that PAGA raises at 1000 x 500 no longer held here). scFates `scfates_ppt` raises
  `IndexError` (in `scf.tl.pseudotime`: the root tip's milestone is missing) at all three sizes. Palantir
  fails only sometimes: `InvariantError: fate probabilities must sum to 1` at 1000 x 500, but the same call
  can pass on a second try, and at the other sizes it ran and gave finite pseudotime across both parts (the
  C1 finding). `paga_dpt` runs at all sizes (cells in the other part get no pseudotime and a
  `clusters_cut_off` notice where PAGA has the cluster). Palantir's own graph (`DM_Kernel`, knn 30) has the
  same two parts as the prep graph.
- **Rule.** `engine/graph.py::require_connected` raises `DisconnectedGraphError` (a `ValueError`) when the
  graph has more than one part. `Palantir.run` checks the kernel it just built (the graph it uses), after the
  diffusion maps and before the Markov chain; `ScFates.run` checks the prep kNN graph before fitting the
  tree. The message gives the number of parts, the cells in the root's part and in the largest other part,
  and says to subset to one part or raise the neighbours. It names counts, never obs names.
- **Why these two, not paga_dpt.** Palantir and scFates either crash internally or return pseudotime across
  parts they cannot order; paga_dpt handles the split by design (unreachable clusters get no pseudotime and
  a notice), so its behaviour and recorded results are unchanged. Nothing in a threshold, control or result
  changed.
- **How it reaches users.** CLI: `run_method` errors that are `ValueError` already become a
  `BadParameter`. API: the run is a job, so it ends `failed` with the error text
  (`lumencell.engine.graph.DisconnectedGraphError: ...`), listed under the run's `error`; no 500.
- **Tests.** `tests/unit/test_disconnected_graph.py` (the helper, both methods on `two_groups` 1000 x 500,
  paga_dpt still runs) and `tests/api/test_disconnected_graph.py` (job failure and message through the API).
  Both failed before the change with the internal errors above.
- **R19 (controller ruling): a stored explicit zero is not an edge.** `require_connected` drops zeros on a
  copy before counting parts and never changes the caller's matrix. Connectivity is weak (`directed=False`):
  a one-way chain is one part. The root's part is left out of "largest other part". Pinned by
  `tests/unit/test_disconnected_graph.py` (zero-join, one-way chain, 4/2/3 blocks, `ValueError` base, a CLI
  run that exits 2 with the message and no traceback, and Palantir reporting the root's own part).
- **What a re-run would change.** Re-running the recorded verdict calibration/verification now leaves out
  every `two_groups` palantir and scfates_ppt run, because they raise `DisconnectedGraphError`. For example
  palantir seed 5, recorded as a supported control in the 28 of 31, would raise. So the observation that
  "palantir gives finite pseudotime across two_groups' parts" (C1) describes behaviour before #40. The
  recorded files stay as they are.

## Phase 3 plan amendment (X-1)

Docs only, 2026-11-19, checked against `main` at 3f5ab9a (`trust/verdict.py`, `engine/tree_support.py`,
`bench/runner.py` `stability_row`, `trust/stability.py`, `api/models.py`). Phase 3's goals, datasets, gates
and download steps are unchanged. Each change in
`docs/superpowers/plans/2026-11-12-phase-3-tcell-data.md` is marked "(amended 2026-11-19, phase 2 X-1)" and
listed in the plan's "Amendments" block.
- **Depends on:** the verdict as shipped (verification failed, `VERDICT_VERIFIED = False`; after H-1 only
  `supported` is capped; all taus None, so C2-C4 are `excluded`; best level without clones is provisional,
  a failed C5 with clones gives weak).
- **Interfaces:** `verdict(..., verified=None)` and `ClonalEvidence` named; `clonal_evidence(c, n_lineages)`
  and `to_api_tree_support(v)`; `cell_cycle_check` returns `models.CellCycleCheck` and matches phase 2
  Task 4 as specified, but `trust/cell_cycle.py` is not on `main` yet, so "as Task 4 lands".
- **Preflight:** asserts `clonal_evidence`'s parameters are `['c', 'n_lineages']`; the `cell_cycle` import
  passes only after Task 4 merges.
- **Out of scope (C5 in benchmark rows):** a benchmark row is `provisional` or `not_supported`.
- **Task 8 Step 6:** Daniel levels are never `supported`; the report says why.
- **Task 12:** the verdict call is `clonal_evidence(agreement, len(result.lineages))`; skipping stability
  does not change the level; the levels a run can have are listed; `render_markdown` gets a fixed verdict
  footnote; the slow test that allowed only `supported`/`weak`/`not_supported` now uses
  `_assert_capped_verdict` (provisional, weak or not_supported, tied to the checks); a fast cap test is
  added (with a patched-True case so the cap is what excludes `supported`); Step 6 reads the levels.
- **Task 13:** the text states the failed verification and the cap; `weak` and `not_supported` still lead
  a paragraph; `provisional` is "not shown to be a tree"; a test checks the statement and that no unit
  is `supported`.

## Phase 3 Task 1 decisions

- **Phase 3 base:** `main` at 35ac459. Branch `pivot/39-registry-sources`.
- **Where:** `bench/download.py` (new) holds `Source`, `Build`, `fetch_source`, `build_dataset`, `build_key`,
  `local_path`, `is_cached`, `fetch`; `bench/datasets.py` re-exports them. `RealDataset` takes
  `sources`, `licence`, `build`. The three cached datasets keep their paths and recorded SHA-256s.
- **Licences (read 2026-11-19, hosting pages only, no data file downloaded):**
  - pancreas: no licence stated on https://github.com/theislab/scvelo_notebooks (no LICENSE file, badge or README
    mention; the `data/Pancreas` folder page was not opened). Recorded as "no licence stated"; cite Bastidas-Ponce 2019.
  - bonemarrow: CC BY 4.0 per the figshare API record of article 14465559 (file 27686835,
    https://figshare.com/articles/dataset/human_cd34_bone_marrow_h5ad/14465559).
  - gastrulation_erythroid: CC BY 4.0 per the figshare API record of article 14465562 (file 27686871,
    https://figshare.com/articles/dataset/Gastrulation_erythroid_lineage/14465562).
  - The article pages themselves returned 403 to the web fetcher; the figshare API (`api.figshare.com/v2/articles/<id>`)
    metadata gave the licence and confirmed each file id and byte size.
- **Decisions on small ambiguities:**
  - `fetch_source` and `build_dataset` use `os.replace` (atomic, overwrites) rather than `Path.rename`.
  - `build_dataset` checks each source for: downloaded, recorded checksum, matching checksum, in that order.
    A changed-bytes source at the same size passes `is_downloaded` and is caught by the checksum.
  - `build_dataset` does not create `datasets_dir()`; sources exist there already when it runs.
  - CLI: with every source present and only the derived file missing, it builds with no question.
    Otherwise it lists only the missing sources and the total, asks once, prints each new file's SHA-256, then builds.
  - `Fetchable` is a Protocol with read-only properties so `RealDataset` (frozen) satisfies it.
  - The CLI test reuses the `TwoFile` helper from `tests/unit/test_datasets.py`.

### Fix round 1

- **P3 redirects:** downloads use `build_opener(HTTPSHandler(context), _HttpsOnlyRedirect)`; any redirect whose target
  is not https raises `FetchError`. Host changes are allowed (figshare to S3). Tested on the handler directly.
- **P4 empty build:** a builder that leaves the temporary output at 0 bytes raises `FetchError("... empty file")`;
  no derived file is left.
- **P5 one directory per download:** `fetch_source` makes a fresh `mkdtemp` directory (`.download-*`) under the datasets
  directory for each source, writes and verifies the file there, `os.replace`s it into place, and always removes the
  directory. Two downloads never share a directory.
- **Unverified files:** `fetch_source` still accepts `sha256=None`; the CLI line then adds "unverified until its
  SHA-256 is recorded in the registry". The total line now also gives bytes.
- **`is_cached` with the derived file present but a source missing:** the plan says every source downloaded AND the
  derived file exists, so it is False; implemented as written. `bench fetch` then re-downloads the missing source only.
- **Network seam:** `download._open(request)` is the one place a connection is made; tests patch it, and a separate test
  pins the TLS context (`CERT_REQUIRED`, `cafile=certifi.where()`).


## Phase 3 Task 2 decisions and result

- **Branch** `pivot/40-daniel-fetch`. Owner OK for this STOP step: Shreyan, in chat on 2026-11-19, "Yes, approve all" (relayed by the parent session): Daniel 2022, GEO GSE188666, nine files; not the 14.7 GB GSE188670_RAW.tar.
- **HEAD sizes** all 200; total 1,168,794,468 bytes (1.169 GB, inside 1.0-1.2 GB). The plan's 574 and 519 "MB" for the two mtx files are MiB (602,240,886 and 544,512,106 bytes).
- **SHA-256 as recorded on first download (2026-11-19):**
  - GSE188666_scRNA_LCMV.mtx.gz: be4a6e2825d9cdffbb4bebb094e43be00ed388faa11a5bac95c2418bdd576ed6
  - GSE188666_scRNA_LCMV_barcodes.tsv.gz: c75bc8d63c793b4c74215d46d40d243118b5fb8f4e1aff248af0451ee15196fe
  - GSE188666_scRNA_LCMV_features.tsv.gz: dee1b30b23a9ceae7117d71b894ec7e4d5441fc934fbcb242493ef57b2e06e2d
  - GSE188666_scRNA_LCMV_metadata.tsv.gz: 0594a0463b0d667307736fa45453d9aae9ed761ce105d53b3a24c30dffa46393
  - GSE188666_scRNA_D5_LCMV.mtx.gz: 14bc8bf47f68d11b27dbabfc8ec6613385f58a83578a0a81b7187817602657b3
  - GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz: 3ec9ca5d7650e5ab428121161c7533d1da42d72ca1187dd9941e8de0207151c8
  - GSE188666_scRNA_D5_LCMV_features.tsv.gz: e3367683bbf603ef45124593bcc9a2868ffe393e422df1a5a035c24e7e76a576
  - GSE188666_scRNA_D5_LCMV_metadata.txt.gz: 9f8943dea624e8f56c894421eaad3b723b9e2698a3f8d148122fa516bab7c4ec
  - GSE188666_RAW.tar: db4932008b1b9e603c2f076aa7cf60302f38d8b66108188f3827a121a676c9a9
- **Tar members: 23, not the plan's 24** (controller ruling P6: accept; GSE188666_RAW.tar is the approved file and what GEO serves; no other source, including GSE188670, is approved). All end in `_tcr_filtered_contig_annotations.csv.gz`: GSM5687816-GSM5687834 and GSM6381666-GSM6381669. Which sample lacks a TCR file is left to Task 3 / 5a. The tar was never extracted. Members:
  - GSM5687816_Liver_Chronic_D21_CD8_TetNeg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687817_Liver_Chronic_D21_CD8_TetPos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687818_Lung_Chronic_D21_CD8_TetNeg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687819_Lung_Chronic_D21_CD8_TetPos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687820_Spleen_Acute_D21_CD8_TetNeg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687821_Spleen_Acute_D21_CD8_TetPos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687822_Spleen_Acute_D8_CD8_TetNeg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687823_Spleen_Acute_D8_CD8_TetPos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687824_Spleen_Chronic_D21_CD8_CX3CR1Neg_Slamf6Pos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687825_Spleen_Chronic_D21_CD8_CX3CR1Neg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687826_Spleen_Chronic_D21_CD8_CX3CR1Pos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687827_Spleen_Chronic_D21_CD8_PD1Pos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687828_Spleen_Chronic_D21_CD8_TetMid_tcr_filtered_contig_annotations.csv.gz
  - GSM5687829_Spleen_Chronic_D21_CD8_TetNeg_rep1_tcr_filtered_contig_annotations.csv.gz
  - GSM5687830_Spleen_Chronic_D21_CD8_TetNeg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687831_Spleen_Chronic_D21_CD8_TetPos_rep1_tcr_filtered_contig_annotations.csv.gz
  - GSM5687832_Spleen_Chronic_D21_CD8_TetPos_tcr_filtered_contig_annotations.csv.gz
  - GSM5687833_Spleen_Chronic_D8_CD8_TetNeg_tcr_filtered_contig_annotations.csv.gz
  - GSM5687834_Spleen_Chronic_D8_CD8_TetPos_tcr_filtered_contig_annotations.csv.gz
  - GSM6381666_Spleen_Chronic_D5_CD8_TetNeg_rep1_tcr_filtered_contig_annotations.csv.gz
  - GSM6381667_Spleen_Chronic_D5_CD8_TetNeg_rep2_tcr_filtered_contig_annotations.csv.gz
  - GSM6381668_Spleen_Chronic_D5_CD8_TetPos_rep1_tcr_filtered_contig_annotations.csv.gz
  - GSM6381669_Spleen_Chronic_D5_CD8_TetPos_rep2_tcr_filtered_contig_annotations.csv.gz
- `gzip -t` passed on all eight .gz files. The plan's tests pin no count of 24, so no test changed.
- **Licence:** the GEO disclaimer page and PMC11225711 both returned a reCAPTCHA page (not bypassed); LICENCE records both as UNVERIFIED.

## Phase 3 Task 3 decisions

- **Branch** `pivot/43-daniel-metadata`. The research note is [research/2026-11-19-daniel2022-metadata.md](research/2026-11-19-daniel2022-metadata.md); the script is `benchmarks/inspection/daniel2022_metadata.py`. Nothing was downloaded; no matrix and no contig file was opened.
- **No sample lacks a TCR file.** GSE188666 is 23 captures, each with one GEX and one TCR GSM (46 = 2 x 23). The "24" of the spec and the dataset note was a miscount (GSM5687816-GSM5687834 is 19 GSMs, not 20). Task 5a's `len(CAPTURES)` is therefore 23, "the number the note records". Controller ruling P6 stands; nothing is missing.
- **Two record inconsistencies are left to Task 5a's join check,** not resolved from names: (1) GSM6381667 and GSM6381668 have titles (D5 TetPos_rep1, TetNeg_rep2) that disagree with their attached file names (TetNeg_rep2, TetPos_rep1); (2) the D21 spleen TetNeg/TetPos and their `_rep1` pairs (GSM5687810-13 with GSM5687829-32) cross sequencers when paired by title. The note pairs by title and flags both.
- **Label column found: `ident`** (11 labels in the main file, 7 in the D5 file; the vocabularies differ). Sample column: `orig.ident`. Task 5a sets `LABEL_COLUMN`; Task 6 applies its rules.
- **The D5 metadata repeats the chronic spleen D8/D21 cells** (43,820 of them, every one also in the main file) with its own labels and embedding. Found with a scratch script outside the repo that read the two files' index, `orig.ident` and `ident` only (`usecols`), beyond the fixed output of the committed script (now `--compare`, see Fix round 1). It bears on Task 6's Days rule ("processed together", "same vocabulary"), which the records answer "not stated" and the vocabularies answer "no".
- **Records read beyond the GSMs:** the 46 linked BioSample records (efetch, text; their only attributes are source name and cell type) and the GSE188666/GSE188670 series records. They are sample records, not data files. No mouse identifier is in any of them.
- **The paper was not read:** PMC is behind a CAPTCHA (Task 2) and Europe PMC has no open full text (`fullTextXML` HTTP 500). The note marks mice per capture and pooling as not stated.
- **The script refuses a withheld column as `--labels` or `--samples`** (exit 2, no values printed), tested; the plan's two tests are as written, plus one for the tag vocabulary.
- **Freeze kept:** `tcr_cdr3s_aa` and `tcr_cdr3s_nt` were listed by name and dtype only. No clone x state, clone x capture or clone x lineage count was computed.
- **Fix round 1 (PR #48 review):**
  - Every withholding rule is pinned by its own parametrised test (`Clonotype`, `ir_foo`, `airr:x`, `raw_clonotype_id`, `my_cdr3`, `Xclon`, `foo_gene`, upper-case `TRB_V` and `CHAIN_x`): a planted value never reaches stdout or stderr and the column is listed by name only. Four mutants (drop the regex, drop `clon`, drop `_gene`, drop `.lower()`) each fail it.
  - A `--labels` or `--samples` column absent from a file now exits non-zero with "is not a column of <file>" (it was a soft "cross-tab skipped", exit 0); both flags are tested.
  - The output shape (value-count order, cross-tab header with samples as rows and labels as columns, the "first index values" line) is pinned.
  - The 43,820-cell overlap and the main-vs-D5 label cross-tab are now reproducible with `--compare --labels ident --samples orig.ident` on the committed script. It reads only the index and those two columns of each file, applies the same withholding guard, and prints cells in both, only in main, only in D5 and the cross-tab of the shared cells. Tested on synthetic files; not run on the real files. The earlier "scratch script" is superseded by it.

## Phase 3 Task 4 decisions

- **Branch** `pivot/44-contig-clones`, from `main` 7d7a212. Tested on synthetic contig tables only; no Daniel file was read.
- **Barcodes without a qualifying contig are kept.** Output has one row per distinct `(capture, barcode)` in the input, with counts 0 and NaN ids, so Task 5 can tell "seen, unlabelled" from "absent". Rows are sorted by `(capture, barcode)`.
- **`raw_clonotype_id` per barcode** is the first non-missing, non-`"None"` value among its contigs in input order, prefixed `capture:`. It ignores contig qualification (it is the 10x label, reported only).
- **`is_cell`** is parsed with `parse_flag` when the column is present; a missing column is not an error.
- **Missing required columns** raise `BuildError` naming the column names (ours, not data values).
- **`read_members` errors** state only a count, never member names (untrusted). A wanted name that is a directory or other non-regular member counts as missing. Duplicate wanted members: the last one wins (tar semantics).
- **`cdr3_nt` regex** is applied with `fullmatch` and is case-sensitive: lowercase bases do not qualify.

Fix round 1:
- **Zero contig rows raise `BuildError("no contigs")`**, not an empty frame. The plan fixes no empty-input contract, and an empty contig table in a real build means a wrong or truncated file, so the builder should stop rather than label every cell as unseen. Task 5b maps a cell with no contig to NaN itself.
- **A blank `raw_clonotype_id`** (empty or whitespace-only) is missing, like `"None"` and NaN.
- **`read_members`:** a file that is not a tar raises `BuildError` naming the path (ours) and the error class. Wanted names never appear in the "other" list; for a repeated wanted name the last member wins, so a later directory hides an earlier file and the name counts as missing.
- **Pinned by tests:** the TRA/TRB count rules (3 TRA + 1 TRB, repeated TRA), first non-missing raw id, `parse_flag` near misses, lowercase `cdr3_nt` ignored and `ANA` accepted, sorted index and `other`, no member names in errors, and `is_cell` NaN, `"None"` and numpy bools. Each fix test was checked against a re-applied mutant.

## Phase 3 Task 5a decisions

- **Branch** `pivot/46-daniel-formats`, from `main` df236f5. Script `benchmarks/inspection/daniel2022_formats.py`; the Formats and Join check are sections 9 and 10 of [research/2026-11-19-daniel2022-metadata.md](research/2026-11-19-daniel2022-metadata.md). Nothing was downloaded; the tar was read in memory, never extracted.
- **Formats:** both matrices are MatrixMarket coordinate integer, genes x cells (31,053 x 96,750 and 31,053 x 77,614); features are one gene-symbol column and the D5 list equals the main one; barcodes are `<tag><16-nt core>-1`; the 23 contig headers are the 10x layout and already match `CONTIG_COLUMNS` (`COLUMN_RENAMES = {}`). `chain` holds Multi, TRA, TRB, TRD, TRG; every flag is `True`/`False`; `high_confidence` and `is_cell` are always True. `LABEL_COLUMN = "ident"`.
- **The D5 matrix repeats the chronic spleen D8/D21 cells:** 43,820 of its barcodes are also main barcodes, matching the metadata (Task 3).
- **Capture map:** a contig barcode `b` is the cell `tag + b`; the tag is the GEX title plus `_`. Every member joins its own tag at 0.58-0.96, every other tag at 0.02 or less (table in the note).
- **Conflict 1 (GSM6381667/68) resolved by the join, toward the file names:** by title 0.01 and 0.01; by file name 0.75 and 0.84. The `Capture` tags and label fields (`sort_gate` gp33- for GSM6381667, gp33+ for GSM6381668) follow the file names, because the labels must describe the cells the barcodes join; the TCR GSMs' own titles and characteristics describe each other's file. `gsm` stays the GSM the file is attached to.
- **Conflict 2 (GSM5687829-32) resolved by the join, toward the title pairing:** own shares 0.95, 0.85, 0.96, 0.87; the sequencer pairing's tags get at most 0.02. The sequencers simply differ.
- **`matrix`:** `d5` for the four D5 captures, `main` for the other 19 (the main matrix holds all of them; D5 lacks 24 of the 11 shared captures' cells). Whether D5's copies are used is Task 6's.
- **Withholding is an allow-list:** of a contig member only `barcode`, `chain`, `productive`, `high_confidence`, `full_length` and `is_cell` are read (`usecols`). Each of the 12 other real columns and an unknown column is pinned by a planted-value test. Map errors give counts only, never member names.
- **Beyond the plan:** the script also prints the orientation it reads from the size line, the barcode suffixes and the D5-in-main barcode count; `test_the_capture_map_joins_on_the_real_files` (`real_data`, `slow`) re-checks the 0.5 and best-other rule from `CAPTURES` on the cached files (passes locally).
- **Freeze kept:** no value of `cdr3`, `cdr3_nt`, `v/d/j/c_gene`, `raw_clonotype_id`, `raw_consensus_id`, `reads`, `umis`, `contig_id` or `length` was read. The shares count barcodes only. No clone x state, clone x capture or clone x lineage count was computed.
- **Fix round 1:**
  - The format script validates every contig member's header (it must contain `barcode`) before printing anything from a member; a headerless member is an error naming a count and position, never the member name or its first line.
  - A member that is not valid gzip or UTF-8 gives the script's error, with no bytes echoed.
  - D5-equals-main is tested for reordered features and for one changed name of equal length (both False).
  - All 23 `CAPTURES` are pinned literally from the note, plus a test that tissue, day, arm and sort gate agree with the member name; `MATRIX_ORIENTATION` and `COLUMN_RENAMES` are pinned.
  - Tests for the union join (a D5-only tag), a member with 0 barcodes (`n/a`), one unmapped member, `neither` orientation and the exactly-100,000-value integer check.
  - The note's section 10 states that the join share counts every contig row's barcode.

## Phase 3 Task 5b decisions

- **Branch** `pivot/49-daniel-builder`, from `main` ddb39c3. `build` and `BUILD_VERSION = 1` in `bench/builders/daniel2022.py`; fixture `tests/unit/daniel_fixture.py`. Freeze kept: the builder never ran on the real files, no real file was opened (the formats came from the note's sections 2, 9 and 10), and no clone x state, clone x capture or clone x lineage count was computed.
- **obs names are `<gsm>_<core>`, the core without `-1`.** A cell of capture `c` is a matrix barcode that fullmatches `c.tag + [ACGT]{16} + "-1"`; the core is the 16 nt. Its contigs are those whose barcode is exactly `<core>-1` (so `c.tag + b` is the matrix barcode); a contig barcode of any other shape joins nothing. Every real barcode ends in `-1` (note section 9), so dropping it loses nothing.
- **Metadata join:** the metadata index is the matrix barcode itself (`<tag><core>-1`, note sections 2 and 9), so cells are matched by exact string. Kept = matrix cells of a listed capture that its own matrix's metadata lists, in matrix order, captures in the listed order. A metadata cell of a listed capture that the matrix lacks is a `BuildError` (the note says each barcode file is its metadata file's cells). Only the index and `LABEL_COLUMN` are read (`usecols`), so `tcr_cdr3s_*` never enter the build; a blank label is NaN. With `LABEL_COLUMN = None` there is no `state` column.
- **D5 handling:** each capture is read from the matrix and metadata its `Capture.matrix` names, and only matrices a listed capture uses are read. The D5 matrix's copies of the 11 chronic spleen D8/D21 captures are never read (those captures are `main`), so nothing is deduplicated: main and D5 cells are different captures by tag. When both matrices are used, their feature lists (and symbols) must be identical, else `BuildError`.
- **obs columns:** `capture`, `arm`, `tissue`, `sort_gate`, `state`, `donor` (categorical), `day` (int64), `clone_id`, `clone_id_trb`, `raw_clonotype_id` (NaN without a call), `receptor_n_tra`, `receptor_n_trb` (int64, 0 without a contig). `var_names` are the feature lines (gene symbols); `var["gene_symbols"]` is added only for a features file with a second column. X is CSR float32 cells x genes, from `scipy.io.mmread` on a gzip stream, transposed for `genes_x_cells`, checked for shape against barcodes and features and for non-negative integer values.
- **`uns["lumencell_build"]`:** exactly the plan's keys (`version`, `params`, `n_cells`, `contig_barcodes`, `join_misses`, `ignored_members`, `labelled_cells`, `labelled_cells_trb`); no extras. `contig_barcodes` counts distinct barcodes over every contig row of the member (Task 5a's denominator). A join miss is a contig barcode with at least one qualifying contig (`receptor_n_tra + receptor_n_trb > 0`) and no kept cell of its own capture; a barcode with only non-qualifying contigs is not a miss.
- **`ignored_members`:** tar members that belong to no `CAPTURES` entry (members of known but unlisted captures are not listed, and not read). A name matching `download.SAFE_NAME` is stored as is; any other (the fixture's `../stray.csv.gz`) is stored as `"<unsafe member name>"`, so the count is kept and the untrusted string never reaches the file. Same key and type as Sheih's builder (#52), which stores names unfiltered.
- **Params are refused, never approximated:** keys other than `captures` and `namespace`, a string or empty list, a repeated or unknown GSM, a namespace that misses a listed capture or names an unlisted one, or a blank donor. Errors name counts and our own constants, never a barcode, label, member name or GSM.
- **Atomic write:** the builder writes to a `mkstemp` file in `out`'s directory and `os.replace`s it onto `out`; on any error the temp file is removed and `out` is not created. Tests check that nothing else is left beside `out` and that a malformed matrix leaves the directory empty.
- **Beyond the plan's list of tests:** the kept-cell order and the `obs`/`uns` values, a D5 capture read from the D5 matrix while the main capture's repeated cell is read from main, mismatched features refused, a two-column features file, a failed build, and `LABEL_COLUMN = None`.
- **For Task 6/7:** the real main matrix has 190,891,292 stored values; `mmread` holds them all before the subset. Seat A's estimate for the code at 6d07e75: about 4.2 GB after `mmread` and about 8 GB peak (the COO arrays plus the CSR plus the float64 check temporaries). After P13 below the expected peak is about 5.5 GB (mmread's 16 bytes a value plus the CSR's 12, 28.7 measured on the synthetic matrix, against 37.7 before); an estimate, unmeasured on real data. Column-subsetting before the CSR (reading only the listed captures' columns) is deferred: Task 8 should do it if the real build fails for memory.
- **Fix round 1** (controller rulings after review; every new test was seen to fail on the old code or under a named mutant):
  - **P12, a listed capture that joins nothing is refused.** A listed capture that keeps 0 cells raises `BuildError` (its GSM from `CAPTURES` and nothing else); so does one with qualifying contig barcodes that all miss (a structural zero: a wrong tag or a changed barcode format). This is not a threshold and does not replace Task 8's 0.5 share check. The old "no cell kept at all" check is gone: the per-capture check covers it. Tests: `test_a_listed_capture_with_no_kept_cell_is_refused`, `test_a_listed_capture_whose_metadata_keeps_nothing_is_refused`, `test_a_capture_whose_qualifying_barcodes_all_miss_is_refused` (old code built an empty or all-NaN capture).
  - **P13, memory of `_read_matrix`.** `del m` right after the CSR is built; `_are_counts` checks an integer dtype with `min() >= 0` only (no float64 temporaries) and a float dtype in 1M-value slices. Test: `test_reading_an_integer_matrix_makes_no_float_temporaries` (tracemalloc peak under 33 bytes a value; old code 37.7, now 28.7).
  - **Join misses.** `test_a_miss_whose_only_qualifying_contig_is_a_tra_is_counted` (mutant: count only TRB), `test_misses_are_counted_under_their_own_capture_and_never_leak_across` (a GSM_B or GSM_D miss counts under its capture; a GSM_B contig on a core only a GSM_A cell has is a GSM_B miss; mutant: match on the core alone), `test_is_cell_false_contigs_give_no_clone_and_no_join_miss` (mutant: drop the `is_cell` filter). The fixture gained a TRA-only miss, a GSM_B and a GSM_D miss and the is_cell=False cells; GSM_B is now Armstrong/lung and GSM_D liver.
  - **Output.** `test_a_failing_write_leaves_no_output_and_no_temp_file` (a partial write that raises: `out` absent, an earlier `out` byte-identical, no temp file; mutant: write to `out` directly), `test_capture_is_the_gsm_prefix_of_every_obs_name`, `test_arm_day_tissue_gate_and_donor_follow_each_capture` (mutants: first capture's value), `test_the_obs_day_is_int64`, `test_a_blank_label_becomes_nan` (mutant: keep the empty string).
  - **Refusals, one test each:** a metadata cell absent from its matrix; non-integer and negative counts; matrix shape against barcodes and features; duplicate matrix barcodes; duplicate feature names; non-unique metadata ids; metadata with no id column; a features file mixing 1 and 2 fields; D5 features in another order and with the same ids but other symbols (mutants: compare sorted, compare ids only); a namespace that is a list, a non-str donor, a whitespace donor; `test_errors_never_carry_donor_values_or_paths` (no donor value and no temp path in a bad-params or corrupt-matrix message).

## Phase 3 Task 6 decisions

- **Branch** `pivot/50-daniel-declaration`, from `main` 30b98c9. The declaration and each rule's evidence:
  [research/2026-11-20-daniel2022-declaration.md](research/2026-11-20-daniel2022-declaration.md). Freeze kept: written
  from Task 3's note (with Task 5a's sections), the GSM records and the literature; the derived file was not built,
  `bench fetch daniel2022` was not run, no contig or metadata file was opened, and no clone count of any kind was
  computed. The `>50 cells` rule used the per-capture `ident` counts already in the note.
- **Literature:** all six DOIs (Hudson, Zander, Beltra, Im, Utzschneider, Daniel) match Crossref; no correction.
  Full texts read as PMC author manuscripts through NCBI E-utilities `efetch db=pmc` (no CAPTCHA); Europe PMC
  `fullTextXML` returned 500/503 or timed out on retries. Utzschneider 2016 is not in PMC: abstract only (PubMed),
  used only as support for the root. Daniel 2022 (Task 3 found it unreadable) is readable this way; only its
  non-clonal marker description, sorts and Methods were used.
- **Subset:** chronic only; days 8 and 21 (D5 excluded: separate Seurat object, a different 7-label vocabulary, and no
  record of joint processing; Daniel's Methods cluster "Chronic including Day 5" apart); spleen only (liver and lung
  only at D21). State gates: the CX3CR1/SLAMF6 sorts GSM5687824-26; tetramer and PD-1 sorts are not.
- **Namespace rule (b):** the records name no mouse, so (a) cannot apply; the 3 state-gated captures are dropped and 8
  remain (GSM5687827-34), namespace = capture (donor = GSM). Assumption in `reason`: one mouse per capture.
- **Root `TEx_prog`** (Daniel's Texprog: Tcf7, Slamf6, Id3; sorted PD-1+SLAMF6+). PD-1 level and TIM-3 are not stated
  by Daniel for the cluster.
- **Judgement call, lineages:** by transfer of sorted subsets, Hudson and Beltra give Tpex -> CX3CR1+ transitory/Texint
  -> terminal (TEx), and Zander gives the CX3CR1+ subset as an end state (no conversion shown, naive hosts). The
  papers' CX3CR1+ subset is two Daniel labels, `TEx_int` and `TEx_KLR` (Daniel's non-clonal markers; Daniel itself
  splits the CX3CR1+ pool citing Hudson and Zander; the CX3CR1+ sort holds 2,805 TEx_KLR and 1,521 TEx_int). A stage
  that maps to two labels cannot be written as one label per stage without inventing an order, so
  `lineages = {}` and `gated = []`; the lineages are in the note and in `notes`.
- **Clonal: declared not assessable** ("fewer than two gated lineages, so no clonal comparison between lineages",
  plus the namespace assumption), so the two clonal gate lines are left out of `[real.daniel2022.gates]`. The
  outcome does not hinge on the judgement call: mapping Zander's subset to `TEx_KLR` alone gives one gated lineage.
- **New types:** `RealClonal` (frozen; `columns()` returns the four scoring-only columns) and
  `RealDataset.clonal` / `report_time_key` as last fields; `__post_init__` refuses an assessable declaration with
  fewer than two gated lineages. Tests pin the captures, namespace, root, lineages, gated, the whole `RealClonal`
  (reason text included), the subset rules against `CAPTURES`, and that the params pass `daniel2022._check_params`.

Fix round 1:

- **BINDING on Task 7 (the unmeasured-gate gap).** With `gated = []` the three base gates score NaN and
  `check_gates` skips NaN, so a Daniel row would print PASS having measured nothing. `nan = "fail"` was refused: the
  gates stay exactly the plan's TOML. Task 7 renders such a row as "lineage gates: not measured (no gated lineage,
  declared 2026-11-20)", never PASS, in the CLI `bench real` line and the README cell, and must add
  `test_a_dataset_with_no_gated_lineage_is_shown_as_not_measured` (both show that text, neither shows PASS/pass for
  that row). This PR's `test_daniel_base_gates_score_nan_and_are_skipped_without_a_gated_lineage` documents the gap.
  The TOML comment and the changelog entry state it.
- **Distinct gated lineages.** `RealDataset.__post_init__` counts `len(set(gated) & set(lineages))`, so `["A", "A"]`
  and a name absent from `lineages` do not make a declaration assessable.
- **For Task 7/8:** after this merge `make bench-real` (`bench real` with no `--datasets`) fails "not downloaded:
  daniel2022" until Task 8 builds the derived file. `bench fetch daniel2022` builds without asking because its
  sources are cached, as Task 8 Step 3 expects.
- The changelog entry is pinned by `test_the_daniel_changelog_entry_is_dated_and_matches_the_registry`; the research
  note's namespace section no longer quotes Daniel's clonal section.

## Phase 3 Task 7 decisions

- **Branch** `pivot/51-bench-real-clonal`, from `main` fb06f05. Freeze kept: tests use the synthetic `fake` fixture
  (1500-cell bifurcation with clones, a built Daniel file's columns); the Daniel derived file was not built and
  nothing ran on real data. Clone, receptor and namespace columns reach `from_prep` and `stability_row` as
  `clonal_columns` (`ds.clonal.columns()`), so they never reach inference; pinned by two tests.
- **P14 (controller ruling, rendered).** `datasets.lineage_gates_note(ds)` returns `None` when `ds.gated` is non-empty
  (pancreas, bone marrow and gastrulation render exactly as before) and otherwise
  `lineage gates: not measured (no gated lineage, declared <ds.clonal.declared>)`. A dataset with no clonal
  declaration has no date to cite, so it reads `lineage gates: not measured (no gated lineage)`. It replaces the PASS
  (CLI) or `pass` (README cell) mark. A real gate failure still prints FAIL or `fail: ...`, never the note. Gates and
  `check_gates` are untouched (no `nan = "fail"`). The README cell and the CLI line also carry the clonal entry, so a
  Daniel cell shows both, e.g. `Spearman –, branch –, terminal – (lineage gates: not measured (no gated lineage,
  declared 2026-11-20)), clonal: not assessable (declared 2026-11-20, <reason>)`. Test:
  `test_a_dataset_with_no_gated_lineage_is_shown_as_not_measured`.
- **Reported fields only when they apply.** `reported_day_spearman` is emitted only when the dataset sets
  `report_time_key` (otherwise a NaN column would appear for pancreas after a re-run); the clone fields only when
  `ds.clonal` is set. Day values are read as numbers (first number in the text, so `D21` and `21` both work).
  `reported_labelled_share_of_raw` is cells with a `clone_id` among cells with a raw clonotype, per the Task 6 TOML
  text. `reported_n_clone_strings` is an integer. The ARI is the Hubert-Arabie index in numpy; with a zero
  denominator it is 1.0 for identical partitions and 0.0 otherwise.
- **`clonal_reason` is always a string** (`""` when none): `test_gates.py` expected `None` and now expects `""`. The old
  `test_a_real_case_row_is_not_applicable` stub dataset gained `clonal=None, report_time_key=None`.
- **README** regenerates unchanged (no result file has `reported_*` or a clonal status other than not_applicable).
- **Fix round 1** (controller rulings, each test seen failing on the old code or under its mutant first):
  - **Day parsing.** A numeric `report_time_key` column is used directly (`pd.to_numeric`, NaN stays NaN); the
    first-number regex is only for string or object columns (`D5`, `21`). The old code read `1e-05` as 1. Tests: the
    numeric, string and all-unparsable day tests.
  - **P16, labelled share.** `reported_labelled_share_of_raw` is the declared text literally: cells with a
    `clone_id` divided by cells with a `raw_clonotype_id`, NaN with no raw cells. The declaration wins over a
    bounded version, so it can exceed 1 when cells carry a `clone_id` but no raw id (5/4 in the test). Report it
    as defined. Test: `test_the_labelled_share_is_clone_cells_over_raw_cells_as_declared_and_can_exceed_one`.
  - **Recurrence.** Cells whose capture is missing are dropped, not counted as a capture called "nan".
  - **ARI.** The zero-denominator branch is reachable only for two single clusters or two all-singletons, which
    are always identical, so the dead 0.0 outcome is removed (it returns 1.0). Tests cover both-single, all
    singletons, one cluster against singletons (0.0), a partial overlap (-0.5), n = 2, and the exhaustive small-n
    check of the claim.
  - **Pinned values.** Tests for `reported_metrics` on hand-built obs (capture-keyed raw ARI, restriction to cells
    with both ids, TRB-only truth uses `trb_only_key` and the donor column, empty recurrence is `(nan, 0)`),
    `CaseResult.row` rounding, `truth_for` columns, the report cell variants and the exact "Reported, not gated"
    header. The P14 test is split so a real failure (`fail:` in the README, FAIL on the CLI line) is asserted
    on its own.

## Open questions for Shreyan

- **H-1 (closed 2026-11-19): should a failed check keep `weak` under the cap?** Shreyan ruled yes: keep
  `weak`; while the verification has failed only `supported` is capped at `provisional`. See "Phase 2 H-1
  decision".

- **depth_only/N2 (closed 2026-11-17).** The failed check stays recorded as failed; Shreyan
  approved a pre-declared N2 amendment for phase 2 (see "Phase 2 rulings").
- **Reordering cells (closed 2026-11-17).** Within Monte Carlo error accepted (see "Phase 2
  rulings").
- **Verdict verification failed (closed 2026-11-17).** 28 of 31 control runs were verdicted supported
  (max 0); the served verdict is capped at provisional. Shreyan ruled option (a): accept a
  provisional-only verdict, keep the cap, state the failure plainly, and record a pre-declared C1/C2
  redesign as follow-up work (see "Phase 2 Task 9 decisions and result"). Tasks 11 and 12 proceed.
- Commit 33c281d and PR #5 were created at 07:57 under your account, not by the orchestrating
  session. Was that you, or a tool?

## Security review (phase 1)

Reviewed by hand along the brief's untrusted-input paths (the `security-review` skill reviews a
pending diff and phase 1 is already committed, so the diff `142a788..HEAD` was read directly).
At review time phase 1 was not yet merged into `main` (PRs #14 to #23 were a stack), so the diff
was `142a788..HEAD` on `pivot/19-phase1-security`. The stack has since merged, and the Done When
commands were re-run on `main` (PR #24).

| # | Severity | Where | Finding | Outcome |
|---|---|---|---|---|
| S1 | Low | `api/models.py` `ClonalSettings` | `clone_key`, `donor_key`, each strata key and the strata list had no length bound, so one PUT could push an unbounded list or name through `check_settings`, its error messages and the state file. | Fixed: 256 characters per name, 8 strata keys. Test: `test_oversized_settings_are_refused_before_the_file_is_looked_at`. `make contract` changed only `contracts/openapi.json`. |
| S2 | Low | `engine/store.py` `trust_path` / `write_trust` / `read_trust` | The trust file name was checked (`TRUST_NAME`, `fullmatch`) but the run id was not; callers validate it today (`_RUN_ID` in `api/runs.py`, computed ids in the CLI), so nothing was reachable. | Fixed as defence in depth: `RUN_ID` fullmatch in `trust_path`. Test: `test_trust_paths_refuse_a_run_id_that_is_not_a_plain_name`. |
| S3 | Info | `bench/calibration.py` / `tests/clonal/test_calibration.py:32` | `LUMENCELL_CALIBRATION_OUT` is written as given, following symlinks, and the parent directory must exist. The variable is set by whoever runs the slow test, not by file or request data. | Recorded, no change. |
| S4 | Info | `engine/store.py` `read_trust`, `prep/housekeeping.py` | Trust files are read without a size cap and through symlinks. Anyone who can plant files in the cache directory already owns the user's account. Housekeeping anchors `_TRUST_TEMP` with `fullmatch`, skips a symlinked `trust/` folder, and `check_prune` refuses symlinked removals. | Recorded, no change. |

Checked and found sound, no change:
- Path traversal: trust file names are `clonal-<16 hex>` from a hash; the column names never reach a path or a job key (`score_clonal:<rid>:<hash>`); `TRUST_NAME` and the housekeeping patterns use `fullmatch`, so a trailing newline does not slip through (tested for names).
- Terminal and log injection: every message that names a column uses `!r` (`check_settings`, `scoring_only_message`, the notes), so an escape or newline in a column name prints escaped. Pinned by `test_problems_escape_control_characters_in_column_names`. The CLI prints one raw column name, the `symbol_source` var column in `cli.py` (a name from the user's own file, shown to that user in their terminal).
- HTML and JSON injection: the web app has no view of the notes or the clonal settings yet (only the generated types), so there is no sink. Re-check when the phase 2 UI renders `notes` (React text nodes escape; keep `dangerouslySetInnerHTML` out). Trust files are written with `allow_nan=False` and read back as a JSON object only.
- Many strata keys and many distinct values: strata keys must exist in `obs`, so the list is bounded by the file's columns, and now by 8.

### Done When (phase 1)

Run on the stack head `pivot/19-phase1-security` first, then re-run on `main` (aac08ed, after #14-#23 merged, 2026-11-14). The `main` run gave the same results: `make ci` exit 0; `make bench-synth` exit 0 with 0 regressions; 79 isolation/labels/scores/properties tests passed; 16 trust API tests passed; 3 constants tests passed; `clonal agreement: ok (own cells, m = 2)`; `clonal-calibration-2026-11-13.json` still `passed=False` (the recorded `depth_only/N2` failure, unchanged).

| Command | Result |
|---|---|
| `make ci; echo $?` | exit 0 (run after the fixes): PASS |
| `make bench-synth`, exit code and `grep -c '^REGRESSION'` | exit 0, 0 regressions: PASS |
| `ls benchmarks/results/real-*.json \| tail -1`; `grep -n 'TR-gene exclusion' benchmarks/thresholds.toml` | `real-2026-11-13.json` (the Task 9 file); changelog entry at `thresholds.toml:169`: PASS |
| isolation, labels, scores, properties tests | 79 passed, exit 0: PASS |
| calibration file `passed` | `clonal-calibration-2026-11-13.json False`: **FAILED**. `depth_only/N2` failed calibration in Task 8 (57/200 vs bound 18) and was recorded as a real failed check; see "Task 8 decisions and result". No data, threshold, seed or definition was changed. |
| `tests/api/test_clonal_trust.py` | 16 passed: PASS |
| `lumencell synth` then `lumencell run` with the clone options (temp cache and state dirs) | output contains `clonal agreement: ok (own cells, m = 2)`: PASS |
| `tests/unit/test_trust_constants.py`; first commit touching thresholds, baseline or results | 3 passed; `6ff4676` "Declare clonal agreement definitions ... (#7)", the Task 1 commit: PASS |
| `git diff 142a788 -- docs/HANDOFF.md` | hunks cover lines 3, 32, 49, 71-73, 89-90, 159-160 and 185: PASS |
| this section | present: PASS |

## Security review (phase 2)

Reviewed by hand (the security-review skill was not invoked) along the brief's untrusted-input paths, on the diff `2217f76..HEAD` (phase 1's last
merge to `main`, c2b3e97 as the head). Tasks 2-4 (Tirosh lists, cell-cycle stage, cell-cycle check) merged afterwards (#39, #42, #46); their items are closed below.

| # | Severity | Where | Finding | Outcome |
|---|---|---|---|---|
| S5 | Low | `trust/stability.py` `to_api_stability`, `not_measured` | A replicate's `error` (`ExcType: message`) and a failed job's error reach the served `notes` and the verdict reasons verbatim. An exception message can carry a name from the user's file, with control characters, and has no length bound. | Fixed: `safe_text` escapes non-printable characters and caps a note at 300 characters. Tests: `test_a_replicate_error_is_one_escaped_line_of_bounded_length`, `test_a_not_measured_note_is_escaped_and_bounded`. |
| S6 | Low | `engine/store.py` `read_trust` | Trust files are written strict (`allow_nan=False`) but were read with Python's lenient `json.loads`, which accepts `NaN` and `Infinity`. A hand-edited or damaged file would validate and then be served. | Fixed: `parse_constant` refuses them, so the file reads as unreadable and is healed like any other. Test: `test_a_trust_file_with_a_non_json_number_is_unreadable`. |
| S7 | Info | `api/runs.py` `_trust_lookup` | A GET resubmits a job that is `done` while its file is missing or invalid, with no count. If a bug made a job write an invalid file every time, each GET after the job finished would queue another hours-long stability run. Failed and cancelled jobs are not resubmitted, so only a deterministic "done but invalid" output reaches this. | Recorded, no change: it needs a code bug, not input, and the missing-file case (a user clearing the cache) is the reason for the rule. A retry counter would be new state in the job manager. |
| S9 | Low | `api/runs.py` `Runs.clonal_agreement` | Found by the review seat. A failed clonal job's `job.error` was served raw as the agreement's note, and it reaches the verdict reasons through C5. Same class as S5. | Fixed: wrapped in `safe_text`. Test: `test_a_failed_score_note_is_escaped_and_bounded`; `test_safe_text_pins_the_cap_the_marker_and_every_control_character` pins the 300 cap, the marker and DEL, C1 and U+2028. |
| S10 | Info | `api/runs.py`, `api/session.py` (job and run status) | The raw `job.error` is also served in job and run status (phase 1 surface). JSON escapes control characters on the wire and only the local user's own file can put text there. | Recorded, no change: it is the job's own diagnostic, not a note composed for the verdict; capping it would hide the traceback summary the user needs. |
| S8 | Info | `tests/unit/test_verdict_calibration.py` | `LUMENCELL_VERDICT_CALIBRATION_OUT` and `LUMENCELL_VERDICT_VERIFICATION_OUT` (and the `..._N2_OUT` one) are written as given, as `LUMENCELL_CALIBRATION_OUT` is (S3). They are set by whoever runs the slow test. | Recorded, no change, as S3. |

Checked and found sound, no change:
- Obs names in replicate roots: the root is `cell:<obs_name>` and is found by exact match on the string index (`_root_index`); a repeated name and an absent name are both refused with `StabilityError`, and the message uses `!r`, so control characters print escaped. Existing tests in `tests/unit/test_stability.py` pin both refusals. The name is never a path, a job key or a regex.
- Lineage names: no lineage name reaches the verdict reasons or the stability notes (the reasons are fixed sentences with numbers). The clonal and stability notes that carry an exception or a job error go through `safe_text` (S5, S9); job and run status carry the raw error (S10).
- Trust file names and job keys: the stability file is `stability-<16 hex>` from `stability_hash()`; `TRUST_NAME` and `RUN_ID` use `fullmatch` in `trust_path`; the job key is `stability:<run id>:<hash>`, where the run id passed `_RUN_ID` at the route. The worker refuses a payload `trust_name` that is not `stability_name()`.
- The stability cancel path: the handler writes nothing until the end and `write_trust` renames a temp file, removing it on any exception; `test_a_cancelled_stability_job_leaves_nothing` (a terminated worker) and `test_trust_refuses_nan_and_leaves_nothing` pin it. A terminated worker's stray temp is swept by housekeeping (phase 1, S4).
- The GET-healing path: a failed or cancelled job is never resubmitted by a GET (`test_a_cancelled_stability_keeps_the_verdict_provisional_and_is_not_resubmitted`, `test_a_failed_stability_job_is_reported_with_its_error_and_not_resubmitted`), and the job manager keys jobs, so concurrent GETs share one job. The one open edge is S7.
- Strict JSON: writes use `allow_nan=False`; undefined stability values become `null` (`_num`); reads are now strict (S6).

Tasks 2-4 items, closed when their PRs merged:
- Gene symbols in cell-cycle matching (case, whitespace, the `var` symbol column): resolved in #42 (Task 3). `match_genes` takes the first matching var symbol, the view's var names go through `make_index_unique`, and duplicate symbols are scored once. The stage is report-only: any exception while checking or scoring becomes `not_scored` with the error named and one stderr line, so it never fails prep (`test_a_scoring_error_is_reported_and_prep_still_completes`, `test_duplicate_symbols_are_scored_once`).
- The package-data TSV parsing of the Tirosh lists: resolved in #39 (Task 2). `_rows` checks the header row and the column count of every row and strips fields, so a missing header or a short row is refused naming the line (fix round 1). The xlsx source was read with stdlib `zipfile` and `ElementTree` only, and every script that read a downloaded source ran with `python -I`.
- The downloaded sources are never committed: resolved in #39. They were kept outside the repo, with SHA-256 recorded in "Phase 2 Task 2 decisions and result"; `git status` is clean of them.
- Lineage names in the cell-cycle check (new, #46, Task 4): the names reach only the served notice text, never a path, a job key or a log line. An unreadable run cell file is reported as `not_scored` with the exception type named and one stderr line (ruling R20), not hidden as `cell_cycle: null`; any non-OSError still propagates.

## Done When (phase 2)

Run on `main` 570fb5e (every phase 2 task, H-1 and the disconnected-kNN error merged) from the worktree root of `pivot/45-phase2-close`,
2026-11-19, exactly as the plan writes them. Nothing was changed to make any of them pass. Two fail, both because the
verdict verification failed (Task 9, ruling F14, owner rulings (a) and H-1); the cell-cycle command now passes (Tasks 2-4 merged).

| Command | Exit | Output | Result |
|---|---|---|---|
| `make ci; echo $?` | `0` | ruff passed; `1098 passed, 36 deselected, 3 warnings in 313.63s`; web `36 passed`; web-e2e `1 passed` | PASS |
| `make bench-synth > "$TMPDIR/synth.log" 2>&1; echo $?` (background) | `0` | the 12-case matrix ran; the log ends with 29 `gate` lines for the experimental methods' recorded gate misses (for example `scfates_ppt: binary2/realistic/seed0: pt_spearman_global = 0.869 (needs >= 0.9)`), which are the known baseline, not regressions | PASS |
| `grep -c '^REGRESSION' "$TMPDIR/synth.log"` | `1` (grep's own status for no match) | `0` | PASS |
| `uv run pytest tests/unit/test_prep.py tests/unit/test_prep_cache.py tests/unit/test_cell_cycle_genes.py tests/api/test_cache_healing.py -q; echo $?` | `0` | `87 passed, 1 warning in 86.93s` (the warning is the non-unique var names from `loader.py`) | PASS (this failed on c2b3e97: `test_cell_cycle_genes.py` did not exist; Tasks 2-4 added it) |
| `uv run pytest tests/api/test_robustness_trust.py -q; echo $?` | `0` | `17 passed in 26.75s` | PASS |
| the calibration/verification `python -c` command | `0` | `benchmarks/results/verdict-calibration-2026-11-17.json True benchmarks/results/verdict-verification-2026-11-17.json False` | **FAIL**: the second value is `False`, not `True` |
| `grep -c '^\[tree_support.calibrated\]' benchmarks/thresholds.toml` | `0` | `1` | PASS |
| `uv run pytest tests/unit/test_trust_constants.py tests/unit/test_verdict_calibration.py -q` | `1` | `1 failed, 35 passed in 708.58s (0:11:48)`; the failure is `test_verification_on_fresh_seeds` | **FAIL** |
| `git merge-base --is-ancestor $(git log --format=%h -1 --grep='Declare the T-cell robustness') $(git log --format=%h --reverse -- 'benchmarks/results/verdict-*.json' | head -1); echo $?` | `0` | c8a50d9 is an ancestor of 0d0608a | PASS |
| this file has "Security review (phase 2)" that records each finding as resolved or recorded | n/a | present: S5 to S10, the checked-and-sound list, and the Tasks 2-4 items, all closed | PASS |

Why the two FAILs: the verification of the verdict on fresh seeds 3-5 failed (28 of 31 control runs supported, max 0). The
recorded `passed` is `False` and the slow test re-measures the verification and asserts `ver.passed`, so it fails for the same
reason (`test_the_recorded_verification_follows_from_its_records` passes: the recorded `False` follows from its records).
Shreyan accepted a provisional-only verdict (ruling (a)): a supported verdict is capped at provisional, and since H-1 a failed
check is served `weak`. Both would pass only with a verdict rule that passes verification on seeds not used before, which
needs the pre-declared C1/C2 redesign recorded as follow-up work; none is planned for phase 3.

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
| `pivot/10-clonal-calibration` | #14 | merged (merge commit, 4d044c4), phase 1 Task 8 |
| `pivot/11-bench-real-rerun` | #15 | merged (merge commit, 37ec0d1), phase 1 Task 9 |
| `pivot/12-clonal-settings` | #16 | merged (merge commit, 1cf292a), phase 1 Task 10 |
| `pivot/13-run-trust-store` | #17 | merged (merge commit, 3c97dbe), phase 1 Task 11 |
| `pivot/14-job-priorities` | #18 | merged (merge commit, 0b04559), phase 1 Task 12 |
| `pivot/15-score-clonal-job` | #19 | merged (merge commit, 61687b2), phase 1 Task 13 |
| `pivot/16-run-trust-api` | #20 | merged (merge commit, 24183d7), phase 1 Task 14 |
| `pivot/17-cli-clonal` | #21 | merged (merge commit, 077f196), phase 1 Task 15 |
| `pivot/18-handoff-corrections` | #22 | merged (merge commit, dd6fc6a), phase 1 Task 16 |
| `pivot/19-phase1-security` | #23 | merged (merge commit, aac08ed), phase 1 Task 17 |
| `pivot/20-done-when-on-main` | #24 | merged (merge commit, 2217f76), phase 1 Done When re-run on `main` |
| `pivot/21-phase2-start` | #25 | merged (squash, e1c63d4), phase 2 start: current state, Shreyan's two rulings |
| `pivot/22-robustness-definitions` | #26 | merged (squash, c8a50d9), phase 2 Task 1 |
| `pivot/23-synth-controls` | #27 | merged (squash, 7257f5e), phase 2 Task 5 |
| `pivot/25-stability-engine` | #28 | merged (squash, 4387948), phase 2 Task 6 |
| `pivot/26-tree-support` | #29 | merged (squash, 14b8979), phase 2 Task 7 (with the PAGA no-edge fix) |
| `pivot/24-n2-amendment` | #30 | phase 2 N2a; merged with a merge commit so the recorded commit 2b15c77 stays on `main` |
| `pivot/28-stability-job` | #31 | merged (squash), phase 2 Task 10 |
| `pivot/27-verdict-calibration` | #32 | merged (merge commit, 67361c8, keeps the recorded commit 0a584d5), phase 2 Task 8 |
| `pivot/29-verdict-verification` | #33 | merged (merge commit, keeps the recorded commit 37197c8), phase 2 Task 9 (FAILED, verdict capped) |
| `pivot/30-robustness-api` | #34 | merged (squash), phase 2 Task 11 (stability and the capped verdict served; cell cycle deferred to Task 4) |
| `pivot/31-bench-real-stability` | #35 | merged (merge commit, 4961426, keeps the measured commit 99c3786), phase 2 Task 12 |
| `pivot/32-phase2-security` | #36 | merged (squash, 1167b2d), phase 2 Task 13 (security review; Done When recorded) |
| `pivot/33-phase2-final-review` | #37 | merged (squash, 35752e9), phase 2 final review recorded; H-1 asked of Shreyan |
| `pivot/34-weak-under-cap` | #38 | merged (squash, 53d0dc1), phase 2 H-1 (weak served as is under the verdict cap) |
| `pivot/35-cell-cycle-genes` | #39 | merged (squash), phase 2 Task 2 (Tirosh 2016 lists, verified mouse ortholog table) |
| `pivot/38-phase3-plan-x1` | #41 | merged (squash, 35eadc1), phase 3 plan amended to the verdict and interfaces phase 2 shipped (X-1) |
| `pivot/36-disconnected-knn-error` | #40 | merged (squash), clear error for a disconnected kNN graph (palantir, scfates_ppt) |
| `pivot/40-daniel-fetch` | #44 | merged (squash, 7d7a212), phase 3 Task 2, Daniel 2022 GEO files recorded (sizes, SHA-256, 23 contig files) |
| `pivot/39-registry-sources` | #43 | merged (squash, 4c05695), phase 3 Task 1, multi-file sources, licences, build step (phase 3 base 35ac459) |
| `pivot/37-cell-cycle-stage` | #42 | merged (squash, 8db251b), phase 2 Task 3 (report-only cell_cycle stage, left out of prep_key) |
| `pivot/42-cell-cycle-check` | #46 | merged (squash, 570fb5e), phase 2 Task 4 (per-run cell-cycle check, served with run and dataset; final-review doc fixes) |
| `pivot/41-sheih-verify` | #45 | merged (squash), phase 3 Task 9, Sheih 2020 TCR tables join to cells (PASS, share 0.99) |
| `pivot/43-daniel-metadata` | #48 | merged (squash), phase 3 Task 3, Daniel 2022 sample records and metadata inspected without touching clones |
| `pivot/44-contig-clones` | #47 | merged (squash, 6932b67), phase 3 Task 4, paired and TRB-only clone calling from 10x contigs |
| `pivot/45-phase2-close` | #49 | merged (squash), closes phase 2: Done When recorded on 570fb5e, final state |
| `pivot/47-sheih-declaration` | #51 | merged (squash, 5ef3151), phase 3 Task 10, case-study declaration (clone definition, namespace, strata, root) before any analysis |
| `pivot/46-daniel-formats` | #50 | merged (squash), phase 3 Task 5a, Daniel 2022 file formats and the capture map (23 members joined) |
| `pivot/48-sheih-builder` | #52 | merged (squash), phase 3 Task 11, case-study fetch and build (62,167 cells, 42,362 paired-labelled; no analysis) |
| `pivot/49-daniel-builder` | #53 | merged (squash, 30b98c9), phase 3 Task 5b, the Daniel 2022 builder (contigs joined to cells by capture, clones called; synthetic tests only) |
| `pivot/50-daniel-declaration` | #54 | merged (squash, fb06f05), phase 3 Task 6, the Daniel 2022 declaration before any run (subset, namespace, clone definition, truth, gates; clonal not assessable) |
| `pivot/51-bench-real-clonal` | #55 | merged (squash, a7d69c4), phase 3 Task 7, clonal scores and reported fields for real datasets (synthetic `fake` fixture only; Daniel not built) |
| `pivot/52-bench-real-daniel` | #56 | open, phase 3 Task 8, the real benchmark with Daniel 2022 rows (`real-2026-11-20.json`); gates not measured, clonal not assessable |

## Phase 3 Task 8 decisions and result

- **Branch** `pivot/52-bench-real-daniel`, from `main` a7d69c4. Step 1: all four datasets were downloaded, so no
  download was asked for. Step 3: `bench fetch daniel2022` built the derived file from cached sources in about
  10 s (peak resident memory about 1.8 GB when sampled, far under the 5.5 GB estimate); no memory problem.
- **Build shares** (`1 - join_misses / contig_barcodes`, all at or above 0.5; `tests/real` `test_the_daniel_build_joins_every_capture` passes):
  GSM5687827 0.904, 828 0.794, 829 0.954, 830 0.855, 831 0.961, 832 0.875, 833 0.582 (lowest, D8 gp33-), 834 0.783.
  `n_cells` 33,188, equal to Task 3's total for the kept captures (per-capture cells 1989, 5875, 2057, 8914, 2031,
  8194, 2526, 1602). 26,714 cells carry a `clone_id`, 29,271 a `clone_id_trb`. The 0.5 was not touched.
- **Run:** `make bench-real` exited 0 (7 real-data tests passed, then `bench real`). The results file
  `benchmarks/results/real-2026-11-20.json` records `commit` a7d69c4, which is `main`, not a branch commit, so no
  merge commit is needed.
- **Daniel rows (all three methods):** `clonal_status == "not_assessable"`, `clonal_reason` starts `declared 2026-11-20`,
  every gated lineage metric NaN and rendered "lineage gates: not measured (no gated lineage, declared
  2026-11-20)", never PASS or FAIL. Stability ok: paga_dpt pt 0.989 / lineage 0.833, palantir 0.989 / 0.667,
  scfates_ppt 0.990 / 1.0. Verdict level `provisional` for all three (tree_likeness 0.632); the level is never
  `supported` because the verdict verification failed on 2026-11-17, and `provisional` reads as "not shown to be a
  tree", not as support. Cluster counts 19 (spurious 6, 6, 2).
- **Reported fields (not gated):** `reported_day_spearman` 0.291 (paga_dpt), 0.341 (palantir), 0.125 (scfates_ppt);
  identical across methods: `reported_clone_recurrence` 0.0612, `reported_labelled_share_of_raw` 0.8979 (26,714 over
  29,750, recomputed independently), `reported_n_clone_strings` 7,242, `reported_raw_clonotype_ari` 1.0,
  `reported_trb_status` not_assessable. The ARI of 1.0 means clone_id and the 10x raw clonotype id partition the
  jointly labelled cells identically within a capture. That is true by construction, not independent evidence:
  the builder derives both ids from the same contigs, and among the 26,714 cells with both they map 1:1 within
  capture. The ARI shows the clone call reproduces the 10x grouping where it labels a cell; it hides the 3,036
  raw-labelled cells (10.2%) that fail the qualifying rule, which the labelled share (0.898) shows. Day Spearman is
  low to moderate and is reported as is.
- **Root sits late in every pseudotime (ungated, unexplained):** `root_pt_median` of the `TEx_prog` root cells is
  0.70 (paga_dpt), 0.80 (palantir) and 0.984 (scfates_ppt), against 0.07-0.58 on the other datasets. With no gated
  lineage `root_pt_excess` is NaN and no gate reads it, so "stability ok, provisional" must not be read as a
  healthy trajectory: the root, though explicit, is not early in any method's ordering. Reported as is.
- **Review (controller, PR #56):** `data:validate-data` on the Daniel rows: values are finite wherever a status is
  ok (stability, tree_likeness, reported fields), shares lie in [0, 1], and the reported fields match their
  declared definitions (labelled share and recurrence recomputed from the built file by the review seat). The real
  build test now also checks the listed captures are the declared ones, every declared capture has cells, the
  label/day/capture/donor columns are populated, some cells carry each clone id, and `n_cells` equals the obs count.
- **Step 7:** pancreas, bone marrow and gastrulation metrics are identical to `real-2026-11-17.json`
  (runtime aside). The only difference is `clonal_reason`, `null` before and `""` now (Task 7, always a string).
- **Not assessable checks:** README `grep -c "clonal: not assessable (declared"` prints 3. Gate failures on the other
  datasets stay failures. Nothing was re-run with other settings.

## Phase 3 Task 9 decisions and result

- **Owner's OK (STOP steps):** Shreyan, in chat on 2026-11-19, "Yes, approve all" (relayed by the parent session). It covered the largest post-infusion Sheih clonotype table, the first line of `GSE125881_raw.expMatrix.csv.gz` (streamed, capped at 16 MB) and the Step 5 whole-matrix fetch (not needed). The Cao 2025 fallback was approved only for its metaData, cells_colnames and aggregated contig files, and only if Sheih failed.
- **Result: PASS.** CLL-2 d12 (GSM3583882), 5,817 of 5,873 distinct table barcodes found among the matrix cell IDs, share 0.99 (rule: at least 0.5). Mapping: `<16-nt core>-1` becomes `<16-nt core>-<aggr ID>`. The Cao fallback was not needed; nothing from GSE207192 was downloaded.
- **Approval wording:** the owner's approval named "the largest post-infusion clonotype table" for Task 9 and `GSE125881_RAW.tar` for Task 10. The tables exist only inside that 2 MB tar, so Task 9 fetched the tar; the controller ruled this within the approval (ledger P8).
- **Listing differs from the plan:** the suppl folder lists only `GSE125881_RAW.tar` (2,048,000 bytes, 16 clonotype tables inside, 6.9 KB to 297 KB each) and the matrix (114,573,707 bytes). The tables exist only inside the tar, so the whole 2 MB tar was fetched (SHA-256 `769fa5b71b12b7efe433f51164d4370e0af00c20bebc647029cc4fbd9a943487`) and read in memory, never extracted. It sits at `/private/tmp/claude-501/-Users-shreyan-lumenCell-V2/453b6f74-ae4a-4614-a97a-fbc583e58d2c/scratchpad/sheih-verify/GSE125881_RAW.tar` with the matrix header beside it (`header.txt`, SHA-256 `23955ade...`), not in the dataset cache; the chosen member's SHA-256 is `93f8ef95...` and the largest table is 305,886 bytes (members run 6,905 to 305,886). Full hashes and a negative control (wrong suffix gives 0.00) are in the research note; Tasks 10 and 11 re-verify the tar's hash when they fetch it again.
- **Matrix is genes x cells** (62,167 cells, aggr IDs 1-16), so no whole-matrix fetch. Tables are one row per clonotype x chain, `;`-separated barcode lists, no productive flag. Details, the GSM table and the licence lines are in [research/2026-11-19-sheih-verification.md](research/2026-11-19-sheih-verification.md); Tasks 10 and 11 read it.
- **Script:** `benchmarks/inspection/barcode_join.py` prints headers, dtypes, counts, barcodes and the share, never a receptor value.

## Phase 3 Task 10 decisions

- **Declared before analysis:** `[case_study.sheih2020]` in `benchmarks/thresholds.toml`, with a 2026-11-19 changelog entry giving each value and its reason. Written from the paper and Task 9's note only; nothing was downloaded and nothing was computed from clonotype, CDR3 or gene values. Pinned by `tests/unit/test_case_study_declaration.py`.
- **Paper read:** yes, Europe PMC's open full text of PMC6954177 (nature.com redirected to a sign-in step, which was not followed; no CAPTCHA was met or bypassed).
- **Root `cluster:timepoint=IP`.** The paper names no differentiation root (its Monocle 3 trajectory is on infusion-product cells only and gave no clear trajectories). It treats the infusion product as the source of the in vivo repertoire (Results, "Contribution of infused CAR-T clusters to in vivo repertoire", Fig. 7), which is also the plan's fallback (the earliest sample). The order comes from the trial design, never from clones. Task 11 must write obs `timepoint` with GEO's "time point" values as written (IP, d12, d21, d28, d29, d38, d83, d89, d102, d112); every patient has `IP`.
- **Chains `paired`** (Daniel's rule, 1-2 of each of TRA and TRB). Judgement call: the GEO tables are labelled "unpaired" (one row per clonotype x chain), but every row names its chain and lists its cells' barcodes, so the chains group by cell barcode, which is the plan's test for `paired`. The paper also used Cell Ranger vdj's paired clonotype calling. Cost: cells with one chain detected are unlabelled.
- **CDR3 `nt`** (`cdr3_nts` is present). **Filters `[]`:** productive, high_confidence and full_length (and is_cell) are all absent from the tables; whether the authors kept only productive chains is UNVERIFIED.
- **Unit:** per patient primary (4 patients x 3 methods = 12 runs with stability, about a quarter of the cells each), patient as the clone namespace, `strata_keys = ["timepoint"]` (N3), pooled run as a labelled diagnostic only.
- **Review fixes (controller):** the changelog's root bullet now says IP as the start is an assumption, not named by the paper; `test_the_declaration_has_its_dated_changelog_entry` ties the entry's date to `declared` and checks it states the assumption and the per-patient cost (a re-dated entry fails it).

## Phase 3 Task 11 decisions

- **Owner's OK (STOP step 1):** Shreyan, in chat on 2026-11-19, "Yes, approve all" (relayed by the parent session), item 3: "Task 10: the rest of GSE125881 (raw.expMatrix.csv.gz and GSE125881_RAW.tar), ~111 MB". Nothing else was downloaded (no Cao files).
- **Downloads** (`fetch_source` only, into the datasets cache, 2026-11-20). HEAD sizes equal Task 9's: matrix 114,573,707 bytes, tar 2,048,000 bytes. SHA-256: matrix `243e6e2dd09158e08dcaecfbf35b23122ac3096c1039355b7590474e76d8c02c`; tar `769fa5b71b12b7efe433f51164d4370e0af00c20bebc647029cc4fbd9a943487`, equal to Task 9's. Both are in `SOURCES` of `bench/builders/sheih2020.py`. The tar is read in memory, never extracted.
- **Build:** `uv run lumencell bench fetch sheih2020` builds `sheih2020-9720dc62e9bc.h5ad` (build version 2 after fix round 1; the version 1 file `78b787b0a59c` was removed; 167 s, 3.4 GB peak, the matrix read 500 genes at a time into CSR). A case study is fetched by its id; `bench fetch all` still covers `DATASETS` only.
- **n_cells = 62,167, equal to the paper's 62,167** (20,223 genes, 116,507,933 non-zeros). Every cell maps to a sample by its aggr suffix; an unknown suffix would be a `BuildError`.
- **Labelled totals (the declared rule):** 42,362 cells with a paired `clone_id` (1-2 distinct TRA and 1-2 distinct TRB CDR3 nt) and 50,903 with a TRB-only `clone_id_trb`. Per patient (paired / TRB-only): CLL-1 12,745 / 14,951; CLL-2 9,310 / 11,521; NHL-6 13,575 / 16,807; NHL-7 6,732 / 7,624. No clone x time point or clone x cluster table was computed (case-study freeze).
- **Joined share per sample** (distinct table barcodes found among the sample's cells; floor 0.5, enforced by the build): minimum 0.974 (GSM3583873, NHL-7 IP); all 16 are between 0.974 and 1.000. CLL-2 d12 is 0.990 (56 misses of 5,873, so 5,817 found), as in Task 9. The distinct-barcode and miss counts per sample are in `uns["lumencell_build"]` (`contig_barcodes`, `join_misses`, `joined_share`, `cells_per_sample`).
- **Clones from the tables.** Each table row (clonotype x chain) is expanded to one contig-style row per listed barcode (`<core>-1` becomes `<core>-<aggr id>`), flags set to "True" (the tables have none, filters = []), and handed to `call_clones`, so the `clone_id` and `clone_id_trb` strings are hashed exactly as Daniel's. `raw_clonotype_id` is dropped: the table's clonotype ids are per chain, so a cell sits under several of them and none is the cell's clonotype.
- **Multi-value CDR3 fields.** The format check (not values) found none: all 41,756 rows hold one `[ACGTN]+` string in `cdr3_nts` and a chain of TRA or TRB. A row whose chain is not TRA/TRB or whose `cdr3_nts` is not one bare `[ACGTN]+` string (a `;`-list, an empty field) refuses the build (fix round 1; it was counted and dropped), because a dropped row could make a 3-TRA cell look like a 2-TRA one and label it. A `;`-list is never split. `cdr3_unusable_rows` stays in `uns["lumencell_build"]` (0 in every sample). The `barcodes` field is the multi-valued one: it is split on `;`; a barcode not of the form `<16 nt>-1` stays in the denominator as a miss (it does not raise by itself).
- **Cells with more than 2 distinct CDR3 on a chain** are unlabelled in the paired call (declared `max_per_chain = 2`); a cell with no table row gets NaN and 0 counts.
- **obs:** `patient` (the clone namespace, passed as a clonal column), `timepoint` (GEO's value as written: IP, d12, d21, d28, d29, d38, d83, d89, d102, d112; a stratum, visible), `day` (float; the sampling day in the label, NaN for IP, which is not a sampling day), `clone_id`, `clone_id_trb`, `receptor_n_tra`, `receptor_n_trb`. The clone and receptor columns match `RECEPTOR_COLUMN`. obs_names are the matrix cell ids. There is no `sample` column: the GSM encodes the patient, so it would leak the namespace to methods in the pooled diagnostic; the GSM appears only in `uns["lumencell_build"]` (`cells_per_sample` and the per-sample counts).
- **Untrusted values:** cell ids, table barcodes and tar member names are validated by regex or compared to the fixed member list, never used as paths, and never put in an error message; params the builder cannot honour (chains other than paired, cdr3 other than nt, any filter, max_per_chain other than 2) raise `BuildError`.
- **Tests:** `tests/unit/test_sheih2020.py` (fixture: 5 genes, 6 cells, 2 samples, a duplicate barcode core across samples, a stray `../` tar member) and `test_case_study_files_match_their_checksums` in `tests/real` (real_data).

- **Fix round 1** (controller rulings P9 to P11, plus the mutation seat's tests):
  - **P9:** `day` is float with NaN for IP. Day 0 would enter a day Spearman and, with root = IP, correlate by construction. **Task 12 must use `timepoint` only as a stratum or root, never its category codes as time** (categories are not time-ordered), and `day` only over non-IP cells.
  - **P10:** `cdr3_unusable_rows > 0` in any sample raises `BuildError` (a TRG row, a `;`-list, an empty or non-ACGTN CDR3). The real data has 0. A CDR3 containing `N` is used (the call-clones alphabet is ACGTN).
  - **P11:** obs `sample` dropped (leaks the patient). Duplicate gene names in the matrix raise `BuildError`. `BUILD_VERSION` is 2 (the output changed), so the derived file name changed.
  - Tests added: all 16 samples pinned literally from Task 9's GSM table (VDJ GSM and member names written out); two-digit aggr ids (`-10`, `-12`); the join-share floor at 0.75 / 0.5 builds and 0.43 refuses; matrix 1.5, -1 and ragged rows; extra params refused; unknown-suffix error without barcode text; missing member; duplicate cell ids; `chunk_rows < 1`; missing required column; zero usable barcodes; a `-2` table barcode counts as a miss; `cells_per_sample`; the two SHA-256 literals, sizes and URLs pinned.
  - Mutants seen to fail: `group(2)[:1]`, `<` to `<=` on the floor, subset instead of equal params, `chunk_rows < 0`, dropped negative or fractional check, no unknown-suffix check.

## Verify the current state

```bash
git fetch origin && git rev-list --left-right --count main...origin/main
make install
make ci
make bench-synth
```
