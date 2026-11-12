# Handoff: T-cell clonotype work

A new session starts with: "Read docs/superpowers/HANDOFF.md and continue."
This file is updated inside every task's PR. `docs/HANDOFF.md` is the milestone 1 note; its
stale lines are corrected once, by phase 1's handoff-corrections task.

## Current state (2026-11-12)

- **Phase:** phase 1 has started. Both approval gates are passed: Shreyan approved the spec
  (PR #5) and the plans (PR #6, merged as 153a7f0) on 2026-11-12. Task 1 (clonal definitions,
  gates and runtime constants) is on branch `pivot/03-clonal-definitions`.
- **Spec:** approved and merged in PR #5 (142a788):
  [specs/2026-11-12-tcell-clonotypes-design.md](specs/2026-11-12-tcell-clonotypes-design.md).
  Branch `pivot/02-phase-plans` amends it in two places (see Decisions).
- **Plans**, merged in PR #6:
  - [phase 1, clonotypes end to end](plans/2026-11-12-phase-1-clonotypes.md): 17 tasks
  - [phase 2, T-cell robustness](plans/2026-11-12-phase-2-robustness.md): 13 tasks
  - [phase 3, real T-cell data](plans/2026-11-12-phase-3-tcell-data.md): 15 tasks
- **Next:** phase 1 Task 2, after Task 1's PR is reviewed and squash-merged. Execute phase 1 with `superpowers:subagent-driven-development`, which Shreyan chose:
  - one task per subagent, on the model the task names
  - TDD throughout
  - a fresh reviewer per task, plus the pr-review-toolkit agents on each PR
  - a squash-merge once CI is green
  - branches `pivot/NN-<slug>`, numbered from 03 in execution order
- **Handoffs:** at the end of every phase, and whenever the session nears 35% of its context.

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
- **Gate failures don't fail the benchmark commands.** `lumencell bench real` and
  `bench synth` exit 0 on gate failures; only regressions fail `bench-synth`.

## Open questions for Shreyan

- Commit 33c281d and PR #5 were created at 07:57 under your account, not by the orchestrating
  session. Was that you, or a tool?

## Branches and PRs

| Branch | PR | State |
|---|---|---|
| `pivot/01-tcell-spec` | #5 | merged (squash, 142a788) |
| `pivot/02-phase-plans` | #6 | merged (squash, 153a7f0) |
| `pivot/03-clonal-definitions` | see the PR list | phase 1 Task 1, open for review |

## Verify the current state

```bash
git fetch origin && git rev-list --left-right --count main...origin/main
make install
make ci
make bench-synth
```
