# Handoff: T-cell clonotype work

A new session starts with: "Read docs/superpowers/HANDOFF.md and continue."
This file is updated inside every task's PR. `docs/HANDOFF.md` is the milestone 1 note; its
stale lines get corrected once, at the end of phase 1.

## Current state (2026-11-12)

- **Phase:** design, specifically approval gate 1, the spec.
- **Spec:** [specs/2026-11-12-tcell-clonotypes-design.md](specs/2026-11-12-tcell-clonotypes-design.md),
  on branch `pivot/01-tcell-spec` with its own PR.
- **Next:** once Shreyan approves the spec PR, squash-merge it. Then run
  `superpowers:writing-plans` and write one plan per phase in `docs/superpowers/plans/`. Each
  plan has Global Constraints, Review Focus, Not in This Plan, a File Map and Done when, with
  tasks small enough for one subagent each. Open the plans as a PR and stop for approval
  gate 2.

## Done

- Reconnaissance by read-only subagents covering loader, bench, prep, API, robustness-2 and
  web research. The research notes are in `docs/superpowers/research/`.
- Baseline at 865ce22: `make install` exit 0. `make ci` exit 0 in 192 s, with 353 Python
  tests passed and 12 deselected, 36 web tests and 1 Playwright test.
- Spec revision 2, after an independent method review and an architecture review (both on
  Opus). Spec §11 lists what each review changed. The two method blockers were receptor genes
  leaking clones in through HVGs, and the need for a lineage-free baseline, because T-cell
  clones are transcriptionally coherent.
- Shreyan's answers, all three the recommended options:
  - **Daniel truth:** the paper's state labels, with the root at Tpex.
  - **Case study:** Sheih 2020, after verifying it has per-cell barcodes; Cao 2025 is the
    fallback.
  - **Stability:** runs automatically as a follow-on job.

## Found along the way

- **Local `main` was 13 commits behind `origin/main`** on 2026-11-12, so it had no CI and no
  `docs/backend-robustness-2.md`. It was fast-forwarded to 865ce22. Check
  `git rev-list --left-right --count main...origin/main` before reading code.
- **`main` has no branch protection.** Keeping work on PRs is discipline, not enforcement.
- **The `need_from_me` section of the brief was missing.** It is taken to mean dataset
  downloads and new runtime dependencies (spec §8, decision 16).
- **Large categorical columns are described one way and served another.** Above 1,000
  categories, `dataset_view.py` lists a column as `string` but serves it as categorical codes.
  Phase 1 fixes this.
- **Gate failures don't fail the benchmark commands.** `lumencell bench real` and
  `bench synth` exit 0 on gate failures; only regressions fail `bench-synth`.

## Decisions

See spec §8 for each decision and what it costs if wrong.

## Open questions for Shreyan

- None beyond the spec review.

## Branches and PRs

| Branch | PR | State |
|---|---|---|
| `pivot/01-tcell-spec` | see `gh pr list` | open, awaiting approval |

## Verify the current state

```bash
git fetch origin && git rev-list --left-right --count main...origin/main
make install
make ci
make bench-synth
```
