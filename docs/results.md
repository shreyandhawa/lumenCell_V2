# Results so far

This page says what the first real-data case study found, and what it does not show.
The case study is on Sheih et al. 2020, a study of CAR-T cells in four patients.

## What was declared before the data

The analysis plan was written down before any real data were scored. It was merged in
[PR #51](https://github.com/shreyandhawa/lumenCell_V2/pull/51), at commit 5ef3151.
The plan was written from the paper and the file headers only.

- A clone is a set of cells that share the same paired T-cell receptor chains.
- Clones are counted within one patient. The same receptor in two patients is two clones.
- The clone null shuffles labels within each time point, so clones that grow at one time point do not count as lineage agreement.
- The trajectory starts from the infusion product cells. This start is an assumption.
- Each patient is run alone, with each of three methods. That is 12 runs.
- A run with all four patients pooled is a diagnostic only, never the headline.
- Clone size 2 is the main test. Clone size 3 is a sensitivity check. There are 1000 shuffles, with seed 0.

## What failed

The verdict rule failed its check on fresh seeds, on 2026-11-17. Its verification ran on
31 negative-control runs. 28 of 31 were rated supported. The rule was not changed to fix this.

Because of that failure:

- No run is shown as supported.
- A supported verdict is capped at provisional ("Not shown to be a tree").
- A failed check is shown as weak ("Clones do not back these lineages") or not_supported ("Not supported: cells left without a pseudotime").

## What was found

Sheih 2020, per patient. There are 12 method runs: 4 patients, 3 methods.

- 8 of 12 are weak ("Clones do not back these lineages").
- 4 of 12 are provisional ("Not shown to be a tree").
- None is supported. None is not_supported.

| Patient | paga_dpt | palantir | scfates_ppt |
|---|---|---|---|
| CLL-1 | provisional | weak | weak |
| CLL-2 | provisional | weak | weak |
| NHL-6 | weak | weak | weak |
| NHL-7 | weak | provisional | provisional |

A run passes the clone check when its headline p is at most 0.01 and its lift over a
k-means clustering baseline (Δη) is above 0. Four of the 12 runs pass it.

## What it does not show

- **It does not show a tree.** Provisional means "not shown to be a tree". It is never support.
- **The verdict is not evidence of a tree.** Its rule cannot tell a tree from these controls. Supported means "connected", not "a tree".
- **Clones give no direction.** Shared clones say two states are related. They never say which came first.
- **The direction comes from an assumed root.** The paper's main text and Methods name no trajectory root (the supplementary figures were not read).
- **It is four patients.** Two have CLL and two have NHL. Nothing here applies beyond them.
- **There is no batch correction.** Pooled runs may mix patient structure into the lineages, so the pooled run is a diagnostic only.
- **Clone scores fall under stricter nulls.** In every per-patient run, the clone score against the time-point null or the pseudotime-depth null is below half its value against the within-patient null.
- **The lineages add little.** Expression clustering matches most of the clone score. Δη is at most 0.050 in absolute value.
- **Each method scores its own cells.** So the three methods are not compared like for like.
- **Blood only.** The infusion product is a cultured, stimulated cell population.

Full details: [the case study](case-studies/sheih-2020.md) and [the handoff](superpowers/HANDOFF.md).
