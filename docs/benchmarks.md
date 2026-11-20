# Benchmarks

Generated from the committed results by `uv run lumencell bench report --readme`; do not edit by hand.

<!-- bench:start -->
**Synthetic** (4 topologies × 3 seeds, `realistic` preset; baseline from commit `857d03d`). Ranges over all cases:

| Method | Pseudotime (global) | Pseudotime (per lineage) | Branch accuracy | Clonal excess ratio | Gates failed |
|---|---|---|---|---|---|
| PAGA + DPT | 0.86–0.99 | 0.95–0.99 | 0.81–0.99 | 0.301–0.542 | 9 of 12 cases |
| Palantir | 0.90–0.99 | 0.91–0.99 | 0.90–0.99 | 0.367–0.564 | 7 of 12 cases |
| scFates | 0.87–0.99 | 0.98–0.99 | 0.97–1.00 | 0.443–0.616 | 7 of 12 cases |

Cases failing each gate (a case can fail several): PAGA + DPT: `branch_accuracy` 3, `clonal_excess_ratio` 6, `pt_spearman_global` 3, `root_pt_excess` 1; Palantir: `clonal_excess_ratio` 7, `root_pt_excess` 2; scFates: `clonal_excess_ratio` 5, `pt_spearman_global` 2.

The clonal gate is not a branch gate: a method that merges two true lineages can score a clonal ratio above 1, and lineage recall catches that.

**Real datasets** (`benchmarks/results/real-2026-11-20.json`, commit `a7d69c4`). Spearman is against stage order (embryo time for gastrulation):

| Method | bonemarrow | daniel2022 | gastrulation_erythroid | pancreas |
|---|---|---|---|---|
| PAGA + DPT | Spearman 0.73, branch 0.62, terminal 0.99 (fail: branch_accuracy, lineage_recall), verdict not_supported | Spearman –, branch –, terminal – (lineage gates: not measured (no gated lineage, declared 2026-11-20)), verdict provisional, clonal: not assessable (declared 2026-11-20, fewer than two gated lineages, so no clonal comparison between lineages; namespace = capture (rule b, state-gated GSM5687824, GSM5687825 and GSM5687826 dropped): one mouse per capture is assumed; the records do not say which mice a capture holds or whether a capture pools mice) | Spearman 0.75, branch –, terminal 0.94 (pass), verdict provisional | Spearman 0.95, branch 0.63, terminal 0.93 (fail: branch_accuracy, lineage_recall), verdict provisional |
| Palantir | Spearman 0.73, branch 0.94, terminal 0.99 (pass), verdict provisional | Spearman –, branch –, terminal – (lineage gates: not measured (no gated lineage, declared 2026-11-20)), verdict provisional, clonal: not assessable (declared 2026-11-20, fewer than two gated lineages, so no clonal comparison between lineages; namespace = capture (rule b, state-gated GSM5687824, GSM5687825 and GSM5687826 dropped): one mouse per capture is assumed; the records do not say which mice a capture holds or whether a capture pools mice) | Spearman 0.72, branch –, terminal 0.58 (fail: terminal_share), verdict provisional | Spearman 0.94, branch 0.84, terminal 0.94 (pass), verdict provisional |
| scFates | Spearman 0.67, branch 0.79, terminal 0.95 (fail: pt_spearman_lineage), verdict provisional | Spearman –, branch –, terminal – (lineage gates: not measured (no gated lineage, declared 2026-11-20)), verdict provisional, clonal: not assessable (declared 2026-11-20, fewer than two gated lineages, so no clonal comparison between lineages; namespace = capture (rule b, state-gated GSM5687824, GSM5687825 and GSM5687826 dropped): one mouse per capture is assumed; the records do not say which mice a capture holds or whether a capture pools mice) | Spearman 0.80, branch –, terminal 0.61 (fail: terminal_share), verdict provisional | Spearman 0.94, branch 0.75, terminal 0.80 (fail: branch_accuracy), verdict provisional |

Reported, not gated (daniel2022):

| Method | clone recurrence | day spearman | labelled share of raw | n clone strings | raw clonotype ari | trb excess ratio | trb p | trb status | cells | effective clones | η at m = 3 | p at m = 3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PAGA + DPT | 0.06 | 0.29 | 0.90 | 7242 | 1.00 | – | – | not_assessable | 0 | – | – | – |
| Palantir | 0.06 | 0.34 | 0.90 | 7242 | 1.00 | – | – | not_assessable | 0 | – | – | – |
| scFates | 0.06 | 0.12 | 0.90 | 7242 | 1.00 | – | – | not_assessable | 0 | – | – | – |

A supported verdict is capped at provisional, because the verdict's verification failed on fresh seeds (2026-11-17): at x1.0 the tree-likeness check (C2) separates blob and sheet from the matrix in the recorded calibration but not loop or two_groups, and the stability checks do not separate trees from the controls. A weak or not_supported verdict is a real failure and is not capped. The verdict is reported, not gated, and is not evidence of a tree.
<!-- bench:end -->
