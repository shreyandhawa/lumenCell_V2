"""Clonal agreement: T-cell receptor clones as evidence for a run's lineages.

Cells that share a receptor descend from one ancestor, so a run whose lineages keep each clone
together agrees with evidence the inference never saw. Clone labels are used only for scoring.

This package holds pure functions on arrays (numpy, scipy and pandas, plus the runtime constants
in ``lumencell.trust.constants``). It imports nothing from the engine, the API, the benchmark or
anndata, so the benchmark (oracle-cells mode) and the worker (own-cells mode) run the same code.

- ``labels``: raw labels and donors to clone, donor and stratum codes. A clone is the pair
  (donor, label).
- ``scores``: normalised entropy and purity, the effect η against three nulls, a k-means
  baseline that uses no lineages, and lineage-pair coupling.

Every definition was declared in ``benchmarks/thresholds.toml`` (``[clonal]``) before any score
was computed, and spec §4.3 states them in full.

Sources:

- STARTRAC: Zhang et al., Nature 2018, doi:10.1038/s41586-018-0694-x. Purity of a clone's
  distribution over states.
- LRT: Xie et al., PLoS Comput Biol 2023, doi:10.1371/journal.pcbi.1011300. Lineage bias of
  clones.
- CoSpar: Wang et al., Nat Biotechnol 2022, doi:10.1038/s41587-022-01209-1. Coupling of states
  through shared clones.
- Vinh, Epps & Bailey, JMLR 11(95):2837-2854, 2010, jmlr.org/papers/v11/vinh10a.html (no DOI).
  Chance-corrected information measures, the model for η.
- Phipson & Smyth, Stat Appl Genet Mol Biol 9(1), 2010, doi:10.2202/1544-6115.1585.
  Permutation p-values that are never zero.
- Sturm et al. (scirpy), Bioinformatics 36(18):4817-4818, 2020, doi:10.1093/bioinformatics/btaa611,
  and Schattgen et al. (CoNGA), Nat Biotechnol 40:54-63, 2022, doi:10.1038/s41587-021-00989-2.
  Clones are transcriptionally coherent, hence the lineage-free baseline.
"""
