# Clonal coherence metric: prior work and candidate definitions

Date: 2026-11-12. Purpose: exact method definitions from primary sources before designing lumenCell's "clonal coherence" metric.

## How to read the verification marks

- VERIFIED = read in the primary source (docs page, source code, or API record) during this research.
- UNVERIFIED = not read from a primary source. It is recalled, inferred, or reported only by a secondary summary.
- Method note: pages were read through a summarising fetch tool (a small model paraphrases the page), so wording is paraphrase, not quotation. Source code was downloaded and read directly for STARTRAC and CoSpar. The Nature, PMC, PubMed and 10x pages blocked automated fetching (bot check, 429 or 403) and were not worked around.
- Several items below are my design recommendations. They are labelled "RECOMMENDATION" and are not from any source.

---

## 1. scirpy

Sources:
- Docs: https://scirpy.scverse.org/
- Data structure: https://scirpy.scverse.org/en/latest/data-structure.html
- define_clonotypes: https://scirpy.scverse.org/en/latest/generated/scirpy.tl.define_clonotypes.html
- define_clonotype_clusters: https://scirpy.scverse.org/en/latest/generated/scirpy.tl.define_clonotype_clusters.html
- chain_qc: https://scirpy.scverse.org/en/latest/generated/scirpy.tl.chain_qc.html
- clonal_expansion: https://scirpy.scverse.org/en/latest/generated/scirpy.tl.clonal_expansion.html
- Tutorial: https://scirpy.scverse.org/en/latest/tutorials/tutorial_3k_tcr.html
- Source: https://github.com/scverse/scirpy (file `src/scirpy/tl/_clonotypes.py`)
- PyPI JSON: https://pypi.org/pypi/scirpy/json
- Paper: Sturm et al., Bioinformatics 2020;36(18):4817-4818, DOI 10.1093/bioinformatics/btaa611 (DOI from a search result; paper text not read; PMC7751015 per search).

### Version, licence, dependencies (VERIFIED from the PyPI JSON and GitHub API)
- Current version 0.26.1, uploaded 2026-11-05. `requires_python >= 3.12`.
- Licence: BSD 3-Clause (PyPI and GitHub API agree).
- Core `requires_dist`: adjusttext, airr>=1.6.1, anndata>=0.9, awkward>=2.1, igraph, joblib, logomaker, mudata>=0.3, networkx, numba>=0.60, numpy, pandas>=2.3, pooch, python-levenshtein, scanpy>=1.9.3, scikit-learn, scipy, scverse-misc, squarify, tqdm (pycairo on Windows only).
- Extras: cupy, dandelion (sc-dandelion), diversity (scikit-bio), parasail, rpack (rectangle-packer).
- Weight as an optional extra: heavy (scanpy, numba, igraph, awkward, mudata, scikit-learn, plus Python 3.12 or newer).
- RECOMMENDATION: lumenCell should not import scirpy. A user's `.h5ad` carries the derived `clone_id` column in `obs` after a scirpy run, so the metric only needs pandas and numpy on the columns. Offer scirpy only as a documented way to produce those columns.

### Data layout
- VERIFIED: since v0.13.0, receptor data is an awkward array in `adata.obsm["airr"]`. Axis 0 is cells (aligned to `obs`), axis 1 is chains per cell (variable length), axis 2 is a record of AIRR Rearrangement fields (for example `locus`, `junction_aa`). Earlier versions expanded these into `obs` columns.
- VERIFIED: `obsm["chain_indices"]` holds VJ and VDJ index lists and a `multichain` flag. It is created by `ir.pp.index_chains()`. Chains are separated into VJ (TRA, TRG, IGK, IGL) and VDJ (TRB, TRD, IGH). Up to two chain pairs per cell are modelled ("dual IR").
- VERIFIED: `ir.get.airr()` and `ir.get.airr_context()` temporarily expose chain fields in `obs`.
- VERIFIED: the recommended container is MuData with modalities `gex` and `airr`. Scirpy functions default to the `airr` modality. Results are written to `mdata.mod["airr"].obs[key]` and to `mdata.obs["airr:key"]`. A plain AnnData with `obsm["airr"]` also works. `mdata.pull_obs()` exposes `modality:column` names.

### Columns written by common calls (VERIFIED)
- `ir.tl.chain_qc` writes `receptor_type`, `receptor_subtype` and `chain_pairing` (`key_added` can change the names).
  - `receptor_type`: TCR | BCR | ambiguous | multichain | no IR.
  - `receptor_subtype`: TRA+TRB | TRG+TRD | IGH | IGH+IGL | IGH+IGK | IGH+IGK/L | multichain | ambiguous | no IR.
  - `chain_pairing`: single pair | orphan VJ | orphan VDJ | extra VJ | extra VDJ | two full chains | ambiguous | multichain | no IR.
- `ir.tl.define_clonotypes(adata, *, key_added="clone_id", ...)` writes `clone_id` and `clone_id_size` (`f"{key_added}_size"`, the number of cells in the cell's clonotype). The distance result goes to `adata.uns[key_added]`.
  - It is an alias of `define_clonotype_clusters` with stringent defaults: nucleotide CDR3 identity, connected components of the distance graph.
  - The tutorial calls it with `receptor_arms="all", dual_ir="primary_only"`.
- `ir.tl.define_clonotype_clusters(...)` writes `cc_{sequence}_{metric}` plus `_size`. Examples from the tutorial: `cc_aa_tcrdist`, `cc_aa_tcrdist_size`, `cc_aa_tcrdist_same_v`.
  - Defaults: `sequence="aa"`, `metric="identity"`, `receptor_arms="all"`, `dual_ir="any"`, `within_group="receptor_type"`, `partitions="connected"`.
  - The tutorial uses `ir.pp.ir_dist(metric="tcrdist", sequence="aa", cutoff=15)` first. The cutoff belongs to `ir_dist`, not to this function.
  - These are looser, similarity-based groups. They are not evidence of a common ancestor.
- `ir.tl.clonal_expansion(adata, target_col="clone_id", expanded_in=None, breakpoints=(1, 2), key_added="clonal_expansion")` writes categories "1", "2" and ">2" (three or more) by default. `expanded_in=<sample column>` computes expansion within each sample.
- Other `obs` columns seen in tutorials: `clonal_expansion`, `is_convergent`, `clonotype_modularity`, `clonotype_modularity_fdr`. `has_ir` is not on the pages read (UNVERIFIED as a current column). `IR_VJ_1_*` and `IR_VDJ_1_*` are legacy (pre-0.13) column names (UNVERIFIED here; the docs say receptor fields moved to `obsm`).
- Typical `obs` after a workflow, as it appears in an exported MuData `.h5ad`: `airr:receptor_type`, `airr:receptor_subtype`, `airr:chain_pairing`, `airr:clone_id`, `airr:clone_id_size`, `airr:clonal_expansion`, plus the user's own sample or donor columns (`gex:sample`, `gex:patient` in the tutorial). In a plain AnnData the names have no `airr:` prefix. `clone_id` values are strings or categories; NaN for cells without a receptor. The exact dtype is UNVERIFIED.

### How scirpy keeps samples and donors apart
- VERIFIED: `define_clonotypes` and `define_clonotype_clusters` do not separate samples by default. The only default constraint is `within_group="receptor_type"` (T and B cells are not mixed).
- VERIFIED in source: `within_group` accepts `str | Sequence[str] | None`. A list such as `["receptor_type", "patient"]` forces clonotypes to stay within one donor. Any `obs` column can be used.
- VERIFIED: identical receptors in two donors get the same `clone_id` unless a donor column is in `within_group`. The tutorial does not use `within_group` for samples and treats cross-patient sharing as "public clonotypes".
- VERIFIED: `clonal_expansion(expanded_in=...)` and `clonotype_imbalance(replicate_col=...)` take sample columns for per-sample statistics.
- RECOMMENDATION: lumenCell should key clones by the pair (donor, clone_id) and never trust a bare `clone_id` across donors.

### 10x Cell Ranger VDJ naming
- Sources (search snippets only; the pages themselves returned 429 or 403): https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-annotations-vdj and https://software.10xgenomics.com/single-cell-vdj/software/pipelines/7.1/output/annotation
- Reported by search snippets (UNVERIFIED against the page): `raw_clonotype_id` is a column of `all_contig_annotations.csv` and `filtered_contig_annotations.csv`; it holds the clonotype the barcode was assigned. In multiplexed runs it carries a sample prefix, for example `sample1_clonotype23`. `clonotypes.csv` has `clonotype_id`, `frequency` (cells with that clonotype), `proportion`, and `cdr3s_aa`. The `cdr3s_aa` content is described differently across Cell Ranger versions.
- Differences from scirpy (reasoning, UNVERIFIED against 10x docs): Cell Ranger's ID is assigned per run, so it is only unique within that run or sample. It is not comparable across runs. Its exact grouping rule (full set of CDR3 nucleotide sequences) was not confirmed from the page. Scirpy's `clone_id` is recomputed from chain sequences by the user's chosen rule and arms, so the two IDs usually do not match one to one. scirpy's `read_10x_vdj` ingests the contig files; the sibling AnnData does not keep `raw_clonotype_id` as the clone key.

---

## 2. STARTRAC (Zhang et al. 2018)

- Paper: Zhang L, Yu X, Zheng L, et al. "Lineage tracking reveals dynamic relationships of T cells in colorectal cancer." Nature 564(7735):268-272 (2018). DOI 10.1038/s41586-018-0694-x. PMID 30479382 (Europe PMC record: not open access, not in PMC).
- URLs: https://doi.org/10.1038/s41586-018-0694-x ; https://pubmed.ncbi.nlm.nih.gov/30479382/
- R package: https://github.com/Japrin/STARTRAC (install: `devtools::install_github("Japrin/STARTRAC")`). DESCRIPTION: Version 0.1.0, authors Xianwen Ren and Liangtao Zheng, "License: file LICENSE"; the LICENSE file is GNU GPL v3. README says GPL-3.0. VERIFIED.
- Licence consequence: GPL-3. Do not copy or translate the code into lumenCell. Reimplement from the formulas below and cite the paper.
- The paper's Methods text was NOT read (paywalled, blocked). Everything below is VERIFIED from the package source (`R/startrac.R`, `R/functions.R`), raw files: https://raw.githubusercontent.com/Japrin/STARTRAC/master/R/startrac.R and .../R/functions.R. Agreement with the paper's wording is UNVERIFIED.

### Inputs
- A data frame `cell.data` with required columns `Cell_Name`, `clone.id`, `patient`, `majorCluster`, `loc` (tissue site). VERIFIED.
- Two tables are built: clone x location counts and clone x cluster counts. `clone.id` is assumed unique across patients (a clone-to-patient map is built from unique (patient, clone.id) pairs).

### Indices (VERIFIED from `Startrac.calIndex` and `Startrac.pIndex`)
Notation: n_ck = cells of clone c in cluster k. N_k = sum_c n_ck. p_c|k = n_ck / N_k. H(x) = -sum f log2 f with f = x / sum(x), zero terms dropped.
- Expansion, per cluster k: `expa_k = 1 - H(n_.k) / log2(K_k)`, where H is the entropy of the clone-size distribution in cluster k and K_k = number of clones with at least one cell in k. Gini is also reported. The code carries a TODO for the case K_k = 1 (log2 1 = 0, undefined). No minimum clone size filter appears in the index code.
- Migration, per clone c: `migr_c = H(n_c.^loc)` = entropy of the clone's distribution over locations. Raw base-2 entropy by default (`normEntropy=FALSE`); with TRUE it is divided by log2(number of locations).
- Transition, per clone c: `tran_c = H(n_c.^cluster)` = entropy of the clone's distribution over clusters. Same normalisation rule (divide by log2 K when `normEntropy=TRUE`).
- Cluster-level migration and transition: `tran_k = sum_c (n_ck / N_k) * tran_c`. The weight is the clone's share of cells in that cluster. Migration is the same with `migr_c`.
- Pairwise transition pTrans(k1, k2), for each pair of distinct clusters: take the two columns k1, k2 of the clone x cluster table. For each clone, compute h_c = H(n_c,k1, n_c,k2) (entropy over just those two clusters; NaN set to 0; maximum 1 bit). Then `pTrans(k1,k2) = sum_c [ (n_c,k1 + n_c,k2) / sum_c (n_c,k1 + n_c,k2) ] * h_c`. The value is symmetric in (k1, k2) by construction, so it carries no direction. Pairwise migration pMigr is the same per cluster over pairs of locations (within a cluster, subtable clones x two locations).
- Interpretation: high tran = clones spread over many states; pTrans = how mixed the clones in a state pair are, weighted by how many cells the clones have in that pair.

### Per-patient handling and aggregation (VERIFIED from `Startrac.run`)
- `Startrac.run(cell.data, proj, cores, n.perm, verbose)` builds one pooled object (all patients) and, when there is more than one patient, one object per patient restricted to patients with more than 30 cells (`patient.size > 30`). Per-patient and pooled results are stored separately. How the paper combines them (for example a mean or a test across patients) is UNVERIFIED.

### Null / permutation (VERIFIED)
- With `n.perm` set, each permutation shuffles the `clone.id` vector among cells within each patient (`clone.id[sample(nrow(.dat))]` per patient). Cluster and location labels stay fixed, so clone sizes and patient membership are preserved. The code keeps a commented-out global-shuffle alternative.
- Empirical p (`getSig`): `p = #(perm >= observed) / #perm`, one-sided "greater", with no +1 correction (so p can be 0).
- Clone size cutoff: none in the code. The paper's clone definition and any minimum size are UNVERIFIED.

---

## 3. LRT (Xie et al. 2023)

- Paper: Xie J, Jeon H, Xin G, Ma Q, Chung D (author list from a search result, UNVERIFIED). "LRT: Integrative analysis of scRNA-seq and scTCR-seq data to investigate clonal differentiation heterogeneity." PLOS Computational Biology, July 2023. DOI 10.1371/journal.pcbi.1011300. PMID 37428794; PMC10358952 (PMC page blocked by a bot check).
- URLs: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011300 (read via summary) ; https://doi.org/10.1371/journal.pcbi.1011300
- Repository: https://github.com/JuanXie19/LRT (R package with two Shiny apps; from search result, repo not opened; licence UNVERIFIED).
- Your guess about the authors and topic was right. The acronym's expansion was not verified.

### What it computes (read from the PLOS article, via summary)
- Trajectories: Slingshot on scRNA-seq only (cluster MST plus simultaneous principal curves). Pseudotime = arc length along the curve at the cell's projection.
- Clonotypes: defined by paired TCR alpha and beta sequences. Filters: CD8 analysis removed clonotypes with fewer than 5 cells; CD4 analysis kept clones with at least 2 cells from a single sample.
- Clonotype clustering: NOT by TCR sequence distance. A Dirichlet multinomial mixture on the clonotype x cell-state frequency matrix groups clonotypes by their distribution over cell states. EM fits it, K is chosen by Laplace approximation, and each clonotype goes to the cluster with maximum posterior.
- Trajectory biasedness, tested per clonotype cluster k by permutation. The null: "randomly shuffle the clonotype cluster membership for the cells while maintaining the size of each clonotype cluster." Statistics:
  1. Lineage bias: R_klm = log2(cells of cluster k in lineage l / cells in lineage m). Two-sided.
  2. Localisation: IQR of the cluster's pseudotime. One-sided (smaller IQR).
  3. Early versus late: median pseudotime against the global median. Two-sided.
  p = proportion of permuted statistics more extreme than observed. The number of permutations is not stated on the page read.
- Clone sharing: STARTRAC expansion, migration and transition indices are reported for clusters (details are in S1 Text, not read). There is no separate between-state clone-sharing score.
- Not stated on the page read: the rule for assigning a clone or cell to a lineage where Slingshot curves overlap (the trunk), and any Slingshot parameters. UNVERIFIED.
- Relevance: LRT shuffles a clustering label (not raw clone labels) and tests bias against the global pattern. Its pseudotime localisation and early/late statistics are the template for candidate C below.

---

## 4. CoSpar and LARRY coupling

### CoSpar
- Paper: Wang SW, Herriges MJ, Hurley K, Kotton DN, Klein AM. "CoSpar identifies early cell fate biases from single-cell transcriptomic and lineage information." Nature Biotechnology 40(7):1066-1074 (2022). DOI 10.1038/s41587-022-01209-1 (from a search result; paper text not read). Note: the title in the brief ("Dynamic inference ...") is the software's tagline, not the paper title.
- Code: https://github.com/ShouWenWang-Lab/cospar (MIT licence). Docs: https://cospar.readthedocs.io/en/stable/ . Source files read: `cospar/tool/_utils.py` (`get_normalized_covariance`), `cospar/tool/_map.py` (`fate_coupling`), `cospar/tool/_clone.py` (`clonal_fate_bias`). All below VERIFIED from code. Agreement with the paper's Methods is UNVERIFIED.
- Inputs: a cells x clones matrix `X_clone`, cell-state labels (`state_info`), time points. Fate coupling between fate clusters l and l' can be computed from the learned transition map or directly from clonal data (`source="X_clone"`).
- Clone-based coupling (`tl.fate_coupling(source="X_clone", method=...)`): first build `coarse_X_clone`, the clone x cluster count matrix, normalised within each cluster and then within each clone per time point if `normalize=True` (this was off in the original publication); `ignore_cell_number=True` binarises it. Then `get_normalized_covariance(coarse_X_clone.T, method)`:
  - `method="SW"` (default): Y_ll' = sum_i X_il X_il'; then Y_ll' / sqrt(Y_ll Y_l'l'). This is the cosine similarity of cluster columns over clones (map-based version uses fate probabilities P_i(C_l) over cells instead of clone counts).
  - `method="Weinreb"`: Y_ll' = Cov_i(X_il, X_il') divided by (mean_l * mean_l'), with np.cov and +0.0001 added to the means. It is the normalised covariance from Weinreb and Klein 2020.
  - `method="Jaccard"`: 1 - Jaccard distance between binarised columns.
- Clone fate bias (`tl.clonal_fate_bias`): Fisher exact test of a clone's enrichment in a fate cluster against the rest, accounting for clone size; Benjamini-Hochberg q-value; the bias is -log(q).

### LARRY / Weinreb
- Weinreb C, Rodriguez-Fraticelli A, Camargo FD, Klein AM. "Lineage tracing on transcriptional landscapes links state to fate during differentiation." Science 367(6479):eaaw3381 (2020). DOI 10.1126/science.aaw3381. Europe PMC: PMID 31974159, PMCID PMC7608074. The full text was not retrievable, so the paper's own "clonal coupling score" definition is UNVERIFIED.
- Weinreb C, Klein AM. "Lineage reconstruction from clonal correlations." PNAS 117(29):17041-17048 (2020). DOI 10.1073/pnas.2000238117. PMC7382289, read via Europe PMC XML. VERIFIED: it defines the "normalized covariance" of clone-size vectors over cell types, shows it can reconstruct a hierarchy (CLiNC), and shows it obeys symmetries that cross-tree transitions break. CoSpar's docstring attributes `method="Weinreb"` to this paper. The exact formula as typeset in the paper was not extracted; the formula in the bullet above comes from the CoSpar code.
- Caveat for lumenCell: these coupling scores assume clones are labelled at one time and fate observed later. They treat cells of one clone as a branching process. TCR clones sampled at one time point do not satisfy that assumption in the same way. Use the coupling formula as a descriptive statistic, not as the full inference model.

---

## 5. Standard practice: minimum clone size and nulls

- scirpy `clonal_expansion` default breakpoints (1, 2): categories 1, 2 and 3 or more, computed overall or per sample with `expanded_in`. "Expanded" is therefore usually 2 or more (or 3 or more in the strict bin). VERIFIED.
- LRT: at least 5 cells per clonotype (CD8 analysis); at least 2 cells in a single sample (CD4 analysis). VERIFIED (summary).
- STARTRAC code: no clone size filter in the index functions; per-patient objects only for patients with more than 30 cells. Paper-level clone filtering UNVERIFIED.
- CoSpar: no explicit minimum in `fate_coupling`; clones present in no selected fate are dropped. Singletons are noted as a bias source in a docstring comment ("biases towards large clusters ... most of them just a few cells").
- Nulls seen in prior work: STARTRAC shuffles the `clone.id` column among cells within a patient (clone sizes kept). LRT shuffles clonotype-cluster membership across cells keeping cluster sizes. Both are permutations of clone labels with cluster/state labels fixed.
- Statistical point (reasoning, not a cited fact): a singleton clone has zero entropy by construction, so including singletons inflates any purity score. Exclude them (minimum size at least 2) and report sensitivity at 2, 3 and 5.

---

## 6. Candidate metric definitions (RECOMMENDATION; derivations from prior work noted)

### Common set-up
- Cells i with donor d(i), clone key c(i) = (donor, clone_id) (NaN clone_id cells excluded), pseudotime t_i, lineage weights w_il for L >= 2 lineages with sum_l w_il = 1.
  - Palantir: w_il = terminal-state fate probabilities.
  - PAGA+DPT, scFates: one-hot on the lineage (leaf branch or milestone path) the cell belongs to.
- Commitment u_i in [0,1] handles the shared trunk. A cell on the shared trunk carries no information about which lineage a clone chose.
  - Hard-label methods: u_i = 0 for cells on non-leaf (shared) segments; u_i = 1 on leaf segments. In scFates with nested branches, cells on an internal segment can instead take a weight spread over descendant leaves; then u_i = 0 for the headline score.
  - Palantir: u_i = 1 - H(w_i) / log L, or a threshold max_l w_il >= tau. The threshold is a hyperparameter that needs a default and sensitivity check.
- Pool P = {i : u_i > 0}. Clone eligibility: at least m_min cells of the clone in P (default m_min = 2; report at 3 and 5).
- Soft counts: n_cl = sum_{i in c, i in P} u_i w_il. N_c = sum_l n_cl. p_cl = n_cl / N_c.
- Null (shared by all candidates): B permutations (default 1000) that shuffle clone keys among the cells of P within each donor (STARTRAC code, LRT). This preserves each donor's clone-size multiset and every cell's w_il, u_i, t_i. Report the observed score, null mean and sd, effect = obs - null mean, and p = (1 + #{null >= obs}) / (1 + B) (the +1 avoids p = 0, which the STARTRAC code allows).
  - Soft weights also make single-cell uncertainty look like clonal spread. The null contains the same effect, so use the null-subtracted effect as the headline, not the raw score.
  - Stricter variant: shuffle within donor x pseudotime bin, to separate clone-lineage agreement from the simple fact that clones are concentrated at similar maturation stages.
- Interpretation limit: shared clones show that two states are related, not which came first. Siblings that split late legitimately share clones, so a low score is not proof of a wrong lineage split. Report per-pair values next to the scalar.

### Candidate A: clone purity across lineages (entropy)
- Derives from: STARTRAC clone-level transition entropy (tran_c) and its cell-share weighting; the LRT lineage-bias idea.
- Inputs: n_cl, m_min, L.
- Formula: H_c = -sum_l p_cl log p_cl / log L (normalised to [0,1]). Run score A = 1 - sum_c N_c H_c / sum_c N_c (mass-weighted mean purity, in [0,1]; weight = the clone's committed cell mass, as STARTRAC weights clones by share of cells).
- Soft weights: used directly through n_cl. Option: Monte-Carlo hardening (sample each cell's lineage from w_i) as a robustness check.
- Trunk cells: removed from P through u_i. They never lower purity.
- Null: as above. Primary output is A_obs - mean(A_null), plus z and p.
- Why primary: one number, bounded, one-sided, directly answers "do clones stay within lineages".

### Candidate B: lineage-pair coupling
- Derives from: CoSpar `fate_coupling` (SW cosine, Weinreb normalised covariance) and STARTRAC pTrans.
- Inputs: matrix X (eligible clones x L) of soft counts n_cl.
- Formulas (pick one, record in the output):
  - Cosine (CoSpar "SW"): Y_ll' = sum_c X_cl X_cl' / sqrt(sum_c X_cl^2 * sum_c X_cl'^2).
  - Normalised covariance (Weinreb): Y_ll' = Cov_c(X_cl, X_cl') / (mean_c X_cl * mean_c X_cl'). Needs many clones; unstable with few.
  - pTrans on soft counts (STARTRAC): pTrans_ll' = sum_c [ (X_cl + X_cl') / S ] * H2(X_cl, X_cl'), S = sum_c (X_cl + X_cl'), H2 = entropy of the two-state share in bits (symmetric; 0 to 1).
- Per pair, compute the null mean and sd of Y_ll' (or pTrans) from the same permutations. Output: the L x L matrix of observed values, null-standardised z, and p.
- Scalar summary: mean over pairs of (null mean - observed) / null sd, i.e. "decoupling". High = clones stay within lineages more than chance. For L = 2 there is a single pair.
- Soft weights and trunk: as for A (n_cl restricted to P). Because it is pairwise, sibling lineages are expected to couple more than distant ones; the matrix should be read against the tree.
- Tree-consistency check (optional): order lineage pairs by tree distance; expect coupling to fall with tree distance. This is a design idea, not from a cited method.

### Candidate C: pseudotime coherence of clones (secondary, diagnostic)
- Derives from: LRT localisation (IQR of pseudotime) and early/late (median) statistics, with LRT's shuffle null.
- Inputs: t_i rank-normalised to [0,1] per run (makes PAGA/DPT, Palantir and scFates scales comparable), clone keys, optional w_il.
- Formula: R2 = 1 - [ sum_c sum_{i in c} (t_i - mean_c(t))^2 ] / [ sum_i (t_i - mean_d(i)(t))^2 ], over eligible clones; the denominator centres on the donor mean. Optional per-lineage version using weights w_il in the sums.
- Trunk: all cells with pseudotime are usable, since the trunk has valid pseudotime.
- Null: as above. Directionless: it measures that clones occupy narrow pseudotime windows, not who is earlier.
- Warning: low coherence is expected for clones that legitimately span naive-to-effector states. Treat C as a diagnostic, and do not combine it with A or B into one headline.

### Minimum implementation notes
- Run the metric per donor, then combine donors: pool the donor-wise sums for A, or average the donor-wise effects weighted by committed cells. Donors with fewer than a configurable number of eligible clones (STARTRAC uses 30 cells per patient) are skipped and reported as skipped, not silently dropped. How to combine donors is a design choice, not taken from a source.
- Reimplement the formulas from scratch. Do not copy STARTRAC (GPL-3). CoSpar is MIT and scirpy BSD-3, so their code could be adapted with attribution, but there is no need to.
- Cite: Zhang et al. 2018 (STARTRAC); Xie et al. 2023 (LRT); Wang et al. 2022 (CoSpar); Weinreb et al. 2020 Science (LARRY) and Weinreb and Klein 2020 PNAS (normalised covariance); Sturm et al. 2020 (scirpy).

---

## Source and DOI index

- scirpy docs: https://scirpy.scverse.org/ ; repo https://github.com/scverse/scirpy ; PyPI https://pypi.org/project/scirpy/
- scirpy paper: https://doi.org/10.1093/bioinformatics/btaa611
- STARTRAC paper: https://doi.org/10.1038/s41586-018-0694-x ; PubMed https://pubmed.ncbi.nlm.nih.gov/30479382/
- STARTRAC code: https://github.com/Japrin/STARTRAC
- LRT paper: https://doi.org/10.1371/journal.pcbi.1011300 ; PLOS URL above ; PMC https://pmc.ncbi.nlm.nih.gov/articles/PMC10358952
- LRT code: https://github.com/JuanXie19/LRT
- CoSpar paper: https://doi.org/10.1038/s41587-022-01209-1 ; code https://github.com/ShouWenWang-Lab/cospar ; docs https://cospar.readthedocs.io/
- LARRY: https://doi.org/10.1126/science.aaw3381
- Normalised covariance: https://doi.org/10.1073/pnas.2000238117 ; PMC https://pmc.ncbi.nlm.nih.gov/articles/PMC7382289
- 10x Cell Ranger VDJ outputs: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-annotations-vdj

## Open items (UNVERIFIED, need a human or a paywall-capable reader)

1. STARTRAC paper Methods text: exact clone definition, any minimum clone size, how per-patient values are combined.
2. LARRY (Science 2020) clonal coupling score as defined in that paper.
3. CoSpar paper Methods text on fate coupling and `normalize` (the docstring says normalisation was off in the original publication).
4. 10x Cell Ranger clonotype definition and `raw_clonotype_id` semantics from the 10x page itself (blocked by 429/403).
5. LRT lineage-assignment rule, permutation count, repository licence.
6. Whether `has_ir` and the exact dtype of `clone_id` exist in current scirpy output.
