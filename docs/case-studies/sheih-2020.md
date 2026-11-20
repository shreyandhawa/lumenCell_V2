# Case study: CD8 CAR-T cells in Sheih et al. 2020

The numbers come from `benchmarks/case_studies/results/sheih2020-2026-11-20.json`, generated on
2026-11-20 at commit 5b4f83f by `benchmarks/case_studies/sheih2020.py`. The tables under Results are
generated from that file by `lumencell.bench.case_study.render_markdown` and are never edited by
hand; `tests/unit/test_case_study_writeup.py` fails if they drift from it. A number in the prose is
written exactly as it appears in those tables.

**The short answer: no patient's cells are shown to fit a tree.** Of the 12 per-patient runs, 8 are
`weak` (the clonal check failed) and 4 are `provisional`, which means "not shown to be a tree", never
support. No run is `supported`, because the verdict rule failed its verification on 2026-11-17.

## Question

Sheih et al. sampled CD8 CAR-T cells from the infusion product and from blood at three later time
points in four patients, with gene expression and T-cell receptor chains for each cell. For each
patient, and for each of lumenCell's three methods (`paga_dpt`, `palantir`, `scfates_ppt`), do the
cells fit a branching trajectory, and do the T-cell receptor clones, which no method sees, agree with
the lineages the method returns? Clone sharing can only say that two states are related. It never
says which came first, and nothing here claims an order from it.

## Data

- **Paper:** Sheih A., Voillet V., Hanafi L.-A. et al. Clonal kinetics and single-cell
  transcriptional profiling of CAR-T cells in patients undergoing CD19 CAR-T immunotherapy. Nat
  Commun 11:219 (2020). doi:10.1038/s41467-019-13880-1, PMC6954177. Read as open full text from PMC
  (NCBI E-utilities) on 2026-11-20; the DOI, title, journal, volume and article number were checked
  against Crossref the same day.
- **Accession:** GEO GSE125881 (the paper's Data availability statement names it for the scRNA-seq
  data). The case study uses two files from it: the expression matrix
  `GSE125881_raw.expMatrix.csv.gz` and `GSE125881_RAW.tar`, which holds one clonotype table per
  sample. Their sizes and SHA-256 are recorded in `src/lumencell/bench/builders/sheih2020.py`.
- **Licence:** the article is CC BY 4.0 (its licence statement, read in the PMC full text). The
  data is a public GEO accession; GEO's own terms were not read, because the GEO HTML pages sit
  behind a CAPTCHA, which is not bypassed (unverified).
- **Patients:** CLL-1 and CLL-2 (chronic lymphocytic leukaemia) and NHL-6 and NHL-7 (non-Hodgkin
  lymphoma). These are the four patients of the paper's scRNA-seq cohort (Results, "Single-cell
  transcriptome analysis of CAR-T cells over time"). Their CD8 CAR-T cells were made from bulk CD8
  T cells (CLL) or central memory CD8 T cells (NHL), and cells were sorted on the EGFRt marker that
  the CAR construct carries (same section, and Methods).
- **Time points:** the infusion product (IP) and three blood samples per patient, with GEO's labels
  as written (Task 9's note, [research/2026-11-19-sheih-verification.md](../superpowers/research/2026-11-19-sheih-verification.md)):

  | Patient | Time points |
  |---|---|
  | CLL-1 | IP, d21, d38, d112 |
  | CLL-2 | IP, d12, d29, d83 |
  | NHL-6 | IP, d12, d29, d102 |
  | NHL-7 | IP, d12, d28, d89 |

  The paper's windows are early (day 7 to 14), late (day 26 to 30) and very late (day 83 to 112)
  (Results, "Single-cell transcriptome analysis of CAR-T cells over time"). CLL-1's d21 and d38
  samples fall outside the first two windows although GEO labels them Early and Contraction; why
  is not stated in the paper (unverified). lumenCell uses the labels only as strata and to pick the
  root, never as a time axis.
- **Cells:** 62,167, the number the paper reports after its quality filters (Methods,
  "Transcriptional profiling of single CD8+ CAR-T cells"), and the number in the derived file. Per
  patient (the results file's dataset block): CLL-1 18,893, CLL-2 12,990, NHL-6 21,762, NHL-7
  8,522. 42,362 cells carry a paired clone label under the declared rule (CLL-1 12,745, CLL-2
  9,310, NHL-6 13,575, NHL-7 6,732; Task 11 in HANDOFF).

## Declaration

Everything below was declared in
[`[case_study.sheih2020]` of `benchmarks/thresholds.toml`](../../benchmarks/thresholds.toml), with
a dated changelog entry (2026-11-19), and merged in PR #51 as commit
[5ef3151](https://github.com/shreyandhawa/lumenCell_V2/commit/5ef31511db66c143abe8217be25887d903be8ee3)
before the case-study data was built and before any score was computed. It was written from the
paper and the file headers only.

- **Clone definition:** `clone_id`, from the paired chains: a cell with 1 or 2 distinct TRA and 1 or
  2 distinct TRB CDR3 nucleotide sequences gets the sorted set as its label (`chains = "paired"`,
  `cdr3 = "nt"`, `max_per_chain = 2`). The GEO tables are labelled "unpaired" (one row per clonotype
  and chain), but every row names its chain and lists its cells' barcodes, so the chains group by
  cell. No filter is applied (`filters = []`): the tables carry no productive, high-confidence or
  full-length flag, and whether the authors kept only productive chains is not stated (unverified).
- **Namespace:** the patient (`donor_key = "patient"`), so the same receptor in two patients is two
  clones.
- **Strata:** `strata_keys = ["timepoint"]`. The third null (N3) shuffles labels within patient and
  time point, so clones that expand at one time point are not read as lineage agreement.
- **Root:** `cluster:timepoint=IP`, the infusion-product cells. The paper's main text and Methods
  name no trajectory root (unverified: supplementary figures not read). The IP is assumed as the
  start because it is the earliest sample and the population the paper compares blood against
  (Fig. 5, Fig. 7). This is an assumption of the declaration. Every
  pseudotime below runs away from this root because the root was declared, not because of anything
  in the data.
- **Unit:** each method runs once per patient, with stability on each (4 x 3 = 12 runs), because
  lumenCell has no batch correction. A pooled run over all four patients is declared as a labelled
  diagnostic only, never the headline.
- **Clonal parameters:** the `[clonal]` defaults: clone size m = 2 primary and m = 3 as the
  sensitivity check, B = 1000 shuffles, seed 0. No threshold was added or changed.

## What the paper reports

- **Bulk TCRB sequencing, another cohort.** In ten other patients, CAR-T clonal diversity is highest
  in the IP and falls after infusion (Fig. 1c), the IP and blood share clonotypes, more so early
  than late (Fig. 1b), and individual clones follow different count patterns over the first 30 days
  (Fig. 3c). These patients are not in GSE125881.
- **Expression over time.** Blood profiles move further from the IP's at each later time point
  (Fig. 5a). IP cells express more glycolysis and AP-1 genes and CD69, which the authors put down to
  stimulation during manufacturing; cytotoxicity genes are higher early and late and fall very late;
  MKI67 falls at the late and very late points (Fig. 5b, 5d). No exhaustion signature rises
  consistently (Fig. 5c, 5e).
- **Four IP clusters.** Clustering finds four clusters in every patient's IP (Fig. 6b, 6c). Cluster 2
  is marked by cytotoxicity genes, and cluster 4 consists of cells in S or G2/M (Supplementary
  Fig. 4a). The authors read cluster 3's mitochondrial genes as a sign of pre-apoptotic cells
  (Discussion).
- **No trajectory.** Monocle 3 pseudotime and RNA velocity on IP cells did not yield clear
  trajectories to "indicate a developmental relationship between clusters" (Results, "Contribution
  of infused CAR-T clusters to in vivo repertoire"; Supplementary Fig. 4b, 4c). That analysis
  covers IP cells only (Methods, "Single-cell trajectory analysis"). The main text and Methods run
  no trajectory across time points and name no root (unverified: supplementary figures not read).
- **Clones and IP clusters (NHL-6 and NHL-7).** Clonotypes are called IRF or DRF when their relative
  frequency is higher or lower at the early time point than in the IP (Fisher's exact test, FDR 5%;
  Methods). Most IRF clones sit in IP clusters 2 and 4 (85.6% in NHL-6 and 84.6% in NHL-7), which
  hold 42.7% and 27.1% of IP cells (Fig. 7c, Fig. 6b). All IRF clonotypes are seen again at the late
  time point, against 22% and 5% of DRF clonotypes (Table 1). The share of persisting clonotypes
  found in IP cluster 2 grows at later time points (Fig. 7d). These are comparisons of clone
  frequencies between sampled time points; the order is the trial's sampling order.
- **Processing.** The authors regressed out UMI count, mitochondrial share and patient before
  clustering (Seurat `ScaleData`; Methods), and called paired TCR clonotypes with Cell Ranger vdj
  (Methods, "TCRB sequencing of single CD8+ CAR-T cells"). lumenCell does no such regression.

## Results

The per-patient tables come first and the pooled diagnostic after them. Columns: L is the number of
lineages the method returned. Cells, Clones and Effective clones are what the clonal score rests on:
cells in exactly one lineage that belong to a clone of at least m = 2 cells in the same patient,
those clones, and their inverse Simpson count. η is the normalised clonal excess against each null
(N1 within patient, N2 within patient and pseudotime decile, N3 within patient and time point); the
headline p is the largest of the three. η_base is the same score for a k-means partition with k = L
on the same PCA, which uses no lineage; Δη = η − η_base, with its paired permutation p. Purity
follows STARTRAC (Zhang et al., Nature 2018, doi:10.1038/s41586-018-0694-x) and LRT (Xie et al., PLoS
Comput Biol 2023, doi:10.1371/journal.pcbi.1011300). The lineage-pair coupling follows CoSpar (Wang
et al., Nat Biotechnol 2022, doi:10.1038/s41587-022-01209-1), on raw scored-cell counts rather than
CoSpar's normalisation. A dash is a value that could not be computed (for example a lineage with no
scored cell).

<!-- case-study:start -->
### Per patient (primary)

| Patient | Method | L | Verdict | Clonal status | Cells | Clones | Effective clones | Strata | η (N1) | η (N2) | η (N3) | Headline p | η_base | Δη | p_Δη | η (m = 3) | p (m = 3) | Purity | Purity without largest |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CLL-1 | paga_dpt | 5 | provisional | ok | 2749 | 652 | 127.706 | 4 | 0.695 | 0.363 | 0.247 | 0.000999 | 0.659 | 0.036 | 0.000999 | 0.626 | 0.000999 | 0.852 | 0.870 |
| CLL-1 | palantir | 5 | weak | ok | 5521 | 786 | 197.902 | 3 | 0.125 | 0.028 | 0.020 | 0.0759 | 0.156 | -0.031 | 0.996 | 0.129 | 0.0769 | 0.911 | 0.914 |
| CLL-1 | scfates_ppt | 5 | weak | ok | 3233 | 631 | 198.641 | 4 | 0.571 | 0.052 | 0.073 | 0.000999 | 0.577 | -0.005 | 0.766 | 0.468 | 0.0509 | 0.877 | 0.880 |
| CLL-2 | paga_dpt | 5 | provisional | ok | 3569 | 1271 | 649.190 | 4 | 0.213 | 0.128 | 0.024 | 0.000999 | 0.188 | 0.024 | 0.002 | 0.135 | 0.154 | 0.713 | 0.719 |
| CLL-2 | palantir | 5 | weak | ok | 3360 | 1177 | 595.694 | 3 | 0.114 | 0.053 | 0.021 | 0.0699 | 0.125 | -0.011 | 0.792 | 0.088 | 0.147 | 0.730 | 0.734 |
| CLL-2 | scfates_ppt | 3 | weak | ok | 2664 | 968 | 510.421 | 3 | 0.143 | 0.124 | 0.011 | 0.155 | 0.125 | 0.018 | 0.153 | 0.112 | 0.364 | 0.608 | 0.621 |
| NHL-6 | paga_dpt | 3 | weak | ok | 5939 | 548 | 58.793 | 4 | 0.415 | 0.103 | 0.089 | 0.000999 | 0.422 | -0.008 | 0.997 | 0.402 | 0.000999 | 0.485 | 0.504 |
| NHL-6 | palantir | 3 | weak | ok | 9539 | 458 | 35.065 | 4 | 0.138 | 0.038 | 0.031 | 0.000999 | 0.159 | -0.022 | 1 | 0.135 | 0.000999 | 0.673 | 0.681 |
| NHL-6 | scfates_ppt | 5 | weak | ok | 3684 | 320 | 44.068 | 4 | 0.313 | 0.006 | 0.011 | 0.026 | 0.363 | -0.050 | 1 | 0.293 | 0.025 | 0.595 | 0.608 |
| NHL-7 | paga_dpt | 3 | weak | ok | 3511 | 299 | 36.363 | 4 | 0.216 | 0.134 | 0.009 | 0.015 | 0.184 | 0.032 | 0.000999 | 0.202 | 0.014 | 0.479 | 0.488 |
| NHL-7 | palantir | 3 | provisional | ok | 3432 | 276 | 34.329 | 4 | 0.064 | 0.049 | 0.028 | 0.005 | 0.055 | 0.010 | 0.019 | 0.059 | 0.00599 | 0.667 | 0.665 |
| NHL-7 | scfates_ppt | 3 | provisional | ok | 1584 | 185 | 35.051 | 3 | 0.070 | 0.027 | 0.034 | 0.000999 | 0.021 | 0.049 | 0.000999 | 0.068 | 0.000999 | 0.351 | 0.358 |

Verdict: the rule failed its verification on 2026-11-17, so supported is never shown; a passed check is shown as provisional. Tree-likeness and both stability checks are left out (no calibrated threshold), so the verdict reads coverage and clonal support only. A failed check is shown as weak or not_supported.

#### Lineage-pair coupling (Benjamini–Hochberg)

| Patient | Method | a | b | Observed | z | p_adj |
|---|---|---|---|---|---|---|
| CLL-1 | paga_dpt | 0 | 1 | 0.455 | -9.160 | 1 |
| CLL-1 | paga_dpt | 0 | 2 | 0.080 | -21.234 | 1 |
| CLL-1 | paga_dpt | 0 | 3 | 0.018 | -11.841 | 1 |
| CLL-1 | paga_dpt | 0 | 4 | – | – | – |
| CLL-1 | paga_dpt | 1 | 2 | 0.118 | -65.719 | 1 |
| CLL-1 | paga_dpt | 1 | 3 | 0.046 | -15.847 | 1 |
| CLL-1 | paga_dpt | 1 | 4 | – | – | – |
| CLL-1 | paga_dpt | 2 | 3 | 0.338 | -9.496 | 1 |
| CLL-1 | paga_dpt | 2 | 4 | – | – | – |
| CLL-1 | paga_dpt | 3 | 4 | – | – | – |
| CLL-1 | palantir | 0 | 1 | – | – | – |
| CLL-1 | palantir | 0 | 2 | – | – | – |
| CLL-1 | palantir | 0 | 3 | – | – | – |
| CLL-1 | palantir | 0 | 4 | – | – | – |
| CLL-1 | palantir | 1 | 2 | 0.192 | -1.652 | 1 |
| CLL-1 | palantir | 1 | 3 | 0.257 | 0.242 | 1 |
| CLL-1 | palantir | 1 | 4 | – | – | – |
| CLL-1 | palantir | 2 | 3 | 0.523 | -7.998 | 1 |
| CLL-1 | palantir | 2 | 4 | – | – | – |
| CLL-1 | palantir | 3 | 4 | – | – | – |
| CLL-1 | scfates_ppt | 0 | 1 | 0.043 | -22.117 | 1 |
| CLL-1 | scfates_ppt | 0 | 2 | 0.056 | -12.005 | 1 |
| CLL-1 | scfates_ppt | 0 | 3 | 0.207 | -0.443 | 1 |
| CLL-1 | scfates_ppt | 0 | 4 | 0.469 | -5.590 | 1 |
| CLL-1 | scfates_ppt | 1 | 2 | 0.617 | 3.086 | 0.02 |
| CLL-1 | scfates_ppt | 1 | 3 | 0.000 | -2.624 | 1 |
| CLL-1 | scfates_ppt | 1 | 4 | 0.037 | -10.637 | 1 |
| CLL-1 | scfates_ppt | 2 | 3 | 0.000 | -2.181 | 1 |
| CLL-1 | scfates_ppt | 2 | 4 | 0.030 | -7.431 | 1 |
| CLL-1 | scfates_ppt | 3 | 4 | 0.109 | -0.915 | 1 |
| CLL-2 | paga_dpt | 0 | 1 | 0.134 | -0.007 | 1 |
| CLL-2 | paga_dpt | 0 | 2 | 0.095 | -1.725 | 1 |
| CLL-2 | paga_dpt | 0 | 3 | 0.044 | -5.630 | 1 |
| CLL-2 | paga_dpt | 0 | 4 | 0.045 | -5.715 | 1 |
| CLL-2 | paga_dpt | 1 | 2 | 0.045 | -1.557 | 1 |
| CLL-2 | paga_dpt | 1 | 3 | 0.052 | -2.894 | 1 |
| CLL-2 | paga_dpt | 1 | 4 | 0.064 | -2.738 | 1 |
| CLL-2 | paga_dpt | 2 | 3 | 0.360 | 0.771 | 1 |
| CLL-2 | paga_dpt | 2 | 4 | 0.224 | -1.657 | 1 |
| CLL-2 | paga_dpt | 3 | 4 | 0.550 | -7.572 | 1 |
| CLL-2 | palantir | 0 | 1 | – | – | – |
| CLL-2 | palantir | 0 | 2 | 0.208 | -1.814 | 1 |
| CLL-2 | palantir | 0 | 3 | – | – | – |
| CLL-2 | palantir | 0 | 4 | 0.353 | 0.227 | 1 |
| CLL-2 | palantir | 1 | 2 | – | – | – |
| CLL-2 | palantir | 1 | 3 | – | – | – |
| CLL-2 | palantir | 1 | 4 | – | – | – |
| CLL-2 | palantir | 2 | 3 | – | – | – |
| CLL-2 | palantir | 2 | 4 | 0.543 | -9.554 | 1 |
| CLL-2 | palantir | 3 | 4 | – | – | – |
| CLL-2 | scfates_ppt | 0 | 1 | 0.521 | -6.660 | 1 |
| CLL-2 | scfates_ppt | 0 | 2 | 0.298 | -2.080 | 1 |
| CLL-2 | scfates_ppt | 1 | 2 | 0.424 | 0.928 | 0.527 |
| NHL-6 | paga_dpt | 0 | 1 | 0.342 | -128.512 | 1 |
| NHL-6 | paga_dpt | 0 | 2 | 0.248 | -117.054 | 1 |
| NHL-6 | paga_dpt | 1 | 2 | 0.613 | -66.632 | 1 |
| NHL-6 | palantir | 0 | 1 | – | – | – |
| NHL-6 | palantir | 0 | 2 | 0.766 | -80.239 | 1 |
| NHL-6 | palantir | 1 | 2 | – | – | – |
| NHL-6 | scfates_ppt | 0 | 1 | 0.132 | -0.793 | 1 |
| NHL-6 | scfates_ppt | 0 | 2 | 0.248 | -58.538 | 1 |
| NHL-6 | scfates_ppt | 0 | 3 | 0.622 | -25.912 | 1 |
| NHL-6 | scfates_ppt | 0 | 4 | 0.612 | -38.313 | 1 |
| NHL-6 | scfates_ppt | 1 | 2 | 0.000 | -1.623 | 1 |
| NHL-6 | scfates_ppt | 1 | 3 | 0.045 | -1.294 | 1 |
| NHL-6 | scfates_ppt | 1 | 4 | 0.037 | -1.409 | 1 |
| NHL-6 | scfates_ppt | 2 | 3 | 0.142 | -40.147 | 1 |
| NHL-6 | scfates_ppt | 2 | 4 | 0.138 | -46.176 | 1 |
| NHL-6 | scfates_ppt | 3 | 4 | 0.954 | 2.230 | 0.02 |
| NHL-7 | paga_dpt | 0 | 1 | 0.111 | -24.311 | 1 |
| NHL-7 | paga_dpt | 0 | 2 | 0.115 | -24.474 | 1 |
| NHL-7 | paga_dpt | 1 | 2 | 0.777 | -46.745 | 1 |
| NHL-7 | palantir | 0 | 1 | 0.826 | -17.956 | 1 |
| NHL-7 | palantir | 0 | 2 | – | – | – |
| NHL-7 | palantir | 1 | 2 | – | – | – |
| NHL-7 | scfates_ppt | 0 | 1 | 0.929 | -2.173 | 1 |
| NHL-7 | scfates_ppt | 0 | 2 | 0.665 | -7.523 | 1 |
| NHL-7 | scfates_ppt | 1 | 2 | 0.645 | -7.387 | 1 |

#### Stability

| Patient | Method | Kind | R | Failed | Share matched median (C3) | Share matched min | Spearman median (C4) | Spearman min | Mean Jaccard median | Mean Jaccard min | Lineages min |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CLL-1 | paga_dpt | reseeds | 2 | 0 | 0.600 | 0.400 | 1.000 | 1.000 | 0.733 | 0.643 | 4 |
| CLL-1 | paga_dpt | subsamples | 3 | 0 | 0.400 | 0.400 | 0.990 | 0.987 | 0.721 | 0.620 | 2 |
| CLL-1 | palantir | reseeds | 2 | 0 | 0.800 | 0.600 | 0.999 | 0.998 | 0.881 | 0.802 | 4 |
| CLL-1 | palantir | subsamples | 3 | 0 | 0.400 | 0.200 | 0.916 | 0.915 | 0.695 | 0.673 | 2 |
| CLL-1 | scfates_ppt | reseeds | 2 | 0 | 0.700 | 0.600 | 0.986 | 0.983 | 0.789 | 0.715 | 4 |
| CLL-1 | scfates_ppt | subsamples | 3 | 0 | 0.800 | 0.800 | 0.973 | 0.970 | 0.806 | 0.660 | 4 |
| CLL-2 | paga_dpt | reseeds | 2 | 0 | 1.000 | 1.000 | 1.000 | 1.000 | 0.938 | 0.933 | 5 |
| CLL-2 | paga_dpt | subsamples | 3 | 0 | 0.600 | 0.200 | 0.947 | 0.943 | 0.772 | 0.733 | 1 |
| CLL-2 | palantir | reseeds | 2 | 0 | 0.900 | 0.800 | 0.999 | 0.999 | 0.894 | 0.832 | 5 |
| CLL-2 | palantir | subsamples | 3 | 1 | 0.400 | 0.000 | 0.964 | 0.000 | 0.583 | 0.000 | 0 |
| CLL-2 | scfates_ppt | reseeds | 2 | 0 | 0.500 | 0.333 | 0.989 | 0.988 | 0.633 | 0.629 | 6 |
| CLL-2 | scfates_ppt | subsamples | 3 | 0 | 1.000 | 0.000 | 0.908 | 0.860 | 0.759 | 0.463 | 5 |
| NHL-6 | paga_dpt | reseeds | 2 | 0 | 0.667 | 0.667 | 1.000 | 1.000 | 0.916 | 0.913 | 2 |
| NHL-6 | paga_dpt | subsamples | 3 | 0 | 1.000 | 0.667 | 0.985 | 0.983 | 0.874 | 0.831 | 2 |
| NHL-6 | palantir | reseeds | 2 | 0 | 0.333 | 0.333 | 0.999 | 0.999 | 0.779 | 0.772 | 2 |
| NHL-6 | palantir | subsamples | 3 | 0 | 1.000 | 0.333 | 0.937 | 0.907 | 0.848 | 0.735 | 2 |
| NHL-6 | scfates_ppt | reseeds | 2 | 0 | 1.000 | 1.000 | 0.996 | 0.996 | 0.927 | 0.896 | 5 |
| NHL-6 | scfates_ppt | subsamples | 3 | 0 | 1.000 | 0.800 | 0.986 | 0.984 | 0.883 | 0.864 | 4 |
| NHL-7 | paga_dpt | reseeds | 2 | 0 | 0.667 | 0.667 | 1.000 | 1.000 | 0.865 | 0.763 | 2 |
| NHL-7 | paga_dpt | subsamples | 3 | 0 | 0.667 | 0.333 | 0.980 | 0.966 | 0.728 | 0.548 | 2 |
| NHL-7 | palantir | reseeds | 2 | 0 | 0.667 | 0.667 | 0.999 | 0.999 | 0.721 | 0.631 | 2 |
| NHL-7 | palantir | subsamples | 3 | 0 | 0.333 | 0.333 | 0.898 | 0.878 | 0.647 | 0.640 | 2 |
| NHL-7 | scfates_ppt | reseeds | 2 | 0 | 1.000 | 1.000 | 0.994 | 0.993 | 0.832 | 0.803 | 6 |
| NHL-7 | scfates_ppt | subsamples | 3 | 0 | 1.000 | 1.000 | 0.991 | 0.991 | 0.856 | 0.834 | 5 |

#### Cell cycle

| Patient | Method | Cell-cycle status | Dataset share cycling | Warnings |
|---|---|---|---|---|
| CLL-1 | paga_dpt | scored | 0.210 | Lineage Cluster 5 is mostly cycling cells; it may be a cell-cycle state rather than a fate. |
| CLL-1 | palantir | scored | 0.210 | – |
| CLL-1 | scfates_ppt | scored | 0.210 | Lineage Cluster 8 is mostly cycling cells; it may be a cell-cycle state rather than a fate. |
| CLL-2 | paga_dpt | scored | 0.183 | Lineage Cluster 8 is mostly cycling cells; it may be a cell-cycle state rather than a fate. |
| CLL-2 | palantir | scored | 0.183 | – |
| CLL-2 | scfates_ppt | scored | 0.183 | – |
| NHL-6 | paga_dpt | scored | 0.192 | – |
| NHL-6 | palantir | scored | 0.192 | – |
| NHL-6 | scfates_ppt | scored | 0.192 | – |
| NHL-7 | paga_dpt | scored | 0.266 | – |
| NHL-7 | palantir | scored | 0.266 | Lineage Cluster 4 is mostly cycling cells; it may be a cell-cycle state rather than a fate. |
| NHL-7 | scfates_ppt | scored | 0.266 | Lineage Cluster 4 is mostly cycling cells; it may be a cell-cycle state rather than a fate. |

### Pooled diagnostic (not the headline: patient structure is not corrected)

| Patient | Method | L | Verdict | Clonal status | Cells | Clones | Effective clones | Donors scored | Strata | η (N1) | η (N2) | η (N3) | Headline p | η_base | Δη | p_Δη | η (m = 3) | p (m = 3) | Purity | Purity without largest |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pooled | paga_dpt | 5 | provisional | ok | 18000 | 3246 | 372.766 | 4 | 16 | 0.518 | 0.302 | 0.013 | 0.002 | 0.484 | 0.034 | 0.000999 | 0.461 | 0.000999 | 0.840 | 0.851 |
| pooled | palantir | 5 | weak | ok | 27759 | 3805 | 226.064 | 4 | 16 | 0.438 | 0.034 | 0.004 | 0.0799 | 0.420 | 0.019 | 0.000999 | 0.403 | 0.031 | 0.769 | 0.772 |
| pooled | scfates_ppt | 5 | weak | ok | 12243 | 1526 | 123.945 | 4 | 13 | 0.074 | 0.000 | 0.002 | 0.475 | 0.156 | -0.082 | 1 | 0.076 | 0.506 | 0.542 | 0.551 |

#### Lineage-pair coupling (Benjamini–Hochberg)

| Patient | Method | a | b | Observed | z | p_adj |
|---|---|---|---|---|---|---|
| pooled | paga_dpt | 0 | 1 | 0.032 | -22.501 | 1 |
| pooled | paga_dpt | 0 | 2 | 0.188 | -61.071 | 1 |
| pooled | paga_dpt | 0 | 3 | 0.130 | -63.006 | 1 |
| pooled | paga_dpt | 0 | 4 | 0.018 | -5.328 | 1 |
| pooled | paga_dpt | 1 | 2 | 0.043 | 0.290 | 1 |
| pooled | paga_dpt | 1 | 3 | 0.097 | -4.617 | 1 |
| pooled | paga_dpt | 1 | 4 | 0.016 | 0.064 | 1 |
| pooled | paga_dpt | 2 | 3 | 0.622 | -55.788 | 1 |
| pooled | paga_dpt | 2 | 4 | 0.043 | 1.411 | 0.979 |
| pooled | paga_dpt | 3 | 4 | 0.032 | -0.253 | 1 |
| pooled | palantir | 0 | 1 | 0.181 | -73.226 | 1 |
| pooled | palantir | 0 | 2 | 0.427 | -1.978 | 1 |
| pooled | palantir | 0 | 3 | 0.654 | -58.292 | 1 |
| pooled | palantir | 0 | 4 | 0.313 | 0.005 | 1 |
| pooled | palantir | 1 | 2 | 0.056 | -38.107 | 1 |
| pooled | palantir | 1 | 3 | 0.114 | -50.430 | 1 |
| pooled | palantir | 1 | 4 | 0.033 | -32.787 | 1 |
| pooled | palantir | 2 | 3 | 0.207 | -5.417 | 1 |
| pooled | palantir | 2 | 4 | 0.869 | 2.969 | 0.00999 |
| pooled | palantir | 3 | 4 | 0.120 | -3.054 | 1 |
| pooled | scfates_ppt | 0 | 1 | 0.632 | -57.474 | 1 |
| pooled | scfates_ppt | 0 | 2 | 0.025 | -0.847 | 1 |
| pooled | scfates_ppt | 0 | 3 | 0.099 | -0.482 | 1 |
| pooled | scfates_ppt | 0 | 4 | 0.608 | -59.876 | 1 |
| pooled | scfates_ppt | 1 | 2 | 0.074 | -0.609 | 1 |
| pooled | scfates_ppt | 1 | 3 | 0.114 | -0.754 | 1 |
| pooled | scfates_ppt | 1 | 4 | 0.966 | -1.008 | 1 |
| pooled | scfates_ppt | 2 | 3 | 0.037 | -0.653 | 1 |
| pooled | scfates_ppt | 2 | 4 | 0.070 | -0.720 | 1 |
| pooled | scfates_ppt | 3 | 4 | 0.139 | -0.591 | 1 |

#### Stability

| Patient | Method | Kind | R | Failed | Share matched median (C3) | Share matched min | Spearman median (C4) | Spearman min | Mean Jaccard median | Mean Jaccard min | Lineages min |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pooled | paga_dpt | reseeds | 2 | 0 | 0.600 | 0.200 | 1.000 | 1.000 | 0.938 | 0.875 | 1 |
| pooled | paga_dpt | subsamples | 3 | 0 | 0.200 | 0.200 | 0.364 | 0.022 | 0.644 | 0.489 | 1 |
| pooled | palantir | reseeds | 2 | 0 | 0.400 | 0.000 | 0.992 | 0.991 | 0.695 | 0.582 | 1 |
| pooled | palantir | subsamples | 3 | 0 | 0.400 | 0.000 | 0.945 | 0.904 | 0.678 | 0.582 | 1 |
| pooled | scfates_ppt | reseeds | 2 | 0 | 0.500 | 0.200 | 0.852 | 0.747 | 0.589 | 0.532 | 5 |
| pooled | scfates_ppt | subsamples | 3 | 0 | 0.400 | 0.000 | 0.758 | 0.755 | 0.636 | 0.508 | 5 |

#### Cell cycle

| Patient | Method | Cell-cycle status | Dataset share cycling | Warnings |
|---|---|---|---|---|
| pooled | paga_dpt | scored | 0.209 | – |
| pooled | palantir | scored | 0.209 | Lineage Cluster 12 is mostly cycling cells; it may be a cell-cycle state rather than a fate. |
| pooled | scfates_ppt | scored | 0.209 | – |
<!-- case-study:end -->

## Reading

**How to read a verdict.** The verdict rule failed its verification on 2026-11-17: 28 of 31
negative-control runs were verdicted supported. So no run is shown as supported. Tree-likeness (C2)
and both stability checks (C3 and C4) are left out, since none has a calibrated threshold, so the
verdict reads coverage (C1) and clonal support (C5) only. C1 passed in every run (coverage is 1.0 in
the results file). C5 passes when the headline p is at most 0.01 and Δη is above 0. A `weak` run is
one whose clonal check failed: the clones did not meet the declared bar for backing its lineages,
and that is the finding. A `provisional` run passed both checks that are read, which means "not shown to be a tree"; it is not
support. Stability is reported as numbers only.

**Pseudotime direction** comes from the declared root (the IP cells), which is an assumption: the
paper's main text and Methods name no root (unverified: supplementary figures not read). Nothing
about direction is read from the clones.

**One bar per run.** The 0.01 bar is per run and was declared in advance. With B = 1000 the smallest
possible p is 0.000999, above the Bonferroni bar of 0.01/12 for 12 runs, and the runs are not
independent (three methods on the same cells of each patient). So no corrected, family-wise claim can
be made. The count below is "4 of 12 runs passed C5", nothing more.

### CLL-1

- **paga_dpt: provisional, not shown to be a tree.** The headline p is 0.000999 and Δη is 0.036
  (p_Δη 0.000999), so C5 passed. η falls from 0.695 against N1 to 0.363 against N2 and 0.247
  against N3, so part of the clonal concentration does not survive the stricter nulls; η_base
  is 0.659, so expression clustering alone matches most of it. The m = 3 check agrees (η 0.626, p
  0.000999). No lineage pair shares more clones than chance (every computed p_adj is 1). Reseed
  Spearman is 1.000; subsample share matched has a median of 0.400. Lineage Cluster 5 is flagged as mostly
  cycling, so clonal agreement there may follow a cycling state rather than a fate.
- **palantir: weak.** The clonal check failed: the headline p is 0.0759 and Δη is -0.031 (p_Δη
  0.996), so the clonal agreement with these lineages meets neither the declared p bar once time
  point is controlled nor the baseline. η is 0.125 against N1 and 0.020 against N3. Only 3 time points
  have scored cells. No coupling pair passes (every computed p_adj is 1).
  Reseed Spearman is 0.999, subsample share matched 0.400 (minimum 0.200). No cell-cycle flag.
- **scfates_ppt: weak.** The clonal check failed on the baseline: the headline p is 0.000999, but
  Δη is -0.005 (p_Δη 0.766), so the lineages do not beat the k-means baseline (η 0.571,
  η_base 0.577). η drops to 0.052 against N2. At m = 3 the p rises to 0.0509. Lineages 1 and 2 share
  more clones than chance (z 3.086, p_adj 0.02), which says they are related, not in what order.
  Reseed Spearman is 0.986, subsample share matched 0.800. Lineage Cluster 8 is flagged as mostly
  cycling.

### CLL-2

- **paga_dpt: provisional, not shown to be a tree.** The headline p is 0.000999 and Δη is 0.024
  (p_Δη 0.002). The effect is small and mostly does not survive the time-point null: η is 0.213 against N1
  and 0.024 against N3. It does not hold at m = 3 (p 0.154), so it is fragile. Its scored cells
  spread over the most effective clones of any run (649.190). No coupling pair passes. Reseeds agree fully (share matched
  1.000, Spearman 1.000), but subsamples do not: share matched has a median of 0.600 and a minimum of
  0.200, and one subsample found a single lineage (lineages minimum 1). Lineage Cluster 8 is flagged
  as mostly cycling.
- **palantir: weak.** The clonal check failed: the headline p is 0.0699 and Δη is -0.011 (p_Δη
  0.792). At m = 3, p is 0.147. No coupling pair passes. One subsample replicate failed (an
  `InvariantError` in Palantir's fate probabilities), which the stability engine counts as no
  agreement; that one failure is why the subsample minimums are 0.000 and the lineages minimum is 0.
  The subsample Spearman median is 0.964. No cell-cycle flag.
- **scfates_ppt: weak.** The clonal check failed on p: the headline p is 0.155, although Δη is 0.018
  (p_Δη 0.153). η is 0.143 against N1 and 0.011 against N3, so almost none of it survives the time-point
  null. At m = 3, p is 0.364. No coupling pair passes (lineages 1 and 2: p_adj 0.527). Subsample share
  matched has a median of 1.000 but a minimum of 0.000. No cell-cycle flag.

### NHL-6

- **paga_dpt: weak.** The clonal check failed on the baseline: the headline p is 0.000999 (also
  0.000999 at m = 3), but Δη is -0.008 (p_Δη 0.997): η 0.415 against η_base 0.422. The clones are
  concentrated (58.793 effective clones from 548), yet the k-means baseline scores them as high as
  the lineages do. Every coupling z is strongly negative (lineage pairs share fewer clones than
  chance) and no pair passes. Reseed Spearman is 1.000, subsample share matched 1.000. No cell-cycle
  flag.
- **palantir: weak.** The clonal check failed on the baseline: the headline p is 0.000999, but Δη is
  -0.022 (p_Δη 1). Only one lineage pair has scored cells on both sides, and it does not pass.
  Reseed share matched is 0.333 with Spearman 0.999. No cell-cycle flag.
- **scfates_ppt: weak.** The clonal check failed on p: the headline p is 0.026, and Δη is -0.050
  (p_Δη 1), the largest shortfall against the baseline of any per-patient run. η is 0.313 against N1
  but 0.006 against N2, so most of the clonal concentration does not survive the pseudotime-depth
  null. Lineages 3 and 4 share more clones than chance (z 2.230, p_adj 0.02): related, in no
  stated order. Stability is high (reseed share matched 1.000, Spearman 0.996), which shows these
  lineages are reproducible, not that they are a tree. No cell-cycle flag.

### NHL-7

- **paga_dpt: weak.** The clonal check failed on p: the headline p is 0.015, although Δη is 0.032
  (p_Δη 0.000999). η is 0.216 against N1 and 0.009 against N3, so almost none of it survives the time-point
  null. At m = 3, p is 0.014. No coupling pair passes. Reseed Spearman is 1.000; subsample share
  matched 0.667 (minimum 0.333). No cell-cycle flag.
- **palantir: provisional, not shown to be a tree.** The headline p is 0.005 and Δη is 0.010 (p_Δη
  0.019), so C5 passed by a small margin; the m = 3 check agrees (p 0.00599). The effect is small (η
  0.064 against N1, 0.028 against N3). No coupling pair passes. Subsample share matched is 0.333 and
  subsample Spearman 0.898, the lowest of the per-patient runs. Lineage Cluster 4 is flagged as mostly
  cycling.
- **scfates_ppt: provisional, not shown to be a tree.** The headline p is 0.000999 and Δη is 0.049
  (p_Δη 0.000999), the largest gain over the baseline of any run, and the m = 3 check agrees (p
  0.000999). It rests on the fewest cells and clones of any run (1584 cells, 185 clones) and the
  effect is small (η 0.070 against N1, η_base 0.021); purity is 0.351. No coupling pair passes.
  Its subsample Spearman, 0.991, is the highest of any run (share matched 1.000). Lineage
  Cluster 4 is flagged as mostly cycling.

### Across patients

- **No tree is shown in any patient.** 8 of 12 runs are weak and 4 are provisional; none is
  supported or not_supported. 7 of 12 runs have a headline p at most 0.01, and 4 of them also beat
  the clustering baseline (Δη above 0), which C5 requires: CLL-1 and CLL-2 `paga_dpt`, NHL-7
  `palantir` and `scfates_ppt`. The 3 runs with p 0.000999 that are weak (CLL-1 `scfates_ppt`, NHL-6
  `paga_dpt` and `palantir`) fail on Δη at or below 0.
- **No method holds across patients.** `paga_dpt` passes C5 in CLL-1 and CLL-2 only, `palantir` and
  `scfates_ppt` in NHL-7 only. Every patient had enough clones to score (the fewest is 185).
- **Expression clustering matches most of the N1 clonal score.** Δη is at most 0.050 in absolute
  value in every per-patient run and is at or below 0 in 6 of 12. Relative to η, the exception is
  NHL-7 `scfates_ppt` (η 0.070, η_base 0.021). T-cell clones tend to be transcriptionally alike,
  which is why the baseline exists (spec §4.3; unverified for this dataset); here the lineages add
  little to what the k-means partition scores.
- **Much of the clonal concentration does not survive the time-point and pseudotime-depth nulls.**
  η against N2 or N3 is below half of η against N1 in every per-patient run (for example 0.415 to
  0.089 in NHL-6 `paga_dpt`, 0.313 to 0.006 in NHL-6 `scfates_ppt`). Against N1 alone every run has
  p 0.000999; every headline p above 0.01 comes from N3 or N2. Two readings fit, and these numbers
  do not separate them: clones may be specific to a time point, and the lineages may themselves
  follow time point, since the root is defined by it and pseudotime runs from it. Either way, η
  against N1 overstates what the lineages add. The paper reports a repertoire that changes after
  infusion in its bulk TCRB cohort of ten other patients (Fig. 1b, 1c) and, in this cohort,
  expression that changes with time after infusion (Fig. 5a); lumenCell does not test a link
  between either and these lineages.
- **The passes are small and partly fragile.** Of the 4 C5 passes, CLL-2 `paga_dpt` fails at m = 3
  (p 0.154) and NHL-7 `palantir` has Δη 0.010 with p_Δη 0.019.
- **Cell cycle.** 5 of 12 runs flag a lineage as mostly cycling, and all 4 C5 passes are among them.
  The clonal agreement in those runs may follow a cycling state rather than a fate. The paper's IP
  cluster 4 is a cycling cluster (Supplementary Fig. 4a); whether lumenCell's flagged lineages match
  it was not checked.
- **The largest clone.** Dropping the largest clone moves purity by at most 0.019 in any
  per-patient run; η was not re-scored without it. Effective clones range from 34.329 (NHL-7
  `palantir`) to 649.190 (CLL-2 `paga_dpt`).
- **Coupling.** 2 per-patient lineage pairs share more clones than chance after Benjamini–Hochberg
  (CLL-1 `scfates_ppt` lineages 1 and 2, NHL-6 `scfates_ppt` lineages 3 and 4, both p_adj 0.02). Most
  z are negative: lineages share fewer clones than chance, as expected when clones sit within one
  lineage. Coupling says two lineages are related, never which came first.
- **Stability.** Reseed Spearman medians are 0.986 to 1.000 (minimums 0.983 to 1.000) and subsample
  Spearman medians 0.898 to 0.991 across the per-patient runs. Subsample share matched minimums fall
  to 0.000 in CLL-2 `palantir` and `scfates_ppt`, against reseed minimums of 0.333 to 1.000;
  share matched medians span 0.333 to 1.000 for both reseeds and subsamples. Stability says a result
  is reproducible, not that it is a tree, and it is not in the verdict.

### Pooled diagnostic

The pooled run is a diagnostic of patient structure, never a finding. The pooled diagnostic puts all
four patients' 62,167 cells into one run with no batch correction. In the pooled diagnostic,
`paga_dpt` is provisional (headline p 0.002, Δη 0.034), `palantir` is weak (p 0.0799) and
`scfates_ppt` is weak (p 0.475, Δη -0.082). Per patient, `scfates_ppt` is provisional only in NHL-7,
and the pooled diagnostic run of `scfates_ppt` is weak. The pooled diagnostic's subsample Spearman
median is 0.364 for `paga_dpt`, against 0.947 to 0.990 per patient, and 0.758 for `scfates_ppt`,
against 0.908 to 0.991 per patient; for `palantir` the pooled diagnostic value, 0.945, is not lower
than the per-patient range (0.898 to 0.964). Pooling without batch correction may mix patient
structure into the lineages; this pooled diagnostic does not test that, and it is a diagnostic
only. One pooled diagnostic lineage pair passes coupling (`palantir` lineages 2 and 4, p_adj
0.00999), and like the rest of the pooled diagnostic it is not read as a result.

## Limits

- **Each per-patient run holds on average a quarter of the cells** (14% to 35%: 8,522 to 21,762
  of 62,167), and the
  clonal scores rest on fewer still (1584 to 9539 scored cells).
- **Each method scores its own cells.** Only cells in exactly one of a run's lineages are scored, so
  the three methods are scored on different cells of the same patient and their numbers are not a
  like-for-like comparison.
- **lumenCell has no batch correction,** which is why pooling is not the headline. The authors
  regressed out patient before clustering; lumenCell runs each patient alone instead.
- **Unpaired chains were not used as such.** The GEO tables are labelled "unpaired"; the declaration
  groups their chains by cell barcode into paired labels. Cells with only one chain detected have no
  label, and whether the tables hold only productive chains is not stated (unverified).
- **Only four patients,** two CLL and two NHL, with CAR-T cells made from different starting CD8
  populations. No result here generalises beyond them.
- **The infusion product is manufactured.** The root is a cultured product, stimulated during
  manufacture (Methods), and the paper reads its expression differences as partly due to that
  stimulation (Fig. 5b). Pseudotime from it measures distance from a cultured state, and the root
  itself is an assumption of the declaration.
- **The verdict cannot show support.** Its rule failed verification; only C1 and C5 are read, so a
  provisional verdict is never support.
- **Blood only.** All post-infusion samples are blood; the paper notes that the CAR-T cells at the
  tumour site were not sampled (Discussion).
- **Time-point labels.** CLL-1's d21 and d38 sit outside the paper's early and late windows
  (unverified why).

## Reproduce

```bash
# Fetch the two GEO files (about 117 MB; sizes and SHA-256 are checked) and build the derived file.
uv run lumencell bench fetch sheih2020
# The case-study numbers: 4 patients x 3 methods with stability, then the pooled diagnostic run.
uv run python -I benchmarks/case_studies/sheih2020.py
# Refill the generated tables of this document from the newest results file.
uv run python -c "import json, pathlib; from lumencell.bench.case_study import render_markdown; from lumencell.bench.report import replace_section; p = sorted(pathlib.Path('benchmarks/case_studies/results').glob('sheih2020-*.json'))[-1]; replace_section(pathlib.Path('docs/case-studies/sheih-2020.md'), render_markdown(json.loads(p.read_text())), '<!-- case-study:start -->', '<!-- case-study:end -->')"
# Check that the tables match the newest results and that the text keeps its rules.
uv run pytest tests/unit/test_case_study_writeup.py -q
```

A re-run writes a new dated results file; its numbers replace these only after the tables are
refilled and this text is checked against them. The 2026-11-20 run took 14 minutes.
