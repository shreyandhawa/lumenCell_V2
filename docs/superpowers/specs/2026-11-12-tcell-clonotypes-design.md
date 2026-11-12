# T-cell clonotypes as lineage evidence: design

Date: 2026-11-12. Status: draft for Shreyan's approval, revision 2. Branch `pivot/01-tcell-spec`.
Baseline: `main` at 865ce22, `make ci` green (353 Python tests, 36 web tests, 1 end-to-end test).
Background notes: [datasets](../research/2026-11-12-tcr-datasets.md),
[clonal metrics](../research/2026-11-12-clonal-metrics.md).

Revision 2 folds in an independent method review and an architecture review of revision 1
(§11 lists what changed and why).

## 1. Purpose

lumenCell infers pseudotime and lineages. Every benchmark so far is developmental, and none
contains T cells. This work points lumenCell at T-cell data and uses T-cell receptor (TCR)
clonotypes as evidence that is independent of the inference: cells that share a receptor
descend from one ancestor. A hiring reviewer will ask two questions, and the work exists to
answer both:

1. *What did you find with this tool?* This is answered by the Daniel 2022 benchmark and the
   Sheih 2020 case study (phase 3).
2. *How do you validate a trajectory with no ground truth?* This is answered by clonal
   agreement against null models and against a lineage-free baseline (phase 1), stability
   under reseeding and subsampling, and a plain verdict when the data doesn't support a tree
   (phase 2).

Success means each phase leaves working, tested software, and every reported number can be
defended. A gate that fails honestly is an acceptable result. A number that passes because
something was tuned on the evaluation data is not.

This applies known ideas. The clonal scores derive from STARTRAC (Zhang et al. 2018), LRT
(Xie et al. 2023) and CoSpar (Wang et al. 2022). The normalised effect size follows
chance-corrected information measures (Vinh et al. 2010). The p-value follows Phipson and
Smyth (2010). Sources are in §10.

## 2. Global constraints

These rules come from Shreyan's brief and bind every phase.

**Science**
- Clone labels are used only for scoring. They never reach inference: not through obs (a
  prep stage, a method, a root, a terminal, `cluster_key`, `label_key`), and not through
  expression, since receptor segment genes are excluded from HVG selection (§4.2).
- Each metric, each null model and the baseline are defined here and in
  `benchmarks/thresholds.toml` before any result is seen. Every score is reported against its
  nulls.
- Clones are keyed by donor. The same receptor string in two donors is two clones.
- The minimum clone size is stated, and every score reports how many cells and clones it rests
  on, and the effective number of clones. Single-cell clones never contribute.
- Clone sharing shows that two states are related, not which came first. No output or text
  claims direction from clones.
- No method default and no threshold is tuned on T-cell data. Both T-cell datasets are for
  evaluation only.
- Thresholds are never lowered to fit results. New gates and definitions go into
  `thresholds.toml` with their reasons before the first run that measures them.
- In the case study, the biology comes from the dataset's paper and is cited. Anything
  unverified is marked. If the data doesn't fit a tree, that is the result.

**Existing decisions that stay**
- Roots are always explicit.
- Inference never uses the UMAP.
- The input file is never modified.
- A method's `experimental` flag matches its synthetic gates
  (`tests/golden/test_synthetic.py`).
- API changes start in `src/lumencell/api/models.py` and go through `make contract`.
- The server process never runs prep or methods; heavy work runs in the worker. A GET that
  finds a missing artefact may submit a job, as cache healing does
  (`docs/backend-robustness-2.md`).
- State-file writes are locked and atomic (`api/state.py`). Run and cache directories are
  written by atomic rename.
- `engine/` never imports `bench/`.

**Process**
- Test-first. One branch and PR per task, named `pivot/NN-<slug>`. Squash-merge only when CI
  is green, after a review by a different agent from the author. Nothing is pushed to `main`
  directly.
- No `Co-Authored-By` trailers and no "Generated with" footers. Gate each commit on the exit
  code of `make ci`.
- No new runtime dependency without asking. scanpy and SciPy are declared runtime
  dependencies, and Hypothesis is a dev dependency. scikit-learn arrives only through scanpy,
  so it is not imported directly.
- Datasets are never committed. Downloads follow the existing fetch pattern: ask first, record
  the SHA-256, keep certificate checks on. Tests that need real data are marked `real_data`
  and skipped by default. CI uses synthetic data and small fixtures only.
- Values read from a user's file (column names, labels) are untrusted. They never reach a file
  path or a shell command. Files derived from them use hashed names.
- In `web/`, only `make contract` output and the smallest fixture updates needed to keep
  `make ci` green: the `dataset` and `run` fixtures in `web/src/test/fixtures.ts`. Web views
  come later, in their own prompt.
- **Runtime thresholds.** `thresholds.toml` is not in the wheel, so the values the worker
  needs at run time (`min_clones`, the effect floors, the verdict τs, the cell-cycle ratios)
  live as constants in `src/lumencell/trust/constants.py`. A test asserts they equal the TOML.

## 3. Terms

- **Clone label**: the per-cell value of the obs column the user names as the clone column.
  Missing, empty, `"nan"` and `"None"` values mean the cell has no receptor, and it is not
  scored.
- **Donor**: the per-cell value of the column the user names as the donor column. Clone
  identity is the pair (donor, clone label). A file without a donor column must say
  explicitly that it holds one donor; nothing is assumed.
- **Clone namespace**: the donor, or the finest unit that the sample records show holds one
  individual. For Daniel 2022 this is decided from the GEO sample records (§6.1).
- **Membership**: the boolean `(L, n)` matrix of lineage masks that every method returns
  (`engine/schema.py`). It overlaps on the shared trunk. Palantir's masks are thresholded fate
  probabilities, and its orphan cells are forced into their most likely fate.
- **Single-lineage cell**: a cell inside exactly one lineage mask.
- **Scored cells**:
  - In **own-cells mode** (API, case study), the run's single-lineage cells that carry a clone
    label whose clone has at least `min_clone_size` such cells.
  - In **oracle-cells mode** (benchmark gates), the same set defined on the *truth* lineages.
    The method's assignment of those cells is then used: a cell in one method lineage counts
    there, a cell in several counts as an equal fraction in each, and a cell in none counts
    equally across all L.
- **Null**: the same score recomputed after shuffling clone identities among the scored cells
  within strata. Three nulls are defined in §4.3.

## 4. Phase 1: Clonotypes end to end

### 4.1 Clonal settings and reading labels

The clone and donor columns are a **dataset-level setting**, like the root. They are stored in
the state file through `StateStore.update`, never guessed, and never part of the run id. The
loader needs no change, because obs is always kept in the server.

- New model `ClonalSettings { clone_key: str, donor_key: str | None, single_donor: bool,
  strata_keys: list[str] }`.
  - Exactly one of `donor_key` or `single_donor = true` must be given.
  - `strata_keys` is optional (see the third null in §4.3).
- Validation runs in the server:
  - Every named column exists in obs, and the clone and donor columns are different columns.
  - The clone column is categorical or string.
  - The donor column has between 1 and 1,000 distinct values.
  - Any failure answers 422 with the existing error body.
- `DatasetInfo` gains `clonal: ClonalSettings | None` and `clonal_candidates: list[str]`.
  Candidates are columns whose names match scirpy or 10x conventions (`clone_id`,
  `raw_clonotype_id`, `clonotype`, with or without an `airr:` prefix). They are suggestions
  only, never selected automatically.
- Routes:
  - `PUT /api/v1/dataset/clonal` sets the setting and returns a summary: cells with a label,
    clones, clones of size 2 or more, donors, and the clone-size quantiles. It only counts
    values, so it may run in the server. The PUT also rebuilds the server's cached inference
    context (`dataset_view.py`).
  - `DELETE /api/v1/dataset/clonal` clears the setting.
- **Fix the high-cardinality mismatch.** Today a column with more than 1,000 categories is
  listed as `string` but served as categorical codes (`dataset_view.py:162` against `:254`).
  Clone columns always hit this case. Listing and serving must apply one shared rule, with a
  test, so a column is never described one way and served another.
- No scirpy dependency. A file prepared with scirpy already carries `clone_id` in obs, and that
  is all lumenCell reads. Building clones from receptor contigs is a benchmark-builder concern
  (§6.1).

### 4.2 Keeping clones out of inference

There are two routes by which clone identity could reach inference. Both are closed.

**Obs.** `InferenceContext.obs` today carries every obs column (`engine/context.py`).
- The clone and donor columns in `ClonalSettings`, plus any column matching the receptor
  patterns. `strata_keys` columns (for example time point) carry no clone information, so they
  stay visible to methods and may define a root (amended at plan time)
  (`clone_id*`, `clonotype*`, `raw_clonotype_id`, `IR_*`, `airr:*`, `receptor_*`,
  `chain_pairing`, `cdr3*`), become **scoring-only**. They are removed from the obs that
  methods receive.
- `label_key` detection runs on the filtered obs.
- A scoring-only column is refused at three points:
  - in `resolve_root` (`engine/roots.py`), for a `cluster:` root
  - in `run_method` (`engine/registry.py`), as `cluster_key`, `label_key` or a terminal label
    column
  - in `Runs.create_run`, which answers 422 with a clear message
- `score_clonal` refuses any run whose recorded provenance used a scoring-only column.

**Expression.** TCR segment genes are in the count matrix, and in T-cell data they can become
highly variable genes (HVGs). That would pull clone-mates together in the PCA and kNN graph.
- A new `PrepParams.exclude_receptor_genes = True` drops genes matching
  `^TR[ABDG]([VDJ]\d|C\d*$)` from HVG selection. The match is case-insensitive, so it covers
  human and mouse symbols. The pattern was tightened at plan time from `^TR[ABDG][VDJC]`, which
  also caught the signalling gene TRADD. This covers
  V, D and J segments and the constant genes, because TRBC1 against TRBC2 usage is clonal.
- The prep info records how many receptor genes would otherwise have been HVGs.
- The change alters the HVG stage key, so every dataset is re-prepared once. Synthetic gene
  names (`G00000`) match nothing, so synthetic results are unchanged. The three real
  benchmarks are re-run, and any change is logged in the `thresholds.toml` changelog.
- Immunoglobulin genes are not excluded in this work, to avoid changing developmental
  benchmarks. This is listed in §9.

**Run records.** `result.json` gains the run's `prep_key` and an `inference_isolation`
version, and `store_version` is bumped. Runs made before this work drop out of listings, as
the store already does for other versions, so no legacy run is ever scored. This costs one
re-run per dataset.

**Isolation test.** Load one `Dataset`. Run each method twice on it: once on the plain
`Dataset`, and once with a clone column added in memory that copies the true lineage, plus TR
genes in the matrix whose expression encodes the clone. The two results must be identical
value for value. Comparing run ids would prove nothing, because obs does not enter the run id.

### 4.3 Clonal scores

New package `src/lumencell/clonal/`, holding pure functions on arrays with no engine, API or
anndata imports:
- `labels.py` turns raw labels and donors into integer clone ids keyed by (donor, label),
  plus donor and stratum indices.
- `scores.py` holds the metrics, the nulls and the baseline.

The benchmark and the worker both call this package.

**Lineage distribution.** For each clone k, the lineage distribution p_k over the L lineages
is built from (possibly fractional) cell counts. Its normalised entropy is
H_k = −Σ_l p_kl log p_kl / log L, and the mass-weighted mean entropy is
H̄ = Σ_k n_k H_k / Σ_k n_k.

**Primary effect: normalised clonal excess.**

    η = 1 − H̄_obs / E_null[H̄]

A value of 0 means "no better than chance" and 1 means every clone sits in one lineage. The
`log L` cancels, so η is comparable across different numbers of lineages, much like adjusted
mutual information. The raw purity P = 1 − H̄_obs (after STARTRAC's transition index and
LRT's lineage-bias test) is reported alongside it.

**Nulls.** Each shuffles clone identities among the scored cells, B = 1,000 times, with a
seeded generator:
- **N1**: within donor.
- **N2**: within donor × decile of pseudotime. This guards against clones that sit at one depth
  across lineages. Each run uses its own pseudotime; the oracle uses true pseudotime.
- **N3**: within donor × `strata_keys`, when declared (for example time point or tissue).

Shuffling among scored cells keeps the lineage mix of expanded cells fixed. The test then asks
only whether cells of one clone co-locate more than chance. STARTRAC shuffles among all
cells, which mixes in expansion bias between lineages.

**Headline p.** For each null, p = (1 + #{H̄_b ≤ H̄_obs}) / (B + 1), one-sided. The headline
p is the largest over the computed nulls. η is reported against each null; its headline value
is taken against N1.

**Lineage-free baseline.** T-cell clones are transcriptionally coherent (CoNGA, scirpy's
clonotype modularity), so any expression partition beats the shuffle null. Every score is
therefore also computed for a lineage-free partition: k-means with k = L on the same PCA
(`scipy.cluster.vq.kmeans2`, seeded), scored on the same cells. The result is reported as η_base, together with
Δη = η − η_base and a paired permutation p for Δη (the same shuffles, applied to both). The
claim "the lineages carry clonal information beyond expression clustering" needs Δη > 0.

**Reported for every score**
- the mode (own cells or oracle cells)
- cells scored, clones scored, and the effective clone count (inverse Simpson index of scored
  mass)
- donors, strata, m, B and the seeds
- P, H̄, the null mean and sd of H̄ for N1, N2 and N3, η against each null, the p for each null
  and the headline p
- η_base, Δη and p_Δη
- the clone-weighted purity (every clone weighted equally) and the purity with the largest
  clone left out, which show how far a few hyperexpanded clones drive the result
- the same block at m = 3, as a sensitivity check

The primary score uses m = 2.

**Secondary (reported, never gated): lineage-pair coupling.** Build the clone × lineage matrix
X of raw scored-cell counts. This is not CoSpar's normalised default, so the normalisation is
stated in the output. For each lineage pair (a, b), compute the cosine similarity of columns
a and b, and compare it with N1, which gives z and p. Benjamini–Hochberg adjustment is applied
across the pairs. A high z says the two lineages share more clones than chance. It says they
are related, not which came first.

**Statuses**
- `ok`
- `not_configured`: no clonal setting
- `single_lineage`: L < 2
- `too_few_clones`: fewer than `min_clones` = 20 clones scored; the values are still reported
  but never count as a pass
- `refused`: the run used a scoring-only column
- `pending`: a job is queued
- `failed`: an error, with its message

### 4.4 Synthetic clones

A new module `bench/synth_clones.py` works on an AnnData that has already been generated,
using only `uns["lumencell_truth"]`. It imports nothing private from `synth.py`, and the
runner applies it. Generator output is unchanged byte for byte; a test pins the hash of X, obs
and uns for a fixed case.

- **RNG.** A separate generator, `np.random.default_rng([seed, CLONE_STREAM])`. The draw order
  and `benchmarks/baseline.json` are untouched.
- **Labelling.** Each cell independently has a receptor with probability 0.85, close to 10x
  capture rates. This fixes the labelled share at every `fate_bias`.
- **Donors.** Each cell gets one of D donors (default 3), at random.
- **Pools.**
  - Within each donor, every labelled cell on a **leaf** segment joins that leaf's pool with
    probability `fate_bias`, which makes its clone founded after its last fork.
  - All other labelled cells (trunk, intermediate segments, and the remaining leaf cells) join
    the donor's trunk pool. A trunk-founded clone may span every lineage.
- **Clone sizes.** Each pool is cut into clones whose sizes are drawn independently from one
  fixed heavy-tailed distribution with a large singleton share. Only the last clone in a pool
  is truncated. This means the size distribution does not depend on `fate_bias`, and a test
  asserts that the labelled share and the size quantiles match across `fate_bias` levels.
- **What `fate_bias` means.**
  - `fate_bias = 1`: every leaf cell's clone is confined to one leaf, on every topology
    including binary2.
  - `fate_bias = 0`: clone grouping is independent of lineage within a donor, so this is the
    null case.
- **Collisions.** Clone label strings restart in each donor, so every string occurs in several
  donors.
- **Controls**, used in tests:
  - *Depth-only*: clone members are drawn from a pseudotime window across lineages.
  - *Donor confound*: donors differ in lineage composition, and clones are random within each
    donor.
- **Outputs.** obs gets `clone_id` and `donor`, and `uns["lumencell_truth"]["clones"]` gets the
  spec. `Truth` gains optional `clone_ids` and `donors` fields, following how `resolvable` and
  `terminals` were added.

### 4.5 Benchmark metrics and gates

- Every branching case in the synthetic matrix gets a clonal layer with `fate_bias = 0.7`. The
  counts are unchanged, so every existing metric is unchanged.
- **Oracle.** The score on the truth lineages, in oracle-cells mode, gives η_oracle. The
  method is scored on the same cells (oracle-cells mode), which gives η_method.
- **Row fields.** `clonal_excess_ratio = η_method / η_oracle`, `clonal_p` (headline),
  `clonal_delta_base` (Δη against the baseline) and `clonal_status`. `CaseResult.row()` always
  emits these fields, with NaN when they do not apply, so `check_regression` never hits a
  missing key.
- **Gates**, written into `thresholds.toml` with their reasons before the first measuring run:
  - `clonal_excess_ratio ≥ 0.5`: a method whose lineages capture less than half of the
    chance-corrected clonal signal of the true lineages is not separating the fates the clones
    show.
  - `clonal_p ≤ 0.01`: with B = 1,000 the smallest p is about 0.001, and its Monte Carlo error
    at the boundary is about 0.003.
- **The oracle must clear a floor.** The gate is assessable only if η_oracle ≥ 0.05 and the
  oracle's headline p ≤ 0.01. Below that, the ratio is dominated by noise, and the case
  records `not_assessable`.
- **NaN policy per gate.** `check_gates` skips NaN today, and the existing gates keep that
  behaviour. The clonal gates declare `nan = "fail"`, except when the truth has fewer than 2
  lineages, which is structurally not applicable. A method that returns L < 2 on a branching
  truth fails. `too_few_clones`, `failed` and `not_assessable` all fail, and the reason is
  recorded.
- **The clonal gate is not a branch gate.** Merging two true lineages can push the ratio above
  1. Existing gates such as `lineage_recall` catch that, and the report says so.
- `clonal_excess_ratio` joins `HIGHER_IS_BETTER`. The baseline is re-measured once, on purpose,
  in its own commit.
- **Experimental flags.** They cannot change: all three methods already fail at least one
  synthetic gate (checked against `baseline.json`). The golden test still runs.

**Tests of the scores**, in CI unless marked slow:
- Hand-computed examples, including fractional counts in oracle-cells mode.
- Hypothesis properties:
  - η, P and the p-values are unchanged when clones, lineages or cells are relabelled or
    reordered.
  - Singleton clones never change a score.
  - 0 ≤ P ≤ 1.
  - Splitting clones across donors never merges them.
  - Ignoring donors changes the score when label strings collide.
- Calibration (slow, with pre-declared seeds):
  - At `fate_bias = 0`, scored on each method's own masks, the share of headline p ≤ 0.05
    across 200 seeds must not exceed the binomial 99% upper bound. The test is one-sided,
    because permutation p-values are discrete and conservative.
  - The donor-confound control is calibrated under N1, and the same data scored with donors
    ignored is not.
  - The depth-only control is calibrated under N2.
- Monotonicity: the oracle's η rises with `fate_bias` over {0, 0.35, 0.7, 1}, for a
  pre-declared seed set. At `fate_bias = 0`, the oracle's headline p is above 0.01 for that
  set.

### 4.6 API: clonal agreement for each run

New models in `api/models.py`:
- `ClonalScore`, with the fields listed under "Reported for every score" (§4.3).
- `LineageCoupling { a, b, observed, null_mean, z, p, p_adj }`.
- `ClonalAgreement { status, settings, primary: ClonalScore | None,
  sensitivity_m3: ClonalScore | None, coupling: list[LineageCoupling], notes: list[str] }`.
- `RunTrust { clonal: ClonalAgreement | None }`, which phase 2 extends.
- `TrajectoryRun.trust: RunTrust | None`.

Extend the `JobKind` literal, and make `job_info` report unknown kinds by their own name
rather than as `prepare` (`session.py:199`).

**Where the score is computed.** In the worker, by a `score_clonal` job, never in the server.
- The job runs in own-cells mode.
- It is submitted when an `infer` job completes, if settings exist.
- It is submitted for every complete run when the settings change.
- It is submitted by `GET` on a run that has neither the file nor a queued job, which covers a
  server restart and settings that changed mid-queue. This is the cache-healing pattern.
- Its dedup key is `score_clonal:{rid}:{h}`.

**Storage.** `RunStore.write_trust(rid, name, payload)` writes
`<run>/trust/clonal-<h>.json` atomically, using a temporary name that
`prep/housekeeping.py` recognises as a leftover. `<h>` is a SHA-256 prefix of the settings,
the metric version and the parameters, so column names never appear in a path.
`RunStore.save` of an existing complete run keeps its `trust/` folder; today it deletes the
folder, and CLI `run` saves again.

**Job order.** The manager keeps two priorities: `prepare` and `infer` run ahead of
`score_clonal` and `stability`. Otherwise a cache-healing prepare would wait behind a long
background job.

**CLI.** `lumencell run` gains `--clone-key`, `--donor-key`, `--single-donor` and `--strata`,
which compute the clonal block in-process and print it.

### 4.7 Phase 1 done when

- `make ci` is green.
- `make bench-synth` shows no regressions against the deliberately re-measured baseline.
- `make bench-real` has been re-run after the TR-gene exclusion, and any change is logged.
- The isolation test, the property tests and the calibration tests pass. The calibration tests
  are slow and are run once and recorded.
- A run on a synthetic file with clones returns `trust.clonal` through the API and the CLI. A
  restart while scoring is pending heals itself.
- `thresholds.toml` holds the clonal definitions, gates and reasons, dated before the first
  measuring run, and `trust/constants.py` matches it.
- `docs/HANDOFF.md` stale lines are corrected once:
  - line 3: the commit
  - line 32: the 100k-cell budgets are measured
  - line 49: the test counts
  - lines 71–73: `keep_matrix`
  - lines 89–90: time limits and worker replacement
  - lines 159–160 and 185: GitHub Actions CI exists
- `/security-review` is run, with its findings resolved or recorded.

## 5. Phase 2: T-cell robustness

### 5.1 Cell-cycle stage

- New prep stage `cell_cycle`, with input `lognorm` and gene symbols from `ds.gene_symbols`.
  - Outputs: `cc_s`, `cc_g2m`, `cc_phase` and `cc_info` (JSON: status, resolved species, genes
    found, source).
  - All outputs are always written. When the stage cannot score, the per-cell arrays are NaN or
    `unscored` and `cc_info.status = "not_scored"` with a reason. This means the CI datasets,
    whose genes are named `G00000`, always take that path.
  - Its params are `PrepParams` fields: `cc_species ∈ {auto, human, mouse}` and
    `cc_gene_set = "tirosh2016"`.
- **Genes.** The S and G2M lists are from Tirosh et al. 2016, shipped as package data with
  their source. Mouse genes map through a fixed human-to-mouse ortholog table, shipped with its
  source and date. Both are verified during implementation.
- **Species.** `auto` resolves from the case of the matching symbols. With fewer than 20 list
  genes found, the stage reports `not_scored`.
- **Scoring** uses scanpy's `score_genes_cell_cycle` with a fixed seed.
- **Report-only.** `Stage` gains `report_only: bool`. Report-only stages count towards prep
  readiness but are left out of `prep_key`, so adding this stage changes no run id.
  `tests/unit/test_prep_cache.py` gains cases for that.
- **Per-run check.** For each lineage, compare the share of its single-lineage cells in S or
  G2M with the dataset's share. If a lineage's share is at least twice the dataset's share and
  at least 0.5, the run gets a warning: "Lineage X is mostly cycling cells; it may be a
  cell-cycle state rather than a fate." The thresholds go into `thresholds.toml` and
  `trust/constants.py`. Proliferating effector cells are real biology, so this is a hint, not
  a verdict.

### 5.2 Stability

- **Automatic** (Shreyan's choice). When an `infer` job completes, a `stability` job is
  submitted for that run at background priority. The run's result is served at once, and the
  stability block reads `pending` until the job finishes. A GET heals a missing job, the same
  way as for clonal scores. The job can be cancelled.
- **Replicates.** R = 5 per run, all seeded.
  - 2 **reseeds**: Leiden, then PAGA, then the method are re-run with new seeds on the full
    data, so the PAGA tree is never stale.
  - 3 **subsamples**: 80% of cells, uniformly at random, always keeping the reference root
    cells. kNN, Leiden, the diffusion maps and PAGA are recomputed through `run_stage` on a
    `Dataset` restricted to the subset. The full-data normalisation, HVG and PCA are reused.
  - Every replicate uses the root `cell:<obs_name>` of the reference run's resolved root, so a
    `cluster:` root is never resolved again on different data.
  - Terminals given as labels of a user column are passed through unchanged. Terminals on
    Leiden labels are mapped to the replicate's Leiden cluster with the largest overlap.
  - Replicate prep lives in the job's temporary directory and is not added to the stage cache.
- **Comparison with the reference**, in `engine/stability.py`, using an extension of
  `engine/compare.py`'s Hungarian matching on Jaccard that takes a shared-cell mask. Subsample
  results are padded to full length, and Jaccard is computed on shared cells only.
  - Pseudotime Spearman on shared cells.
  - The share of reference lineages matched at Jaccard ≥ 0.6, and the mean Jaccard of the
    matches.
  - The number of lineages in each replicate.
- **Reported**, separately for reseeds and subsamples: the median and minimum across
  replicates, the per-replicate rows, R, the fraction and the seeds.
- **Stored** through `RunStore.write_trust` as `stability-<h>.json`. The `stability` job kind
  gets an explicit time limit.
- **The benchmark** computes stability for every method on the real datasets, and for the
  calibration below.

### 5.3 Tree-support verdict

A plain signal with three levels: `supported`, `weak` or `not_supported`, plus `provisional`
while stability is pending. Each comes with its reasons in plain words.

| Check | Applies | Rule |
|---|---|---|
| C1 connectivity | always | coverage (share of cells with finite pseudotime) ≥ 0.95 |
| C2 tree-likeness | always (dataset level) | weight on PAGA's maximum spanning tree ÷ total PAGA edge weight, counting each edge once (upper triangle), ≥ τ_tree |
| C3 lineage stability | once stability is done | median over replicates of the share of reference lineages matched ≥ τ_lin |
| C4 pseudotime stability | once stability is done | median over replicates of the Spearman ≥ τ_pt |
| C5 clonal support | when clones are configured and L ≥ 2 | status `ok`, headline p ≤ 0.01, and Δη > 0 against the lineage-free baseline |

The verdict is `not_supported` if C1 or C2 fails, `weak` if any applicable check among C3–C5
fails, and `supported` otherwise.

**Thresholds come from synthetic data only, by a rule fixed in advance.**
- **Negative controls** live in a new module `bench/synth_controls.py`, which builds datasets
  directly rather than through the tree generator:
  - `blob`: one population with no expression change
  - `two_groups`: two unconnected populations
  - `loop`: a cyclic programme
  - `sheet`: two independent gradients, a continuum with no tree
  They sit in `[synthetic.controls]`, outside the gated matrix, so they can never affect the
  `experimental` flags.
- **Calibration** uses seeds 0–2 of the 12-case matrix and of the controls, all three methods
  pooled. For each check, τ = the midpoint between the 10th percentile of the matrix values and
  the 90th percentile of the control values. If the 10th percentile is not above the 90th, the
  check is recorded as non-discriminating and left out of the verdict.
- **C2 robustness.** C2 is also calibrated at Leiden resolution × 0.5 and × 2. If it doesn't
  discriminate at all three resolutions, it is left out.
- **Verification** on fresh seeds 3–5:
  - No control run may be verdicted `supported`.
  - The share of matrix runs verdicted `supported` is reported, not gated.
  - DPT on a blob is likely stable, so C4 may well prove non-discriminating. That is an
    acceptable, reported outcome.
- These slow runs are one-off. Their results and the derived τs are committed before the
  verdict is used on real data.

### 5.4 API

- `RunTrust` gains `stability: Stability | None`, `tree_support: TreeSupport | None` and
  `cell_cycle: CellCycleCheck | None`, with the models defined in the phase 2 plan.
- `DatasetInfo` gains the dataset-level C2 value and the cell-cycle info.
- The web change is again only `make contract` and the smallest fixture update.

### 5.5 Phase 2 done when

- `make ci` is green, with no `bench-synth` regressions.
- The cell-cycle stage is cached, tested like the other stages, and changes no run id.
- Stability is queued automatically at background priority, served, healed after a restart,
  and cancellable.
- The verdict τs are derived by the pre-declared rule on seeds 0–2 and verified on seeds 3–5,
  with the results committed.
- `/security-review` is run.

## 6. Phase 3: Real T-cell data and the case study

### 6.1 Benchmark: Daniel et al. 2022 (GSE188670, subseries GSE188666)

Mouse LCMV, acute (Armstrong) and chronic (Clone 13) infection, days 5, 8 and 21, spleen,
liver and lung; 10x 5′ with TCR. It is open access.

1. **Fetch.** With Shreyan's OK at the time, download the GSE188666 expression files and the
   24 per-sample `tcr_filtered_contig_annotations.csv.gz` files, about 1.1 GB. The ATAC tar is
   not fetched.
   - `RealDataset` gains multi-file `sources`, a `licence` field, and a build step
     (`bench/builders/daniel2022.py`) that writes one derived `.h5ad` into the cache. The
     downloads are never modified.
   - The SHA-256 of every file is recorded in the registry.
2. **Declare from metadata alone, before any clone is tabulated.** No clone × state, clone ×
   capture or clone × lineage count is computed before this PR merges. The oracle could be
   computed without running any method, so the freeze has to come first. In one PR:
   - **Subset.** The chronic arm only, because the Tpex root holds there. The days and tissues
     are chosen from the sample records, with reasons.
   - **Namespace.** The finest grouping that the GSM sample records show to hold one mouse. If
     captures are sort-gated by state and the mouse cannot be identified across them, a
     per-capture namespace would hide exactly the sharing between states that matters. The
     clonal result on Daniel is then declared `not_assessable` in advance, rather than computed
     on a misleading namespace. Whether captures pool mice is an assumption, stated as such.
   - **Clone definition.** A cell needs at least one productive, high-confidence, full-length
     TRA and at least one such TRB, with at most two of each; other cells are unlabelled. A
     clone is the sorted set of their CDR3 nucleotide sequences. Requiring both chains stops
     state-dependent chain dropout from splitting clones. A TRB-only definition is reported as
     a sensitivity check.
   - **Truth.** The paper's state labels give the lineages as state lists. The lists are
     justified from non-clonal literature on exhaustion subsets (Hudson et al. 2019, Zander et
     al. 2019, Beltra et al. 2020), not from Daniel's own clonal analysis. The root is at
     progenitor-exhausted (Tpex) cells (Im et al. 2016, Utzschneider et al. 2016). The exact
     label names come from the metadata. If GEO's metadata has no state labels, the fallback
     is day as the only order truth, with no lineage truth and no clonal gate, recorded as
     such.
   - **Gates.** `[real.daniel2022.gates]` reuses the existing real-data gate values, plus the
     two clonal gates in oracle-cells mode against the truth lineages, with the oracle floor.
3. **Run** all three methods with `make bench-real`, and regenerate the README tables.
   - Reported, not gated: pseudotime Spearman against day; clone recurrence across captures as
     an estimate of convergence; agreement of the clone definition with 10x `raw_clonotype_id`;
     and the TRB-only sensitivity.
   - Gate failures are reported as they are.

### 6.2 Case study: Sheih et al. 2020 (GSE125881)

Human CD8 CAR T cells from 4 patients at 4 time points. Licence CC BY 4.0.

1. **Verify first.** With Shreyan's OK, download one clonotype table (about 7 MB) and check
   that it carries per-cell barcodes that join to the expression matrix. If it does not,
   switch to Cao et al. 2025 (GSE207192). The GEO page says 3 patients and the paper says 7;
   that discrepancy is resolved before designing around it.
2. **Declare before analysing.** In one PR, written from the paper and the file headers only:
   - the clone definition (the GEO tables are labelled "unpaired", so the chains used are
     stated)
   - the namespace: each patient
   - `strata_keys = [time point]`, so the headline p includes N3 and clones that expand at one
     time point are not read as lineage agreement
   - the root, taken from the paper's description
3. **Analyse.** Run all three methods, and report for each:
   - the clonal block in own-cells mode, with the m = 3 sensitivity check, all three nulls and
     the lineage-free baseline
   - lineage-pair coupling
   - stability
   - the tree verdict
   - the cell-cycle check
4. **Write it up.** The numbers come from a script committed under `benchmarks/case_studies/`,
   so they are reproducible. The write-up goes to `docs/case-studies/sheih-2020.md`.
   - The biology comes from the paper and is cited. Anything unverified is marked.
   - No direction is claimed from clones.
   - If the data doesn't fit a tree, that is the finding.
5. **Check the statistics.** `data:statistical-analysis` and `data:validate-data` are run on
   the calibration tests and the case-study numbers.

### 6.3 Phase 3 done when

- Daniel 2022 is in the registry with its SHA-256 values, subset, namespace, clone definition,
  truth and gates, all merged before its first run.
- `make bench-real` has produced a results file, and the README tables are regenerated.
- The case-study declaration has merged before the analysis.
- The write-up is built from reproducible numbers, with citations.
- `/security-review` is run.

## 7. Interfaces and data flow

```
.h5ad ─loader─▶ Dataset (obs kept; matrix dropped in server)
   └─ dataset_view: one cardinality rule; ClonalSettings validated ─▶ state file; context rebuilt
infer job (worker, foreground):
   prep cache (HVG without TR genes) ─▶ InferenceContext (obs minus scoring-only columns)
   ─▶ run_method (refuses scoring-only keys) ─▶ RunStore.save (keeps trust/)
   └─▶ manager submits background jobs: score_clonal ─▶ trust/clonal-<h>.json
                                         stability   ─▶ trust/stability-<h>.json
GET run (server): result.json + trust/*.json + verdict ─▶ TrajectoryRun.trust
                  (missing file and no job ─▶ submit, as cache healing does)
bench: synth ─▶ synth_clones (separate RNG) ─▶ Truth (+clone_ids, donors) ─▶ methods
       ─▶ clonal.scores (oracle-cells mode) ─▶ row (+clonal_*) ─▶ gates (per-gate NaN policy)
```

Each unit has one purpose:
- `clonal/labels.py`: labels to clone ids.
- `clonal/scores.py`: scores, nulls and the baseline.
- `bench/synth_clones.py`: clone simulation.
- `bench/synth_controls.py`: the negative controls.
- `engine/stability.py`: replicates and comparison.
- `engine/tree_support.py`: the verdict.
- `trust/constants.py`: runtime thresholds.
- `prep/stages.py`: the receptor-gene exclusion and the report-only `cell_cycle` stage.
- `bench/builders/daniel2022.py`: contigs to the derived `.h5ad`.

## 8. Decisions and what they cost if wrong

| # | Decision | Cost if wrong |
|---|---|---|
| 1 | Clone identity is (donor, label). A donor column, or an explicit single-donor flag, is required. | One extra setting for single-donor files. |
| 2 | Clonal settings are dataset-level state, not part of the run id. | Settings would have to be repeated if a user wants different clone columns per run. |
| 3 | Scoring-only columns are removed from obs, and refused as root, `cluster_key`, `label_key` or terminal column. | A legitimate column matching a receptor pattern can't be used by a method. |
| 4 | TR segment genes are excluded from HVG selection for every dataset. | One re-prep everywhere. A real-data benchmark result could move, which is logged. |
| 5 | `store_version` is bumped, so legacy runs are never scored. | Users re-run once. |
| 6 | No scirpy dependency. | AIRR-format input needs a later optional extra. |
| 7 | The primary effect is η = 1 − H̄/E_null[H̄] over single-lineage cells, with gates in oracle-cells mode. | It discards cells on intermediate segments in multi-fork trees. A changed definition must be logged, never swapped silently. |
| 8 | Three nulls (donor, donor × pseudotime decile, donor × declared strata). The headline p is the largest. | More conservative, so weak real signals may not clear p ≤ 0.01. That is acceptable. |
| 9 | A k-means baseline with k = L; C5 needs Δη > 0. | k-means is one choice of lineage-free partition. Leiden at matched L could differ. |
| 10 | m = 2 primary, m = 3 always reported. | None beyond reading two numbers. |
| 11 | The synthetic gate is an excess ratio ≥ 0.5 with p ≤ 0.01, and an oracle floor of η ≥ 0.05. | If 0.5 proves strict, methods fail it and stay experimental. The bar is not lowered. |
| 12 | Synthetic clones come from a non-depleting pool model, after the counts, on a separate stream. | None to existing results. |
| 13 | Cell cycle is report-only, left out of `prep_key`. | Users wanting regression wait for later work. |
| 14 | Stability runs automatically at background priority: R = 5, with 2 reseeds and 3 subsamples at 80%, PCA reused. | The worker is busy about 5 times longer per run. PCA reuse slightly understates subsampling variance. |
| 15 | Verdict τs come from a quantile midpoint on seeds 0–2, verified on seeds 3–5, against four controls. | A check that doesn't discriminate drops out, which is reported, not hidden. |
| 16 | The Daniel namespace comes from the GSM records. If sort-gated captures hide the mouse, Daniel's clonal result is declared not assessable in advance. | Daniel may give no clonal number at all. That is an honest outcome, and Sheih still provides one. |
| 17 | Daniel clones require paired chains; TRB-only is reported as a sensitivity check. | Fewer labelled cells. |
| 18 | The Daniel truth uses the paper's state labels, justified from non-clonal literature, with the root at Tpex and the chronic arm only. | The labels come from expression, like the other three benchmarks. The clones and the baseline remain the independent checks. |
| 19 | Sheih is verified before use, with Cao 2025 as the fallback; the time-point stratum is declared in advance. | One small download. |
| 20 | The missing `need_from_me` section of the brief is taken to mean dataset downloads and new runtime dependencies. | Something else Shreyan meant to approve gets decided without them. |

## 9. Not in this work

- web views
- the multi-user version
- new trajectory methods
- RNA velocity
- receptor design
- packaging
- the 100k-cell performance budgets
- cell-cycle regression
- excluding immunoglobulin genes
- AIRR or scirpy import in the app
- similarity-based clonotype clusters (only exact clones are used)
- any direction inferred from clones
- per-clone pseudotime coherence

## 10. References

The researchers checked the first group (see the notes). The second group is cited from
memory and is checked during implementation; until then it is marked as unverified.

**Checked**
- Zhang L. et al. STARTRAC. Nature 2018. doi:10.1038/s41586-018-0694-x
- Xie J. et al. LRT. PLoS Comput Biol 2023. doi:10.1371/journal.pcbi.1011300
- Wang S.-W. et al. CoSpar. Nat Biotechnol 2022. doi:10.1038/s41587-022-01209-1
- Weinreb C. et al. LARRY. Science 2020. doi:10.1126/science.aaw3381
- Weinreb C., Klein A.M. PNAS 2020. doi:10.1073/pnas.2000238117
- Daniel B. et al. Nat Immunol 2022. doi:10.1038/s41590-022-01337-5
- Sheih A. et al. Nat Commun 2020. doi:10.1038/s41467-019-13880-1
- Cao et al. Nat Commun 2025. doi:10.1038/s41467-025-59298-w

**Unverified**
- Phipson B., Smyth G.K. Permutation p-values should never be zero. Stat Appl Genet Mol Biol
  2010.
- Vinh N.X., Epps J., Bailey J. Information theoretic measures for clusterings comparison.
  JMLR 2010.
- Sturm G. et al. scirpy. Bioinformatics 2020.
- Schattgen S.A. et al. CoNGA. Nat Biotechnol 2022.
- Suo C. et al. Dandelion. Nat Biotechnol 2024.
- Saelens W. et al. A comparison of single-cell trajectory inference methods. Nat Biotechnol
  2019.
- Wolf F.A. et al. PAGA. Genome Biol 2019.
- Tirosh I. et al. Science 2016.
- Im S.J. et al. Nature 2016.
- Utzschneider D.T. et al. Immunity 2016.
- Hudson W.H. et al. Immunity 2019.
- Zander R. et al. Immunity 2019.
- Beltra J.-C. et al. Immunity 2020.

## 11. Changes in revision 2

**From the method review**
- Receptor genes are excluded from HVG selection (blocker).
- A lineage-free baseline is added, and C5 requires Δη > 0 (blocker).
- η replaces the raw excess, with gates in oracle-cells mode, an oracle floor, and L < 2 on a
  branching truth counted as a failure.
- The synthetic clone model is non-depleting and leaf-confined, with depth-only and
  donor-confound controls and one-sided calibration.
- The effective clone count and hyperexpansion diagnostics are reported.
- The stratified nulls feed the headline p.
- Daniel gets a metadata-first freeze, a namespace from the sample records, and paired chains.
- Verdict calibration is split from verification, with a quantile rule and four controls.
- Citations are added, along with BH adjustment and raw-count coupling.

**From the architecture review**
- Every refusal point for scoring-only columns is named, and `label_key` is filtered.
- `store_version` is bumped, so legacy runs are never scored.
- A GET heals a missing trust job, and background jobs get their own priority.
- `RunStore.write_trust` is added, and `save` keeps `trust/`.
- Stability uses `engine/compare.py` with a shared-cell mask, a `cell:` root, re-run Leiden and
  PAGA, and terminal mapping.
- The cell-cycle stage always emits its outputs, is report-only, and is left out of
  `prep_key`.
- Runtime thresholds live in `trust/constants.py`.
- The per-gate NaN policy and an always-emitted row are specified.
- Synthetic clones live in their own module, which avoids a circular import.
- The `dataset` fixture is listed.
