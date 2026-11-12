# T-cell clonotypes as lineage evidence: design

Date: 2026-11-12. Status: draft for Shreyan's approval. Branch `pivot/01-tcell-spec`.
Baseline: `main` at 865ce22, `make ci` green (353 Python tests, 36 web tests, 1 end-to-end test).
Background notes: [datasets](../research/2026-11-12-tcr-datasets.md),
[clonal metrics](../research/2026-11-12-clonal-metrics.md).

## 1. Purpose

lumenCell infers pseudotime and lineages. Every benchmark so far is developmental, and none
contains T cells. This work points lumenCell at T-cell data and uses T-cell receptor (TCR)
clonotypes as evidence that is independent of the inference: cells that share a receptor
descend from one ancestor. A hiring reviewer will ask two questions, and the work exists to
answer both:

1. *What did you find with this tool?* This is answered by the Daniel 2022 benchmark and the
   Sheih 2020 case study (phase 3).
2. *How do you validate a trajectory with no ground truth?* This is answered by clonal
   agreement against a null model (phase 1), stability under reseeding and subsampling, and a
   plain verdict when the data doesn't support a tree (phase 2).

Success means each phase leaves working, tested software, and every reported number can be
defended. A gate that fails honestly is an acceptable result. A number that passes because
something was tuned on the evaluation data is not.

This applies known ideas. The clonal scores derive from STARTRAC (Zhang et al. 2018), LRT
(Xie et al. 2023) and CoSpar (Wang et al. 2022), and are cited as such (§10).

## 2. Global constraints

These rules come from Shreyan's brief and bind every phase.

**Science**
- Clone labels are used only for scoring. They never reach inference: not a prep stage, a
  method, a root or a terminal (§4.2).
- Each metric and its null model are defined here and in `benchmarks/thresholds.toml` before
  any result is seen. Every score is reported against its null.
- Clones are keyed by donor. The same receptor string in two donors is two clones.
- The minimum clone size is stated, and every score reports how many cells and clones it rests
  on. Single-cell clones never contribute.
- Clone sharing shows that two states are related, not which came first. No output or text
  claims direction from clones.
- No method default is tuned on T-cell data. Both T-cell datasets are for evaluation only.
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
- The server process never runs prep or methods. Heavy work runs in the worker (the
  cache-healing invariant from `docs/backend-robustness-2.md`).
- State-file writes are locked and atomic (`api/state.py`). Run and cache directories are
  written by atomic rename.

**Process**
- Test-first. One branch and PR per task, named `pivot/NN-<slug>`. Squash-merge only when CI
  is green, after a review by a different agent from the author. Nothing is pushed to `main`
  directly.
- No `Co-Authored-By` trailers and no "Generated with" footers. Gate each commit on the exit
  code of `make ci`.
- No new runtime dependency without asking. Hypothesis is already a dev dependency.
- Datasets are never committed. Downloads follow the existing fetch pattern: ask first, record
  the SHA-256, keep certificate checks on. Tests that need real data are marked `real_data`
  and skipped by default. CI uses synthetic data and small fixtures only.
- Values read from a user's file (column names, labels) are untrusted. They never reach a file
  path or a shell command. Files derived from them use hashed names.
- In `web/`, only `make contract` output and the smallest fixture updates needed to keep
  `make ci` green. Web views come later, in their own prompt.

## 3. Terms

- **Clone label**: the per-cell value of the obs column the user names as the clone column.
  Missing, empty, `"nan"` and `"None"` values mean the cell has no receptor, and it is not
  scored.
- **Donor**: the per-cell value of the column the user names as the donor column. Clone
  identity is the pair (donor, clone label). A file without a donor column must say
  explicitly that it holds one donor; nothing is assumed.
- **Clone namespace**: the donor, or the narrowest unit known to hold one individual. For
  Daniel 2022 this is the GEO capture (§6.1).
- **Single-lineage cell**: a cell that falls inside exactly one inferred lineage's membership
  mask. Cells on a shared trunk are in several masks and say nothing about fate.
- **Scored cell**: a single-lineage cell with a clone label whose clone has at least
  `min_clone_size` single-lineage cells.
- **Null**: the same score recomputed after shuffling clone identities among the scored cells
  of each donor.

## 4. Phase 1: Clonotypes end to end

### 4.1 Reading clone labels

The clone and donor columns are a **dataset-level setting**, like the root. They are stored in
the state file through `StateStore.update`, never guessed, and never part of the run id.

- New model `ClonalSettings { clone_key: str, donor_key: str | None, single_donor: bool }`.
  Exactly one of `donor_key` or `single_donor = true` must be given.
- Validation runs in the server:
  - Both columns exist in obs, and they are different columns.
  - The clone column is categorical or string.
  - The donor column has between 1 and 1,000 distinct values.
  - Any failure answers 422 with the existing error body.
- The loader reads both columns from obs inside `load()`, before the matrix is dropped
  (`keep_matrix=False`), the same way `genes_detected` is handled today.
- `DatasetInfo` gains `clonal: ClonalSettings | None` and `clonal_candidates: list[str]`.
  Candidates are columns whose names match scirpy or 10x conventions (`clone_id`,
  `raw_clonotype_id`, `clonotype`, with or without an `airr:` prefix). They are suggestions
  only, never selected automatically.
- Routes: `PUT /api/v1/dataset/clonal` sets the setting and returns a summary: cells with a
  label, clones, clones of size 2 or more, donors, and the clone-size quantiles.
  `DELETE /api/v1/dataset/clonal` clears it. The summary only counts values, so it may run in
  the server.
- **Fix the high-cardinality mismatch.** Today a column with more than 1,000 categories is
  listed as `string` but served as categorical codes (`dataset_view.py`, `obs_payload`). Clone
  columns always hit this case. Listing and serving must apply one shared rule, with a test.
  Above the cap, the column is served as a string-kind payload of per-cell codes plus a
  category table, or refused with a clear error. The plan picks one. Either way, a column is
  never described one way and served another.
- No scirpy dependency. A file prepared with scirpy already carries `clone_id` in obs, and
  that is all lumenCell reads. Building clones from receptor contigs is a benchmark-builder
  concern (§6.1), not an app feature.

### 4.2 Keeping clones out of inference

`InferenceContext.obs` currently carries every obs column to the methods
(`engine/context.py`). Phase 1 changes that:

- The columns in the dataset's `ClonalSettings`, plus any column matching the receptor
  patterns (`clone_id*`, `raw_clonotype_id`, `IR_*`, `airr:*`, `receptor_*`,
  `chain_pairing`, `cdr3*`), are removed from the obs that methods receive.
- A root or terminal defined as `cluster:<key>=<value>` on such a column is refused with a
  clear error.
- **Isolation test.** On a synthetic case, add a clone column that is a perfect copy of the
  true lineage. Each method's result must then be identical, value for value, with and without
  that column, and the run id must be unchanged. This is the test that the check is not
  circular.

### 4.3 Clonal scores

New package `src/lumencell/clonal/`, holding pure functions on arrays with no engine, API or
anndata imports:
- `labels.py` turns raw labels and donors into integer clone ids keyed by (donor, label),
  plus a donor index.
- `scores.py` holds the metrics and the null.

The benchmark and the worker both call this package.

**Inputs**
- The run's lineage membership masks M (n cells × L lineages), taken from the result schema
  the methods already return.
- The clone ids κ and donors d from `labels.py`.
- `min_clone_size` m, the number of permutations B, and a seed.

**Primary metric: clone purity**, after STARTRAC's transition index and LRT's lineage-bias
test. For each clone k among the scored cells, take its lineage distribution p_k over the L
lineages and its normalised entropy H_k = −Σ_l p_kl log p_kl / log L. The score is

    P = Σ_k n_k (1 − H_k) / Σ_k n_k        (weighted by cell mass, 0 ≤ P ≤ 1)

**Null.** Shuffle κ among the scored cells within each donor, B = 1,000 times with a seeded
generator, and recompute P each time. Shuffling only among scored cells keeps the lineage mix
of expanded cells fixed. The null therefore tests only whether cells of one clone sit in the
same lineage more often than chance. Expansion bias between lineages is not counted as
agreement. STARTRAC shuffles among all cells, which mixes the two effects.

An optional `strata` argument shuffles within donor × further obs keys, such as time point or
tissue. The primary null is always donor only. A stratified null, when used, is reported
alongside it and never replaces it.

**Reported for every score**
- observed P, null mean μ, null sd σ, excess Δ = P − μ, and z = Δ / σ (null when σ = 0)
- p = (1 + #{P_b ≥ P}) / (B + 1)
- cells scored, clones scored, donors, m, B and the seed

The primary score uses m = 2. The same score at m = 3 is always reported alongside it as a
sensitivity check.

**Secondary metric (reported, not gated): lineage-pair coupling.** Build the clone × lineage
count matrix X over the scored cells. For each lineage pair (a, b), compute the cosine
similarity of columns a and b (CoSpar's clonal coupling). Compare it with the same
within-donor null and report the observed value, null mean, z and p. A high z means the two
lineages share more clones than chance, which says they are related, not which came first.

**Statuses.** `ok`; `not_configured` (no clonal setting); `single_lineage` (L < 2, where
purity is undefined); `too_few_clones` (fewer than `min_clones` = 20 clones scored; the values
are still reported but never count as a pass); `failed` (an error, with its message).

### 4.4 Synthetic clones

The generator gains an optional `clones` spec. With `clones=None`, its output is unchanged
byte for byte. A test pins the hash of X, obs and uns for a fixed case against the current
generator.

- Clones are assigned to the cells after the counts are drawn, using a separate generator,
  `np.random.default_rng([seed, CLONE_STREAM])`. The existing draw order and
  `benchmarks/baseline.json` are untouched.
- **Donors.** Each cell gets one of D donors (default 3), at random.
- **Founders.** Within each donor, clones are founded until a target share of cells is labelled
  (default 0.85, close to 10x TCR capture rates). With probability `fate_bias`, a founder sits
  at a random point on a segment after the first fork. Otherwise it sits on the trunk. The
  clone's size comes from a heavy-tailed distribution with a large singleton share, and its
  members are drawn from unlabelled cells of the same donor that descend from the founder
  point (`true_segment` plus depth, via `_ancestors`).
  - `fate_bias = 0` means no clonal signal: trunk-founded clones spread across lineages in
    proportion to the lineages' sizes, as the null does.
  - `fate_bias = 1` means every clone is confined to one fate.
- **Collisions.** Clone label strings restart in each donor (`c0`, `c1`, …), so every string
  occurs in several donors. A scorer that ignores donors gives a measurably different answer,
  and a test checks this.
- **Random control.** In `mode = "random"`, labels are assigned with no regard to the tree.
  This is the exact null case and is used to check calibration.
- Clones are written to obs as `clone_id` and `donor`, and the spec goes to
  `uns["lumencell_truth"]["clones"]`.
- `Truth` gains optional `clone_ids` and `donors` fields, following how `resolvable` and
  `terminals` were added.

### 4.5 Benchmark metric and gate

- The **oracle** score is clone purity computed on the true lineage masks.
- For each method, `clonal_excess_ratio = Δ_method / Δ_oracle`. It is reported with `clonal_p`
  and with `clonal_scored_fraction`, the method's scored cells divided by the oracle's.
- Every branching case in the synthetic matrix gets a clonal layer with `fate_bias = 0.7`. The
  counts are unchanged, so every existing metric is unchanged.
- **Gates**, written into `thresholds.toml` with their reasons before the first measuring run:
  - `clonal_excess_ratio ≥ 0.5`. A method whose lineages capture less than half of the clonal
    signal present in the true lineages is not separating the fates the clones show.
  - `clonal_p ≤ 0.01`. With B = 1,000, the smallest attainable p is about 0.001.
- `clonal_excess_ratio` joins `HIGHER_IS_BETTER` for the regression guard. The baseline is
  re-measured once, on purpose, in its own commit.
- **A clonal gate never passes by default.** `check_gates` skips NaN today. For the clonal
  gates, a NaN caused by a structural condition (a linear topology, L < 2 in the truth) counts
  as not applicable. A NaN caused by missing evidence (`too_few_clones`, `failed`, or an
  oracle Δ that is not significant) counts as a failure, and the reason is recorded. Existing
  gates keep their current NaN handling.
- The `experimental` flags cannot change: all three methods already fail at least one
  synthetic gate, and a new gate only adds failures. The golden test still runs.

**Tests of the metric itself**, in CI unless marked slow:
- Hand-computed examples.
- Hypothesis properties:
  - P is unchanged when clones, lineages or cells are relabelled or reordered.
  - Singleton clones never change any score.
  - 0 ≤ P ≤ 1.
  - Splitting one donor's clones across donors never merges clones.
  - Ignoring donors changes the score when label strings collide.
- Calibration (slow): on `mode = "random"` across 200 seeds, the p-values are consistent with a
  uniform distribution (Kolmogorov–Smirnov at α = 0.01), and the share of p ≤ 0.05 lies within
  a binomial 99% interval around 0.05.
- Monotonicity: the oracle's Δ rises with `fate_bias` over {0, 0.35, 0.7, 1}, and at
  `fate_bias = 0` it is within the null's spread.

### 4.6 API: clonal agreement for each run

New models in `api/models.py`:

- `ClonalScore`, with the fields listed under "Reported for every score" (§4.3).
- `LineageCoupling { a, b, observed, null_mean, z, p }`.
- `ClonalAgreement { status, settings, purity: ClonalScore | None,
  purity_min3: ClonalScore | None, coupling: list[LineageCoupling], notes: list[str] }`.
- `RunTrust { clonal: ClonalAgreement | None }`, which phase 2 extends.
- `TrajectoryRun.trust: RunTrust | None`.

Because `json_schema_serialization_defaults_required` is set, the new key is required in the
generated TypeScript. The only web change is therefore `make contract` plus adding
`trust: null` to the `run` fixture in `web/src/test/fixtures.ts`.

**Where the score is computed.** It runs in the worker, never in the server.
- At the end of an `infer` job, if clonal settings exist. The settings travel in the payload
  but are not hashed into the run id.
- In a new `score_clonal` job, queued for every complete run of the dataset when the settings
  change.

Results are written atomically to `<run>/trust/clonal-<h>.json`, where `<h>` is a SHA-256
prefix of the settings, the metric version and the parameters. Column names never appear in a
path. `Runs.get_run` serves the file matching the current settings, and reports `pending`
while a job for it is queued.

The CLI command `lumencell run` gains `--clone-key`, `--donor-key` and `--single-donor`, which
print the clonal block.

### 4.7 Phase 1 done when

- `make ci` is green, and `make bench-synth` shows no regressions against the deliberately
  re-measured baseline.
- The isolation test, the property tests and the calibration test pass. The calibration test
  is slow and is run once and recorded.
- A run on a synthetic file with clones returns a `trust.clonal` block through the API and the
  CLI.
- `thresholds.toml` holds the clonal definitions, gates and reasons, dated before the first
  measuring run.
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

- New prep stage `cell_cycle`, with input `lognorm`.
  - Outputs: `cc_s`, `cc_g2m` and `cc_phase`.
  - It is declared and cached like the other stages. Its params are `PrepParams` fields:
    `cc_species ∈ {auto, human, mouse}` and `cc_gene_set = "tirosh2016"`.
- **Genes.** The S and G2M lists are from Tirosh et al. 2016, shipped as package data with
  their source. Mouse genes map through a fixed human-to-mouse ortholog table, shipped with its
  source and date. Both are verified during implementation.
- **Species.** `auto` resolves from the case of the symbols that match the list. The resolved
  species and the number of list genes found go into the manifest. With fewer than 20 genes
  found, the stage records "not scored" and its reason instead of failing.
- **Scoring** uses scanpy's `score_genes_cell_cycle` with a fixed seed. Scanpy is already
  used in prep.
- **Report-only.** No stage reads its outputs and nothing is regressed out. Adding the stage
  changes `prep_key` and therefore every run id, once. Existing cached runs become orphans and
  are recomputed on demand. This is accepted and recorded in the changelog.
- **Per-run check.** For each lineage, compare the share of its single-lineage cells in S or
  G2M with the dataset's share. If a lineage's share is at least twice the dataset's share and
  at least 0.5, the run gets a warning: "Lineage X is mostly cycling cells; it may be a
  cell-cycle state rather than a fate." The thresholds go into `thresholds.toml`. Proliferating
  effector cells are real biology, so this is a hint, not a verdict.

### 5.2 Stability

- **Automatic** (Shreyan's choice). Every completed `infer` job queues a `stability` job for
  that run. The run's result is served at once, the stability block reads `pending` until the
  job finishes, and other jobs queue behind it. The job can be cancelled like any other.
- **Replicates.** R = 5 per run, all seeded:
  - 2 **reseeds**: Leiden and the method are re-run with new seeds on the full data.
  - 3 **subsamples**: 80% of the cells, uniformly at random, always keeping the root cells. For
    each subsample, the kNN graph, Leiden, the diffusion maps and PAGA are recomputed. The
    full-data normalisation, HVG and PCA are reused. The method is then run with the same root
    and terminals.
  - Replicate prep lives in the job's temporary directory and is not added to the stage cache.
- **Comparison with the reference run, on shared cells.**
  - Pseudotime Spearman.
  - Lineages, matched with the existing `match_lineages` (Hungarian assignment on Jaccard,
    `MIN_MATCH_JACCARD = 0.6`): the share of reference lineages matched and the mean Jaccard of
    the matches.
  - The number of lineages in each replicate.
- **Reported** separately for reseeds and subsamples: the median and minimum of each value,
  the per-replicate rows, R, the fraction and the seeds.
- **Stored** atomically in `<run>/trust/stability-<h>.json`. The `stability` job kind gets an
  explicit time limit in `jobs/manager.py`.
- The benchmark computes stability for every method on the real datasets, and on one seed per
  topology in the synthetic matrix, which keeps `make bench-synth` time bounded.

### 5.3 Tree-support verdict

A plain signal with three levels: `supported`, `weak` or `not_supported`, plus `provisional`
while stability is pending. Each comes with its reasons in plain words. The checks are
declared in `thresholds.toml` before they are measured.

| Check | Applies | Rule |
|---|---|---|
| C1 connectivity | always | coverage (share of cells with finite pseudotime) ≥ 0.95 |
| C2 tree-likeness | always (dataset level) | weight on PAGA's maximum spanning tree ÷ total PAGA edge weight ≥ τ_tree |
| C3 lineage stability | once stability is done | median share of reference lineages matched ≥ τ_lin |
| C4 pseudotime stability | once stability is done | median Spearman ≥ τ_pt |
| C5 clonal support | when clones are configured and L ≥ 2 | clonal `ok` with p ≤ 0.01 |

The verdict is `not_supported` if C1 or C2 fails, `weak` if any of C3–C5 fails, and
`supported` otherwise.

- **Where the thresholds come from.** τ_tree, τ_lin and τ_pt are set from synthetic data only,
  by a rule fixed in advance: the midpoint between the lowest value on the 12-case matrix and
  the highest value on the negative controls. If the two ranges overlap, the check is reported
  as non-discriminating and left out of the verdict. It is not tuned.
- **Negative controls.** New synthetic topologies `blob` (no trajectory: one population with no
  expression change) and `two_groups` (two unconnected populations). They live in a
  `[synthetic.controls]` section, outside the gated matrix, so they can never affect the
  `experimental` flags. The test expectation is that no method's verdict on a control is
  `supported`.

### 5.4 API

`RunTrust` gains `stability: Stability | None`, `tree_support: TreeSupport | None` and
`cell_cycle: CellCycleCheck | None`. The models are defined in the phase 2 plan. The web change
is again only `make contract` and the smallest fixture update.

### 5.5 Phase 2 done when

- `make ci` is green, with no `bench-synth` regressions.
- The cell-cycle stage is cached, and its cache behaviour is tested as in
  `tests/unit/test_prep_cache.py`.
- Stability is queued automatically, served, and cancellable.
- The verdict thresholds are derived by the pre-declared rule and recorded.
- No control case is verdicted `supported`.
- `/security-review` is run.

## 6. Phase 3: Real T-cell data and the case study

### 6.1 Benchmark: Daniel et al. 2022 (GSE188670, subseries GSE188666)

Mouse LCMV, acute (Armstrong) and chronic (Clone 13) infection, days 5, 8 and 21, spleen,
liver and lung; 10x 5′ with TCR. It is open access.

1. **Fetch.** With Shreyan's OK at the time, download the GSE188666 expression files and the
   24 per-sample `tcr_filtered_contig_annotations.csv.gz` files, about 1.1 GB. The ATAC tar is
   not fetched.
   - `RealDataset` gains multi-file sources, a `licence` field, and a build step that writes
     one derived `.h5ad` into the cache. The downloads are never modified.
   - The SHA-256 of every file is recorded in the registry.
2. **Clones.** The builder keeps productive, high-confidence, full-length TRA and TRB contigs
   for each barcode. A clone is the sorted set of their CDR3 nucleotide sequences, namespaced
   by GEO capture (GSM).
   - The capture is used because the paper doesn't identify mice. This never merges two mice,
     but it hides sharing between tissues of one mouse, so cross-tissue coupling is
     understated.
   - As a sanity check, agreement with 10x `raw_clonotype_id` is reported.
3. **Inspect, then declare.** Read the metadata without running any method. Then, in one PR
   and before any method runs:
   - declare the subset (infection arm, days, tissues), with reasons
   - declare the truth from the paper's state labels: lineages as state lists, with the root at
     progenitor-exhausted (Tpex) cells as described in the paper and the earlier literature (Im
     et al. 2016, Utzschneider et al. 2016); the exact label names come from the metadata
   - declare the gates in `[real.daniel2022.gates]`, reusing the existing real-data gate values
     plus the two clonal gates, where the truth lineages play the oracle role
   - If the GEO metadata has no state labels, fall back to day and infection arm as the only
     order truth, with no lineage truth. That fallback is recorded as such.
4. **Run** all three methods with `make bench-real`, and regenerate the README tables. Day
   after infection is reported as a second order check (pseudotime Spearman against day), not
   gated, because progenitor and terminal states coexist on day 21.

### 6.2 Case study: Sheih et al. 2020 (GSE125881)

Human CD8 CAR T cells from 4 patients at 4 time points. Licence CC BY 4.0.

1. **Verify first.** With Shreyan's OK, download one clonotype table (about 7 MB) and check
   that it carries per-cell barcodes that join to the expression matrix. If it does not,
   switch to Cao et al. 2025 (GSE207192). The GEO page says 3 patients and the paper says 7;
   that discrepancy is resolved before designing around it.
2. **Analyse.** Each patient is a donor. Run all three methods with an explicit root taken
   from the paper's description, and report:
   - clonal purity, with its m = 3 sensitivity check, against both the donor-only null and a
     donor × time-point null, so that clones expanding at one time point are not read as
     lineage agreement
   - lineage-pair coupling
   - stability
   - the tree verdict
   - the cell-cycle check
3. **Write it up.** The numbers come from a script committed under `benchmarks/case_studies/`,
   so they are reproducible. The write-up goes to `docs/case-studies/sheih-2020.md`.
   - The biology comes from the paper and is cited. Anything unverified is marked.
   - No direction is claimed from clones.
   - If the data doesn't fit a tree, that is the finding.
4. **Check the statistics.** `data:statistical-analysis` and `data:validate-data` are run on
   the null-model test and the case-study numbers.

### 6.3 Phase 3 done when

- Daniel 2022 is in the registry with its SHA-256 values, truth, subset and gates, all
  committed before its first run.
- `make bench-real` has produced a results file, and the README tables are regenerated. Gate
  failures are reported as they are.
- The case study is written from reproducible numbers, with citations.
- `/security-review` is run.

## 7. Interfaces and data flow

```
.h5ad ──loader──▶ Dataset (obs incl. clone/donor columns; matrix dropped in server)
          │
          ├─ dataset_view: lists columns (one cardinality rule); validates ClonalSettings ─▶ state file
          │
infer job (worker): prep cache ─▶ InferenceContext (obs minus scoring-only columns) ─▶ method
          ─▶ run_method/validate ─▶ RunStore.save ─▶ clonal.scores (if settings) ─▶ trust/clonal-<h>.json
          └─▶ queue stability job ─▶ trust/stability-<h>.json ─▶ tree verdict
GET run (server): result.json + trust/*.json ─▶ TrajectoryRun.trust
bench: synth (+clones, separate RNG) ─▶ Truth (+clone_ids, donors) ─▶ methods ─▶ Scores (+clonal) ─▶ gates
```

Each unit has one purpose:
- `clonal/labels.py`: labels to clone ids.
- `clonal/scores.py`: the metrics and the null.
- `bench/synth_clones.py`: clone simulation.
- `engine/stability.py`: replicates and comparison.
- `engine/tree_support.py`: the verdict.
- `prep/stages.py`: one new stage.
- `bench/builders/daniel2022.py`: contigs to the derived `.h5ad`.

## 8. Decisions and what they cost if wrong

| # | Decision | Cost if wrong |
|---|---|---|
| 1 | Clone identity is (donor, label). A donor column, or an explicit single-donor flag, is required. | One extra setting for single-donor files. |
| 2 | Clonal settings are dataset-level state, not part of the run id. | Settings would have to be repeated if a user wants different clone columns per run. |
| 3 | Scoring-only columns are removed from `InferenceContext.obs`, and receptor-pattern columns too. | A legitimate column matching a pattern can't name lineages or define a root. |
| 4 | No scirpy dependency. Clone labels are read from obs, and the contig parser is benchmark-only. | AIRR-format input needs a later optional extra. |
| 5 | Clone purity over single-lineage cells is the primary metric. Shared-trunk cells are excluded. | It discards cells on intermediate segments in multi-fork trees. A changed definition must be logged, never swapped silently. |
| 6 | The null shuffles among scored cells within each donor. | It ignores confounding by time point or tissue within a donor. For Sheih, a time-stratified null is reported alongside it. |
| 7 | m = 2 primary, m = 3 always reported. | None beyond reading two numbers. |
| 8 | The synthetic gate is an excess ratio ≥ 0.5 against the oracle, with p ≤ 0.01. | If 0.5 proves strict, methods fail it and stay experimental. The bar is not lowered. |
| 9 | Synthetic clones are drawn after the counts on a separate RNG stream. | None to existing results. |
| 10 | Cell cycle is report-only, with a warning. | Users wanting regression wait for later work. The `prep_key` changes once. |
| 11 | Stability runs automatically: R = 5, with 2 reseeds and 3 subsamples at 80%, and the subsample prep recomputed from kNN onward. | Each run occupies the worker about 5 times longer. PCA reuse slightly understates subsampling variance. |
| 12 | Verdict thresholds come from synthetic data by a fixed midpoint rule, with negative controls outside the gated matrix. | A check that doesn't discriminate drops out, which is reported, not hidden. |
| 13 | Daniel clones are namespaced per capture. | Cross-tissue sharing within a mouse is invisible. |
| 14 | The Daniel truth uses the paper's state labels, with the root at Tpex. Day is a reported check. | The labels come from expression, like the other three benchmarks. The clones remain the independent check. |
| 15 | Sheih is verified before use, with Cao 2025 as the fallback. | One small download. |
| 16 | The missing `need_from_me` section of the brief is taken to mean dataset downloads and new runtime dependencies. | Something else Shreyan meant to approve gets decided without them. |

## 9. Not in this work

- web views
- the multi-user version
- new trajectory methods
- RNA velocity
- receptor design
- packaging
- the 100k-cell performance budgets
- cell-cycle regression
- AIRR or scirpy import in the app
- similarity-based clonotype clusters (only exact clones are used)
- any direction inferred from clones
- per-clone pseudotime coherence

## 10. References

Every DOI below except Tirosh 2016, Im 2016 and Utzschneider 2016 was checked by the
researchers (see the notes). Those three are checked during implementation.

- Zhang L. et al. STARTRAC. Nature 2018. doi:10.1038/s41586-018-0694-x
- Xie J. et al. LRT. PLoS Comput Biol 2023. doi:10.1371/journal.pcbi.1011300
- Wang S.-W. et al. CoSpar. Nat Biotechnol 2022. doi:10.1038/s41587-022-01209-1
- Weinreb C. et al. LARRY. Science 2020. doi:10.1126/science.aaw3381
- Weinreb C., Klein A.M. PNAS 2020. doi:10.1073/pnas.2000238117
- Daniel B. et al. Nat Immunol 2022. doi:10.1038/s41590-022-01337-5
- Sheih A. et al. Nat Commun 2020. doi:10.1038/s41467-019-13880-1
- Cao et al. Nat Commun 2025. doi:10.1038/s41467-025-59298-w
- Tirosh I. et al. Science 2016. doi:10.1126/science.aad0501
- Im S.J. et al. Nature 2016. doi:10.1038/nature19330
- Utzschneider D.T. et al. Immunity 2016. doi:10.1016/j.immuni.2016.07.021
