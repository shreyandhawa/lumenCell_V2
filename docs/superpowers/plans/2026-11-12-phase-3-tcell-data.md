# Phase 3: Real T-cell Data and the Case Study — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Benchmark the three methods on Daniel 2022 (chronic LCMV, paired TCR), with the subset, namespace, clone definition, truth and gates frozen from metadata before any clone is counted or any method runs. Then publish a reproducible Sheih 2020 CAR-T case study, scored against clones, the three nulls, the baseline, stability, the tree verdict and the cell-cycle check.

**Architecture:** The real-dataset registry gains multi-file sources, licences and a build step. The build step turns verified downloads into one derived `.h5ad` in the cache, named by a hash of its inputs. Per-dataset builders under `bench/builders/` join expression to TCR contigs and call clones through one shared, pure function. Declarations merge as data (registry entries plus `thresholds.toml`) before anything is measured. The real runner and a case-study module then call the phase 1 and phase 2 code unchanged.

**Tech Stack:** Python 3.13 via uv; numpy, scipy, pandas, anndata, scanpy; stdlib `tarfile`, `gzip`, and `urllib` with certifi; pytest with typer's `CliRunner`. No new dependency.

**Spec:** `docs/superpowers/specs/2026-11-12-tcell-clonotypes-design.md`. Read §2 (constraints), §3 (terms), §6.1–6.3 (this phase), §8 decisions 16–19 and §4.3 (what every clonal score reports). The spec wins over this plan on any conflict. Report the conflict instead of choosing.

**Depends on:** phases 1 and 2, merged. Phase 3 consumes their interfaces by name and never redefines them:
- **Phase 1:** see its plan's Interfaces, `docs/superpowers/plans/2026-11-12-phase-1-clonotypes.md`. Phase 3 uses `clonal/` (`encode_clones`, `score_clonal`, `isolation.RECEPTOR_COLUMN`, `receptor_gene_mask`, `scoring_only_columns`), `InferenceContext.from_prep(..., clonal_columns=…)`, `ClonalRow`, `ClonalRowStatus`, `CaseResult.clonal`, `clonal_row`, the `clonal_*` row fields (including the reported-only `clonal_n_cells`, `clonal_effective_clones`, `clonal_eta_m3` and `clonal_p_m3`) and their gate policy in `check_gates`, `trust.constants` (`ORACLE_ETA_MIN`, `ORACLE_P_MAX`), `api.models` (`ClonalSettings`, `ClonalAgreement`), and `trust.clonal` (`check_settings`, `settings_columns`, `agreement_for_run`).
- **Phase 2:** see its plan's Interfaces, `docs/superpowers/plans/2026-11-12-phase-2-robustness.md`. Task numbers in brackets are phase 2's.
  - `engine/stability.py`: `run_stability(ds, prep, params, reference, spec, config=StabilityConfig(), progress=…) -> StabilityResult`, `RunSpec(method, params, terminals, label_key, clonal_columns=())` and `StabilityError` [T6].
  - `engine/tree_support.py`: `tree_likeness(paga_connectivities, paga_tree)`, `coverage(pseudotime)`, `verdict(...)`, `StabilityEvidence` and `calibrated_thresholds()` [T7, T8].
  - `trust/stability.py`: `to_api_stability`, and `not_measured(status, note)` [T10].
  - `trust/verdict.py`: `stability_evidence`, `clonal_evidence` and `to_api_tree_support` [T11].
  - `trust/cell_cycle.py`: `cell_cycle_check(primary_lineage, names, cc_phase, cc_info) -> CellCycleCheck` [T4], on the prep artifacts `cc_phase` and `cc_info` [T3].
  - `bench/runner.py`: `stability_row(data, prep, params, result, label_key, clonal_columns: tuple[str, ...] = ()) -> StabilityRow`, `CaseResult.stability`, and `run_real_case(dataset_id, method, params, stability=False)`. `bench real` runs stability by default [T12].

**Preflight** (orchestrator, once, on `main` before Task 1):
- `git fetch origin && git rev-list --left-right --count main...origin/main` prints `0	0`.
- This command prints `ok`:

```bash
uv run python -c "from lumencell.clonal.labels import encode_clones; from lumencell.clonal.scores import score_clonal; from lumencell.clonal.isolation import RECEPTOR_COLUMN, receptor_gene_mask, scoring_only_columns; from lumencell.bench.runner import ClonalRow, ClonalRowStatus, CaseResult, clonal_row; from lumencell.bench.synth_clones import add_clones, DEFAULT_LAYER; from lumencell.trust.clonal import agreement_for_run, settings_columns, check_settings; from lumencell.trust.constants import ORACLE_ETA_MIN, ORACLE_P_MAX; from lumencell.api.models import ClonalSettings, ClonalAgreement; from lumencell.engine.stability import run_stability, RunSpec, StabilityError; from lumencell.engine.tree_support import tree_likeness, coverage, verdict, StabilityEvidence, calibrated_thresholds; from lumencell.trust.stability import to_api_stability, not_measured; from lumencell.trust.verdict import stability_evidence, clonal_evidence, to_api_tree_support; from lumencell.trust.cell_cycle import cell_cycle_check; from lumencell.bench.runner import stability_row; import inspect; assert 'clonal_columns' in inspect.signature(stability_row).parameters; print('ok')"
```

If an import fails, stop and report. Phase 3 does not patch phases 1 or 2.

## Global Constraints

Every task's requirements include these (spec §2, plus the phase 3 rules of §6).
- **Scoring only.** Clone labels are used only for scoring. They never reach inference through obs (a prep stage, a method, a root, a terminal, `cluster_key`, `label_key`), and never through expression: genes matching `^TR[ABDG]([VDJ]\d|C\d*$)`, case-insensitive, are excluded from HVG selection. Every clone, receptor or namespace column of a built dataset either matches `RECEPTOR_COLUMN` or is passed as `clonal_columns`.
- **Strata stay visible.** `strata_keys` columns (Sheih's time point) are not scoring-only. They stay visible to methods and may define a root (spec §4.2, amended at plan time).
- **Defined before measured.** Every metric, null and baseline is defined in `benchmarks/thresholds.toml` before any result is seen. Every score is reported against its nulls.
- **Clones are keyed by donor.** Here the donor is the clone namespace, and the same receptor string in two namespaces is two clones.
- **Clone size and counts.** m = 2 is primary and m = 3 is the sensitivity check. Every score reports cells, clones and the effective clone count. Single-cell clones never contribute.
- **No direction from clones.** Clone sharing shows that two states are related, never which came first. No output or text claims direction from clones.
- **No tuning on T-cell data.** No method default and no threshold is tuned on T-cell data. Daniel 2022 and the case study are for evaluation only.
- **Thresholds.** They are never lowered. New gates and definitions go into `thresholds.toml`, with dated reasons, before the first run that measures them.
- **Daniel freeze.** No clone × state, clone × capture or clone × lineage count is computed before Task 6 merges. The oracle needs no method, so the freeze comes first.
- **Case-study freeze.** The declaration (Task 10) merges before the analysis (Task 12). It is written from the paper and the file headers only.
- **Case-study writing.** The biology comes from the dataset's paper and is cited. Anything unverified is marked. If the data doesn't fit a tree, that is the result.
- **Case-study unit.** Each method runs once per patient (4 × 3 runs, stability on each), because lumenCell has no batch correction. The clone namespace is the patient. A pooled run is only a labelled diagnostic, never the headline.
- **Existing decisions.**
  - Roots are always explicit, and inference never uses the UMAP.
  - Neither the input file nor any download is ever modified.
  - A method's `experimental` flag matches its synthetic gates. Real-data results never change it.
  - `engine/` never imports `bench/`.
  - Heavy work never runs in the server. The case study runs in-process, from a script.
- **No API change is planned.** If one becomes necessary, it starts in `src/lumencell/api/models.py` and goes through `make contract`. `web/` then changes only by contract output and the `dataset` and `run` fixtures in `web/src/test/fixtures.ts`.
- **Atomic writes.** Run, cache, derived-dataset and results files are written by atomic rename.
- **Dependencies.** No new runtime dependency. scikit-learn is never imported directly, so the adjusted Rand index is written with numpy.
- **Thresholds at run time.** Phase 3 adds no runtime threshold. It reads `ORACLE_ETA_MIN = 0.05` and `ORACLE_P_MAX = 0.01` from `trust/constants.py`.
- **Downloads.**
  - Datasets are never committed.
  - Every download stops for the owner's OK first, checks size and SHA-256, keeps certificate checks on (certifi) and renames atomically.
  - Tests needing real data are marked `real_data`; slow ones are marked `slow`. CI uses synthetic data and small fixtures only.
- **Untrusted values.**
  - Values read from a user's or a downloaded file (column names, labels, barcodes, tar member names, CSV headers) never reach a file path or a shell command.
  - Derived files use hashed names.
  - Tar archives are read in memory, never extracted.
  - Each download goes into its own directory. Scripts that read downloads run as `uv run python -I <script> <paths>`, from the repo root.
- **Process.**
  - Work test-first, with one branch and one PR per task, named `pivot/NN-<slug>`.
  - Squash-merge only on green CI, after review by a different agent. Nothing is pushed to `main` directly.
  - No `Co-Authored-By` trailers and no "Generated with" footers.
  - Gate every commit on `make ci; echo $?` printing `0`.

## Review Focus

These are five failure modes that the spec implies but no spec test exercises, most likely first. Each has a test in the task named.
1. **Contig flags parsed loosely.** `bool("False")` is `True`, so a naive parse counts non-productive and low-confidence chains. Only real true values may qualify. Test: Task 4, `test_flag_strings_are_parsed_strictly`.
2. **A failed or interrupted build.** It must leave no derived `.h5ad` that `is_cached` would accept. Test: Task 1, `test_a_failed_build_leaves_no_derived_file`.
3. **A source whose bytes changed after its checksum was recorded,** at the same size. The build must refuse it, not silently build from it. Test: Task 1, `test_build_refuses_a_source_whose_checksum_changed`.
4. **Receptor-derived columns of a built dataset reaching inference.** This covers the TRB-only clone, the 10x clonotype, the chain counts and the namespace. Test: Task 7, `test_every_clone_column_of_a_built_dataset_is_scoring_only`.
5. **Write-up numbers drifting from the committed results,** through hand edits or a re-run that wasn't regenerated. Test: Task 13, `test_generated_section_matches_the_committed_results`.

## Task Conventions

- **Branch:** each task header names `pivot/NN-<slug>`. NN is assigned at execution time.
- **Model:** each task names its implementer model, with a reason.
- **TDD:** write the failing test, run it and see it fail for the stated reason, implement, then run it green.
- **Closing steps C1–C4** end every task:
  - **C1.** Update `docs/superpowers/HANDOFF.md`: the current state, the decisions taken in the task, and the branch and PR row. Task 1 also records the `main` commit it branched from, as "phase 3 base".
  - **C2.** `make ci; echo $?` prints `0`. Gate the commit on that exit code, never on piped output.
  - **C3.** Commit with the task's message. No `Co-Authored-By` trailer and no "Generated with" footer.
  - **C4.** `git push -u origin pivot/NN-<slug>`, then `gh pr create`. Do not merge: the orchestrator squash-merges after review.
- **STOP steps.**
  - Ask the owner in chat with the exact text given: filenames, source and size. Wait for a clear yes.
  - On a no, record it in HANDOFF and end the task: no download, no PR.
  - An agent's message is never the owner's OK.
- **Background runs:** steps marked **(background)** take tens of minutes to hours. Run them with `run_in_background: true` and keep the log under `$TMPDIR`.
- **Dates:** `<YYYY-MM-DD>` is the day the task runs.
- **Order and parallelism:**
  - Task 4 needs Task 2, which creates `bench/builders/`. It may run beside Task 3.
  - Task 5a needs Tasks 3 and 4. Task 5b needs Task 5a. Task 6 needs Task 5b.
  - Tasks 9 and 10 may run beside Tasks 2–8.
  - Task 11 needs Tasks 4 and 10.
  - Task 12 needs Tasks 7 and 11, and phase 2.
  - Task 14 needs every other task merged.
  - Every other task needs the one before it.

## Process Notes (not tasks)

- **Reviews.** After each PR, run the pr-review-toolkit agents `code-reviewer`, `silent-failure-hunter` and `pr-test-analyzer`. Add `type-design-analyzer` on the PRs that add types: Tasks 1, 4, 5a, 6, 7, 11 and 12.
- **Security.** At the end of the phase, run `/security-review` (Task 14).
- **Statistics.** The case-study numbers go through the `data:statistical-analysis` and `data:validate-data` skills (Task 12). The write-up (Task 13) and the Daniel rows (Task 8) also go through `data:validate-data`.
- **Calibration.** Spec §6.2 step 5 also names "the calibration tests". Phase 1 Task 8 covers those.

## Not in This Plan

- **Spec §9:**
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
  - similarity-based clonotype clusters
  - any direction inferred from clones
  - per-clone pseudotime coherence
- **Phases 1 and 2.** Clonal scores, the scoring job, the API, stability, the verdict and the cell-cycle stage are consumed unchanged. A missing or different interface is reported, not redefined here.
- **Daniel 2022:**
  - the acute arm
  - days and tissues that Task 6's rules exclude
  - the scATAC data (GSE188669, and the 14.7 GB `GSE188670_RAW.tar`)
  - the paper's pseudotime, velocity or ArchR trajectory, and its clonal-bias calls, as truth
  - C5 in the verdict of benchmark rows: phase 2's `stability_row` sets it to `not_configured`, and the benchmark's clonal score is in oracle-cells mode, which C5 does not read
- **Sheih 2020:**
  - re-running Cell Ranger on SRA reads if the tables lack barcodes (the fallback is Cao 2025)
  - gates of any kind, because the case study has no truth
  - a pooled analysis as the headline: the case study runs per patient, and a pooled run is only a labelled diagnostic
  - batch correction between patients
- **Other datasets.** Cao 2025 comes in only if Task 9 fails Sheih. GSE206230 and the other alternatives in the dataset notes are out.
- **Existing benchmarks.** No existing threshold or developmental benchmark changes.

## File Map

| File | Responsibility | Task |
|---|---|---|
| `src/lumencell/bench/download.py` | sources, builds, fetching, checksums, derived paths | 1 |
| `src/lumencell/bench/datasets.py` | registry: `RealDataset` (sources, licence, build, clonal, report time), `DATASETS["daniel2022"]`, `CaseStudyDataset`, `CASE_STUDIES`, `truth_for` | 1, 6, 7, 11 |
| `src/lumencell/cli.py` | `bench fetch` lists sources, asks once, builds, and knows case studies; `bench real` prints clonal status | 1, 7, 11 |
| `src/lumencell/bench/builders/__init__.py` | package docstring | 2 |
| `src/lumencell/bench/builders/common.py` | strict flags, in-memory tar members, the shared clone call | 4 |
| `src/lumencell/bench/builders/daniel2022.py` | Daniel sources, licence, capture table, label column, builder | 2, 5a, 5b |
| `src/lumencell/bench/builders/sheih2020.py` (or `cao2025.py`) | case-study sources, sample table, builder | 11 |
| `benchmarks/inspection/daniel2022_metadata.py` | metadata-only inspection that withholds receptor values | 3 |
| `benchmarks/inspection/daniel2022_formats.py` | file formats and the barcode join check, withholding receptor values | 5a |
| `benchmarks/inspection/barcode_join.py` | the case-study barcode join check | 9 |
| `benchmarks/thresholds.toml` | `[real.daniel2022.*]`, `[case_study.<id>]`, changelog | 6, 8, 10 |
| `src/lumencell/bench/real_clonal.py` | clonal row and reported fields for real datasets; ARI; recurrence | 7 |
| `src/lumencell/bench/runner.py` | `CaseResult.reported`; `clonal_reason`; `run_real_case` scores clones and passes them to `stability_row` | 7 |
| `src/lumencell/bench/report.py` | clonal status and the "Reported, not gated" table; `replace_section` | 7, 12 |
| `benchmarks/results/real-<YYYY-MM-DD>.json` | first results file with Daniel rows | 8 |
| `README.md` | regenerated benchmark tables | 8 |
| `src/lumencell/bench/case_study.py` | the case-study run, its strict JSON report and its markdown tables | 12 |
| `benchmarks/case_studies/sheih2020.py` | thin script that writes the case-study numbers | 12 |
| `benchmarks/case_studies/results/sheih2020-<YYYY-MM-DD>.json` | committed case-study numbers | 12 |
| `docs/case-studies/sheih-2020.md` (or `cao-2025.md`) | the write-up | 13 |
| `docs/superpowers/research/<YYYY-MM-DD>-daniel2022-metadata.md` | records, labels, gates, identifiability, then formats and join | 3, 5a |
| `docs/superpowers/research/<YYYY-MM-DD>-daniel2022-declaration.md` | literature checks and each declaration rule's outcome | 6 |
| `docs/superpowers/research/<YYYY-MM-DD>-sheih-verification.md` | join verdict, GSM table, formats (or the Cao resolution) | 9 |
| `docs/superpowers/HANDOFF.md` | state, decisions, PRs, security review | every task |
| `tests/unit/test_datasets.py`, `tests/unit/test_cli.py`, `tests/real/test_real_datasets.py` | registry, fetch, build, declarations, real checksums, the Daniel build check | 1, 6, 8, 11 |
| `tests/unit/test_daniel2022.py`, `tests/unit/daniel_fixture.py` | Daniel sources, capture facts, and the builder on a tiny GEO-layout fixture | 2, 5a, 5b |
| `tests/unit/test_daniel_inspection.py` | the inspection scripts never print receptor values | 3, 5a |
| `tests/unit/test_builders_common.py` | clone call, flags, tar reader | 4 |
| `tests/unit/test_real_clonal.py` | real-data clonal path and reported fields | 7 |
| `tests/unit/test_barcode_join.py` | the join check | 9 |
| `tests/unit/test_case_study_declaration.py` | the case-study declaration | 10 |
| `tests/unit/test_sheih2020.py` | case-study builder | 11 |
| `tests/unit/test_case_study.py`, `tests/unit/test_case_study_writeup.py` | case-study numbers and write-up | 12, 13 |

## Done When

These are spec §6.3, as commands. Run them from the repo root on `main` after Task 14. For the Cao branch, read `cao2025` and `cao-2025` for `sheih2020` and `sheih-2020`.
- `make ci; echo $?` → `0`.
- **Daniel is fully declared:** `uv run pytest tests/unit/test_datasets.py tests/unit/test_daniel2022.py -q; echo $?` → `0`. This covers the registry with SHA-256, subset, namespace, clone definition, truth and gates.
- **The declaration merged before the first run.** This prints `0`:

```bash
decl=$(git log --format=%H --reverse -S '[real.daniel2022.gates]' -- benchmarks/thresholds.toml | head -1)
first=$(git log --format=%H --reverse -S daniel2022 -- benchmarks/results | head -1)
git merge-base --is-ancestor "$decl" "$first" && [ "$decl" != "$first" ]; echo $?
```

- **The results file has Daniel rows:** `uv run python -c "import json,glob; f=sorted(glob.glob('benchmarks/results/real-*.json'))[-1]; m=json.load(open(f))['methods']; print(f, sorted(k for k in m if 'daniel2022' in m[k]))"` → the Task 8 file and `['paga_dpt', 'palantir', 'scfates_ppt']`.
- **The README tables are current:** `uv run lumencell bench report --readme && git diff --exit-code README.md; echo $?` → `0`.
- **A declared not-assessable Daniel result is visible:** if `DATASETS["daniel2022"].clonal.assessable` is False, `grep -c "clonal: not assessable (declared" README.md` → `3`, and `uv run pytest tests/unit/test_real_clonal.py -q -k not_assessable; echo $?` → `0`.
- **The case-study declaration merged before the analysis.** This prints `0`:

```bash
decl=$(git log --format=%H --reverse -S '[case_study.sheih2020]' -- benchmarks/thresholds.toml | head -1)
first=$(git log --format=%H --reverse -- benchmarks/case_studies/results | head -1)
git merge-base --is-ancestor "$decl" "$first" && [ "$decl" != "$first" ]; echo $?
```

- **The write-up is built from reproducible numbers, with citations:** `uv run pytest tests/unit/test_case_study.py tests/unit/test_case_study_writeup.py -q; echo $?` → `0`.
- **The security review ran:** `grep -n "Security review (phase 3)" docs/superpowers/HANDOFF.md` → one line, and the section lists every finding as fixed or recorded.

---

### Task 1: Multi-file sources, licences and a build step

**Branch:** `pivot/NN-registry-sources` · **Model:** sonnet — a contained refactor, with every behaviour pinned by a test below.

**Files:**
- Create: `src/lumencell/bench/download.py`
- Modify: `src/lumencell/bench/datasets.py`, `src/lumencell/cli.py` (`bench fetch`), `tests/unit/test_datasets.py`, `tests/unit/test_cli.py`, `tests/real/test_real_datasets.py`

**Interfaces:**
- **Produces, `bench/download.py`.** Builders import from here, so `datasets.py` can import builders without a cycle.
  - **Moved unchanged** from `datasets.py`: `USER_AGENT`, `FetchError`, `datasets_dir()` and `sha256_of(path)`. `datasets.py` re-exports them.
  - **`SAFE_NAME`** `= re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,200}")`.
  - **`Source`:** `@dataclass(frozen=True) class Source: url: str; filename: str; size: int; sha256: str | None`.
    - `__post_init__` raises `ValueError` unless `url` starts with `https://`, `filename` fully matches `SAFE_NAME`, `size > 0`, and `sha256` is `None` or 64 lowercase hex characters.
  - **`Builder`** `= Callable[[Mapping[str, Path], Path, Mapping[str, Any]], None]`. It receives every source's local path keyed by `Source.filename`, the temporary path to write the `.h5ad` to, and the declared params.
  - **`Build`:** `@dataclass(frozen=True) class Build: builder: Builder; version: int; params: Mapping[str, Any]`. The params are JSON-serialisable.
  - **`Fetchable(Protocol)`:** read-only `id: str`, `sources: tuple[Source, ...]` and `build: Build | None`.
  - **`source_path(src: Source) -> Path`** is `datasets_dir() / src.filename`. **`is_downloaded(src) -> bool`** means the file exists and its size equals `src.size`.
  - **`fetch_source(src, progress=None) -> Path`** is today's `fetch` body, for one file:
    - a certifi SSL context and a 60 s timeout
    - size and checksum checks
    - a temporary file, renamed into place
    - the temporary file removed on any error
    - a file already downloaded is returned untouched
  - **`build_key(item: Fetchable) -> str`** is the first 12 hex characters of SHA-256 over `json.dumps({"id": item.id, "version": b.version, "params": b.params, "sources": [s.sha256 for s in item.sources]}, sort_keys=True)`.
  - **`local_path(item) -> Path`:**
    - with no build, `source_path(item.sources[0])`, so single-file datasets keep their current paths
    - with a build, `datasets_dir() / f"{item.id}-{build_key(item)}.h5ad"`
  - **`is_cached(item) -> bool`:** every source is downloaded and `local_path` exists.
  - **`build_dataset(item) -> Path`:**
    - Raises `FetchError` naming the file for any source that is not downloaded, has no recorded `sha256` (the message says `no recorded checksum`), or whose `sha256_of` differs (the message says `checksum`).
    - Otherwise opens `tempfile.mkstemp(prefix=f".{item.id}-", suffix=".tmp.h5ad", dir=datasets_dir())`, calls the builder, and `os.replace`s the result onto `local_path`.
    - The temporary file is removed on any exception. It never writes to a source.
  - **`fetch(item, progress=None) -> Path`:** `fetch_source` for each source not yet downloaded, then `build_dataset` when there is a build and the derived file is missing. Returns `local_path`.
- **Produces, `bench/datasets.py`.**
  - `RealDataset` fields become, in order: `id, title, citation, sources: tuple[Source, ...], licence: str, label_key, root_label, lineages, gated, build: Build | None = None, time_key = None, notes = "", extra = …`. The old `url`, `filename`, `size` and `sha256` move into `Source`.
  - `__post_init__` requires exactly one source when `build is None`, unique source filenames and a non-empty `licence`.
  - `fetch`, `is_cached`, `local_path` and `source_path` are re-exported.
- **Produces, CLI `bench fetch <id>`.**
  - Prints the title, citation and licence.
  - Lists every source not yet downloaded with its URL and size, plus the total.
  - Asks once; `--yes` skips the question.
  - Downloads, prints each new file's SHA-256, then builds, printing `building <path>` and then `built`.
  - When every source is present and only the derived file is missing, it builds without asking, because nothing is downloaded.

- [ ] **Step 1: Write the failing tests** in `tests/unit/test_datasets.py`:

```python
def test_registry_is_complete() -> None:
    for ds in DATASETS.values():
        assert ds.licence.strip() and (ds.build is not None or len(ds.sources) == 1)
        for s in ds.sources:
            assert s.sha256 and len(s.sha256) == 64 and s.size > 0 and s.url.startswith("https://")


def test_single_file_datasets_keep_their_cache_paths(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    assert local_path(DATASETS["pancreas"]) == tmp_path / "pancreas_endocrinogenesis_day15.h5ad"


@pytest.mark.parametrize(
    "bad", [{"filename": "../x"}, {"filename": "a/b"}, {"url": "http://h/x"}, {"size": 0}, {"sha256": "ABC"}]
)
def test_source_fields_are_checked(bad: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        Source(**({"url": "https://h/x", "filename": "x", "size": 1, "sha256": None} | bad))
```

  A fixture `two_file` sets `LUMENCELL_DATASETS_DIR` to `tmp_path`. It builds a `RealDataset` with two `Source`s, whose payloads are `b"alpha"` and `b"beta" * 10`, served by a fake `urlopen` that is keyed by URL and counts its calls. Its builder writes `ad.AnnData(np.ones((2, 2), np.float32))` with `uns["seen"] = sorted(paths)`, and counts its calls. Tests:
  - `test_fetch_downloads_every_source_then_builds`:
    - `fetch(ds) == local_path(ds)`
    - the derived file loads with `uns["seen"]` equal to the sorted filenames
    - each source's bytes equal its payload
    - no `*.tmp*` file remains in `tmp_path`
  - `test_downloaded_sources_are_not_fetched_again`: a second `fetch` makes no `urlopen` call and no builder call.
  - `test_a_failed_build_leaves_no_derived_file` (Review Focus 2):
    - the builder writes bytes to its path, then raises `RuntimeError`, and `fetch` re-raises it
    - `local_path(ds)` does not exist
    - no `*.tmp.h5ad` remains
    - `is_cached(ds)` is False
  - `test_build_refuses_a_source_whose_checksum_changed` (Review Focus 3):
    - after a good fetch, overwrite one source with other bytes of the same length, and delete the derived file
    - `fetch` raises `FetchError` matching `checksum` and the filename
    - the builder count is unchanged
  - `test_build_needs_recorded_checksums`: a source with `sha256=None` makes `build_dataset` raise `FetchError` matching `no recorded checksum`.
  - `test_the_derived_path_follows_the_build_inputs`:
    - changing `params`, `version` or one source's `sha256` changes `local_path`
    - the same inputs give the same path
    - the name fully matches `rf"{ds.id}-[0-9a-f]{{12}}\.h5ad"`
  - Rewrite `test_fetch_checks_the_checksum_and_cleans_up` and `test_fetch_checks_the_size` on `Source` and `fetch_source`.

  In `tests/unit/test_cli.py`, add `test_bench_fetch_lists_every_source_and_asks_once`:
  - `monkeypatch.setitem(dsmod.DATASETS, "twofile", ds)`
  - `invoke(app, ["bench", "fetch", "twofile"], input="n\n")` exits 0, its output contains both URLs, the licence and `skipped`, and no `urlopen` call is made
  - with `--yes`, the output has two `sha256` lines and `built`

- [ ] **Step 2: Run** `uv run pytest tests/unit/test_datasets.py tests/unit/test_cli.py -q`. Expected: FAIL (no `lumencell.bench.download`).
- [ ] **Step 3: Implement** as in Interfaces, then move the three existing entries onto `sources=(Source(...),)`.
  - Fill each entry's `licence` with what the page hosting the file states:
    - the `theislab/scvelo_notebooks` repository, for pancreas
    - the figshare item behind each `ndownloader/files/<id>` URL, found through the figshare article page or API, never by downloading the file
  - If a page states none, write `"not stated at <url> (checked <YYYY-MM-DD>); cite <citation>"`.
  - Record in HANDOFF which pages were checked.
- [ ] **Step 4: Update** `tests/real/test_real_datasets.py`. The checksum test asserts `sha256_of(source_path(s)) == s.sha256` for every source, and is skipped unless `is_cached(ds)`.
- [ ] **Step 5: Run** `uv run pytest tests/unit/test_datasets.py tests/unit/test_cli.py tests/real -q`. Expected: PASS. Real tests skip or pass.
- [ ] **Step 6: C1–C4.** Commit message: `Let real datasets come from several verified files and a build step, with licences`.

---

### Task 2: Download the Daniel 2022 files

**Branch:** `pivot/NN-daniel-fetch` · **Model:** sonnet — a guarded download and a recorded manifest.

**Decision:** the 24 contig files come as `GSE188666_RAW.tar`: one verified file of about 12 MB, which the builder reads in memory. This avoids 24 separate GSM downloads, and it is not the 14.7 GB `GSE188670_RAW.tar`, which also holds the ATAC fragments.

**Files:**
- Create: `src/lumencell/bench/builders/__init__.py` (docstring: one module per derived dataset; a builder writes only the path it is given), `src/lumencell/bench/builders/daniel2022.py` (this task adds `GEO_SUPPL`, `SOURCES` and `LICENCE` only), `tests/unit/test_daniel2022.py`

**Interfaces:**
- Consumes: `Source`, `fetch_source`, `sha256_of` and `datasets_dir` (Task 1).
- Produces:
  - `GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188666/suppl/"`.
  - `SOURCES: tuple[Source, ...]`, with exactly these nine files in this order:
    1. `GSE188666_scRNA_LCMV.mtx.gz`
    2. `GSE188666_scRNA_LCMV_barcodes.tsv.gz`
    3. `GSE188666_scRNA_LCMV_features.tsv.gz`
    4. `GSE188666_scRNA_LCMV_metadata.tsv.gz`
    5. `GSE188666_scRNA_D5_LCMV.mtx.gz`
    6. `GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz`
    7. `GSE188666_scRNA_D5_LCMV_features.tsv.gz`
    8. `GSE188666_scRNA_D5_LCMV_metadata.txt.gz`
    9. `GSE188666_RAW.tar`
  - `LICENCE: str`.

- [ ] **Step 1: STOP for the owner's OK.** Ask: "Phase 3 needs the Daniel 2022 scRNA and TCR files from GEO GSE188666 (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188666/suppl/). That is about 1.11 GB in nine files, into the lumenCell datasets cache (`~/Library/Caches/lumencell/datasets`, or `$LUMENCELL_DATASETS_DIR`): the two expression matrices (574 MB and 519 MB), their barcodes, features and metadata (about 10 MB together), and GSE188666_RAW.tar (12 MB, the 24 TCR contig files). It does not include the 14.7 GB GSE188670_RAW.tar with the ATAC data. May I download them?"
- [ ] **Step 2: Read the exact sizes,** with no data transferred: `curl -sI "<GEO_SUPPL><name>"` for each file, reading `Content-Length`. If a file answers anything but 200, or the total is outside 1.0–1.2 GB, stop and report.
- [ ] **Step 3: Write the failing test:**

```python
EXPECTED = (
    "GSE188666_scRNA_LCMV.mtx.gz",
    "GSE188666_scRNA_LCMV_barcodes.tsv.gz",
    "GSE188666_scRNA_LCMV_features.tsv.gz",
    "GSE188666_scRNA_LCMV_metadata.tsv.gz",
    "GSE188666_scRNA_D5_LCMV.mtx.gz",
    "GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz",
    "GSE188666_scRNA_D5_LCMV_features.tsv.gz",
    "GSE188666_scRNA_D5_LCMV_metadata.txt.gz",
    "GSE188666_RAW.tar",
)


def test_sources_are_the_nine_geo_files() -> None:
    assert tuple(s.filename for s in SOURCES) == EXPECTED
    assert all(s.url == GEO_SUPPL + s.filename and s.sha256 and len(s.sha256) == 64 for s in SOURCES)
    assert 1.0e9 < sum(s.size for s in SOURCES) < 1.2e9 and not any("GSE188670" in s.url for s in SOURCES)


def test_the_licence_is_recorded() -> None:
    assert "GSE188666" in LICENCE and "10.1038/s41590-022-01337-5" in LICENCE
```

- [ ] **Step 4: Run** `uv run pytest tests/unit/test_daniel2022.py -q`. Expected: FAIL (module missing).
- [ ] **Step 5: Write `SOURCES`** with Step 2's sizes and `sha256=None`.
- [ ] **Step 6 (background): Download.** This must print `0`:

```bash
uv run python -c "from lumencell.bench.builders.daniel2022 import SOURCES; from lumencell.bench.download import fetch_source, sha256_of; [print(s.filename, sha256_of(fetch_source(s))) for s in SOURCES]" > "$TMPDIR/daniel-fetch.log" 2>&1; echo $?
```

- [ ] **Step 7: Check integrity,** by names and CRCs only:
  - `gzip -t` on each `.gz` exits 0.
  - `tar -tf GSE188666_RAW.tar` lists 24 names ending in `_tcr_filtered_contig_annotations.csv.gz`. Save the list for Task 5a.
  - If either check fails, stop and report. Never extract.
- [ ] **Step 8: Fill each `sha256`** from the log.
- [ ] **Step 9: Write `LICENCE`** in this form:
  - "Data: GEO GSE188666 (SuperSeries GSE188670); NCBI terms <verified at URL, or 'UNVERIFIED: the page is behind a CAPTCHA, not bypassed'>."
  - "Paper: Daniel et al. 2022, Nat Immunol 23:1614-1627, doi:10.1038/s41590-022-01337-5; <the licence PMC11225711 states, or 'no licence stated on PMC11225711'>."
  - Check https://www.ncbi.nlm.nih.gov/home/about/policies/ and https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html. Never bypass a CAPTCHA.
- [ ] **Step 10: Run** the test. Expected: PASS.
- [ ] **Step 11: C1–C4.** HANDOFF records the sizes, the SHA-256 values and the tar member list. Commit message: `Record the Daniel 2022 GEO files with their sizes and SHA-256`.

---

### Task 3: Inspect the Daniel sample records and metadata, without touching clones

**Branch:** `pivot/NN-daniel-metadata` · **Model:** opus — the namespace and subset follow from what this task finds, and the freeze depends on what it must not compute.

**Freeze rule (binds this task):**
- Read the GSM sample records and the two metadata files only: no matrix and no contig file.
- Compute no clone × state, clone × capture or clone × lineage count.
- A metadata column is *withheld* if its name matches `RECEPTOR_COLUMN`, or contains, case-insensitively, `clon`, `cdr3`, `tcr`, `tra`, `trb`, `chain`, `vdj` or `_gene`. Withheld columns are reported by name and dtype only. No value is printed or summarised, and none is used in any other computation.

**Files:**
- Create: `benchmarks/inspection/daniel2022_metadata.py`, `tests/unit/test_daniel_inspection.py`, `docs/superpowers/research/<YYYY-MM-DD>-daniel2022-metadata.md`

**Interfaces:**
- Consumes: `RECEPTOR_COLUMN` (phase 1), `SOURCES` (Task 2), `source_path` (Task 1).
- **Produces, the script** `uv run python -I benchmarks/inspection/daniel2022_metadata.py --dir <datasets dir> [--labels COL --samples COL]`. For each metadata file it prints only:
  - the row count and the first three index values
  - the *tag vocabulary* of the index, meaning the text around the 16-nt core `[ACGT]{16}`, with counts
  - each column's name, dtype and withheld flag
  - for each non-withheld column with at most 200 distinct values, its value counts
  - for each non-withheld numeric column, its minimum, median and maximum
  - with `--labels`/`--samples`, the cross-tab of those two non-withheld columns
- **Produces, the research note.** Tasks 5 and 6 read it. Its sections, in order:
  1. **Samples:** one row per GSM, with the GSM, its library type (GEX or VDJ), the matching GSM of the other type, title, arm (as written), day, tissue and sort gate (as written). Add the mice in the capture (in the record's words, or `not stated`) and any mouse identifier shared with other GSMs (or `none`). Each value names the record field it came from.
  2. **Metadata:** each file's columns (name, dtype, withheld flag), its index tag vocabulary, and the column that names the sample, with its values mapped to GSMs.
  3. **State labels:** the label column, or `none found`. Give every label with its cell count by arm and day, and say whether the D5 metadata uses the same vocabulary.
  4. **Sort gates:** each distinct gate, its markers, and the GSMs that use it.
  5. **Per-mouse identifiability:** answer each question yes, no or not stated, with the evidence.
     - Is a mouse recorded for each capture?
     - Does any mouse span several captures?
     - Does any capture pool mice?
  6. **Processing:** do the records say the main and D5 matrices were processed together (one Cell Ranger run, or integrated)?
  7. **Withheld:** every withheld column, and the sentence "No clone × state, clone × capture or clone × lineage count was computed."
  8. **Still unverified.**

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_daniel_inspection.py`).
  - **Fixture:** in `tmp_path`, write both metadata files in the GEO names, 6 cells each. Give them a label column `celltype` with `Texprog_fixture` and `Texterm_fixture`, a sample column `orig.ident`, and the withheld columns `clonotype_id`, `cdr3s_aa` and `TRB_v_gene`, holding `PLANTEDCLONE`, `CASSPLANTEDF` and `TRBV-PLANTED`.
  - **Helper:** `run(d)` calls `subprocess.run([sys.executable, "-I", SCRIPT, "--dir", str(d), "--labels", "celltype", "--samples", "orig.ident"], capture_output=True, text=True, check=True)`.

```python
def test_the_inspection_never_prints_receptor_values(fixture_dir: Path) -> None:
    out = run(fixture_dir).stdout
    assert "PLANTED" not in out
    assert all(c in out for c in ("clonotype_id", "cdr3s_aa", "TRB_v_gene"))  # names are listed
    assert "Texprog_fixture" in out and "orig.ident" in out  # labels and samples are listed


def test_the_inspection_writes_nothing(fixture_dir: Path) -> None:
    before = sorted(fixture_dir.iterdir())
    run(fixture_dir)
    assert sorted(fixture_dir.iterdir()) == before
```

- [ ] **Step 2: Run** `uv run pytest tests/unit/test_daniel_inspection.py -q`. Expected: FAIL (no script).
- [ ] **Step 3: Implement** the script. It reads the files with pandas from the paths it is given and writes nothing.
- [ ] **Step 4: Run** the tests. Expected: PASS.
- [ ] **Step 5: Run the script** on the real files: `--dir "$(uv run python -c 'from lumencell.bench.download import datasets_dir; print(datasets_dir())')"`, writing to `$TMPDIR/daniel-metadata.txt`. Then run it again with the label and sample columns it found.
- [ ] **Step 6: Read the 46 GSM records of GSE188666.** Use the GEO text form `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSM>&targ=self&form=text&view=full`, which is how the research notes read them. These are record pages, not data files.
  - If a CAPTCHA appears, do not bypass it. Ask the owner instead before downloading the series SOFT file, because that is a file download.
  - Extract each record's title, source, characteristics (strain, infection, day, tissue, sort, mouse or pool), library strategy and relations.
- [ ] **Step 7: Write the research note,** section by section, as listed in Interfaces. Do not decide the subset or the namespace; Task 6 does.
- [ ] **Step 8: C1–C4.** Commit message: `Inspect the Daniel 2022 sample records and metadata without computing anything from clones`.

---

### Task 4: Clone calling from 10x contigs

**Branch:** `pivot/NN-contig-clones` · **Model:** sonnet — a pure function, with every edge case written below.

**Files:**
- Create: `src/lumencell/bench/builders/common.py`, `tests/unit/test_builders_common.py`

**Interfaces:**
- **Produces, `bench/builders/common.py`.** Tasks 5 and 11 use it.
  - `class BuildError(RuntimeError)`.
  - `CONTIG_COLUMNS = ("barcode", "chain", "productive", "high_confidence", "full_length", "cdr3_nt", "raw_clonotype_id")`. These are the 10x `filtered_contig_annotations` columns. `is_cell` is used when present.
  - `parse_flag(values: pd.Series) -> npt.NDArray[np.bool_]`. It is True only for bool `True` and for the strings `"True"`, `"true"` or `"TRUE"` after stripping. Everything else is False, including `"False"`, `"None"`, `""`, `"1"` and NaN.
  - `read_members(tar_path: Path, wanted: Collection[str]) -> tuple[dict[str, bytes], list[str]]`:
    - opens with `tarfile.open(tar_path, "r:")`
    - returns the bytes of each regular-file member whose name is in `wanted`, read with `extractfile`, and the sorted names of every other member
    - never extracts to disk, and member names are only compared, never joined to a path
    - raises `BuildError` for a wanted name that is missing or is not a regular file
  - `MAX_PER_CHAIN = 2`.
  - `call_clones(contigs: pd.DataFrame) -> pd.DataFrame`.
    - **Input:** one row per contig, with `capture` and `CONTIG_COLUMNS`.
    - **Output:** indexed by `(capture, barcode)`, one row per barcode. Columns are `clone_id`, `clone_id_trb` and `raw_clonotype_id` (object, NaN when absent), and `receptor_n_tra` and `receptor_n_trb` (int).
    - **A contig qualifies** when `productive`, `high_confidence` and `full_length` (and `is_cell`, if present) all parse True, `chain` is `TRA` or `TRB`, and `cdr3_nt` fully matches `[ACGTN]+`.
    - **`receptor_n_*`** counts distinct qualifying `cdr3_nt` per chain.
    - **`clone_id`** is set when `1 <= n_tra <= 2` and `1 <= n_trb <= 2`. Its value is `"c" + sha256(";".join(sorted(f"TRA:{nt}" ... | f"TRB:{nt}" ...)))[:16]`, over the distinct qualifying sequences.
    - **`clone_id_trb`** is set when `1 <= n_trb <= 2`. Its value is `"b" +` the same hash over the TRB elements only.
    - **`raw_clonotype_id`** is `f"{capture}:{value}"`, or NaN when the value is missing or `"None"`.
    - **No namespace in the strings:** the same receptor in two captures gives the same string.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_builders_common.py`). The helper `contig(capture, barcode, chain, nt, **flags)` makes one qualifying row unless a flag is overridden.

```python
def test_flag_strings_are_parsed_strictly() -> None:  # Review Focus 1
    s = pd.Series([True, "True", "true", " TRUE ", False, "False", "None", "", np.nan, "1"], dtype=object)
    assert parse_flag(s).tolist() == [True] * 4 + [False] * 6


def test_a_cell_with_three_trbs_is_unlabelled() -> None:
    out = call_clones(
        frame(
            contig("A", "b1", "TRA", "AAA"), *(contig("A", "b1", "TRB", nt) for nt in ("CCC", "GGG", "TTT"))
        )
    )
    row = out.loc[("A", "b1")]
    assert pd.isna(row.clone_id) and pd.isna(row.clone_id_trb) and row.receptor_n_trb == 3
```

  Also:
  - `test_two_of_each_chain_are_allowed_in_any_order`: cell b1 has TRA `AAA`, `ACA` and TRB `CCC`, `CGC`. Cell b2 has the same rows shuffled. They get the same `clone_id`. A cell b3 sharing only the TRBs has a different `clone_id` and the same `clone_id_trb`.
  - `test_non_qualifying_contigs_are_ignored`, parametrized over `productive="False"`, `high_confidence="False"`, `full_length="False"`, `cdr3_nt="None"` and `chain="IGH"` on the TRA row. `clone_id` is NaN and `clone_id_trb` is set.
  - `test_duplicate_contig_rows_count_once`: two identical TRB rows give `receptor_n_trb == 1`.
  - `test_strings_carry_no_namespace`: the same chains in captures A and B give the same `clone_id`. `raw_clonotype_id` `clonotype1` becomes `"A:clonotype1"` and `"B:clonotype1"`.
  - `test_tar_members_are_read_in_memory`:
    - the tar holds the wanted `x.csv.gz`, a member named `../stray.csv.gz` and a directory `d`
    - `read_members` returns `x`'s bytes, with ignored `["../stray.csv.gz", "d"]`
    - the directory holding the tar lists the same files before and after
    - `read_members(tar, ["missing"])` raises `BuildError`
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_builders_common.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** `common.py` as in Interfaces.
- [ ] **Step 4: Run** the same tests. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Call paired and TRB-only clones from 10x contigs, with strict flags and an in-memory tar reader`.

---

### Task 5a: Inspect the Daniel file formats and map captures to cells

**Branch:** `pivot/NN-daniel-formats` · **Model:** opus — the capture map decides which cells every contig joins, and a mistake here is silent downstream.

**Freeze rule:** read file headers, barcodes, feature rows and the contig flag columns only. Nothing is computed or printed from `cdr3`, `cdr3_nt`, `v_gene`, `d_gene`, `j_gene`, `c_gene`, `raw_clonotype_id`, `raw_consensus_id`, `reads` or `umis`.

**Files:**
- Create: `benchmarks/inspection/daniel2022_formats.py`
- Modify: `src/lumencell/bench/builders/daniel2022.py` (facts only), `tests/unit/test_daniel_inspection.py`, `tests/unit/test_daniel2022.py`, the Task 3 research note (adds its Formats and Join check sections)

**Interfaces:**
- Consumes: `CONTIG_COLUMNS` (Task 4); `SOURCES` (Task 2); Task 3's note.
- **Produces, the format script** `uv run python -I benchmarks/inspection/daniel2022_formats.py --dir <datasets dir> [--capture-map <json: tar member -> barcode tag>]`. It prints:
  - for each matrix, the MatrixMarket banner and size line
  - the barcode and feature line counts, and the first three of each
  - whether the first 100,000 stored values are integers
  - the tag vocabulary of the merged barcodes, and whether the D5 features equal the main features
  - the tar member names and each member's header, read in memory
  - the distinct values of `chain`, `productive`, `high_confidence`, `full_length` and `is_cell`, pooled over members
  - the first member's first three barcodes
  - with `--capture-map`, for each member: the share of its distinct barcodes found in the merged barcodes under its own tag, and the largest share under any other tag
- **Produces, `bench/builders/daniel2022.py`.** These are facts from the note; Task 5b and Task 6 read them.
  - `@dataclass(frozen=True) class Capture: gsm: str; member: str; tag: str; matrix: Literal["main", "d5"]; arm: str; day: int; tissue: str; sort_gate: str`
  - `CAPTURES: tuple[Capture, ...]`, one per tar member, with values as written in the note
  - `LABEL_COLUMN: str | None`, Task 3's label column
  - `COLUMN_RENAMES: dict[str, str]`, the real contig spellings mapped to `CONTIG_COLUMNS`; empty if they already match
  - `MATRIX_ORIENTATION: dict[str, Literal["genes_x_cells", "cells_x_genes"]]`, keyed `"main"` and `"d5"`

- [ ] **Step 1: Write the failing tests.**
  - Extend `tests/unit/test_daniel_inspection.py`: plant `CASSPLANTEDF` in `cdr3`, `PLANTEDNT` in `cdr3_nt`, `TRBV-PLANTED` in `v_gene` and `clonotypePLANTED` in `raw_clonotype_id` of two fixture contig members.
    - `test_the_format_script_never_prints_receptor_values`: no `PLANTED` in the output, and the `chain` vocabulary is printed.
    - `test_the_join_share_is_reported_per_member`: a member whose barcodes all sit under its mapped tag shows `1.00`.
  - In `tests/unit/test_daniel2022.py`, add `test_captures_match_the_inspection`:
    - `len(CAPTURES) == 24`, or the number the note records
    - GSMs and members are unique, and every `member` ends in `_tcr_filtered_contig_annotations.csv.gz`
    - `set(MATRIX_ORIENTATION) == {"main", "d5"}`
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_daniel_inspection.py tests/unit/test_daniel2022.py -q`. Expected: FAIL.
- [ ] **Step 3: Write** the script.
- [ ] **Step 4: Run it on the real files** and write the note's **Formats** section:
  - the orientation of each matrix, from the size line against the barcode and feature counts
  - whether it holds integers
  - the feature columns
  - the barcode and tag format
  - the contig header spellings and flag vocabularies
- [ ] **Step 5: Build the capture map** (tar member → tag) from the GSM records and the tag vocabularies, then re-run with `--capture-map`.
  - Each member's own share must be ≥ 0.5 and above every other tag's share. Otherwise the map is wrong: fix it, or stop and report.
  - Write the shares in the note's **Join check** section, and the facts as the constants in Interfaces.
- [ ] **Step 6: Run** the Step 2 command. Expected: PASS.
- [ ] **Step 7: C1–C4.** Commit message: `Record the Daniel 2022 file formats and map each TCR capture to its cells`.

---

### Task 5b: The Daniel builder

**Branch:** `pivot/NN-daniel-builder` · **Model:** opus — the barcode join decides every clonal number after it.

**Freeze rule:** this task does not build the real derived file. The subset and namespace are declared in Task 6. CI uses the synthetic fixture only.

**Files:**
- Create: `tests/unit/daniel_fixture.py`
- Modify: `src/lumencell/bench/builders/daniel2022.py`, `tests/unit/test_daniel2022.py`

**Interfaces:**
- Consumes: `read_members`, `call_clones`, `BuildError` and `CONTIG_COLUMNS` (Task 4); `Builder` (Task 1); `CAPTURES`, `LABEL_COLUMN`, `COLUMN_RENAMES` and `MATRIX_ORIENTATION` (Task 5a).
- **Produces, `bench/builders/daniel2022.py`:**
  - `BUILD_VERSION = 1`
  - `build(paths: Mapping[str, Path], out: Path, params: Mapping[str, Any]) -> None`, a `Builder`:
    - **Params** are `{"captures": [gsm, ...], "namespace": {gsm: str}}`. An empty list, a GSM not in `CAPTURES`, or a namespace that does not cover every listed capture raises `BuildError`.
    - **Cells kept:** those in the metadata and in a listed capture.
    - **obs:**
      - `obs_names` are `f"{gsm}_{core}"`
      - `capture` (the GSM), `arm`, `day`, `tissue` and `sort_gate`, all from `CAPTURES`
      - `state`, from `LABEL_COLUMN`, only when it is not None
      - `donor`, the namespace value
      - the `call_clones` columns; a cell with no contig gets NaN, and 0 for the counts
    - **X** holds the counts, as CSR float32, cells × genes. The mtx is read with `scipy.io.mmread` from a gzip stream, oriented by `MATRIX_ORIENTATION`.
    - **var:** names are the feature IDs. Add `var["gene_symbols"]` when the features carry symbols in their own column.
    - **Two matrices:** when listed captures come from both, their feature lists must be identical, or `BuildError`.
    - **`uns["lumencell_build"]`** is `{"version", "params", "n_cells", "contig_barcodes": {gsm: int}, "join_misses": {gsm: int}, "ignored_members": [...], "labelled_cells": int, "labelled_cells_trb": int}`. A join miss is a contig barcode with at least one qualifying contig that matches no kept cell of its own capture.

- [ ] **Step 1: Write the failing tests.** `tests/unit/daniel_fixture.py` writes a tiny GEO layout in the formats Task 5a recorded:
  - one main matrix of 8 cells × 6 genes, including `Trbc2` and `Trav1` (or their real symbol spellings)
  - metadata with the label column
  - a `GSE188666_RAW.tar` with contig members for captures `GSM_A` and `GSM_B`, the barcode core `AAACCTGAGAAGGCCT` in both, plus a stray `../stray.csv.gz`

  The fixture monkeypatches `CAPTURES` and `LABEL_COLUMN` to match, and builds with `params = {"captures": ["GSM_A", "GSM_B"], "namespace": {"GSM_A": "m1", "GSM_B": "m2"}}`. In `tests/unit/test_daniel2022.py`:

```python
def test_duplicate_barcodes_across_captures_stay_apart(
    built: ad.AnnData, counts: dict[tuple[str, str], np.ndarray]
) -> None:
    a, b = "GSM_A_AAACCTGAGAAGGCCT", "GSM_B_AAACCTGAGAAGGCCT"
    assert built.obs.loc[a, "clone_id"] != built.obs.loc[b, "clone_id"]
    np.testing.assert_array_equal(built[a].X.toarray().ravel(), counts[("GSM_A", "AAACCTGAGAAGGCCT")])
    np.testing.assert_array_equal(built[b].X.toarray().ravel(), counts[("GSM_B", "AAACCTGAGAAGGCCT")])


def test_a_barcode_join_miss_is_counted_and_dropped(built: ad.AnnData) -> None:
    assert built.uns["lumencell_build"]["join_misses"] == {"GSM_A": 1, "GSM_B": 0}
    assert not any(MISSED_CORE in n for n in built.obs_names)


def test_clone_columns_are_scoring_only(built: ad.AnnData) -> None:
    so = scoring_only_columns(built.obs.columns, ("donor",))
    assert {"clone_id", "clone_id_trb", "raw_clonotype_id", "receptor_n_tra", "receptor_n_trb", "donor"} <= so
    assert not {"state", "day", "tissue", "capture", "sort_gate", "arm"} & so
```

  Also:
  - `test_cells_without_contigs_are_unlabelled`
  - `test_the_clone_call_runs_through_the_build`: the fixture also holds a `GSM_A` cell with three distinct qualifying TRBs, and one whose only TRA is `productive="False"`. Both have `clone_id` NaN. The second keeps its `clone_id_trb`. The first has `receptor_n_trb == 3`.
  - `test_only_listed_captures_are_kept_and_namespaced`: with `{"captures": ["GSM_A"], "namespace": {"GSM_A": "m1"}}`, every `capture` is `GSM_A` and every `donor` is `m1`
  - `test_bad_params_are_refused`, parametrized over an empty list, an unknown GSM and a namespace missing a capture; each raises `BuildError`
  - `test_the_stray_member_is_ignored_and_nothing_is_written_beside_out`
  - `test_the_output_loads_as_counts_with_receptor_symbols`: `load(out)` has `x_kind == "counts"` and `receptor_gene_mask(ds.gene_symbols).sum() == 2`
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_daniel2022.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `build` and `BUILD_VERSION` as in Interfaces.
- [ ] **Step 4: Run** the same tests. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Build the Daniel 2022 dataset: join contigs to cells by capture and call clones`.

---

### Task 6: Declare Daniel 2022 before any run (the freeze)

**Branch:** `pivot/NN-daniel-declaration` · **Model:** opus — scientific declarations made by fixed rules from records and literature, which cannot be revised after measuring.

**Freeze rule:**
- This PR merges before the derived Daniel file is built and before any method runs on Daniel.
- It reads only Task 3's note (with Task 5a's sections), the GSM records and the literature.
- Task 3's freeze rule applies. The PR body says: "No clone was counted and no method was run on Daniel 2022 before this PR."

**Files:**
- Modify: `src/lumencell/bench/datasets.py`, `benchmarks/thresholds.toml`, `tests/unit/test_datasets.py`
- Create: `docs/superpowers/research/<YYYY-MM-DD>-daniel2022-declaration.md`

**Interfaces:**
- Consumes: `SOURCES`, `LICENCE`, `CAPTURES`, `LABEL_COLUMN`, `BUILD_VERSION` and `build` (Tasks 2, 5a and 5b); `Build` (Task 1); `ORACLE_ETA_MIN` and `ORACLE_P_MAX` (phase 1).
- **Produces, `bench/datasets.py`:**
  - `RealClonal`: `@dataclass(frozen=True) class RealClonal: clone_key: str; donor_key: str; trb_only_key: str; raw_clonotype_key: str; capture_key: str; assessable: bool; reason: str; declared: str` (`declared` is the ISO date of this declaration). Its `columns(self) -> tuple[str, ...]` returns `(clone_key, donor_key, trb_only_key, raw_clonotype_key)`, the columns named as scoring-only.
  - `RealDataset` gains `clonal: RealClonal | None = None` and `report_time_key: str | None = None`, as its last fields. `__post_init__` raises `ValueError` when `clonal.assessable` is True and fewer than two lineages are gated.
  - `DATASETS["daniel2022"]`:
    - `sources=daniel2022.SOURCES` and `licence=daniel2022.LICENCE`
    - `build=Build(daniel2022.build, daniel2022.BUILD_VERSION, {"captures": [...], "namespace": {...}})`
    - `label_key="state"`, `root_label`, `lineages` and `gated`
    - `clonal=RealClonal("clone_id", "donor", "clone_id_trb", "raw_clonotype_id", "capture", <assessable>, <reason>, "<YYYY-MM-DD>")`
    - `report_time_key="day"`
    - `notes`, giving the declaration date and the note's path
- **Produces, `benchmarks/thresholds.toml`:** `[real.daniel2022.gates]`, `[real.daniel2022.clonal]`, `[real.daniel2022.reported]` and a dated changelog entry.

**Decision rules.** Apply them in order, and write each outcome with its evidence in the declaration note.

| Item | Rule |
|---|---|
| Arm | Chronic (Clone 13) captures only, because the Tpex root holds there (spec §6.1). |
| Days | Every chronic day in the main matrix. Day 5 (the D5 matrix) is kept only if the records say both matrices were processed together and its labels share the main vocabulary. Otherwise it is excluded, with the reason "a matrix processed apart would put a batch on the day axis". |
| Tissues | Only tissues sampled at every kept day, so tissue is not confounded with day. |
| State gate | A sort on a marker the cited papers use to separate the declared states: Slamf6/Ly108, CX3CR1, TIM-3, CD101, KLRG1, CXCR5 or TCF-1. Tetramer and PD-1 gates are not state gates. |
| Namespace | (a) The records name each capture's mouse, and some mouse spans two or more candidate captures: namespace = mouse, keep every candidate, assessable. (b) Otherwise, drop the state-gated candidates. If any remain: namespace = capture, assessable, with the dropped captures listed. (c) Otherwise, when every candidate is state-gated: keep them all, namespace = capture, and `assessable = False`, reason "captures are sorted by state and the mouse cannot be identified across them (spec §6.1)". Pooled mice, or silent records, go into `reason` as a stated assumption, for example "one mouse per capture is assumed; the records do not say". |
| Clone definition | As `common.call_clones` implements it (spec §6.1). TRB-only (`clone_id_trb`) is the sensitivity check. |
| Root | The label whose described markers match Tpex: TCF-1+, Slamf6/Ly108+, PD-1 intermediate, TIM-3− (Im et al. 2016; Utzschneider et al. 2016). |
| Lineages | One lineage per end state that at least one of Hudson 2019, Zander 2019 or Beltra 2020 describes in chronic LCMV as arising from Tpex. That paper must not show it converting into another declared state, and must show this by a method other than TCR clone sharing (transfer of sorted subsets, fate mapping, time course). An intermediate sits on every lineage those papers place it on. Labels are mapped to the papers' subsets by the papers' markers and by Daniel's non-clonal marker description of each cluster, never by Daniel's clonal analysis. A label that no rule places is on no lineage. A link that rests only on an unreadable source is UNVERIFIED and is not used. **Gated** = each lineage whose root, intermediates and end each map to one label with more than 50 cells in the kept captures (Task 3's counts). |
| Two gated lineages | With fewer than two gated lineages there is no clonal comparison between lineages: `clonal_row` would return `not_applicable` and the gates would be skipped without a word. So this case is declared not assessable: `assessable = False`, reason "fewer than two gated lineages, so no clonal comparison between lineages", whatever the namespace rule chose. If both apply, the reason names both. It is shown like the namespace case (Task 7). |
| Fallback | Applies when `LABEL_COLUMN is None`: `lineages = {}` and `gated = []`. Keep every candidate capture: no clonal score depends on the namespace. Set `time_key = "day"`. The root is the Slamf6/Ly108+ gate's cells when a kept capture uses that gate (`label_key = "sort_gate"`), else the earliest kept day (`label_key = "day"`). `clonal.assessable = False`, with reason "no state labels: no lineage truth, so no clonal gate (spec §6.1)". The only gate is `pt_spearman_global = { min = 0.60 }`, the existing value for real time. |

**Literature.**
- Verify each DOI through Crossref, `https://api.crossref.org/works/<doi>`: the title, journal and year must match.
- Read the abstract and the cited results where they are open (PMC, Europe PMC).
- These candidates are cited from memory and are UNVERIFIED until checked:
  - Hudson W.H. et al., Immunity 2019, 10.1016/j.immuni.2019.11.002
  - Zander R. et al., Immunity 2019, 10.1016/j.immuni.2019.10.009
  - Beltra J.-C. et al., Immunity 2020, 10.1016/j.immuni.2020.04.014
  - Im S.J. et al., Nature 2016, 10.1038/nature19330
  - Utzschneider D.T. et al., Immunity 2016, 10.1016/j.immuni.2016.07.021
- If a DOI does not match, find the right one by a Crossref `query.bibliographic` title search, and note the correction.

- [ ] **Step 1: Write the declaration note.** Cover the literature table (paper, verified DOI, the claim used, the method behind it, and where it was read), then each rule's outcome, then the resulting lists.
- [ ] **Step 2: Write the failing tests** in `tests/unit/test_datasets.py`. `CHRONIC` is the arm string as `CAPTURES` writes it.

```python
def test_daniel_is_declared_from_its_builder() -> None:
    ds = DATASETS["daniel2022"]
    assert ds.sources == daniel2022.SOURCES and ds.licence == daniel2022.LICENCE
    assert (
        ds.build is not None
        and ds.build.builder is daniel2022.build
        and ds.build.version == daniel2022.BUILD_VERSION
    )
    caps = {c.gsm: c for c in daniel2022.CAPTURES}
    chosen = ds.build.params["captures"]
    assert (
        chosen
        and {caps[g].arm for g in chosen} == {CHRONIC}
        and set(ds.build.params["namespace"]) == set(chosen)
    )
    assert ds.clonal is not None and ds.clonal.columns() == (
        "clone_id",
        "donor",
        "clone_id_trb",
        "raw_clonotype_id",
    )
    assert ds.report_time_key == "day" and ds.clonal.reason.strip() and date.fromisoformat(ds.clonal.declared)
    assert all(stages[0] == ds.root_label for stages in ds.lineages.values()) and set(ds.gated) <= set(
        ds.lineages
    )


def test_an_assessable_clonal_declaration_needs_two_gated_lineages() -> None:
    ds = DATASETS["daniel2022"]
    assert len(ds.gated) >= 2 or (
        not ds.clonal.assessable and "fewer than two gated lineages" in ds.clonal.reason
    )
    with pytest.raises(ValueError, match="two gated lineages"):
        replace(ds, gated=ds.gated[:1], clonal=replace(ds.clonal, assessable=True))


def test_daniel_gates_are_the_declared_ones() -> None:
    cfg, ds = load_thresholds()["real"]["daniel2022"], DATASETS["daniel2022"]
    base = {
        "pt_spearman_lineage": {"min": 0.70},
        "branch_accuracy": {"min": 0.80},
        "lineage_recall": {"min": 1.0},
    }
    clonal = {
        "clonal_excess_ratio": {"min": 0.5, "nan": "fail", "status": "clonal_status"},
        "clonal_p": {"max": 0.01, "nan": "fail", "status": "clonal_status"},
    }
    assert cfg["gates"] == base | (clonal if ds.clonal.assessable else {})
    assert cfg["clonal"] == {
        "mode": "oracle_cells",
        "oracle_eta_min": ORACLE_ETA_MIN,
        "oracle_p_max": ORACLE_P_MAX,
    }
    assert set(cfg["reported"]) == {
        "day_spearman",
        "clone_recurrence",
        "raw_clonotype_ari",
        "labelled_share_of_raw",
        "trb_only",
    }
```

  In the fallback, `cfg["gates"] == {"pt_spearman_global": {"min": 0.60}}`, and the lineage assertions become `ds.lineages == {}`. Write the branch that applies.
- [ ] **Step 3: Run** `uv run pytest tests/unit/test_datasets.py -q`. Expected: FAIL.
- [ ] **Step 4: Add the registry entry and the TOML.** The gates are exact. When `assessable` is False, the two clonal lines are left out: a result declared not assessable is not gated. It is never hidden either: Task 7 writes it into every Daniel row of the results JSON and every Daniel cell of the README as `clonal: not assessable (declared <date>, <reason>)`.

```toml
# ---- Daniel 2022 (spec §6.1). Declared <YYYY-MM-DD> from the GEO sample records and metadata alone,
# before any clone was counted and before any method ran on it. See the changelog entry of that date.
[real.daniel2022.gates]
pt_spearman_lineage = { min = 0.70 }   # existing real-data value; against state order along each gated lineage
branch_accuracy = { min = 0.80 }       # the stricter existing value (pancreas 0.80, bone marrow 0.75)
lineage_recall = { min = 1.0 }
clonal_excess_ratio = { min = 0.5, nan = "fail", status = "clonal_status" }   # oracle-cells mode, as on synthetic data
clonal_p = { max = 0.01, nan = "fail", status = "clonal_status" }

[real.daniel2022.clonal]
mode = "oracle_cells"    # scored cells come from the truth lineages, weighted by each method's masks (spec §3)
oracle_eta_min = 0.05    # below either floor the case is not_assessable, which fails both clonal gates
oracle_p_max = 0.01

[real.daniel2022.reported]   # reported, never gated
day_spearman = "Spearman of the method's pseudotime against day post infection, over cells finite in both"
clone_recurrence = "share of distinct clone_id strings found in two or more captures, namespace ignored on purpose: an estimate of convergent receptors if captures hold different mice"
raw_clonotype_ari = "adjusted Rand index between clone_id and 10x raw_clonotype_id, both keyed by capture, over cells labelled by both"
labelled_share_of_raw = "cells with a clone_id divided by cells with a raw_clonotype_id"
trb_only = "the clonal block again with clone_id_trb (1 or 2 qualifying TRB, sorted TRB CDR3 nt set): same mode, nulls and floor"
```

  The changelog entry, dated, states:
  - the kept captures, days and tissues, and why
  - the namespace and the rule branch that chose it, with its assumption
  - the clone definition
  - the root, and each lineage with its papers
  - the gates and their reasons: "existing real-data values; where pancreas and bone marrow differ (branch_accuracy 0.80 against 0.75) the stricter is taken, so the choice cannot be read as lenient; clonal gates and floor as on synthetic data, in oracle-cells mode", or "not assessable, declared before any clone was counted: <reason>"
- [ ] **Step 5: Run** `uv run pytest tests/unit/test_datasets.py tests/unit/test_gates.py -q`. Expected: PASS.
- [ ] **Step 6: C1–C4.** Commit message: `Declare the Daniel 2022 subset, namespace, clone definition, truth and gates before any run`.

---

### Task 7: Clonal scores and reported fields for real datasets

**Branch:** `pivot/NN-bench-real-clonal` · **Model:** sonnet — glue over phase 1's `clonal_row`, plus two small metrics checked by hand.

**Files:**
- Create: `src/lumencell/bench/real_clonal.py`, `tests/unit/test_real_clonal.py`
- Modify: `src/lumencell/bench/runner.py`, `src/lumencell/bench/datasets.py` (`truth_for`), `src/lumencell/bench/report.py`, `src/lumencell/cli.py` (`bench real` line)

**Interfaces:**
- Consumes: `clonal_row`, `ClonalRow`, `ClonalRowStatus`, `CaseResult.clonal` and `score_clonal` (phase 1); `from_prep(..., clonal_columns=…)` (phase 1); `stability_row(..., clonal_columns=…)` and `run_real_case(..., stability=False)` (phase 2); `RealClonal` and `report_time_key` (Task 6).
- **Produces:**
  - **`truth_for(ds, obs)`:** when `ds.clonal` is set and `assessable`, the truth carries `clone_ids = obs[clone_key].to_numpy(dtype=object)` and `donors = obs[donor_key].to_numpy(dtype=object)`.
  - **`CaseResult.reported: dict[str, float | str] = field(default_factory=dict)`,** a new last field, merged into `row()` after the clonal fields. Floats are rounded to 4, and NaN stays NaN.
  - **`row()` also always emits `clonal_reason: str`,** which is `ClonalRow.reason`, or `""` when that is None. This extends phase 1's always-emitted fields, so a declared `not_assessable` and its reason reach the results JSON.
  - **`bench/real_clonal.py`:**
    - `adjusted_rand_index(a: npt.NDArray[np.int64], b: npt.NDArray[np.int64]) -> float`, the Hubert–Arabie ARI from the contingency table. It is NaN when n < 2, and 1.0 when the partitions are identical up to relabelling, including when the expected index equals its maximum.
    - `clone_recurrence(clones: npt.NDArray[np.object_], captures: npt.NDArray[np.object_]) -> tuple[float, int]`. Over cells with a non-missing clone string (`clonal.labels.is_missing`), it returns the share of distinct strings found in two or more captures, and the count of distinct strings. It returns `(nan, 0)` when there are none.
    - `real_clonal_row(ds, truth, result, pca, neighbors_pcs: int) -> ClonalRow`:
      - `ClonalRow("not_applicable")` without `ds.clonal`
      - `ClonalRow("not_assessable", reason=f"declared {ds.clonal.declared}, {ds.clonal.reason}")` when not assessable, with nothing scored
      - otherwise `clonal_row(truth, result, pca, neighbors_pcs)`
    - `reported_metrics(ds, obs, result, truth, pca, neighbors_pcs) -> dict[str, float | str]`. Its keys:
      - `reported_day_spearman`: `scipy.stats.spearmanr` of `result.pseudotime` against `obs[report_time_key]` as numbers, over cells finite in both. It is NaN with fewer than two distinct days.
      - When `ds.clonal` is set, `reported_clone_recurrence` and `reported_n_clone_strings`.
      - `reported_raw_clonotype_ari`, over cells with both `clone_key` and `raw_clonotype_key`, using codes of `(capture, clone_id)` and of `raw_clonotype_id`.
      - `reported_labelled_share_of_raw`.
      - `reported_trb_status`, `reported_trb_excess_ratio` and `reported_trb_p`, from `real_clonal_row` on `replace(truth, clone_ids=obs[trb_only_key].to_numpy(dtype=object))`.
  - **`run_real_case`** keeps phase 2's signature, `(dataset_id, method="paga_dpt", params=None, stability=False)`. It passes `clonal_columns=ds.clonal.columns() if ds.clonal else ()` to `from_prep`, and fills `clonal` and `reported`.
  - **`run_real_case` passes the same columns to phase 2's `stability_row(..., clonal_columns=…)`,** so stability replicates also keep clone columns out of inference. Phase 2 code is not changed. The verdict's C5 stays `not_configured` in benchmark rows, as phase 2 set it.
  - **`report.render()`:**
    - A real-table cell adds one of these, and adds nothing for `not_applicable`:
      - `clonal: <ratio>` when `clonal_status == "ok"`
      - `clonal: not assessable (<clonal_reason>)` for `not_assessable`, which for a declared result reads `clonal: not assessable (declared <date>, <reason>)`
      - `clonal: <status>` for any other status
    - A Daniel row is never shown without its clonal entry.
    - Under the real table, for each dataset whose rows carry `reported_*` fields, a table headed "Reported, not gated" with one row per method. Its columns are the `reported_*` fields, plus phase 1's reported-only clonal fields whenever `clonal_status` is not `not_applicable`: `clonal_n_cells`, `clonal_effective_clones`, `clonal_eta_m3` (headed "η at m = 3") and `clonal_p_m3` (headed "p at m = 3").
  - **`bench real`** adds `clonal=<status>` to each printed line.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_real_clonal.py`).
  - **The `fake` fixture** takes a 1500-cell synthetic bifurcation with `add_clones(DEFAULT_LAYER, seed=0)` and gives it a built Daniel file's columns:
    - `state = true_label` and `capture = donor`
    - `clone_id_trb = clone_id`, with every fifth label blanked
    - `raw_clonotype_id = capture + ":" + clone_id`
    - `receptor_n_tra = receptor_n_trb = 1`
    - `day` set to 5, 8 and 21 by true-pseudotime tertile
  - It writes the file under `LUMENCELL_DATASETS_DIR` and registers it with `monkeypatch.setitem(DATASETS, "fake_tcell", RealDataset(...))`, using:
    - one `Source(...)` of the file's real size and SHA-256
    - `licence="test"`, `label_key="state"`, `root_label="trunk_early"`
    - lineages from the synthetic labels
    - `clonal=RealClonal("clone_id", "donor", "clone_id_trb", "raw_clonotype_id", "capture", True, "test", "2026-01-02")`
    - `report_time_key="day"`

```python
def test_a_real_dataset_with_clones_is_scored_in_oracle_cells_mode(fake: str) -> None:
    row = run_real_case(fake, "paga_dpt").row()
    assert row["clonal_status"] in get_args(ClonalRowStatus) and row["clonal_status"] != "not_applicable"
    assert math.isfinite(row["clonal_eta_oracle"])
    assert {
        "reported_day_spearman",
        "reported_clone_recurrence",
        "reported_raw_clonotype_ari",
        "reported_labelled_share_of_raw",
        "reported_trb_status",
    } <= row.keys()
    assert {
        "clonal_n_cells",
        "clonal_effective_clones",
        "clonal_eta_m3",
        "clonal_p_m3",
    } <= row.keys()  # phase 1, reported only


def test_every_clone_column_of_a_built_dataset_is_scoring_only(
    fake: str, monkeypatch: pytest.MonkeyPatch
) -> None:  # Review Focus 4
    seen: dict[str, InferenceContext] = {}
    real = runner.run_method
    monkeypatch.setattr(
        runner, "run_method", lambda ctx, *a, **k: (seen.setdefault("ctx", ctx), real(ctx, *a, **k))[1]
    )
    run_real_case(fake, "paga_dpt")
    cols = {"clone_id", "clone_id_trb", "raw_clonotype_id", "receptor_n_tra", "receptor_n_trb", "donor"}
    assert cols <= seen["ctx"].scoring_only and not cols & set(seen["ctx"].obs.columns)
    assert {"state", "day", "capture"} <= set(seen["ctx"].obs.columns)


@pytest.mark.slow
def test_stability_replicates_keep_clone_columns_out_too(fake: str, monkeypatch: pytest.MonkeyPatch) -> None:
    ctxs: list[InferenceContext] = []
    real = stability_mod.run_method
    monkeypatch.setattr(
        stability_mod, "run_method", lambda ctx, *a, **k: (ctxs.append(ctx), real(ctx, *a, **k))[1]
    )
    run_real_case(fake, "paga_dpt", stability=True)
    assert ctxs and all(
        {"clone_id_trb", "donor"} <= c.scoring_only and "donor" not in c.obs.columns for c in ctxs
    )


def test_adjusted_rand_index_by_hand() -> None:
    assert adjusted_rand_index(np.array([0, 0, 1, 1]), np.array([0, 0, 1, 2])) == pytest.approx(4 / 7)
    assert adjusted_rand_index(np.array([0, 0, 1, 1]), np.array([5, 5, 7, 7])) == 1.0
    assert math.isnan(adjusted_rand_index(np.array([0]), np.array([0])))


def test_clone_recurrence_by_hand() -> None:
    share, n = clone_recurrence(
        np.array(["X", "X", "Y", "Z", "Z", None], object), np.array(["A", "B", "A", "B", "B", "A"], object)
    )
    assert share == pytest.approx(1 / 3) and n == 3
```

```python
def test_a_declared_not_assessable_result_is_shown_not_omitted(
    fake: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    ds = DATASETS[fake]
    monkeypatch.setitem(
        DATASETS,
        fake,
        replace(ds, clonal=replace(ds.clonal, assessable=False, reason="r", declared="2026-01-02")),
    )
    monkeypatch.setattr(real_clonal_mod, "clonal_row", _raise)  # nothing may be scored
    row = run_real_case(fake, "paga_dpt").row()
    assert (row["clonal_status"], row["clonal_reason"]) == ("not_assessable", "declared 2026-01-02, r")
    assert row["reported_trb_status"] == "not_assessable"
    assert '"clonal_reason": "declared 2026-01-02, r"' in json.dumps(row)  # the results JSON
    assert not [f for f in check_gates(row, BASE_GATES) if f.metric.startswith("clonal")]  # not gated
    monkeypatch.setattr(
        report,
        "latest_real",
        lambda: (RESULTS_PATH, {"commit": "x", "methods": {"paga_dpt": {"daniel2022": row}}}),
    )
    assert "clonal: not assessable (declared 2026-01-02, r)" in report.render()  # the README table
```

  `real_clonal_mod` is `lumencell.bench.real_clonal`, which calls `clonal_row` through its own module namespace. `_raise` raises `AssertionError`. `BASE_GATES` holds the three base Daniel gates. `RESULTS_PATH` is `REPO / "benchmarks/results/real-test.json"`, which is never written.

  Also:
  - `test_day_spearman_is_nan_with_one_day`.
  - `test_the_readme_shows_reported_fields`: `report.latest_real` is monkeypatched to return a results dict with a `daniel2022` row carrying `reported_*`, `clonal_status="ok"` and `clonal_excess_ratio=0.8`. `render()` contains `Reported, not gated`, `clonal: 0.80` and `η at m = 3`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_real_clonal.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** as in Interfaces.
- [ ] **Step 4: Run** `uv run pytest tests/unit/test_real_clonal.py tests/unit/test_datasets.py tests/unit/test_gates.py -q`, then `uv run pytest -m slow tests/unit/test_real_clonal.py -q` (`stability_mod` is `lumencell.engine.stability`). Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Score clonal agreement and report convergence, 10x agreement and TRB-only results on real datasets`.

---

### Task 8: Run the real benchmark with Daniel 2022

**Branch:** `pivot/NN-bench-real-daniel` · **Model:** sonnet — runs, checks and a logged report.

**Files:**
- Create: `benchmarks/results/real-<YYYY-MM-DD>.json`
- Modify: `README.md`, `benchmarks/thresholds.toml` (changelog only), `tests/real/test_real_datasets.py`

**Interfaces:**
- Consumes: Tasks 1–7 merged.
- Produces: the first results file with `daniel2022` rows.

- [ ] **Step 1: Check what is downloaded:**

```bash
uv run python -c "from lumencell.bench.datasets import DATASETS; from lumencell.bench.download import is_downloaded; print({k: all(map(is_downloaded, d.sources)) for k, d in DATASETS.items()})"
```

  If any value is `False`, **stop and ask the owner** before `uv run lumencell bench fetch <id>`, stating the sizes the command prints.
- [ ] **Step 2: Write the build check as a test.** In `tests/real/test_real_datasets.py`, add `test_the_daniel_build_joins_every_capture` (`real_data`, skipped unless `is_cached`). Read `uns["lumencell_build"]` with `anndata.read_h5ad(path, backed="r")`, then assert:
  - every listed capture's joined share, `1 - join_misses / contig_barcodes`, is ≥ 0.5 (a lower share means the capture table is wrong)
  - obs has `state` (unless the fallback applies), `day`, `capture`, `donor`, `clone_id` and `clone_id_trb`
  - the root label has cells

  These are barcode and cell counts, not clone counts. Run `uv run pytest -m real_data tests/real -q -k daniel`. Expected: SKIPPED (not built yet).
- [ ] **Step 3: Build the derived file:** `uv run lumencell bench fetch daniel2022; echo $?` → `0`. This downloads nothing; it verifies every SHA-256 and builds.
- [ ] **Step 4: Run** the Step 2 test. Expected: PASS. If a capture fails, stop and report; never lower the 0.5. Record each share, and `n_cells` against Task 3's label counts for the kept captures, in HANDOFF.
- [ ] **Step 5 (background):** `make bench-real > "$TMPDIR/real.log" 2>&1; echo $?` → `0`. Gate failures do not fail the command.
- [ ] **Step 6: Report.** For each method on Daniel, give:
  - every gated metric with its pass or fail, as it is
  - the clonal status, ratio, p, `delta_base` and oracle η, with the cells, effective clones, and η and p at m = 3
  - the stability fields and the verdict level that phase 2 adds to each row
  - when the declaration says not assessable: every Daniel row in the results file has `clonal_status == "not_assessable"` and a `clonal_reason` starting `declared `, and after Step 9, `grep -c "clonal: not assessable (declared" README.md` prints 3 (one per method)
  - the reported fields
  - Run `data:validate-data` on the Daniel rows: values are finite wherever the status is `ok`, ranges are plausible, and the reported fields match their declared definitions.
  - Gate failures stay failures. Nothing is re-run with other settings.
- [ ] **Step 7: Compare** pancreas, bone marrow and gastrulation with the previous results file. No metric may move, because their code path did not change. If one moved, stop and report.
- [ ] **Step 8: Add a changelog line,** dated: "Daniel 2022 first run, benchmarks/results/real-<date>.json: <one line per method>. Thresholds unchanged."
- [ ] **Step 9: Regenerate the tables:** `uv run lumencell bench report --readme`.
- [ ] **Step 10: C1–C4.** Commit message: `Run the real benchmark on Daniel 2022 and report its gates and clonal agreement as they are`.

---

### Task 9: Verify that the Sheih 2020 TCR tables join to cells

**Branch:** `pivot/NN-sheih-verify` · **Model:** sonnet — a bounded check with a pass rule fixed in advance.

**Pass rule** (fixed before looking): the clonotype table carries per-cell barcodes, either as a column or as barcode lists per clonotype, and at least 0.5 of its distinct barcodes are among the matrix's cell IDs under one stated mapping. Otherwise the case study moves to Cao et al. 2025 (GSE207192).

**Freeze rule:** headers, barcodes and counts only. No CDR3, gene or clonotype value is printed or tabulated.

**Files:**
- Create: `benchmarks/inspection/barcode_join.py`, `tests/unit/test_barcode_join.py`, `docs/superpowers/research/<YYYY-MM-DD>-sheih-verification.md`

**Interfaces:**
- **Produces, the script:** `uv run python -I benchmarks/inspection/barcode_join.py --cells <one cell ID per line> --table <clonotype table> --suffix <text>`.
  - It prints the table's column names, dtypes and row count.
  - It prints the barcode column: one whose values match `[ACGT]{16}(-\d+)?`, or `;`/`,`-separated lists of such values. It also prints three example barcodes.
  - It prints the share of distinct table barcodes found among the cell IDs after mapping (the 16-nt core plus `--suffix`).
- **Produces, the note:**
  - the verdict and the mapping rule
  - the GSM table: GSM, patient, time point label as written, day, and the aggr sample ID ("Sample ID in raw.expmatrix.csv")
  - the table and matrix formats, including the matrix orientation
  - the licence lines
  - Tasks 10 and 11 read it.

- [ ] **Step 1: STOP for the owner's OK.**
  - First read the GEO listing at `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125881/suppl/`. That is names and sizes only, no data.
  - Pick the largest post-infusion clonotype table.
  - Ask: "To check Sheih 2020 (GSE125881) before using it, may I download <name> (<listed size>) from <URL>, and read the first line of GSE125881_raw.expMatrix.csv.gz (streamed, stopping at the first newline or after 16 MB of its 109 MB)?"
  - Note in the message that the spec's "about 7 MB" does not match the listing (6.9 KB to 297 KB per table).
- [ ] **Step 2: Write the failing tests** (`tests/unit/test_barcode_join.py`):
  - `test_join_share_by_suffix`: cells `AAAC…-2`, `AAAG…-2`, `AAAC…-1` and table barcodes `AAAC…-1`, `AAAG…-1`, `CCCC…-1`, with `--suffix -2`. The output shows a share of `0.67`.
  - `test_barcode_lists_in_one_column_are_split`.
  - `test_no_receptor_values_are_printed`: a planted `CASSPLANTEDF` and `TRBV-PLANTED` are absent from the output.
- [ ] **Step 3: Run** them (FAIL), implement the script, then run them (PASS).
- [ ] **Step 4: Download into a fresh directory.**
  - Fetch the table with `LUMENCELL_DATASETS_DIR="$TMPDIR/sheih-verify"` and `fetch_source(Source(url, name, size, None))`, and record its SHA-256.
  - Stream the matrix's first line with `urllib` and the certifi context through `gzip.GzipFile`, capped at 16 MB compressed, into `$TMPDIR/sheih-verify/header.txt`.
- [ ] **Step 5: Work out the orientation.** If the first line lists cell barcodes, write them one per line for `--cells`. If it lists genes (cells are rows), **stop and ask the owner** for the whole matrix (its listed size), then take its first column.
- [ ] **Step 6: Read the 32 GSM records** in their text form (as in Task 3) into the GSM table.
- [ ] **Step 7: Run the script** with the table's aggr sample ID as the suffix, and apply the pass rule.
- [ ] **Step 8: On a pass,** write the note. **On a fail,** take the Cao branch:
  1. **STOP for OK.** Ask: "Sheih's tables do not join (<share>). May I download from GSE207192 the two metaData files (319 KB and 134 KB), the two cells_colnames files (255 KB and 104 KB) and the aggregated filtered_contig_annotations.csv.gz (16 MB)?"
  2. **Resolve 3 versus 7 patients first.** Compare the patient IDs in the metadata, the GSM records' characteristics and Table 2 of PMC12055983. The case study covers the patients present in both the metadata and the TCR file. The note records that count and why the GEO text and the paper differ.
  3. **Run the join check** with the same rule: the aggr suffix against cells_colnames.
  4. **Decide.** If Cao passes, the note names Cao 2025 and every later task takes its Cao branch. If Cao fails, stop and report to the owner; no case study is designed.
- [ ] **Step 9: Write the licence lines:**
  - Sheih: article CC BY 4.0 (PMC6954177), with GEO terms as in Task 2.
  - Cao, if used: CC BY-NC-ND 4.0 (PMC12055983).
- [ ] **Step 10: C1–C4.** Commit message: `Verify that the case-study TCR tables join to expression cells`.

---

### Task 10: Declare the case study before analysing it

**Branch:** `pivot/NN-sheih-declaration` · **Model:** opus — the declaration fixes what the case study may claim.

**Freeze rule:** this declaration is written from the paper and Task 9's note only. The full data is not yet downloaded, and nothing is computed from clonotype values. It merges before Task 11.

**Files:**
- Modify: `benchmarks/thresholds.toml`
- Create: `tests/unit/test_case_study_declaration.py`

**Interfaces:**
- **Produces** `[case_study.sheih2020]` (or `[case_study.cao2025]`), which Tasks 11 and 12 read:

```toml
[case_study.sheih2020]   # declared before any analysis; see the changelog entry of the same date
declared = "<YYYY-MM-DD>"
accession = "GSE125881"
paper = "Sheih et al. 2020, Nat Commun 11:219, doi:10.1038/s41467-019-13880-1"
methods = ["paga_dpt", "palantir", "scfates_ppt"]
root = "cluster:timepoint=<value as written>"   # the paper's starting population: <section or figure>
clone_key = "clone_id"
donor_key = "patient"          # namespace: each patient
single_donor = false
strata_keys = ["timepoint"]    # N3: clones that expand at one time point are not read as lineage agreement
chains = "<paired|trb>"
cdr3 = "<nt|aa>"
filters = ["<each contig flag present and applied>"]
max_per_chain = 2
unit_key = "patient"           # primary analysis: each method runs once per patient (4 x 3 runs, stability on each)
pooled_diagnostic = true       # one pooled run per method, reported only as a labelled diagnostic, never the headline
# Clonal parameters are the [clonal] defaults: m = 2 primary, m = 3 sensitivity, B = 1000, seed 0.
# lumenCell has no batch correction. Pooled patients would let patient structure dominate, and C1/C2
# would then call a batch artefact "no tree". Cost: 12 runs with stability instead of 3, and about a
# quarter of the cells in each run.
```

- **Rules:**
  - **Chains.** Use `paired` (Daniel's rule: 1 or 2 of each of TRA and TRB, the sorted set) if the tables let both chains be grouped by cell barcode. Otherwise use `trb` (1 or 2 TRB), with the reason "the GEO tables are unpaired clonotype tables".
  - **CDR3.** Nucleotide when present; otherwise amino acid, noting that amino acid merges convergent clones.
  - **Filters.** Each of productive, high-confidence and full-length is applied when the tables carry it. Missing ones are listed in the changelog.
  - **Root.** The population the paper describes as the start of the cells' differentiation, cited by section or figure. If the paper names none, use the earliest sample (the infusion product), stated as an assumption. The time point is a stratum, so it stays visible to methods and may define the root (spec §4.2, amended).
  - **Unit of analysis.** Per patient is primary. Each patient's run holds one donor, so its nulls N1 to N3 are within that patient, and the root's time point must exist in every patient. If a patient has too few clones, its runs report `too_few_clones`; patients are never merged to rescue a score. The pooled run is a diagnostic only.
  - **Cao branch.** `accession = "GSE207192"`, `paper = "Cao et al. 2025, Nat Commun 16:4205, doi:10.1038/s41467-025-59298-w"`, and the patients from Task 9, per patient in the same way.
- **Produces** a dated changelog entry giving each choice and its reason. It records the per-patient decision with its cost: 12 runs with stability instead of 3, and fewer cells per run.

- [ ] **Step 1: Write the failing test:**

```python
def test_case_study_is_declared() -> None:
    cs = load_thresholds()["case_study"]["sheih2020"]
    s = ClonalSettings(
        clone_key=cs["clone_key"],
        donor_key=cs["donor_key"],
        single_donor=cs["single_donor"],
        strata_keys=cs["strata_keys"],
    )
    assert (
        s.donor_key == "patient" and s.strata_keys == ["timepoint"] and "timepoint" not in settings_columns(s)
    )
    root = parse_root(cs["root"])
    assert isinstance(root, ClusterRoot) and root.key == "timepoint"
    assert cs["methods"] == ["paga_dpt", "palantir", "scfates_ppt"] and cs["max_per_chain"] == 2
    assert (
        cs["chains"] in ("paired", "trb")
        and cs["cdr3"] in ("nt", "aa")
        and date.fromisoformat(cs["declared"])
    )
    assert cs["unit_key"] == "patient" == cs["donor_key"] and cs["pooled_diagnostic"] is True
```

  If the paper's root falls on another column the builder will write, assert that column instead, and say why in the PR.
- [ ] **Step 2: Run** it. Expected: FAIL (no section).
- [ ] **Step 3: Write** the section and the changelog entry from the paper and Task 9's note.
- [ ] **Step 4: Run** it. Expected: PASS.
- [ ] **Step 5: C1–C4.** Commit message: `Declare the case-study clone definition, namespace, strata and root before any analysis`.

---

### Task 11: Fetch and build the case-study dataset

**Branch:** `pivot/NN-sheih-builder` · **Model:** sonnet — a second builder on the Task 5b pattern, reusing `common.py`.

**Files:**
- Create: `src/lumencell/bench/builders/sheih2020.py` (or `cao2025.py`), `tests/unit/test_sheih2020.py`
- Modify: `src/lumencell/bench/datasets.py`, `src/lumencell/cli.py` (`bench fetch` also looks up case studies), `tests/real/test_real_datasets.py`

**Interfaces:**
- Consumes: `Source`, `Build`, `fetch_source` and `local_path` (Task 1); `read_members`, `parse_flag`, `call_clones` and `BuildError` (Task 4); `[case_study.sheih2020]` (Task 10); Task 9's note.
- **Produces, `bench/datasets.py`:**
  - `@dataclass(frozen=True) class CaseStudyDataset: id: str; title: str; citation: str; sources: tuple[Source, ...]; licence: str; build: Build; notes: str = ""`, which is a `Fetchable`.
  - `CASE_STUDIES: dict[str, CaseStudyDataset]`, holding `"sheih2020"`.
  - `bench fetch all` still covers only `DATASETS`. A case study is fetched by its ID.
- **Produces, `bench/builders/sheih2020.py`:**
  - `BUILD_VERSION = 1`
  - `SOURCES`: the matrix and `GSE125881_RAW.tar` (the 16 tables), with sizes from `curl -sI` and SHA-256 after download
  - `LICENCE`
  - `@dataclass(frozen=True) class Sample: gsm: str; member: str; patient: str; timepoint: str; day: int; aggr_id: int` and `SAMPLES`, from Task 9's GSM table
  - `build(paths, out, params)`, whose `params` equal the declaration's `chains`, `cdr3`, `filters` and `max_per_chain`:
    - It reads the expression csv in chunks of `chunk_rows` rows (default 500, with the orientation from the note) into CSR, never as one dense array.
    - It maps each cell to its sample by the aggr suffix.
    - It writes obs `patient`, `timepoint`, `day`, `clone_id` and the `receptor_n_*` counts, and `uns["lumencell_build"]` as in Task 5b.
    - With 10x-style contig rows, the clones come from `call_clones`. Otherwise it reads the table's per-cell chains under the declared rule, hashed the same way.

- [ ] **Step 1: STOP for the owner's OK.**
  - **Sheih:** "May I download the rest of GSE125881 into the datasets cache: GSE125881_raw.expMatrix.csv.gz (<size>) and GSE125881_RAW.tar (<size>, the 16 clonotype tables), about 111 MB?"
  - **Cao branch:** the post- and pre-infusion RNA mtx, features, cells_colnames, metaData and the aggregated contig file, with their sizes, less what Task 9 fetched. The ADT and CD4 files are not fetched.
- [ ] **Step 2: Write the failing tests** (`tests/unit/test_sheih2020.py`). The fixture is a genes × cells csv.gz of 5 genes and 6 cells over two samples, plus a tar with two tables in Task 9's format. It uses `chunk_rows=2`.
  - `test_the_matrix_is_read_in_chunks_without_changing_values`: X equals the fixture counts, transposed.
  - `test_cells_map_to_patient_and_time_point_by_suffix`.
  - `test_the_same_receptor_in_two_patients_is_two_clones`: the same `clone_id` string and different `patient`, so `encode_clones(clone_id, patient)` gives two IDs.
  - `test_the_declared_chains_are_used`.
  - `test_a_join_miss_is_counted`.
  - `test_params_mirror_the_declaration`: `CASE_STUDIES["sheih2020"].build.params == {k: cs[k] for k in ("chains", "cdr3", "filters", "max_per_chain")}`.
  - `test_case_study_sources_are_complete`: every source has its SHA-256, an https URL and a positive size.
- [ ] **Step 3: Run** `uv run pytest tests/unit/test_sheih2020.py -q`. Expected: FAIL.
- [ ] **Step 4: Implement** the builder, the registry entry with `sha256=None`, and the CLI lookup.
- [ ] **Step 5 (background): Download** each source with `fetch_source` (as in Task 2 Step 6), then fill the SHA-256 values.
- [ ] **Step 6: Build:** `uv run lumencell bench fetch sheih2020; echo $?` → `0`. Every sample's joined share must be ≥ 0.5, or stop. Record `n_cells` against the paper's 62,167 and the labelled total in HANDOFF.
- [ ] **Step 7: Add** `test_case_study_files_match_their_checksums` to `tests/real/test_real_datasets.py` (`real_data`). Run `uv run pytest tests/unit/test_sheih2020.py tests/real -q`. Expected: PASS.
- [ ] **Step 8: C1–C4.** Commit message: `Fetch and build the case-study dataset under its declared clone definition`.

---

### Task 12: The case-study numbers

**Branch:** `pivot/NN-case-study-numbers` · **Model:** opus — statistics-heavy, and it joins phase 1 and phase 2 outputs.

**Files:**
- Create: `src/lumencell/bench/case_study.py`, `benchmarks/case_studies/sheih2020.py`, `benchmarks/case_studies/results/sheih2020-<YYYY-MM-DD>.json`, `tests/unit/test_case_study.py`
- Modify: `src/lumencell/bench/report.py` (`replace_section`)

**Interfaces:**
- Consumes:
  - **Phase 1:** `ClonalSettings`, `check_settings`, `settings_columns`, `agreement_for_run`, `from_prep(..., clonal_columns=…)`, `run_method` and `TrajectoryResult.obs_keys`.
  - **Phase 2:**
    - `run_stability`, `RunSpec` and `StabilityError` [T6]
    - `tree_likeness`, `coverage`, `verdict`, `StabilityEvidence` and `calibrated_thresholds` [T7, T8]
    - `to_api_stability` and `not_measured` [T10]
    - `stability_evidence`, `clonal_evidence` and `to_api_tree_support` [T11]
    - `cell_cycle_check` and the `cc_phase` and `cc_info` artifacts [T3, T4]
    - The case study calls these exactly as the worker and `GET` do. Nothing is re-implemented.
  - Task 10's declaration, and Task 11's `CASE_STUDIES` with `local_path`.
- **Produces, `bench/case_study.py`:**
  - `@dataclass(frozen=True) class CaseStudyDeclaration: settings: ClonalSettings; root: RootSpec; methods: list[str]; unit_key: str; pooled_diagnostic: bool`.
  - `load_declaration(name: str, cfg: dict[str, Any] | None = None) -> CaseStudyDeclaration`, reading `[case_study.<name>]`.
  - `split_units(path: Path, unit_key: str, workdir: Path) -> dict[str, Path]` writes one `.h5ad` per distinct non-missing value of `unit_key`, holding that unit's cells with X, var and obs unchanged. Each file is named `unit-<first 12 hex of sha256(value)>.h5ad` inside `workdir`, so no value reaches a path, and the input file is not modified.
  - `run_units(path: Path, decl: CaseStudyDeclaration, workdir: Path, include_stability: bool = True) -> dict[str, Any]` returns `{"declaration": <the TOML section>, "units": {value: run_case_study(<unit file>, …)}, "pooled_diagnostic": run_case_study(path, …) | None}`. The pooled run happens only when `decl.pooled_diagnostic` is set. Prep caches go under `workdir`.
  - `run_case_study(path: Path, settings: ClonalSettings, root: RootSpec, methods: Sequence[str], include_stability: bool = True) -> dict[str, Any]`:
    - It loads the file, checks `check_settings` (any problem raises `ValueError`) and prepares through `PrepCache`.
    - It builds `InferenceContext.from_prep(ds, prep, params, clonal_columns=settings_columns(settings))` once, then calls `run_method(ctx, method, root)` for each method.
    - It then calls `agreement_for_run(result, obs, settings, pca, neighbors_pcs, result.obs_keys)`. This is own-cells mode with m = 2 and m = 3, N1–N3, the baseline and coupling.
    - **Stability:** `to_api_stability(run_stability(ds, prep, params, result, RunSpec(method, result.params, None, None, settings_columns(settings))))`.
      - It passes the reference run's `result.params`, as the worker passes `graph["params"]`. With `{}`, `method_seed` would be `None`, and the Palantir and scFates reseeds would never reseed, which inflates C3 and C4.
      - A `StabilityError` becomes `not_measured("failed", str(exc))`. With `include_stability=False`, the block is `None`.
    - **Verdict:** `to_api_tree_support(verdict(coverage(result.pseudotime), tree_likeness(prep.artifacts["paga_connectivities"], prep.artifacts["paga_tree"]), <evidence>, clonal_evidence(agreement), calibrated_thresholds()))`. The evidence is `stability_evidence(<stability block>)`, or `StabilityEvidence("pending")` when stability is skipped, which gives `provisional`.
    - **Cell cycle:** `cell_cycle_check(result.primary_lineage, [l.name for l in result.lineages], prep.artifacts["cc_phase"], prep.artifacts["cc_info"])`.
    - **Dataset level:** the tree-likeness value (C2) goes into `report["dataset"]["tree_likeness"]`. A NaN (no PAGA edges) is written as `null`, with `report["dataset"]["tree_likeness_note"] = "no PAGA edges, so tree-likeness is undefined and C2 fails"`. Every other float in the report is finite or already `null`, so `allow_nan=False` never raises on a valid run.
    - It returns `{"dataset": {"name", "n_cells", "n_genes", "fingerprint", "tree_likeness", "tree_likeness_note"}, "settings", "root", "methods": {m: {"n_lineages", "coverage", "runtime_seconds", "clonal", "stability", "tree_support", "cell_cycle"}}, "versions": {"lumencell", "scanpy", "numpy"}}`. Model blocks are `model_dump(mode="json")`.
  - `write_report(report: dict[str, Any], path: Path, commit: str) -> None` adds `generated` and `commit`, writes `json.dumps(..., indent=2, sort_keys=True, allow_nan=False)`, and replaces the file atomically.
  - `render_markdown(report: dict[str, Any]) -> str` takes a `run_units` report. It renders tables only, with no prose, and prints a missing value as `–`.
    - First comes a heading `Per patient (primary)`, holding:
      - one row per patient and method: L, verdict, clonal status, cells, clones, effective clones, strata, η against N1, N2 and N3, headline p, η_base, Δη, p_Δη, η and p at m = 3, purity, and purity without the largest clone
      - a coupling table per patient and method: a, b, observed, z and p_adj
      - a stability table: median and minimum for each check, with reseeds and subsamples apart
      - the cell-cycle warnings
    - Then, only when present, a heading `Pooled diagnostic (not the headline: patient structure is not corrected)`. It holds the same row table for the pooled runs, with donors scored added.
  - `report.replace_section(path: Path, text: str, start: str, end: str) -> None`. `update_readme` now calls it with the bench markers. The case-study markers are `<!-- case-study:start -->` and `<!-- case-study:end -->`.
- **Produces, `benchmarks/case_studies/sheih2020.py`:** `uv run python benchmarks/case_studies/sheih2020.py [--out PATH] [--workdir DIR]`.
  - It loads the declaration and `local_path(CASE_STUDIES["sheih2020"])`.
  - It calls `run_units` (4 patients × 3 methods with stability on each, plus the pooled diagnostic) in `--workdir`, which defaults to a temporary directory.
  - It writes `benchmarks/case_studies/results/sheih2020-<date>.json`.
  - It holds no logic beyond that wiring.

- [ ] **Step 1: Write the failing tests** (`tests/unit/test_case_study.py`).
  - **The `report` fixture** uses a 600-cell bifurcation with `add_clones(DEFAULT_LAYER, seed=0)` and a `timepoint` column of `t0`, `t1`, `t2` by true-pseudotime tertile.
  - It uses settings `clone_id` / `donor` / `["timepoint"]`, the root `cluster:timepoint=t0`, methods `["paga_dpt"]` and `include_stability=False`.

```python
def test_every_block_is_reported(report: dict[str, Any]) -> None:
    m = report["methods"]["paga_dpt"]
    assert [n["null"] for n in m["clonal"]["primary"]["nulls"]] == ["N1", "N2", "N3"]
    assert m["clonal"]["sensitivity_m3"]["min_clone_size"] == 3 and "eta_base" in m["clonal"]["primary"]
    assert {"coupling", "notes"} <= m["clonal"].keys() and {"tree_support", "cell_cycle"} <= m.keys()


def test_the_time_point_can_be_the_root(report: dict[str, Any]) -> None:
    assert report["methods"]["paga_dpt"]["clonal"]["status"] != "refused"
```

  The `units` fixture runs `run_units` on a 900-cell version of the same data, with `unit_key="donor"` (three units: `d0`, `d1`, `d2`), `pooled_diagnostic=True`, methods `["paga_dpt"]` and `include_stability=False`:

```python
def test_each_patient_runs_alone_and_the_pooled_run_is_a_labelled_diagnostic(units: dict[str, Any]) -> None:
    assert set(units["units"]) == {"d0", "d1", "d2"}
    assert all(u["methods"]["paga_dpt"]["clonal"]["primary"]["n_donors"] == 1 for u in units["units"].values())
    assert units["pooled_diagnostic"]["methods"]["paga_dpt"]["clonal"]["primary"]["n_donors"] == 3
    text = render_markdown(units)
    assert text.index("Per patient (primary)") < text.index("Pooled diagnostic (not the headline")

def test_unit_files_are_named_by_hash(tmp_path: Path) -> None:
    # a unit value of "../x" writes unit-<12 hex>.h5ad inside workdir, nothing outside it, and leaves the input's bytes unchanged
```

  Also:
  - `test_the_report_is_strict_json_and_reproducible`: two runs, through `write_report`, are equal after dropping `generated` and `runtime_seconds`. `json.loads(text, parse_constant=_raise)` parses.
  - `test_render_has_one_row_per_patient_and_method_and_no_nan`.
  - `test_stability_gets_the_reference_params`: monkeypatch `case_study.run_stability` to record its `RunSpec`, and run with `include_stability=True`. `spec.params == result.params` and `spec.clonal_columns == settings_columns(settings)`.
  - `test_a_nan_tree_likeness_is_written_as_null_with_a_reason`: monkeypatch `case_study.tree_likeness` to return NaN. `write_report` succeeds, the file has `"tree_likeness": null`, and `tree_likeness_note` is non-empty.
  - `@pytest.mark.slow test_reseeds_use_new_method_seeds`: with Palantir on the 600-cell data and `include_stability=True`, every reseed replicate's `method_seed` is not `None` and differs from `result.params["seed"]`.
  - `test_the_declaration_is_read_from_the_thresholds_file`: `load_declaration("sheih2020")` has `unit_key == "patient"` and `pooled_diagnostic is True`.
  - `@pytest.mark.slow test_stability_and_the_verdict_are_reported`: with `include_stability=True`, `m["stability"]` is set and `m["tree_support"]["level"]` is in `{"supported", "weak", "not_supported"}`.
- [ ] **Step 2: Run** `uv run pytest tests/unit/test_case_study.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the module, `replace_section` and the script.
- [ ] **Step 4: Run** `uv run pytest tests/unit/test_case_study.py tests/unit/test_real_clonal.py -q` and `uv run pytest -m slow tests/unit/test_case_study.py -q` (Palantir is installed by `make install`). Expected: PASS.
- [ ] **Step 5 (background, hours):** `uv run python benchmarks/case_studies/sheih2020.py > "$TMPDIR/case.log" 2>&1; echo $?` → `0`.
- [ ] **Step 6: Check the statistics.** Run `data:statistical-analysis`, then `data:validate-data`, on the results JSON. They must address:
  - whether every headline p includes N3
  - 12 per-patient runs (4 patients × 3 methods) read against one 0.01 bar, and how to state that without claiming more
  - agreement or heterogeneity across patients, and patients with `too_few_clones`
  - where the pooled diagnostic disagrees with the per-patient runs: read as a sign of patient structure, never as a finding
  - Δη and p_Δη against the baseline
  - concentration: effective clones, purity without the largest clone, and donors scored
  - agreement between m = 2 and m = 3
  - coupling after Benjamini–Hochberg, read as relatedness only
  - stability medians against minimums
  - cell-cycle warnings

  Record each conclusion in HANDOFF with its action: none, reported in the write-up, or a reporting bug fixed in this PR (then re-run Step 5). Never change the declaration, the seeds or a threshold because of a result.
- [ ] **Step 7: C1–C4.** Commit message: `Compute the case-study numbers per patient for all three methods, with a pooled diagnostic, and record them`.

---

### Task 13: Write up the case study

**Branch:** `pivot/NN-case-study-writeup` · **Model:** opus — interpretation where overclaiming is the main risk.

**Files:**
- Create: `docs/case-studies/sheih-2020.md` (or `cao-2025.md`), `tests/unit/test_case_study_writeup.py`

**Interfaces:**
- Consumes: the results JSON, `render_markdown` and `replace_section` (Task 12); Task 9's note; Task 10's declaration.

**Rules for the text:**
- The biology comes from Sheih et al. 2020 and is cited by section or figure. Anything not read in the paper is marked `(unverified)`.
- The clonal scores are credited to the work they derive from: STARTRAC (Zhang et al. 2018), LRT (Xie et al. 2023) and CoSpar (Wang et al. 2022), with the DOIs from spec §10.
- No sentence claims direction from clones. Sharing says that states are related, never which came first. Pseudotime direction comes from the declared root, and the text says so.
- If a method's verdict is `weak` or `not_supported`, that leads its paragraph: "doesn't fit a tree" is a finding.
- A number in prose appears exactly as it does in the generated section.
- The sections are:
  - Question
  - Data: paper, accession, licence, patients, time points and cells
  - Declaration: a link to the thresholds section and its merge commit; the clone definition, namespace, strata and root
  - What the paper reports
  - Results: the generated section, with the per-patient tables first and the pooled diagnostic after them
  - Reading:
    - per patient and method: the verdict, clonal agreement against the nulls and the baseline, coupling, stability and cell cycle
    - then what holds across patients and what does not
    - the pooled run appears only as a diagnostic of patient structure, and every prose sentence about it calls it a diagnostic
  - Limits:
    - per-patient runs hold about a quarter of the cells each
    - lumenCell has no batch correction, which is why pooling is not the headline
    - unpaired chains, if the declaration used them
    - only four patients
    - the manufactured infusion product
  - Reproduce: the commands

- [ ] **Step 1: Write the failing tests:**

```python
DOC, START, END = (
    REPO / "docs/case-studies/sheih-2020.md",
    "<!-- case-study:start -->",
    "<!-- case-study:end -->",
)
DIRECTION = (
    "give rise",
    "gives rise",
    "gave rise",
    "differentiat",
    "precursor",
    "derived from",
    "become",
    "became",
    "turn into",
    "turns into",
    "->",
    "→",
)


def test_generated_section_matches_the_committed_results() -> None:  # Review Focus 5
    latest = sorted((REPO / "benchmarks/case_studies/results").glob("sheih2020-*.json"))[-1]
    section = DOC.read_text(encoding="utf-8").split(START, 1)[1].split(END, 1)[0]
    assert section.strip() == render_markdown(json.loads(latest.read_text(encoding="utf-8"))).strip()


CLONAL_WORDS = ("clon", "tcr", "receptor", "shar", "coupl")


def test_no_sentence_about_clones_claims_direction() -> None:
    head, rest = DOC.read_text(encoding="utf-8").split(START, 1)
    for sentence in re.split(r"(?<=[.!?])\s+", head + rest.split(END, 1)[1]):
        if any(w in sentence.lower() for w in CLONAL_WORDS):
            assert not any(p in sentence.lower() for p in DIRECTION), sentence


def test_the_paper_the_data_and_the_methods_are_cited() -> None:
    text = DOC.read_text(encoding="utf-8")
    assert "10.1038/s41467-019-13880-1" in text and "GSE125881" in text  # Sheih et al. 2020 and its data
    for name, doi in (
        ("STARTRAC", "10.1038/s41586-018-0694-x"),
        ("LRT", "10.1371/journal.pcbi.1011300"),
        ("CoSpar", "10.1038/s41587-022-01209-1"),
    ):  # spec §10, checked
        assert name in text and doi in text, name


def test_per_patient_results_lead_and_pooled_is_only_a_diagnostic() -> None:
    text = DOC.read_text(encoding="utf-8")
    assert text.index("Per patient (primary)") < text.index("Pooled diagnostic")
    head, rest = text.split(START, 1)
    for sentence in re.split(r"(?<=[.!?])\s+", head + rest.split(END, 1)[1]):
        if "pooled" in sentence.lower():
            assert "diagnostic" in sentence.lower(), sentence
```

- [ ] **Step 2: Run** `uv run pytest tests/unit/test_case_study_writeup.py -q`. Expected: FAIL (no document).
- [ ] **Step 3: Write** the document, and fill the generated section with `replace_section`.
- [ ] **Step 4: Run** `data:validate-data` on the document. The conclusions must follow from the numbers, the limits must be stated, and the claims must be marked as verified or not. Fix what it finds.
- [ ] **Step 5: Run** the tests. Expected: PASS.
- [ ] **Step 6: C1–C4.** Commit message: `Write up the case study from its committed numbers`.

---

### Task 14: Security review of phase 3

**Branch:** `pivot/NN-phase3-security` · **Model:** sonnet — runs the review and fixes small findings test-first.

**Files:**
- Modify: whatever the findings require, and `docs/superpowers/HANDOFF.md`

- [ ] **Step 1: Invoke the `security-review` skill** on the diff from the "phase 3 base" commit (recorded by Task 1) to `main`. Point it at:
  - downloaded content reaching paths or a shell: tar member names, barcodes, metadata column names and values, CSV headers
  - the in-memory tar reader
  - temporary and derived files in the datasets directory
  - certifi and checksum enforcement
  - `LUMENCELL_DATASETS_DIR`
  - the `bench fetch` confirmation
  - the inspection scripts' `-I` runs
  - the case-study script's output path
- [ ] **Step 2: Fix each finding.** Write a failing test that shows it, fix it, and see the test pass. A finding too large for this PR gets its own task in HANDOFF instead.
- [ ] **Step 3: Add a "Security review (phase 3)" section to HANDOFF.** It lists every finding and says whether it was fixed (with its test) or recorded, and why.
- [ ] **Step 4 (background): Run the Done When commands** at the top of this plan, and record each result in HANDOFF.
- [ ] **Step 5: C1–C4.** Commit message: `Resolve or record the phase 3 security review findings`.
