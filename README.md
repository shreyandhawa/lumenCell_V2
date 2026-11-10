# lumenCell

Interactive explorer for pseudotime and lineage trajectories in single-cell RNA-seq data.

`lumencell launch data.h5ad` prepares the data in the background, asks where the
trajectory starts, runs a trajectory method, and lets you inspect lineages, gene
trends along pseudotime and how far to trust the result, in the browser.

> Early development (M1). Usable on real data, but every method is still labelled
> experimental: see the benchmark below.

## Principles

- **The root is never guessed.** Pseudotime is measured from a start you choose
  (`--root cluster:<column>=<value>`, a cell, or a root stored in the file).
  Suggestions come with their evidence and are never applied for you.
- **Inference runs in PCA / diffusion space, never on the UMAP.** Curves drawn on a
  2-D embedding are projections, labelled as such, with a fidelity check.
- **Accuracy is measured, not assumed.** Every method is benchmarked against known
  trajectories, and a method that misses any gate is marked experimental in the app.
- **Your file is never modified.** Results live in a cache; `lumencell export`
  writes them into a copy.

## Quick start

Not on PyPI yet. From a clone, with Python 3.13, [uv](https://docs.astral.sh/uv/) and Node 22+:

```bash
make install     # Python (all methods) and web dependencies
make web-dist    # build the web app into the package
uv run lumencell launch your-data.h5ad --root "cluster:clusters=Ductal"
```

`X` (or `layers["counts"]` / `.raw`) should hold raw counts, or log-normalised values.

Other commands: `lumencell prepare` (precompute), `lumencell roots` (suggestions),
`lumencell run` (run from the terminal), `lumencell export` (write results into a
copy of the file), `lumencell synth` (synthetic data with a known trajectory).

## Methods

| Method | Install | What it is |
|---|---|---|
| PAGA + DPT | core | Diffusion pseudotime with lineages along the PAGA cluster tree (scanpy) |
| Palantir | `[palantir]` | Waypoint pseudotime and fate probabilities from a Markov chain (Setty et al. 2019) |
| scFates | `[scfates]` | A principal tree fitted in diffusion space; the drawn curve is the fitted tree |

The `palantir` and `scfates` extras pull in GPL-licensed dependencies (mellon,
elpigraph); lumenCell itself is MIT and only imports them.

## Benchmark

Synthetic data (Splatter-style counts from a known tree, no embedding given) runs
in CI; real datasets are a release gate (`lumencell bench fetch all`, then
`make bench-real`). Gates and every change to them are in
[`benchmarks/thresholds.toml`](benchmarks/thresholds.toml).

<!-- bench:start -->
**Synthetic** (4 topologies × 3 seeds, `realistic` preset; baseline from commit `2071de6`). Ranges over all cases:

| Method | Pseudotime (global) | Pseudotime (per lineage) | Branch accuracy | Gates failed |
|---|---|---|---|---|
| PAGA + DPT | 0.86–0.99 | 0.95–0.99 | 0.81–0.99 | 7 of 12 cases |
| Palantir | 0.90–0.99 | 0.91–0.99 | 0.90–0.99 | 2 of 12 cases |
| scFates | 0.87–0.99 | 0.98–0.99 | 0.97–1.00 | 2 of 12 cases |

**Real datasets** (`benchmarks/results/real-2026-11-10.json`, commit `abe3f73`). Spearman is against stage order (embryo time for gastrulation):

| Method | bonemarrow | gastrulation_erythroid | pancreas |
|---|---|---|---|
| PAGA + DPT | Spearman 0.73, branch 0.92, terminal 0.95 (pass) | Spearman 0.75, branch –, terminal 0.94 (pass) | Spearman 0.95, branch 0.63, terminal 0.93 (fail: branch_accuracy, lineage_recall) |
| Palantir | Spearman 0.75, branch 0.64, terminal 0.84 (fail: branch_accuracy) | Spearman 0.72, branch –, terminal 0.59 (fail: terminal_share) | Spearman 0.94, branch 0.84, terminal 0.94 (pass) |
| scFates | Spearman 0.78, branch 0.94, terminal 0.96 (pass) | Spearman 0.80, branch –, terminal 0.61 (fail: terminal_share) | Spearman 0.94, branch 0.75, terminal 0.80 (fail: branch_accuracy) |
<!-- bench:end -->

## Development

```bash
make install   # uv sync --all-extras + npm ci (the palantir/scfates extras are optional)
make ci        # lint, types, tests, API contract check, web build, end-to-end test
make dev-server   # API on :8765 with a synthetic dataset; then `make dev-web` for the app on :5173
```
