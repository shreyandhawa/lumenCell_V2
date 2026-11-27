# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Two audiences, weighted equally:

- **Researchers** (computational biologists) who run `lumencell launch` on their own `.h5ad` and stay in the app for long analysis sessions: picking a root, running methods, reading gene trends and the evidence.
- **First-time viewers** (PIs, reviewers, portfolio visitors) who open `lumencell demo` or the README screenshots and need to understand what they are looking at within a minute.

## Product Purpose

Explore pseudotime and lineage trajectories in single-cell RNA-seq data, and show how far to trust the result. Success is a researcher who can say which lineages the data supports, along what order, and with what evidence.

## Positioning

- The root is never guessed. Pseudotime is measured from a start the user chooses; suggestions come with their evidence and are never applied for them.
- Inference runs in PCA or diffusion space, never on the UMAP. Curves drawn on a 2-D embedding are labelled as projections.
- Every method is benchmarked against known trajectories, and a method that misses any gate is marked experimental in the app.
- Clonal evidence and null baselines are shown beside the result, in the Evidence sheet.
- The input file is never modified.

## Operating Context

The workflow is Prepare, then Pick a root, then Run (PAGA + DPT, Palantir or scFates), then Inspect: colour by observation or gene, follow gene trends along pseudotime, focus a lineage, and read the Evidence. The address bar always names the view on screen, so a link reopens it. Export PNG saves the plot with its legend and provenance.

## Capabilities and Constraints

- 2-D, 3-D and split views on deck.gl, at 100k to 500k cells; frame times are measured in `benchmarks/results/frame-times-2026-11-25.json`.
- Canvas colours are set from TypeScript, so every theme value must reach both CSS and the canvas.
- Runs locally against the `lumencell` Python server; no accounts or hosting.

## Evidence on Hand

- Screenshots: `docs/images/*.png`.
- Real datasets: Sheih 2020 (patient NHL-6), Daniel 2022, bone marrow.
- Benchmark and calibration results: `benchmarks/results/`.
- No testimonials, users or adoption figures exist; none are to be invented.

## Product Principles

1. Never decide for the researcher: choices like the root are theirs, with evidence beside them.
2. Show trust next to the result, failures included.
3. The next step is always obvious, for a first-time viewer as much as an expert.
4. A view on screen is a view you can share.

## Accessibility & Inclusion

Colour-blind-safe data colours, WCAG AA text contrast in both light and dark modes, keyboard access to every control, and reduced motion respected.
