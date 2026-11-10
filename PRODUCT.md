# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Computational biologists working with single-cell RNA-seq: comfortable with `.h5ad`, scanpy and the
command line. They launch lumenCell locally on their own data (`lumencell launch data.h5ad`), at a
desk, often for long sessions, and judge every result critically before it goes into a figure or a
paper. They know what pseudotime, PAGA, diffusion maps and UMAP are; they do not need these
explained, and they distrust tools that hide how a result was made.

## Product Purpose

An interactive explorer for pseudotime and lineage trajectories. It prepares the data in the
background, asks where the trajectory starts, runs a trajectory method, and lets the user inspect
lineages, gene trends along pseudotime and how far to trust the result, in the browser. Success is a
researcher who can tell, from the app alone, what the trajectory is, which genes change along it,
and whether the result is sound enough to use.

## Positioning

lumenCell refuses the shortcuts that made its predecessor wrong:
- **The root is never guessed.** Pseudotime is measured from a start the user chooses; suggestions
  come with their evidence and are never applied for them.
- **Inference never runs on the UMAP.** It runs in PCA / diffusion space; curves on any 2-D
  embedding are projections, labelled as such, with a per-lineage fidelity check.
- **Accuracy is measured, not assumed.** Every method is benchmarked against known trajectories, and
  a method that misses any gate is marked experimental in the app.
- **The user's file is never modified.** Results live in a cache; `lumencell export` writes a copy.

## Operating Context

- A local Python server (`lumencell launch`) serves the web app on 127.0.0.1; one dataset per
  server. The app must work offline (no CDN fonts or assets).
- The workflow is ordered: prepare (background stages) → choose a root → run a method (PAGA + DPT,
  Palantir, scFates) → inspect (lineages, pseudotime window, trust panel, gene trends and the top
  temporal genes heatmap, run comparison).
- Datasets range from a few thousand to about 100k cells; per-cell arrays arrive as binary (LCB1)
  and are drawn with deck.gl on the GPU.
- Desktop browsers on laptops and external monitors; not a phone product.

## Capabilities and Constraints

- Embeddings: whatever the file stores (UMAP, t-SNE, FA, ...), a computed UMAP, and the diffusion
  components (DC1/DC2, and DC1–DC3 in 3-D), the space the trajectory is inferred in.
- Colour cells by annotation, gene expression, pseudotime, lineage, or disagreement between runs;
  lasso selection with a summary; focus a lineage; filter by a pseudotime window; rename lineages.
- Every method is currently flagged experimental; the app must say so wherever a method is chosen
  or a run is judged.
- Terminal names and method parameters cannot yet be set in the app (the API supports them).
- Reference benchmark scores are not yet shown in the trust panel.
- Terminology: "root", "lineage", "pseudotime", "milestone", "shared trunk", "unassigned",
  "fidelity", "experimental". Colour is spelt the British way in the UI copy.

## Brand Commitments

- Name: lumenCell (lowercase l, capital C).
- Voice: plain, precise and honest about uncertainty; states what was done and how sure it is,
  never sells.

## Evidence on Hand

- Benchmark results in the README (synthetic and real datasets), generated from
  `benchmarks/results/`. The app must not invent scores, datasets or claims beyond these.
- No logo, illustrations or photography exist.

## Product Principles

1. Never decide for the user what only the user can know (the root, terminal states).
2. Show provenance and doubt next to the result, not behind it.
3. The plot is the product: everything else frames it and gets out of its way.
4. Fast and precise over decorative: motion only where it explains a change in the data or the view.

## Accessibility & Inclusion

- Colour-blind-aware categorical palette (Okabe-Ito first) and a perceptually uniform sequential
  scale (viridis) for data; these are data encodings, not decoration.
- Every control keyboard-reachable and named; reduced-motion respected.
