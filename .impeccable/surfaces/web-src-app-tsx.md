---
version: 1
slug: "web-src-app-tsx"
primary_target: "web/src/App.tsx"
related_targets: []
---

# lumenCell explorer (web app)

## Scope and mode

The whole single-page explorer served by `lumencell launch` (`web/src`): toolbar, cell field, Cells
panel, workflow inspector, timeline and gene sheet. Mode: **Operate**. Built code-led.

## Audience, job and constraints

Computational biologists on their own `.h5ad`, at a desk, in long sessions. Job: prepare → choose a
root → run a method → inspect lineages, gene trends and how far to trust the result. Up to ~100k
cells drawn on the GPU; must work offline; every method is experimental and the UI says so; the root
is never guessed; curves on 2-D embeddings are projections with a fidelity check.

## Chosen direction

Live imaging (fluorescence time-lapse microscopy). Memorable moment: pressing play sweeps pseudotime
across the field like a time-lapse, cells lighting up as the window passes through them, with a
mono timestamp in the field's corner. Second moment: switching embeddings morphs the same cells
into their new positions, and DC 3D lifts the diffusion map into its third component.

## Direction contract

THESIS: The field is an imaging stage and pseudotime is its clock: the result is watched developing,
not read off a static plot. Refuses the category default of a grey dashboard with boxed side panels
and one accent colour.

OWN-WORLD: Darkfield black field (#06080B) under two panel plies (#0C0F14, then #151A21 for popovers,
the only ply with a shadow); cool-tinted ink. Channel colours each own a job: GFP green for live,
active and primary; DAPI blue for selection, focus and links; mCherry for errors; amber for warnings
and experimental. Data keeps Okabe-Ito and viridis. Atkinson Hyperlegible Next for the UI (gene
symbols like Il1b, Cd34 and 0/O, 1/l/I stay unambiguous); its Mono for measurements. A brightfield
(light) variant is offered as a toggle. Raised: numbers change in place inside fixed columns (from
split-flap); a change stays marked until seen (boarding pass); no box per section, hierarchy from
the type scale (metro tiles); three fixed plies rank everything (acetate); every pseudotime graphic
shares one to-scale 0–1 axis (labanotation).

STORY: The researcher sees their cells first, sets a root and runs a method in the inspector; the
curves draw from the root to the tips; they press play and watch differentiation sweep the field;
trust signals sit beside the result, and hovering any lineage, cluster or trend lights its cells.

FIRST VIEWPORT: A 48px toolbar (wordmark, dataset and size, embedding switch UMAP · DC 2D · DC 3D,
colour-by, panel toggles). Below it the full-bleed black field with the cells; a 260px Cells panel
floats on the left (legend, selection) and a 320px workflow inspector on the right (stepper:
Prepare, Root, Run, Inspect), where the primary action lives. The timeline docks bottom-centre of
the field with the mono pseudotime stamp top-left; the gene sheet rests as a 44px bar under it.

FORM: Live imaging (fluorescence time-lapse), position 4 of 7 on the ordered list; seed:65393df55f628888

FINISH: unreviewed and undocumented is unfinished; this build ends with the finish review, the verdict, DESIGN.md, and every shipping raster carrying its provenance

## Unresolved

- Whether a channel-style single-hue LUT should be offered beside viridis for gene expression (not
  in this build; viridis stays the encoding).
