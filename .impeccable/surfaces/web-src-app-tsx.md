---
version: 1
slug: "web-src-app-tsx"
primary_target: "web/src/App.tsx"
related_targets: ["web/src/components/Sidebar.tsx","web/src/components/Scatter.tsx","web/src/components/trajectory/TrustPanel.tsx"]
---

## Scope

The lumenCell explorer: the single-page app in `web/` (sidebar of steps, the 2-D/3-D stage, pseudotime bar and gene trends, the Evidence sheet). Mode: Operate.

## Audience and job

Researchers analysing their own data and first-time viewers, equally. Job: get from a prepared dataset to a trusted reading of its lineages. Pain points to fix: what to do next is unclear, the Evidence sheet is hard to read, view controls are hard to find, the screen is too dark and dense.

## Direction contract

THESIS: A trajectory reads like a metro map: you are at one station, and lines branch where lineages split. It refuses the category default of a black stage with one neon accent.

OWN-WORLD: A deep-blue sidebar shell (#0B3D91 light, #16305E dark) with white type; a near-white (#FAFAF7) or night-navy (#0E1116) stage; white or #171B22 panels; one signal-red "you are here" accent (#D7263D / #FF4D5E) for the root pin, the focused lineage and the current station only; primary buttons are ink (white on the blue line). Data keeps colour-blind-safe lines. Atkinson Hyperlegible Next and Mono, as on wayfinding signs.

STORY: A first-time viewer sees which station they are at and what the next stop is. A researcher reads the verdict before the details and finds every view control by its name.

FIRST VIEWPORT: Left, the blue sidebar: title, dataset, then the four steps as stations on one vertical line, the current one ringed in red with its primary action under it. Centre, the stage. Bottom, the pseudotime bar. Bottom right, one labelled view toolbar. Right, Evidence: verdict and warnings first, details folded.

FORM: Lineage lines, the user's pick, first of four shown (Lineage lines, Fate landscape, Sun print, Stain atlas, plus a standard); a pinned direction beats the roll, which dealt two wayfinding worlds (a guide map, a split-flap concourse) as challengers. seed:2f826f6659e7ea09

FINISH: unreviewed and undocumented is unfinished; this build ends with the finish review, the verdict, DESIGN.md, and every shipping raster carrying its provenance

## Memorable moment

The step rail as a metro line with a red "you are here" station.

## Unresolved

Whether the top of viridis and the yellow categorical colour hold up on the light stage; checked in the browser during the build.
