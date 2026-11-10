---
name: lumenCell
description: A darkfield imaging stage for pseudotime and lineage trajectories in single-cell data.
colors:
  field: "#06080b"
  panel: "#0c0f14"
  raised: "#151a21"
  control: "#141920"
  control-hover: "#1b212a"
  line: "#1d242e"
  line-strong: "#2b3440"
  ink: "#e6eaf0"
  ink-2: "#a3acbb"
  ink-3: "#7f8899"
  live: "#3ddc84"
  live-ink: "#03140a"
  focus: "#5b8cff"
  danger: "#ff5a76"
  warn: "#f2b544"
  brightfield-field: "#e6eaf0"
  brightfield-panel: "#ffffff"
  brightfield-ink: "#10141a"
  brightfield-live: "#12804a"
  brightfield-focus: "#2f5fd9"
typography:
  stamp:
    fontFamily: "Atkinson Hyperlegible Mono Variable, ui-monospace, monospace"
    fontSize: "1.75rem"
    fontWeight: 600
    lineHeight: 1
    letterSpacing: "-0.03em"
  title:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.9375rem"
    fontWeight: 600
    lineHeight: 1.375rem
    letterSpacing: "-0.01em"
  body:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.8125rem"
    fontWeight: 400
    lineHeight: 1.25rem
  label:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.75rem"
    fontWeight: 500
    lineHeight: 1.125rem
  data:
    fontFamily: "Atkinson Hyperlegible Mono Variable, ui-monospace, monospace"
    fontSize: "0.75rem"
    fontWeight: 400
    lineHeight: 1.125rem
    letterSpacing: "-0.01em"
    fontFeature: "tnum"
rounded:
  control: "8px"
  inner: "6px"
  panel: "16px"
  pill: "9999px"
spacing:
  gap: "12px"
  panel-x: "16px"
  section-top: "20px"
components:
  button-primary:
    backgroundColor: "{colors.live}"
    textColor: "{colors.live-ink}"
    rounded: "{rounded.control}"
    height: "32px"
    padding: "0 12px"
  button-secondary:
    backgroundColor: "{colors.control}"
    textColor: "{colors.ink}"
    rounded: "{rounded.control}"
    height: "32px"
    padding: "0 12px"
  button-secondary-hover:
    backgroundColor: "{colors.control-hover}"
  button-ghost:
    textColor: "{colors.ink-2}"
    rounded: "{rounded.control}"
    height: "28px"
    padding: "0 10px"
  input:
    backgroundColor: "{colors.control}"
    textColor: "{colors.ink}"
    rounded: "{rounded.control}"
    height: "32px"
  panel:
    backgroundColor: "{colors.panel}"
    rounded: "{rounded.panel}"
  popover:
    backgroundColor: "{colors.raised}"
    textColor: "{colors.ink}"
    rounded: "12px"
  tag-warn:
    backgroundColor: "{colors.warn}"
    textColor: "{colors.warn}"
    rounded: "{rounded.pill}"
    height: "20px"
---

# Design System: lumenCell

## Overview

**Creative North Star: "The Darkfield Time-Lapse"**

lumenCell is an imaging stage. The plot is a black darkfield, the cells are the light on it, and pseudotime is the stage's clock. Everything else (toolbar, panels, sheets) frames the field and steps out of its way. The world is fluorescence time-lapse microscopy, read for its discipline rather than its gear: channel colours that each own one job, an acquisition stamp in the corner, and a result you watch develop rather than read off a static chart.

It is a working tool for computational biologists, so density is high and every control is a standard one (Base UI underneath). Character comes from four things only: the darkfield palette, Atkinson Hyperlegible type, a tight scale with tabular numbers, and the signature move, which is pseudotime played as a time-lapse across the field. Motion explains changes in the data or the view and nothing else.

**Key Characteristics:**
- Full-bleed darkfield stage with floating panels; the camera centres the cells in the area the panels leave free.
- Channel colours with fixed jobs; data keeps Okabe-Ito and viridis.
- One shared 0–1 pseudotime axis across the timeline, lineage bars, trends and heatmap.
- Numbers in Atkinson Hyperlegible Mono with tabular figures; they change in place.
- Brightfield is the light variant, chosen from the toolbar.

## Colors

A near-black field under two cool panel plies, lit by four channel colours that each mean one thing.

### Primary
- **GFP Green** (live): the one filled action colour. Play, the primary action of a step, the playback window, finished steps and running progress.
- **DAPI Blue** (focus): selection, keyboard focus rings, the lasso, the unseen-run marker, highlighted list items.

### Tertiary
- **mCherry** (danger): errors only.
- **Amber** (warn): warnings and the experimental flag.

### Neutral
- **Darkfield** (field): the stage under the cells. Brightfield replaces it with a cool grey (brightfield-field), a clear step below its white panels.
- **Panel Slate** (panel) and **Raised Slate** (raised): the first and second plies.
- **Control** and **Control Hover**: input, select and secondary-button fills.
- **Hairline** (line, line-strong): panel borders and dividers.
- **Ink, Ink 2, Ink 3**: primary text, secondary text, and quiet labels and captions (each meets 4.5:1 on its ground).

### Named Rules
**The Channel Rule.** Each channel colour owns one job: green is live and primary, blue is selection and focus, amber is caution, mCherry is error. A colour never borrows another's job, and none is used as decoration.

**The One Green Button Rule.** At most one filled green control is in view at a time. Once a run is on the field, re-running is a secondary button and Play is the only green fill.

**The Data Keeps Its Encoding Rule.** Cells and charts use Okabe-Ito for categories and viridis for values in both themes. The UI palette never tints the data.

## Typography

**Body Font:** Atkinson Hyperlegible Next (self-hosted variable, with system sans fallback)
**Mono Font:** Atkinson Hyperlegible Mono (self-hosted variable)

**Character:** Chosen because gene symbols (Il1b, Cd34) and identifiers (0/O, 1/l/I) must never be misread; the mono shares the same disambiguated letterforms for measurements.

### Hierarchy
- **Stamp** (600, 1.75rem, line-height 1, mono): the pseudotime acquisition stamp in the field's top-left corner. Ink-3 at rest, ink while a window is set or playing.
- **Title** (600, 0.9375rem): panel section headings and the wordmark.
- **Body** (400, 0.8125rem): the default for panels and controls.
- **Label** (500, 0.75rem): control labels, tabs, segmented options, hints.
- **Data** (400, 0.75rem, mono, tabular figures): counts, pseudotime values, fidelity, IDs and parameters.

### Named Rules
**The Mono Means Measurement Rule.** Mono is only for numbers, IDs, parameters and the stamp, never as a "technical" costume for prose or headings.

**The No Eyebrow Rule.** Section titles carry their own weight. There are no uppercase tracked labels above headings.

## Layout

The stage fills the window under a 48px toolbar. Two panels float over it: the Cells panel on the left (264px) and the workflow inspector on the right (336px). The gene sheet sits between them at the bottom, at one of three heights (a 44px bar, a peek capped at 40% of the stage, or about 55% when full). The timeline (52px) docks above the sheet. All of these sit 12px from the stage edges and from each other.

The camera's view padding follows the open panels, so the cells centre in the free area and glide to recentre when a panel moves. Below 1200px of stage width the Cells panel starts closed. Below 900px both panels start closed and cover the field when opened, and the toolbar wraps onto two rows. Inside panels, sections are separated by space (20px above each heading) rather than boxes.

## Elevation & Depth

Depth is three fixed plies: the field, the panels on it, and popovers (selects, comboboxes, tooltips) above them. Panels separate from the field by tone and a 1px hairline only; popovers alone cast a shadow. In the 3-D view, depth is carried by the cells themselves: a shader dims and thins cells with distance from the camera (the furthest keep 30% of their light), and a DC1/DC2/DC3 triad turns with the orbit.

### Shadow Vocabulary
- **Ply shadow** (`0 12px 32px -8px rgb(0 0 0 / 0.55), 0 2px 6px rgb(0 0 0 / 0.35)` in darkfield): popovers and tooltips only.

### Named Rules
**The Three Plies Rule.** Field, panel, popover. Nothing else is raised, and only the top ply has a shadow.

## Shapes

Controls have gently rounded 8px corners and their inner pieces 6px. Floating panels and the sheet use 16px, and popovers 12px. Pills are reserved for small chips and tags. Lines are 1px hairlines; there are no coloured side borders and no nested cards.

## Components

### Buttons
- **Shape:** gently rounded (8px), 32px tall; 28px in dense rows.
- **Primary:** GFP green fill with dark ink text, one per view.
- **Secondary:** control fill, ink text; hover lifts to control-hover.
- **Ghost:** text only in ink-2, control fill on hover.
- **Press:** every pressable scales to 0.97 in 120ms on a strong ease-out. Hover effects only apply to fine pointers.

### Segmented control (embedding switch)
- One-of-n choice on a control-filled track. The selected segment's raised background slides between options on a critically damped spring.

### Inputs and selects
- **Style:** control fill, 8px radius, 32px tall, no border. The focus ring is a 2px DAPI blue outline.
- **Popups:** the raised ply with the ply shadow. They scale from their trigger (0.97 to 1, 180ms ease-out) and leave faster than they arrive (120ms).

### Panels
- Panel slate, 16px radius, 1px hairline, no shadow. They slide in and out on a spring; a keyboard toggle (`[` and `]`) moves them instantly.

### The field (signature component)
- **Morph:** switching embeddings moves every cell to its new position on the GPU (700ms ease-in-out), with the camera flying alongside. DC 3D lifts the diffusion map into its third component.
- **Time-lapse:** a 0.14-wide pseudotime window sweeps from the root to the tips over 7s. Cells inside are fully lit and those it has passed hold an opaque afterglow; in darkfield, a bright exposure front leads the window. Cells ahead stay dark.
- **Draw-on:** when a run first lands, its curves draw from the root to the tips (900ms) and the root pulses once.
- **Highlight:** pointing at a legend entry, a lineage, a root suggestion or a time on a trend chart lights those cells and dims the rest through a 220ms colour transition.
- **Labels:** lineage labels sit on a field-coloured outline halo, never on a box.

### Workflow stepper
- Prepare, Root, Run, Inspect. Done steps show a green check and the current step a blue number. Steps change with a short rise and fade (240ms).

## Do's and Don'ts

### Do:
- **Do** keep the field full-bleed and the panels floating; centre the cells with the view padding, never by shrinking the canvas.
- **Do** put every pseudotime graphic on the one shared 0–1 axis.
- **Do** keep motion to what explains a change: the morph, the time-lapse, the draw-on and highlights. UI transitions stay at or under 240ms, and exits are faster than entrances.
- **Do** keep a change marked until it has been seen; a finished run keeps its blue dot until it is opened.
- **Do** say a method is experimental wherever it is chosen or judged.

### Don't:
- **Don't** box each panel section or put uppercase tracked labels above headings.
- **Don't** use a channel colour outside its job, or tint the data with UI colours.
- **Don't** animate keyboard-initiated panel toggles.
- **Don't** let semi-transparent layers stack into brightness on the field; dimmed cells are drawn opaque and pre-mixed toward the field.
- **Don't** load fonts or assets from a CDN; the app runs offline.
