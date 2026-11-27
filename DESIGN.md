---
name: lumenCell
description: Lineage lines, a trajectory explorer read like a metro map
colors:
  line-blue: "#0b3d91"
  line-blue-raised: "#1a4fa6"
  line-blue-rule: "#3a68b8"
  night-line: "#16305e"
  night-line-raised: "#1f3d72"
  night-line-rule: "#34528a"
  on-line: "#ffffff"
  on-line-muted: "#c9d6ee"
  on-line-faint: "#8fa6d1"
  night-on-line: "#e8edf7"
  night-on-line-muted: "#9fb0ce"
  night-on-line-faint: "#6f83a8"
  you-are-here: "#d7263d"
  you-are-here-night: "#ff4d5e"
  platform-paper: "#fafaf7"
  sheet-white: "#ffffff"
  sheet-raised: "#f1f3f6"
  sheet-rule: "#d9dce3"
  ink: "#111418"
  ink-muted: "#5b6270"
  ink-faint: "#8c93a0"
  night-stage: "#0e1116"
  night-sheet: "#171b22"
  night-sheet-raised: "#1f2530"
  night-rule: "#2a303b"
  night-ink: "#eef1f6"
  night-ink-muted: "#9aa3b2"
  night-ink-faint: "#6b7382"
  caution-amber: "#b07a00"
  caution-amber-night: "#f2a900"
typography:
  headline:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "1.25rem"
    fontWeight: 600
    letterSpacing: "-0.025em"
  title:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "1rem"
    fontWeight: 600
  body:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.43
  small:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.8125rem"
    fontWeight: 400
  label:
    fontFamily: "Atkinson Hyperlegible Next Variable, ui-sans-serif, system-ui, sans-serif"
    fontSize: "0.75rem"
    fontWeight: 600
    letterSpacing: "0.06em"
  data:
    fontFamily: "Atkinson Hyperlegible Mono Variable, ui-monospace, SFMono-Regular, Menlo, monospace"
    fontSize: "0.8125rem"
    fontWeight: 400
rounded:
  control: "4px"
  field: "6px"
  sheet: "12px"
  station: "9999px"
spacing:
  inset: "12px"
  sheet: "16px"
  step: "10px"
  section: "16px"
components:
  button-primary:
    backgroundColor: "{colors.ink}"
    textColor: "{colors.sheet-white}"
    rounded: "{rounded.control}"
    height: "32px"
    padding: "0 12px"
  button-secondary:
    backgroundColor: "{colors.sheet-raised}"
    textColor: "{colors.ink}"
    rounded: "{rounded.control}"
    height: "32px"
    padding: "0 12px"
  view-tool:
    textColor: "{colors.ink}"
    rounded: "{rounded.field}"
    height: "32px"
    padding: "0 10px 0 8px"
  view-tool-on:
    backgroundColor: "{colors.sheet-raised}"
    textColor: "{colors.ink}"
  station-current:
    backgroundColor: "{colors.you-are-here}"
    textColor: "{colors.sheet-white}"
    rounded: "{rounded.station}"
    size: "22px"
  station-done:
    backgroundColor: "{colors.on-line}"
    textColor: "{colors.line-blue}"
    rounded: "{rounded.station}"
    size: "22px"
  field:
    backgroundColor: "{colors.sheet-raised}"
    textColor: "{colors.ink}"
    rounded: "{rounded.field}"
    height: "36px"
  sheet:
    backgroundColor: "{colors.sheet-white}"
    rounded: "{rounded.sheet}"
    padding: "16px"
---

# Design System: lumenCell

## Overview

**Creative North Star: "Lineage lines"**

A trajectory reads like a metro map. You are at one station, the line runs on to the next, and lines branch where lineages split. The sidebar is the line itself: a deep-blue band the four steps (Prepare, Pick a root, Run, Inspect) sit on as stations, with a red "you are here" on the step in hand. Everything else is a working screen in the transit authority's register: plain white or night-navy sheets, one legible sans, data set in a mono, and nothing decorative.

The world lends type, palette, density and one signature move. It never supplies the layout or the controls: buttons, selects, folds and toggles are standard web controls. The plot belongs to the data, so the canvas stays quiet: a near-white or night-navy stage that the scientific palettes sit on unchanged.

Both themes are first-class. Dark is the default; a saved choice wins, and on a first visit the system setting decides.

**Key Characteristics:**
- A coloured shell: the deep-blue sidebar carries the steps as stations on one line.
- One red, "you are here", and only for where you are.
- White or night-navy sheets on a quiet stage; depth from tone, not shadow.
- Atkinson Hyperlegible throughout, mono for numbers, ids and data.
- Scientific data colours, never restyled to match the world.

## Colors

A committed blue shell, a quiet neutral ground and one signal red, with amber kept for caution.

### Primary
- **Line Blue** (light) / **Night Line** (dark): the sidebar shell, the band the steps run along. Its raised step holds fields and buttons inside the shell; its rule step is the stretch of line still ahead.
- **On-line white and its muted and faint steps**: text, icons, done stations and the travelled stretch of line on the shell.

### Secondary
- **You-are-here Red** (light) / its brighter night step: the root pin, the focused lineage (in the list and its curve on the plot) and the current station.

### Neutral
- **Platform Paper** / **Night Stage**: the stage behind the plot and the window ground.
- **Sheet White** / **Night Sheet**, with their raised steps: the Evidence sheet, the pseudotime bar, the view tools, legends and pop-overs; raised for fields, pressed tools and hovered rows.
- **Ink, Ink Muted, Ink Faint** (and night steps): body text, supporting text, and disabled controls only.
- **Sheet Rule** / **Night Rule**: section rules and borders.
- **Caution Amber** / its night step: experimental, provisional and weak results, always beside a label and an icon. On the blue shell it takes the night step to stay visible.

### Named Rules
**The You-Are-Here Rule.** Red marks where you are and nothing else: the root, the focused lineage, the current step. It is never a button, a data colour, an error colour or decoration.

**The Line Rule.** Inside the sidebar the usual tokens point at the shell's blues (`.shell` in `web/src/index.css`), so any control placed there takes the line's colours without its own styling. Never hand-colour a sidebar control.

**The Data Stays Data Rule.** Cells are coloured by viridis (continuous) and an Okabe-Ito-based categorical set, colour-blind-safe. The theme never restyles them; near-reds are ordered last so no category reads as the red.

## Typography

**Body Font:** Atkinson Hyperlegible Next (with ui-sans-serif, system-ui)
**Data Font:** Atkinson Hyperlegible Mono (with ui-monospace, Menlo)

**Character:** a typeface made for legibility, as a wayfinding sign is; the mono sets numbers, ids, gates and cell counts so they line up and read as measurements.

### Hierarchy
- **Headline** (600, 1.25rem, tight tracking): the product name at the top of the shell.
- **Title** (600, 1rem): the current step's name and sheet titles.
- **Body** (400, 0.875rem): the default size of the whole interface, including the verdict's explanation.
- **Small** (400, 0.8125rem): supporting lines under checks, warnings, captions.
- **Label** (600, 0.75rem, 0.06em, uppercase): section heads (Verdict, Warnings, Checks, Lineages).
- **Data** (mono, 0.8125rem): values, thresholds, ids, counts, provenance.

### Named Rules
**The Mono Is Measurement Rule.** Mono is for numbers, ids and data only, never as a "technical" costume for prose.

## Layout

A full-window stage with three floating sheets inset 12px from the window: the 280px shell on the left, the 360px Evidence sheet on the right, and the pseudotime bar along the bottom between them, which takes the full width when Evidence is closed. The plot is centred in the space the sheets leave. Inside sheets: 16px padding, 16px between ruled sections, 10px above and below each step. Desktop is the target; the fixed sheet widths are not built for phones.

## Elevation & Depth

One elevation: the sheets that float over the stage (Evidence, the pseudotime bar and gene trends, the view tools and toggle, legends, pop-overs, the selection bar) carry a single soft shadow, because white on the near-white stage is too close in tone to separate on its own. Inside sheets depth is tonal: raised steps mark fields, pressed tools and hovered rows. The blue shell needs no shadow; its colour separates it.

### Shadow Vocabulary
- **Float** (light `box-shadow: 0 1px 2px rgb(17 20 24 / 0.06), 0 6px 20px rgb(17 20 24 / 0.08)`; dark `0 2px 10px rgb(0 0 0 / 0.35)`): every sheet floating over the stage, via `.lift`.

### Named Rules
**The One Lift Rule.** There is one shadow and it means "floats over the plot". Nothing inside a sheet gets its own.

## Shapes

Gently rounded sheets (12px), fields and tools (6px), buttons (4px). Stations and markers are full circles; the line is a 3px stroke with rounded branch stubs where lineages leave the trunk.

## Components

### Buttons
- **Shape:** slightly rounded (4px), 32px high.
- **Primary:** ink, as on a station sign: white with blue letters on the line, dark with white letters on a sheet; one per step. Disabled, it drops to the raised tone with faint text, so it never reads as ready.
- **Secondary:** the raised sheet tone with ink text.
- **Hover / Focus:** slight fade on hover; a 2px ink focus outline offset 2px. Every pressable thing scales to 0.97 on press.

### Stations (signature component)
- **Done:** a filled on-line disc with a check.
- **Current:** a You-are-here Red disc with the step number, ringed in the line's ink, and `aria-current="step"`.
- **Ahead:** a hollow disc with a faint outline and the step number.
- **Line:** a 3px stroke through the stations, muted where travelled and faint where ahead; it ends at Inspect.
- **Folding:** a station behind you folds to its one-line status; its header reopens it. Only the current station is open, so its action is always in reach.
- **Aside:** each station's short status (done, the chosen root, the method, the lineage count) sits right of its name in muted sans and wraps rather than cuts off.

### Lineage list
- One trunk with a short branch per lineage ending in a small dot; the focused lineage's dot is the red station, matching its curve on the plot.

### View tools
- One labelled group at the plot's lower right: an icon and its word per tool (Lasso, Display, Export, Tour, Rotate, Reset), filled while on or open.

### Evidence sheet
- Verdict first (the amber mark and its explanation as body text), then warnings, then the checks the verdict reads, with left-out checks on one line. Clones, reference results, curve fidelity and provenance are folded; a fold may carry a one-line summary when its content would change a reader's trust (reference results show how many gates fail).

### Inputs / Fields
- **Style:** the raised tone, 6px radius, 36px high, no border; a 2px ink focus outline.

## Do's and Don'ts

### Do:
- **Do** put the red only where the user is: root, focused lineage, current step.
- **Do** add sidebar controls with the ordinary tokens and let `.shell` re-point them.
- **Do** keep numbers, ids and thresholds in the mono.
- **Do** keep a failure visible when folding detail (state the count in the fold's summary).
- **Do** check every new text pair for WCAG AA in both themes.

### Don't:
- **Don't** use red for errors or failed checks; they keep their icon and words in ink.
- **Don't** recolour cells to match the theme; data colours are scientific choices.
- **Don't** use Ink Faint for text anyone has to read; it is for disabled controls.
- **Don't** add a second shadow; inside sheets, use the raised tone.
- **Don't** make a button red; primary buttons are ink.
