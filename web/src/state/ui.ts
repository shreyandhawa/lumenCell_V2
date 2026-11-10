/** UI-only state. Server data lives in TanStack Query, never here. */
import { create } from "zustand";

export type ColorBy =
  | { kind: "obs"; column: string }
  | { kind: "gene"; gene: string }
  | { kind: "pseudotime" }
  | { kind: "lineage" }
  | { kind: "disagreement"; a: string; b: string };

/** A transient emphasis on the field while the pointer rests on something that names cells. */
export type Highlight =
  | { kind: "category"; column: string; value: string }
  | { kind: "lineage"; id: number }
  | { kind: "pseudotime"; at: number };

export type SheetDetent = "bar" | "peek" | "full";
export type Theme = "dark" | "light";

export interface HoverPoint {
  index: number;
  x: number;
  y: number;
}

export interface UiState {
  embedding: string | null;
  colorBy: ColorBy | null;
  runId: string | null;
  /** Inclusive pseudotime window; [0, 1] means no filter. */
  ptRange: [number, number];
  focusLineage: number | null;
  /** 1 = selected. null = no selection. */
  selection: Uint8Array | null;
  selectionCount: number;
  lassoActive: boolean;
  hovered: HoverPoint | null;
  highlight: Highlight | null;
  /** Genes shown in the trends panel. */
  trendGenes: string[];
  playing: boolean;
  panels: { cells: boolean; inspector: boolean };
  /** False when the last panel toggle came from the keyboard: it should not animate. */
  panelsAnimate: boolean;
  sheet: SheetDetent;
  /** Live drag offset of the gene sheet, in px (negative = taller); null when not dragging. */
  sheetDrag: number | null;
  theme: Theme;
  /** Runs the user has opened; a run that finishes elsewhere stays marked until then. */
  seenRuns: string[];

  setEmbedding: (name: string) => void;
  setColorBy: (c: ColorBy | null) => void;
  setRun: (id: string | null) => void;
  setPtRange: (r: [number, number]) => void;
  setFocusLineage: (id: number | null) => void;
  setSelection: (mask: Uint8Array | null) => void;
  setLasso: (on: boolean) => void;
  setHovered: (p: HoverPoint | null) => void;
  setHighlight: (h: Highlight | null) => void;
  addTrendGene: (gene: string) => void;
  removeTrendGene: (gene: string) => void;
  setPlaying: (on: boolean) => void;
  togglePanel: (name: "cells" | "inspector", animate?: boolean) => void;
  setSheet: (d: SheetDetent) => void;
  setSheetDrag: (px: number | null) => void;
  setTheme: (t: Theme) => void;
  markRunSeen: (id: string) => void;
}

export const MAX_TREND_GENES = 6;

const THEME_KEY = "lumencell.theme";

function savedTheme(): Theme {
  try {
    return globalThis.localStorage?.getItem(THEME_KEY) === "light" ? "light" : "dark";
  } catch {
    return "dark";
  }
}

export const useUi = create<UiState>((set) => ({
  embedding: null,
  colorBy: null,
  runId: null,
  ptRange: [0, 1],
  focusLineage: null,
  selection: null,
  selectionCount: 0,
  lassoActive: false,
  hovered: null,
  highlight: null,
  trendGenes: [],
  playing: false,
  panels: { cells: true, inspector: true },
  panelsAnimate: true,
  sheet: "bar",
  sheetDrag: null,
  theme: savedTheme(),
  seenRuns: [],

  setEmbedding: (embedding) => set({ embedding }),
  setColorBy: (colorBy) => set({ colorBy }),
  setRun: (runId) =>
    set((s) => ({
      runId,
      focusLineage: null,
      ptRange: [0, 1],
      playing: false,
      seenRuns: runId && !s.seenRuns.includes(runId) ? [...s.seenRuns, runId] : s.seenRuns,
    })),
  setPtRange: (ptRange) => set({ ptRange }),
  setFocusLineage: (focusLineage) => set({ focusLineage }),
  setSelection: (mask) => {
    let count = 0;
    if (mask) for (let i = 0; i < mask.length; i++) count += mask[i]!;
    set({ selection: mask && count > 0 ? mask : null, selectionCount: count });
  },
  setLasso: (lassoActive) => set({ lassoActive }),
  setHovered: (hovered) => set({ hovered }),
  setHighlight: (highlight) => set({ highlight }),
  addTrendGene: (gene) =>
    set((s) => (s.trendGenes.includes(gene) ? s : { trendGenes: [...s.trendGenes, gene].slice(-MAX_TREND_GENES) })),
  removeTrendGene: (gene) => set((s) => ({ trendGenes: s.trendGenes.filter((g) => g !== gene) })),
  setPlaying: (playing) => set({ playing }),
  togglePanel: (name, animate = true) =>
    set((s) => ({ panels: { ...s.panels, [name]: !s.panels[name] }, panelsAnimate: animate })),
  setSheet: (sheet) => set({ sheet, sheetDrag: null }),
  setSheetDrag: (sheetDrag) => set({ sheetDrag }),
  setTheme: (theme) => {
    try {
      globalThis.localStorage?.setItem(THEME_KEY, theme);
    } catch {
      // storage unavailable (private window): the choice lasts for this page only
    }
    set({ theme });
  },
  markRunSeen: (id) => set((s) => (s.seenRuns.includes(id) ? s : { seenRuns: [...s.seenRuns, id] })),
}));
