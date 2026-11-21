/** UI-only state. Server data lives in TanStack Query, never here. */
import { create } from "zustand";

export type ColorBy =
  | { kind: "obs"; column: string }
  | { kind: "gene"; gene: string }
  | { kind: "pseudotime" }
  | { kind: "lineage" };

export interface UiState {
  embedding: string | null;
  colorBy: ColorBy | null;
  runId: string | null;
  focusLineage: number | null;
  /** Genes the user picked, shown first in the gene trends. */
  trendGenes: string[];
  theme: "dark" | "light";
  view: "2d" | "3d";
  /** The pseudotime cursor shared by the bar and the gene trends; null until the user plays or scrubs. */
  cursor: number | null;
  playing: boolean;
  sheetOpen: boolean;
  /** Point size as a multiple of the default for this many cells. */
  pointSize: number;
  /** Far cells fade toward the stage in 3-D. */
  fog: boolean;

  setEmbedding: (name: string) => void;
  setColorBy: (c: ColorBy | null) => void;
  setRun: (id: string | null) => void;
  setFocusLineage: (id: number | null) => void;
  addTrendGene: (gene: string) => void;
  setTheme: (t: "dark" | "light") => void;
  setView: (v: "2d" | "3d") => void;
  setCursor: (t: number | null) => void;
  setPlaying: (on: boolean) => void;
  setSheetOpen: (open: boolean) => void;
  setPointSize: (k: number) => void;
  setFog: (on: boolean) => void;
}

export const MAX_TREND_GENES = 4;

export const useUi = create<UiState>((set) => ({
  embedding: null,
  colorBy: null,
  runId: null,
  focusLineage: null,
  trendGenes: [],
  theme: "dark",
  view: "2d",
  cursor: null,
  playing: false,
  sheetOpen: true,
  pointSize: 1,
  fog: true,

  setEmbedding: (embedding) => set({ embedding }),
  setColorBy: (colorBy) => set({ colorBy }),
  setRun: (runId) => set({ runId, focusLineage: null, cursor: null, playing: false }),
  setFocusLineage: (focusLineage) => set({ focusLineage }),
  addTrendGene: (gene) =>
    set((s) => (s.trendGenes.includes(gene) ? s : { trendGenes: [gene, ...s.trendGenes].slice(0, MAX_TREND_GENES) })),
  setTheme: (theme) => set({ theme }),
  setView: (view) => set({ view }),
  setCursor: (cursor) => set({ cursor }),
  setPlaying: (playing) => set({ playing }),
  setSheetOpen: (sheetOpen) => set({ sheetOpen }),
  setPointSize: (pointSize) => set({ pointSize }),
  setFog: (fog) => set({ fog }),
}));
