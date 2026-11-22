/** UI-only state. Server data lives in TanStack Query, never here. */
import { create } from "zustand";
import { type Axes, DEFAULT_AXES } from "@/lib/axes";

export type ColorBy =
  | { kind: "obs"; column: string }
  | { kind: "gene"; gene: string }
  | { kind: "pseudotime" }
  | { kind: "lineage" }
  /** Palantir's per-cell fate entropy; only runs that serve `extras.entropy` offer it. */
  | { kind: "entropy" };

/** One string per colouring, as the Colour by menu and a shared link write it. */
export function colorByKey(c: ColorBy | null): string {
  if (!c) return "";
  if (c.kind === "obs") return `obs:${c.column}`;
  if (c.kind === "gene") return `gene:${c.gene}`;
  return c.kind;
}

export function parseColorBy(v: string): ColorBy | null {
  if (v.startsWith("obs:") && v.length > 4) return { kind: "obs", column: v.slice(4) };
  if (v.startsWith("gene:") && v.length > 5) return { kind: "gene", gene: v.slice(5) };
  if (v === "pseudotime" || v === "lineage" || v === "entropy") return { kind: v };
  return null;
}

export type View = "2d" | "3d" | "split";

export interface UiState {
  embedding: string | null;
  colorBy: ColorBy | null;
  runId: string | null;
  focusLineage: number | null;
  /** Genes the user picked, shown first in the gene trends. */
  trendGenes: string[];
  theme: "dark" | "light";
  /** One stage that morphs between the 2-D embedding and DC1-DC3, or both side by side. */
  view: View;
  /** The pseudotime cursor shared by the bar and the gene trends; null until the user plays or scrubs. */
  cursor: number | null;
  playing: boolean;
  sheetOpen: boolean;
  /** Point size as a multiple of the default for this many cells. */
  pointSize: number;
  /** Far cells fade toward the stage in 3-D. */
  fog: boolean;
  /** The diffusion components on the 3-D view's x, y and z. */
  axes: Axes;
  /** Graphs drawn over the plot. */
  overlays: { paga: boolean; knn: boolean };

  setEmbedding: (name: string) => void;
  setColorBy: (c: ColorBy | null) => void;
  setRun: (id: string | null) => void;
  setFocusLineage: (id: number | null) => void;
  addTrendGene: (gene: string) => void;
  setTheme: (t: "dark" | "light") => void;
  setView: (v: View) => void;
  setCursor: (t: number | null) => void;
  setPlaying: (on: boolean) => void;
  setSheetOpen: (open: boolean) => void;
  setPointSize: (k: number) => void;
  setFog: (on: boolean) => void;
  setAxes: (a: Axes) => void;
  setOverlay: (name: "paga" | "knn", on: boolean) => void;
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
  axes: DEFAULT_AXES,
  overlays: { paga: false, knn: false },

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
  setAxes: (axes) => set({ axes }),
  setOverlay: (name, on) => set((s) => ({ overlays: { ...s.overlays, [name]: on } })),
}));
