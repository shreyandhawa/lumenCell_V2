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
  /** Inclusive pseudotime window; [0, 1] means no filter. */
  ptRange: [number, number];
  focusLineage: number | null;
  /** 1 = selected. null = no selection. */
  selection: Uint8Array | null;
  selectionCount: number;
  lassoActive: boolean;
  hovered: number | null;

  setEmbedding: (name: string) => void;
  setColorBy: (c: ColorBy | null) => void;
  setRun: (id: string | null) => void;
  setPtRange: (r: [number, number]) => void;
  setFocusLineage: (id: number | null) => void;
  setSelection: (mask: Uint8Array | null) => void;
  setLasso: (on: boolean) => void;
  setHovered: (i: number | null) => void;
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

  setEmbedding: (embedding) => set({ embedding }),
  setColorBy: (colorBy) => set({ colorBy }),
  setRun: (runId) => set({ runId, focusLineage: null, ptRange: [0, 1] }),
  setPtRange: (ptRange) => set({ ptRange }),
  setFocusLineage: (focusLineage) => set({ focusLineage }),
  setSelection: (mask) => {
    let count = 0;
    if (mask) for (let i = 0; i < mask.length; i++) count += mask[i]!;
    set({ selection: mask && count > 0 ? mask : null, selectionCount: count });
  },
  setLasso: (lassoActive) => set({ lassoActive }),
  setHovered: (hovered) => set({ hovered }),
}));
