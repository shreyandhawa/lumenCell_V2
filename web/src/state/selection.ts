/** The cell selection, shared by the plot (dimming, lasso) and the sidebar (composition). */
import { create } from "zustand";
import { countSelected, type SelectOp } from "@/lib/selection";

interface SelectionState {
  /** One byte per cell, 1 = selected; null = no selection (nothing dimmed). */
  mask: Uint8Array | null;
  count: number;
  /** How the next lasso or Alt-click combines with the selection. */
  op: SelectOp;
  setMask: (mask: Uint8Array | null) => void;
  setOp: (op: SelectOp) => void;
}

export const useSelection = create<SelectionState>((set) => ({
  mask: null,
  count: 0,
  op: "new",
  setMask: (mask) => {
    const count = mask ? countSelected(mask) : 0;
    // An empty selection is no selection: nothing should dim.
    set({ mask: count ? mask : null, count });
  },
  setOp: (op) => set({ op }),
}));
