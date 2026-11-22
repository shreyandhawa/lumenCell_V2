/** Camera per deck.gl view id. Null until the user or an animation moves it (the view then fits). */
import type { OrbitViewState } from "@deck.gl/core";
import { create } from "zustand";
import type { ViewCamera } from "@/lib/viewUrl";

interface CameraState {
  cameras: Record<string, OrbitViewState | null>;
  /** Where the next morph into 3-D ends instead of the default tilt: a shared link's camera. */
  arrival: ViewCamera | null;
  setCamera: (id: string, vs: OrbitViewState | null) => void;
  setArrival: (c: ViewCamera | null) => void;
}

export const useCamera = create<CameraState>((set) => ({
  cameras: {},
  arrival: null,
  setCamera: (id, vs) => set((s) => ({ cameras: { ...s.cameras, [id]: vs } })),
  setArrival: (arrival) => set({ arrival }),
}));
