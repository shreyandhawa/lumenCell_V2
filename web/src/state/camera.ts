/** Camera per deck.gl view id. Null until the user or an animation moves it (the view then fits). */
import type { OrbitViewState } from "@deck.gl/core";
import { create } from "zustand";

interface CameraState {
  cameras: Record<string, OrbitViewState | null>;
  setCamera: (id: string, vs: OrbitViewState | null) => void;
}

export const useCamera = create<CameraState>((set) => ({
  cameras: {},
  setCamera: (id, vs) => set((s) => ({ cameras: { ...s.cameras, [id]: vs } })),
}));
