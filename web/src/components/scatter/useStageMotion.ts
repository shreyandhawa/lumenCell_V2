import type { OrbitViewState } from "@deck.gl/core";
import { animate, cubicBezier } from "motion";
import { useEffect, useRef, useState } from "react";
import { type Camera, FLAT, lerp, lerpCam, TILTED } from "@/lib/camera";
import type { Normalised, Vec3 } from "@/lib/embedding";
import { centroid } from "@/lib/trajectory";
import { useCamera } from "@/state/camera";
import type { View } from "@/state/ui";

/** The stage's camera id, one view filling the window. */
export const VIEW_ID = "stage";
/** The split view: the 2-D embedding on the left, DC1-DC3 on the right, each with its own camera. */
export const PANES = ["flat", "deep"] as const;
export type Pane = (typeof PANES)[number];
export const EASE_OUT = [0.23, 1, 0.32, 1] as const;
const EASE_IN_OUT = [0.77, 0, 0.175, 1] as const;
const MORPH_SECONDS = 0.6;
/** The cells' GPU transition. It runs on the same curve as the camera, so the root pin and curves stay on them. */
export const MORPH = { duration: MORPH_SECONDS * 1000, easing: cubicBezier(...EASE_IN_OUT) };
/** Auto-rotate speed, degrees per second. */
const SPIN = 12;

/**
 * How the stage moves on its own: the 2-D <-> 3-D morph (``t``, 0 = the 2-D embedding, 1 = DC1-DC3),
 * auto-rotate, reset, and the turn toward a focused lineage. The cameras live in the camera store.
 */
export function useStageMotion({
  view,
  deep,
  viewState,
  paneStates,
  fitZoom,
  paneFit,
  focus,
  focusPath,
  reduce,
}: {
  view: View;
  /** DC1-DC3; until they load, the morph into 3-D waits. */
  deep: Normalised | null;
  /** The stage's camera, null until the first fit. */
  viewState: OrbitViewState | null;
  /** Side by side, each pane's camera; null otherwise. */
  paneStates: Record<Pane, OrbitViewState> | null;
  fitZoom: number;
  paneFit: number;
  focus: number | null;
  /** The focused lineage's curve in DC1-DC3. */
  focusPath: readonly Vec3[] | undefined;
  reduce: boolean;
}) {
  const setCamera = useCamera((s) => s.setCamera);
  const [t, setT] = useState(0);

  // Live values for the animations below, which must not restart when these change.
  const live = useRef({ viewState, t, fitZoom, paneStates });
  useEffect(() => {
    live.current = { viewState, t, fitZoom, paneStates };
  });

  /** Move a camera (the stage's unless a pane is named) and the morph toward toT together. */
  const tween = (
    to: Camera,
    toT: number,
    seconds: number,
    ease: readonly number[],
    id: string = VIEW_ID,
    from: OrbitViewState | null = live.current.viewState,
  ) => {
    const t0 = live.current.t;
    if (!from) return () => {};
    if (reduce) {
      setT(toT);
      setCamera(id, { ...from, ...to });
      return () => {};
    }
    const c = animate(0, 1, {
      duration: seconds,
      ease: [...ease] as [number, number, number, number],
      onUpdate: (p) => {
        setT(lerp(t0, toT, p));
        setCamera(id, { ...from, ...lerpCam(from as Camera, to, p) });
      },
    });
    return () => c.stop();
  };

  // The 2-D <-> 3-D morph: positions and camera move together. It waits for the stage's first
  // camera and, into 3-D, for the diffusion components; it runs again when they arrive.
  const staged = viewState !== null;
  useEffect(() => {
    if (view === "split") return; // the panes do not morph; the stage morphs when it is shown again
    const toT = view === "3d" ? 1 : 0;
    if (live.current.t === toT || (toT === 1 && !deep) || !staged) return;
    // A shared link's 3-D camera, once; otherwise the default tilt.
    const { arrival, setArrival } = useCamera.getState();
    if (arrival && toT === 1) setArrival(null);
    const cam = (toT === 1 && arrival) || { target: [0, 0, 0] as Vec3, zoom: live.current.fitZoom, ...(toT ? TILTED : FLAT) };
    return tween(cam, toT, MORPH_SECONDS, EASE_IN_OUT);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- tween reads live values; only these start a morph
  }, [view, deep, staged]);

  // Auto-rotate, in 3-D only (the right pane side by side); any drag, key or reset stops it.
  const [spinning, setSpinning] = useState(false);
  const spin = spinning && (paneStates !== null || t === 1) && !reduce;
  useEffect(() => {
    if (!spin) return;
    let last = performance.now();
    let id = requestAnimationFrame(function step(now) {
      const panes = live.current.paneStates;
      const vs = panes ? panes.deep : live.current.viewState;
      if (vs) setCamera(panes ? "deep" : VIEW_ID, { ...vs, rotationOrbit: (vs.rotationOrbit ?? 0) + ((now - last) / 1000) * SPIN });
      last = now;
      id = requestAnimationFrame(step);
    });
    return () => cancelAnimationFrame(id);
  }, [spin, setCamera]);

  // Focusing a lineage in 3-D (the right pane side by side) turns the camera around it.
  useEffect(() => {
    const centre = focusPath && centroid(focusPath);
    const pane = paneStates?.deep;
    const vs = pane ?? live.current.viewState;
    if ((!pane && live.current.t !== 1) || !vs || !centre) return;
    const to = { target: centre, zoom: vs.zoom, rotationX: vs.rotationX, rotationOrbit: vs.rotationOrbit };
    return pane ? tween(to, live.current.t, 0.45, EASE_OUT, "deep", pane) : tween(to, 1, 0.45, EASE_OUT);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- only a focus change moves the camera
  }, [focus]);

  const reset = () => {
    setSpinning(false);
    if (paneStates) {
      for (const id of PANES)
        tween({ target: [0, 0, 0], zoom: paneFit, ...(id === "deep" ? TILTED : FLAT) }, t, 0.3, EASE_OUT, id, paneStates[id]);
      return;
    }
    tween({ target: [0, 0, 0], zoom: fitZoom, ...(view === "3d" ? TILTED : FLAT) }, t, 0.3, EASE_OUT);
  };

  return { t, spin, setSpinning, reset };
}
