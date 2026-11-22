import type { OrbitViewState } from "@deck.gl/core";
import { useEffect, useRef } from "react";
import { useDataset, useRuns } from "@/api/queries";
import { decodeView, encodeView, type SharedView, type ViewCamera } from "@/lib/viewUrl";
import { useCamera } from "@/state/camera";
import { useUi } from "@/state/ui";

/**
 * Opens the view a link names (once, when the dataset and its runs have loaded), then keeps the
 * address bar naming the view on screen, so copying the address shares the picture.
 */
export function useSharedView() {
  const fingerprint = useDataset().data?.fingerprint;
  const runs = useRuns().data?.runs;
  const opened = useRef(false);

  useEffect(() => {
    if (!fingerprint || !runs) return;
    if (!opened.current) {
      opened.current = true;
      const shared = decodeView(window.location.hash, fingerprint);
      if (shared) open(shared, new Set(runs.map((r) => r.run_id)));
    }
    let timer: number | undefined;
    // At most every 300 ms: cameras change every frame while they move.
    const write = () => {
      timer ??= window.setTimeout(() => {
        timer = undefined;
        window.history.replaceState(null, "", `#${encodeView(shown(), fingerprint)}`);
      }, 300);
    };
    write();
    const offUi = useUi.subscribe(write);
    const offCamera = useCamera.subscribe(write);
    return () => {
      offUi();
      offCamera();
      window.clearTimeout(timer);
    };
  }, [fingerprint, runs]);
}

function open(v: SharedView, runIds: Set<string>) {
  const ui = useUi.getState();
  if (v.run && runIds.has(v.run)) ui.setRun(v.run); // clears focus and cursor, so those come after it
  if (v.colorBy) ui.setColorBy(v.colorBy);
  ui.setAxes(v.axes);
  ui.setPointSize(v.pointSize);
  if (v.focus !== null) ui.setFocusLineage(v.focus);
  if (v.cursor !== null) ui.setCursor(v.cursor);
  const camera = useCamera.getState();
  for (const [id, c] of Object.entries(v.cameras)) {
    // 3-D opens with the morph, which then flies to the shared camera instead of the default tilt.
    if (id === "stage" && v.view === "3d") camera.setArrival(c);
    else camera.setCamera(id, c as OrbitViewState);
  }
  ui.setView(v.view);
}

/** The view on screen, as a link stores it: only the cameras of the views shown. */
function shown(): SharedView {
  const ui = useUi.getState();
  const { cameras } = useCamera.getState();
  const kept: Record<string, ViewCamera> = {};
  for (const id of ui.view === "split" ? ["flat", "deep"] : ["stage"]) {
    const c = cameras[id];
    if (c) {
      const [x = 0, y = 0, z = 0] = c.target;
      kept[id] = { target: [x, y, z], zoom: Number(c.zoom), rotationX: c.rotationX ?? 0, rotationOrbit: c.rotationOrbit ?? 0 };
    }
  }
  return {
    run: ui.runId,
    colorBy: ui.colorBy,
    view: ui.view,
    cameras: kept,
    focus: ui.focusLineage,
    cursor: ui.cursor,
    axes: ui.axes,
    pointSize: ui.pointSize,
  };
}
