/** Geometry of the shell, shared by the panels and the camera (which keeps the cells out from under them). */
import type { Insets } from "@/lib/geometry";
import type { SheetDetent } from "@/state/ui";

export const GAP = 12;
export const CELLS_W = 264;
export const INSPECTOR_W = 336;
export const TIMELINE_H = 52;

export function sheetHeight(detent: SheetDetent, stageHeight: number): number {
  if (detent === "bar") return 44;
  if (detent === "peek") return 276;
  return Math.max(276, Math.round(stageHeight * 0.55));
}

export function stageInsets(
  panels: { cells: boolean; inspector: boolean },
  detent: SheetDetent,
  stageHeight: number,
  hasTimeline: boolean,
): Insets {
  return {
    top: GAP,
    left: panels.cells ? GAP + CELLS_W + GAP : GAP,
    right: panels.inspector ? GAP + INSPECTOR_W + GAP : GAP,
    bottom: GAP + sheetHeight(detent, stageHeight) + (hasTimeline ? GAP + TIMELINE_H : 0) + GAP,
  };
}
