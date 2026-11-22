/**
 * The view as a URL hash, so a link opens the same picture: run, colouring, view, cameras, focused
 * lineage, pseudotime cursor, 3-D axes and point size. It names the dataset's fingerprint and is
 * ignored on any other dataset. The selection is never stored.
 *
 * A link is input from anyone: every field is checked, and a bad one is left out, not guessed.
 */
import { type Axes, DEFAULT_AXES } from "@/lib/axes";
import { type ColorBy, colorByKey, parseColorBy, type View } from "@/state/ui";

export interface ViewCamera {
  target: [number, number, number];
  zoom: number;
  rotationX: number;
  rotationOrbit: number;
}

export interface SharedView {
  run: string | null;
  colorBy: ColorBy | null;
  view: View;
  /** By deck.gl view id: "stage", or "flat" and "deep" side by side. */
  cameras: Record<string, ViewCamera>;
  focus: number | null;
  cursor: number | null;
  axes: Axes;
  pointSize: number;
}

const VERSION = "1";
const VIEWS: readonly View[] = ["2d", "3d", "split"];
const CAMERA_IDS = ["stage", "flat", "deep"];
const round = (x: number) => Number(x.toFixed(4));

export function encodeView(v: SharedView, fingerprint: string): string {
  const p = new URLSearchParams({ v: VERSION, data: fingerprint, view: v.view });
  if (v.run) p.set("run", v.run);
  if (v.colorBy) p.set("color", colorByKey(v.colorBy));
  for (const [id, c] of Object.entries(v.cameras))
    p.set(`cam.${id}`, [...c.target, c.zoom, c.rotationX, c.rotationOrbit].map(round).join(","));
  if (v.focus !== null) p.set("focus", String(v.focus));
  if (v.cursor !== null) p.set("t", round(v.cursor).toString());
  if (v.axes.join() !== DEFAULT_AXES.join()) p.set("axes", v.axes.join(","));
  if (v.pointSize !== 1) p.set("size", String(v.pointSize));
  return p.toString();
}

/** The view a hash names, or null when it names none or another dataset. */
export function decodeView(hash: string, fingerprint: string): SharedView | null {
  const p = new URLSearchParams(hash.replace(/^#/, ""));
  if (p.get("v") !== VERSION || p.get("data") !== fingerprint) return null;
  const num = (key: string, ok: (x: number) => boolean) => {
    const raw = p.get(key);
    const x = raw === null || raw === "" ? NaN : Number(raw);
    return Number.isFinite(x) && ok(x) ? x : null;
  };
  const cameras: Record<string, ViewCamera> = {};
  for (const id of CAMERA_IDS) {
    const parts = p.get(`cam.${id}`)?.split(",").map(Number);
    if (parts?.length === 6 && parts.every(Number.isFinite)) {
      const [x, y, z, zoom, rotationX, rotationOrbit] = parts as [number, number, number, number, number, number];
      cameras[id] = { target: [x, y, z], zoom, rotationX, rotationOrbit };
    }
  }
  const axes = p.get("axes")?.split(",").map(Number);
  const view = p.get("view") as View;
  const run = p.get("run");
  return {
    run: run && /^[A-Za-z0-9_-]{1,64}$/.test(run) ? run : null,
    colorBy: parseColorBy(p.get("color") ?? ""),
    view: VIEWS.includes(view) ? view : "2d",
    cameras,
    focus: num("focus", (x) => Number.isInteger(x) && x >= 0),
    cursor: num("t", (x) => x >= 0 && x <= 1),
    axes:
      axes?.length === 3 && axes.every((a) => Number.isInteger(a) && a >= 1) && new Set(axes).size === 3 ? (axes as Axes) : DEFAULT_AXES,
    pointSize: num("size", (x) => x >= 0.5 && x <= 3) ?? 1,
  };
}
