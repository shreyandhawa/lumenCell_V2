export type Point = readonly [number, number];

/** Even-odd ray casting. */
export function pointInPolygon(x: number, y: number, poly: readonly Point[]): boolean {
  let inside = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const [xi, yi] = poly[i]!;
    const [xj, yj] = poly[j]!;
    if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

/**
 * Selection mask (1 = inside) for interleaved [x, y, x, y, ...] points. The lasso passes
 * cells' *screen* positions, so selection works the same in 2-D and in an orbited 3-D view.
 */
export function selectInPolygon(points: Float32Array, poly: readonly Point[]): Uint8Array {
  const n = points.length / 2;
  const mask = new Uint8Array(n);
  if (poly.length < 3) return mask;
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  for (const [x, y] of poly) {
    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
    minY = Math.min(minY, y);
    maxY = Math.max(maxY, y);
  }
  for (let i = 0; i < n; i++) {
    const x = points[i * 2]!;
    const y = points[i * 2 + 1]!;
    if (x < minX || x > maxX || y < minY || y > maxY) continue;
    if (pointInPolygon(x, y, poly)) mask[i] = 1;
  }
  return mask;
}

export interface Bounds {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
}

export function bounds(positions: Float32Array): Bounds {
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  for (let i = 0; i < positions.length; i += 2) {
    const x = positions[i]!;
    const y = positions[i + 1]!;
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
  }
  if (!Number.isFinite(minX)) return { minX: 0, minY: 0, maxX: 1, maxY: 1 };
  return { minX, minY, maxX, maxY };
}

/** Orthographic view state that fits the bounds into a width x height canvas. */
export function fitView(b: Bounds, width: number, height: number, padding = 0.08) {
  const w = Math.max(b.maxX - b.minX, 1e-9);
  const h = Math.max(b.maxY - b.minY, 1e-9);
  const scale = Math.min(width / w, height / h) * (1 - 2 * padding);
  return {
    target: [(b.minX + b.maxX) / 2, (b.minY + b.maxY) / 2, 0] as [number, number, number],
    zoom: Math.log2(Math.max(scale, 1e-9)),
  };
}

/* ---- The stage: every embedding drawn in one shared frame ------------------------------- */

/** Maps an embedding's own coordinates into the stage frame (centred, x and y within [-1, 1]). */
export interface StageTransform {
  center: [number, number, number];
  scale: number;
}

export interface Stage {
  /** n x 3 interleaved positions in the stage frame; z is 0 for 2-D embeddings. */
  positions: Float32Array;
  dims: 2 | 3;
  transform: StageTransform;
  /** Largest distance of a cell from the centre, for fitting a rotated 3-D view. */
  radius: number;
}

/**
 * Put an embedding on the stage. Scale comes from x and y only, so the same diffusion
 * components share one frame in 2-D and 3-D, and switching between them only lifts z.
 */
export function toStage(coords: Float32Array, dims: number): Stage {
  const d = dims >= 3 ? 3 : 2;
  const n = Math.floor(coords.length / dims);
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < d; k++) {
      const v = coords[i * dims + k]!;
      if (!Number.isFinite(v)) continue;
      if (v < lo[k]!) lo[k] = v;
      if (v > hi[k]!) hi[k] = v;
    }
  }
  const center: [number, number, number] = [0, 0, 0];
  for (let k = 0; k < d; k++) center[k] = Number.isFinite(lo[k]!) ? (lo[k]! + hi[k]!) / 2 : 0;
  const half = Math.max((hi[0]! - lo[0]!) / 2, (hi[1]! - lo[1]!) / 2, 1e-12);
  const scale = Number.isFinite(half) ? 1 / half : 1;
  const positions = new Float32Array(n * 3);
  let radius = 0;
  for (let i = 0; i < n; i++) {
    let r2 = 0;
    for (let k = 0; k < 3; k++) {
      const v = k < d ? (coords[i * dims + k]! - center[k]!) * scale : 0;
      positions[i * 3 + k] = v;
      if (Number.isFinite(v)) r2 += v * v;
    }
    if (r2 > radius) radius = r2;
  }
  return { positions, dims: d as 2 | 3, transform: { center, scale }, radius: Math.sqrt(radius) };
}

/** Map one point (2 or 3 coordinates) from embedding space into the stage frame. */
export function onStage(p: readonly number[], t: StageTransform): [number, number, number] {
  return [
    (p[0]! - t.center[0]) * t.scale,
    (p[1]! - t.center[1]) * t.scale,
    p.length > 2 ? (p[2]! - t.center[2]) * t.scale : 0,
  ];
}

export interface Insets {
  top: number;
  right: number;
  bottom: number;
  left: number;
}

/**
 * Zoom and target that fit the stage into the part of the canvas not covered by panels.
 * The target is offset so the stage centres in that free area, not under a panel.
 */
export function fitStage(radius: number, width: number, height: number, insets: Insets, dims: 2 | 3) {
  const freeW = Math.max(width - insets.left - insets.right, 80);
  const freeH = Math.max(height - insets.top - insets.bottom, 80);
  // 2-D stages span [-1, 1] in x and y; a rotating 3-D stage needs its whole sphere.
  const extent = dims === 3 ? Math.max(radius, 1) : 1;
  const scale = (Math.min(freeW, freeH) / (2 * extent)) * 0.86;
  const zoom = Math.log2(Math.max(scale, 1e-9));
  return { zoom, offset: centreOffset(insets, zoom) };
}

/** World-space shift that moves the view's centre into the middle of the free area. */
export function centreOffset(insets: Insets, zoom: number): [number, number] {
  const scale = 2 ** zoom;
  return [(insets.right - insets.left) / 2 / scale, (insets.bottom - insets.top) / 2 / scale];
}
