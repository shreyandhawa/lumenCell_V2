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

/** Selection mask (1 = inside) for interleaved [x, y, x, y, ...] positions. */
export function selectInPolygon(positions: Float32Array, poly: readonly Point[]): Uint8Array {
  const n = positions.length / 2;
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
    const x = positions[i * 2]!;
    const y = positions[i * 2 + 1]!;
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
