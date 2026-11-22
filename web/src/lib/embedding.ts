export type Vec3 = [number, number, number];

export interface Normalised {
  /** n x 3 positions centred on 0, longest side 1 (z = 0 for a 2-D embedding). */
  xyz: Float32Array;
  /** The same transform for one point (curves, milestones). */
  map: (p: ArrayLike<number>) => Vec3;
  /** The box that holds every cell, after the transform: the depth fog spans it. */
  bounds: [Vec3, Vec3];
}

/** Centre an (n x dims) embedding on 0 and scale its longest side to 1, as xyz. */
export function normalise(src: Float32Array, dims: 2 | 3): Normalised {
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < src.length; i++) {
    const v = src[i]!;
    const k = i % dims;
    if (Number.isFinite(v)) {
      if (v < lo[k]!) lo[k] = v;
      if (v > hi[k]!) hi[k] = v;
    }
  }
  const mid = [0, 0, 0];
  let span = 0;
  for (let k = 0; k < dims; k++) {
    if (hi[k]! < lo[k]!) continue; // no finite value on this axis
    mid[k] = (lo[k]! + hi[k]!) / 2;
    span = Math.max(span, hi[k]! - lo[k]!);
  }
  span ||= 1;
  const n = src.length / dims;
  const xyz = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < dims; k++) xyz[i * 3 + k] = (src[i * dims + k]! - mid[k]!) / span;
  }
  const map = (p: ArrayLike<number>): Vec3 => [0, 1, 2].map((k) => (k < dims ? (p[k]! - mid[k]!) / span : 0)) as Vec3;
  const edge = (k: number, v: number) => (k < dims && hi[k]! >= lo[k]! ? (v - mid[k]!) / span : 0);
  const bounds: [Vec3, Vec3] = [
    [0, 1, 2].map((k) => edge(k, lo[k]!)) as Vec3,
    [0, 1, 2].map((k) => edge(k, hi[k]!)) as Vec3,
  ];
  return { xyz, map, bounds };
}
