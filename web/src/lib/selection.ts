/** Cell selection: a mask with one byte per cell (1 = selected), built from lassos and clicks. */

export type SelectOp = "new" | "add" | "subtract" | "intersect";

/**
 * Screen pixels (top-left origin) of every cell, as interleaved [x, y, ...], in one pass with the
 * viewport's pixel projection matrix (column-major), the same maths as deck.gl's `viewport.project`.
 * Cells behind the camera or without a position get NaN, which no lasso contains.
 */
export function projectAll(xyz: Float32Array, m: ArrayLike<number>): Float32Array {
  const n = xyz.length / 3;
  const out = new Float32Array(n * 2);
  for (let i = 0; i < n; i++) {
    const x = xyz[i * 3]!;
    const y = xyz[i * 3 + 1]!;
    const z = xyz[i * 3 + 2]!;
    const w = m[3]! * x + m[7]! * y + m[11]! * z + m[15]!;
    if (!(w > 0)) {
      out[i * 2] = out[i * 2 + 1] = Number.NaN;
      continue;
    }
    out[i * 2] = (m[0]! * x + m[4]! * y + m[8]! * z + m[12]!) / w;
    out[i * 2 + 1] = (m[1]! * x + m[5]! * y + m[9]! * z + m[13]!) / w;
  }
  return out;
}

/** Apply a new selection to the current one. */
export function combine(current: Uint8Array | null, next: Uint8Array, op: SelectOp): Uint8Array {
  if (op === "new" || !current) return op === "subtract" ? new Uint8Array(next.length) : next;
  const out = new Uint8Array(next.length);
  for (let i = 0; i < out.length; i++) {
    const a = current[i]!;
    const b = next[i]!;
    out[i] = op === "add" ? a | b : op === "subtract" ? a & (b ^ 1) : a & b;
  }
  return out;
}

export function countSelected(mask: Uint8Array): number {
  let n = 0;
  for (let i = 0; i < mask.length; i++) n += mask[i]!;
  return n;
}

/** Cells per category, in the selection and overall (codes below 0 are missing and not counted). */
export function composition(mask: Uint8Array, codes: ArrayLike<number>, nCategories: number) {
  const selected = new Array<number>(nCategories).fill(0);
  const all = new Array<number>(nCategories).fill(0);
  for (let i = 0; i < codes.length; i++) {
    const c = codes[i]!;
    if (c < 0 || c >= nCategories) continue;
    all[c]!++;
    if (mask[i]) selected[c]!++;
  }
  return { selected, all };
}
