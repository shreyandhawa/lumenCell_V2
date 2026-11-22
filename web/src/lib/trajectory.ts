/** Geometry for drawing a run's trajectory over the cells: milestones and the pseudotime playhead. */
import type { Vec3 } from "@/lib/embedding";

export type MilestoneKind = "root" | "fork" | "tip" | "cluster";

export interface PlacedMilestone {
  id: number;
  kind: MilestoneKind;
  label: string;
  pseudotime: number;
  position: Vec3;
}

type Lerp = (a: Vec3, b: Vec3) => Vec3;

/**
 * Milestones of the run (or PAGA's clusters: anything with an id) at their projected positions: in
 * the 2-D view, the 3-D view, or between them while morphing (`mix`). One missing from a view keeps
 * the other's position.
 */
export function placeMilestones<M extends { id: number }>(
  milestones: readonly M[],
  flat: ReadonlyMap<number, Vec3> | null,
  deep: ReadonlyMap<number, Vec3> | null,
  mix: Lerp,
): (M & { position: Vec3 })[] {
  const out: (M & { position: Vec3 })[] = [];
  for (const m of milestones) {
    const a = flat?.get(m.id);
    const b = deep?.get(m.id);
    const position = a && b ? mix(a, b) : (a ?? b);
    if (position) out.push({ ...m, position });
  }
  return out;
}

/**
 * Where a lineage stands at pseudotime `cursor`: the mean position of its cells within `half` of
 * the cursor. Null when no cell of the lineage is that close (a gap in pseudotime).
 */
export function cursorPoint(
  pseudotime: ArrayLike<number>,
  member: ArrayLike<number>,
  xyz: Float32Array,
  cursor: number,
  half = 0.02,
): Vec3 | null {
  let n = 0;
  const sum = [0, 0, 0];
  for (let i = 0; i < pseudotime.length; i++) {
    const pt = pseudotime[i]!;
    if (!member[i] || !(Math.abs(pt - cursor) <= half)) continue;
    const x = xyz[i * 3]!;
    const y = xyz[i * 3 + 1]!;
    const z = xyz[i * 3 + 2]!;
    if (!Number.isFinite(x + y + z)) continue;
    sum[0]! += x;
    sum[1]! += y;
    sum[2]! += z;
    n++;
  }
  return n ? (sum.map((s) => s / n) as Vec3) : null;
}

/** The curve from its start to the vertex nearest `p`, then on to `p`: the stretch already travelled. */
export function prefixTo(path: readonly Vec3[], p: Vec3): Vec3[] {
  let best = 0;
  let bestD = Infinity;
  path.forEach((q, k) => {
    const d = (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2;
    if (d < bestD) [best, bestD] = [k, d];
  });
  return [...path.slice(0, best + 1), p];
}

/** The mean of a curve's points: where the camera turns around when the lineage is focused. */
export function centroid(path: readonly Vec3[]): Vec3 | null {
  if (!path.length) return null;
  const s = [0, 0, 0];
  for (const q of path) for (let k = 0; k < 3; k++) s[k]! += q[k]!;
  return s.map((v) => v / path.length) as Vec3;
}
