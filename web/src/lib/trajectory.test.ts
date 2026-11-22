import { describe, expect, it } from "vitest";
import type { Vec3 } from "@/lib/embedding";
import { centroid, cursorPoint, placeMilestones, prefixTo } from "@/lib/trajectory";

const half: (a: Vec3, b: Vec3) => Vec3 = (a, b) => a.map((v, k) => (v + b[k]!) / 2) as Vec3;

describe("placeMilestones", () => {
  const ms = [
    { id: 0, kind: "root" as const, label: "IP", pseudotime: 0 },
    { id: 1, kind: "fork" as const, label: "3", pseudotime: 0.4 },
    { id: 2, kind: "tip" as const, label: "7", pseudotime: 1 },
  ];
  it("mixes the two views and keeps the run's kind and label", () => {
    const flat = new Map<number, Vec3>([
      [0, [0, 0, 0]],
      [1, [2, 0, 0]],
    ]);
    const deep = new Map<number, Vec3>([
      [0, [0, 0, 4]],
      [1, [2, 2, 0]],
    ]);
    expect(placeMilestones(ms, flat, deep, half)).toEqual([
      { ...ms[0], position: [0, 0, 2] },
      { ...ms[1], position: [2, 1, 0] },
    ]);
  });
  it("falls back to the one view that places a milestone, and drops unplaced ones", () => {
    const deep = new Map<number, Vec3>([[2, [1, 1, 1]]]);
    expect(placeMilestones(ms, null, deep, half)).toEqual([{ ...ms[2], position: [1, 1, 1] }]);
  });
});

describe("cursorPoint", () => {
  // four cells: two on the lineage near pt 0.5, one off the lineage, one far in pseudotime
  const pt = new Float32Array([0.49, 0.51, 0.5, 0.9]);
  const member = new Uint8Array([1, 1, 0, 1]);
  const xyz = new Float32Array([0, 0, 0, 2, 4, 6, 100, 100, 100, 50, 50, 50]);
  it("averages the lineage's cells near the cursor only", () => {
    expect(cursorPoint(pt, member, xyz, 0.5)).toEqual([1, 2, 3]);
  });
  it("is null where the lineage has no cells near the cursor", () => {
    expect(cursorPoint(pt, member, xyz, 0.7)).toBeNull();
  });
  it("skips cells without a pseudotime", () => {
    expect(cursorPoint(new Float32Array([Number.NaN]), new Uint8Array([1]), new Float32Array(3), Number.NaN)).toBeNull();
  });
});

describe("prefixTo", () => {
  const path: Vec3[] = [
    [0, 0, 0],
    [1, 0, 0],
    [2, 0, 0],
    [3, 0, 0],
  ];
  it("keeps the curve up to the nearest vertex, then ends on the point", () => {
    expect(prefixTo(path, [1.2, 0.5, 0])).toEqual([[0, 0, 0], [1, 0, 0], [1.2, 0.5, 0]]);
  });
  it("is just the point at the start of the curve", () => {
    expect(prefixTo(path, [-1, 0, 0])).toEqual([[0, 0, 0], [-1, 0, 0]]);
  });
});

describe("centroid", () => {
  it("is the mean of the points, and null for an empty curve", () => {
    expect(centroid([[0, 0, 0], [2, 4, 6]])).toEqual([1, 2, 3]);
    expect(centroid([])).toBeNull();
  });
});
