import { describe, expect, it } from "vitest";
import { normalise } from "@/lib/embedding";

const row = (a: Float32Array, i: number) => Array.from(a.subarray(i * 3, i * 3 + 3));

describe("normalise", () => {
  it("centres a 2-D embedding on 0, scales its longest side to 1 and sets z to 0", () => {
    // x spans 0..4, y spans 10..12: x is the longest side
    const { xyz, map, bounds } = normalise(new Float32Array([0, 10, 4, 12, 2, 11]), 2);
    expect(row(xyz, 0)).toEqual([-0.5, -0.25, 0]);
    expect(row(xyz, 1)).toEqual([0.5, 0.25, 0]);
    expect(row(xyz, 2)).toEqual([0, 0, 0]);
    expect(map([4, 12])).toEqual([0.5, 0.25, 0]);
    expect(bounds).toEqual([
      [-0.5, -0.25, 0],
      [0.5, 0.25, 0],
    ]);
  });

  it("keeps depth for a 3-D embedding", () => {
    const { xyz, bounds } = normalise(new Float32Array([0, 0, -1, 1, 1, 1]), 3);
    expect(row(xyz, 0)).toEqual([-0.25, -0.25, -0.5]);
    expect(row(xyz, 1)).toEqual([0.25, 0.25, 0.5]);
    expect(bounds).toEqual([
      [-0.25, -0.25, -0.5],
      [0.25, 0.25, 0.5],
    ]);
  });

  it("ignores missing values for the bounds and keeps them missing", () => {
    const { xyz } = normalise(new Float32Array([0, 0, Number.NaN, 2, 2, 2]), 2);
    expect(row(xyz, 0)).toEqual([-0.5, -0.5, 0]);
    expect(Number.isNaN(xyz[3]!)).toBe(true);
    expect(row(xyz, 2)).toEqual([0.5, 0.5, 0]);
  });

  it("does not divide by zero when every cell sits on one point", () => {
    const { xyz, bounds } = normalise(new Float32Array([3, 3, 3, 3]), 2);
    expect(Array.from(xyz)).toEqual([0, 0, 0, 0, 0, 0]);
    expect(bounds).toEqual([
      [0, 0, 0],
      [0, 0, 0],
    ]);
  });
});
