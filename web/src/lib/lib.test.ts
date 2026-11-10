import { describe, expect, it } from "vitest";
import {
  categoryPalette,
  colorsForCodes,
  colorsForLineages,
  colorsForValues,
  hexToRgb,
  MISSING,
  robustDomain,
  sequential,
  SHARED,
} from "@/lib/colors";
import { bounds, fitView, pointInPolygon, selectInPolygon } from "@/lib/geometry";

const rgba = (buf: Uint8Array, i: number) => Array.from(buf.subarray(i * 4, i * 4 + 4));

describe("colors", () => {
  it("parses hex", () => {
    expect(hexToRgb("#0072B2")).toEqual([0, 114, 178]);
    expect(hexToRgb("#fff")).toEqual([255, 255, 255]);
  });

  it("maps the sequential ramp ends to viridis' ends and clamps", () => {
    expect(sequential(0).slice(0, 3)).toEqual([68, 1, 84]);
    expect(sequential(1).slice(0, 3)).toEqual([253, 231, 37]);
    expect(sequential(-5)).toEqual(sequential(0));
  });

  it("colours codes, with missing codes greyed out", () => {
    const palette = categoryPalette(2);
    const buf = colorsForCodes(new Int16Array([0, 1, -1, 7]), palette);
    expect(rgba(buf, 0)).toEqual([...palette[0]!]);
    expect(rgba(buf, 2)).toEqual([...MISSING]);
    expect(rgba(buf, 3)).toEqual([...MISSING]);
  });

  it("uses the file's colours only when there is one per category", () => {
    expect(categoryPalette(2, ["#ff0000", "#00ff00"])[1]!.slice(0, 3)).toEqual([0, 255, 0]);
    expect(categoryPalette(3, ["#ff0000"])[0]!.slice(0, 3)).toEqual([0, 114, 178]);
  });

  it("uses a robust domain so one outlier cannot wash out the scale", () => {
    const values = new Float32Array(1000).map((_, i) => i / 999);
    values[0] = 1e6;
    const [lo, hi] = robustDomain(values);
    expect(lo).toBeLessThan(0.02);
    expect(hi).toBeLessThan(1.01);
  });

  it("greys out NaN values", () => {
    const buf = colorsForValues(new Float32Array([0, Number.NaN, 1]), [0, 1]);
    expect(rgba(buf, 1)).toEqual([...MISSING]);
  });

  it("colours lineages, the shared trunk and unassigned cells differently", () => {
    const buf = colorsForLineages(new Int16Array([0, -1, -2]));
    expect(rgba(buf, 1)).toEqual([...SHARED]);
    expect(rgba(buf, 2)).toEqual([...MISSING]);
    expect(rgba(buf, 0)).not.toEqual(rgba(buf, 1));
  });
});

describe("geometry", () => {
  const square = [
    [0, 0],
    [2, 0],
    [2, 2],
    [0, 2],
  ] as const;

  it("tests points against a polygon", () => {
    expect(pointInPolygon(1, 1, square)).toBe(true);
    expect(pointInPolygon(3, 1, square)).toBe(false);
  });

  it("builds a selection mask from interleaved positions", () => {
    const positions = new Float32Array([1, 1, 3, 3, 0.5, 1.5]);
    expect(Array.from(selectInPolygon(positions, square))).toEqual([1, 0, 1]);
    expect(Array.from(selectInPolygon(positions, square.slice(0, 2)))).toEqual([0, 0, 0]);
  });

  it("computes bounds, skipping non-finite points", () => {
    expect(bounds(new Float32Array([0, 0, Number.NaN, 5, 4, 2]))).toEqual({ minX: 0, minY: 0, maxX: 4, maxY: 2 });
  });

  it("fits a view to the canvas", () => {
    const v = fitView({ minX: 0, minY: 0, maxX: 10, maxY: 10 }, 1000, 500, 0);
    expect(v.target).toEqual([5, 5, 0]);
    expect(v.zoom).toBeCloseTo(Math.log2(50));
  });
});
