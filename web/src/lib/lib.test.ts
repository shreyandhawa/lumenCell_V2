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
import { bounds, fitView, onStage, pointInPolygon, selectInPolygon, toStage } from "@/lib/geometry";
import { playWindow, progressFor, softRange } from "@/lib/playback";

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

describe("stage", () => {
  it("centres an embedding and scales x and y into [-1, 1], padding 2-D with z = 0", () => {
    const s = toStage(new Float32Array([0, 0, 10, 4, 4, 2]), 2);
    expect(s.dims).toBe(2);
    const expected = [-1, -0.4, 0, 1, 0.4, 0, -0.2, 0, 0];
    expected.forEach((v, i) => expect(s.positions[i]).toBeCloseTo(v, 6));
  });

  it("puts 2-D and 3-D versions of the same components in one frame, so a lift only moves z", () => {
    const flat = toStage(new Float32Array([0, 0, 2, 2]), 2);
    const deep = toStage(new Float32Array([0, 0, -1, 2, 2, 5]), 3);
    expect(deep.dims).toBe(3);
    for (const i of [0, 1]) {
      expect(deep.positions[i * 3]).toBeCloseTo(flat.positions[i * 3]!);
      expect(deep.positions[i * 3 + 1]).toBeCloseTo(flat.positions[i * 3 + 1]!);
    }
    expect(deep.positions[5]).toBeGreaterThan(deep.positions[2]!);
  });

  it("maps curve points into the same frame as the cells", () => {
    const s = toStage(new Float32Array([0, 0, 10, 4]), 2);
    const [x, y, z] = onStage([10, 4], s.transform);
    expect([x, y, z].map((v) => Math.round(v * 1e6) / 1e6)).toEqual([1, 0.4, 0]);
  });

  it("selects on screen positions, which is what the lasso passes in 2-D and 3-D", () => {
    const screen = new Float32Array([5, 5, 50, 50]);
    expect(Array.from(selectInPolygon(screen, [[0, 0], [10, 0], [10, 10], [0, 10]]))).toEqual([1, 0]);
  });
});

describe("playback", () => {
  it("sweeps a window from the root to the tips", () => {
    expect(playWindow(0)).toEqual([0, 0.14]);
    expect(playWindow(1)).toEqual([0.86, 1]);
    expect(playWindow(2)).toEqual([0.86, 1]);
  });

  it("resumes from the window on screen, and from the start when there is none", () => {
    expect(progressFor([0, 1])).toBe(0);
    expect(progressFor(playWindow(0.5))).toBeCloseTo(0.5, 2);
  });

  it("keeps soft edges inside the window", () => {
    const [a, b] = softRange([0.4, 0.6]);
    expect(a).toBeGreaterThan(0.4);
    expect(b).toBeLessThan(0.6);
    expect(a).toBeLessThan(b);
  });
});
