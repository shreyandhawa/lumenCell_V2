import { OrbitView } from "@deck.gl/core";
import { describe, expect, it } from "vitest";
import { combine, composition, countSelected, projectAll } from "@/lib/selection";

const mask = (...v: number[]) => new Uint8Array(v);

describe("projectAll", () => {
  it("matches deck.gl's viewport.project for a tilted orbit camera", () => {
    const view = new OrbitView({ id: "t", orbitAxis: "Y", padding: { left: 300, right: 0, top: 0, bottom: 80 } });
    const viewport = view.makeViewport({
      width: 1440,
      height: 900,
      viewState: { target: [0.1, -0.05, 0], zoom: 9, rotationX: 18, rotationOrbit: -32 },
    })!;
    const xyz = new Float32Array([0, 0, 0, 0.4, -0.2, 0.3, -0.5, 0.5, -0.25]);
    const px = projectAll(xyz, viewport.pixelProjectionMatrix as unknown as number[]);
    for (let i = 0; i < 3; i++) {
      const [x, y] = viewport.project([xyz[i * 3]!, xyz[i * 3 + 1]!, xyz[i * 3 + 2]!]);
      expect(px[i * 2]).toBeCloseTo(x!, 2);
      expect(px[i * 2 + 1]).toBeCloseTo(y!, 2);
    }
  });

  it("gives NaN for a cell without a position", () => {
    const identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
    const px = projectAll(new Float32Array([Number.NaN, 0, 0, 2, 3, 0]), identity);
    expect(Number.isNaN(px[0]!)).toBe(true);
    expect([px[2], px[3]]).toEqual([2, 3]);
  });
});

describe("combine", () => {
  const current = mask(1, 1, 0, 0);
  const next = mask(0, 1, 1, 0);
  it("replaces, adds, subtracts and intersects", () => {
    expect(Array.from(combine(current, next, "new"))).toEqual([0, 1, 1, 0]);
    expect(Array.from(combine(current, next, "add"))).toEqual([1, 1, 1, 0]);
    expect(Array.from(combine(current, next, "subtract"))).toEqual([1, 0, 0, 0]);
    expect(Array.from(combine(current, next, "intersect"))).toEqual([0, 1, 0, 0]);
  });
  it("starts from the new lasso when nothing is selected, except a subtraction", () => {
    expect(Array.from(combine(null, next, "add"))).toEqual([0, 1, 1, 0]);
    expect(Array.from(combine(null, next, "intersect"))).toEqual([0, 1, 1, 0]);
    expect(Array.from(combine(null, next, "subtract"))).toEqual([0, 0, 0, 0]);
  });
});

describe("counts", () => {
  it("counts selected cells and the categories inside and outside the selection", () => {
    const m = mask(1, 0, 1, 1, 0);
    expect(countSelected(m)).toBe(3);
    // codes: two of category 0, two of category 1, one missing
    expect(composition(m, new Int16Array([0, 0, 1, 1, -1]), 2)).toEqual({ selected: [1, 2], all: [2, 2] });
  });
});
