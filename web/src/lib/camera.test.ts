import { describe, expect, it } from "vitest";
import { depthRange, FLAT, fitZoom, lerpCam, TILTED } from "@/lib/camera";

describe("camera", () => {
  it("fits the unit-span data to 80% of the shorter side, with a floor for tiny stages", () => {
    expect(fitZoom(1000, 500)).toBeCloseTo(Math.log2(400), 9);
    expect(fitZoom(50, 50)).toBeCloseTo(Math.log2(96), 9);
  });

  it("interpolates every camera field and lands exactly on both ends", () => {
    const a = { target: [0, 0, 0] as [number, number, number], zoom: 8, ...FLAT };
    const b = { target: [1, -2, 4] as [number, number, number], zoom: 10, ...TILTED };
    expect(lerpCam(a, b, 0)).toEqual(a);
    expect(lerpCam(a, b, 1)).toEqual(b);
    expect(lerpCam(a, b, 0.5)).toEqual({ target: [0.5, -1, 2], zoom: 9, rotationX: 9, rotationOrbit: -16 });
  });

  it("bounds the fog by the box of cells: nearest point to farthest corner", () => {
    const box = [
      [-1, -1, -2],
      [1, 1, 2],
    ] as [[number, number, number], [number, number, number]];
    // straight on a face: the face is nearest, a back corner farthest
    const [near, far] = depthRange([0, 0, 10], box);
    expect(near).toBeCloseTo(8, 9);
    expect(far).toBeCloseTo(Math.hypot(1, 1, 12), 9);
    // off to one side: the nearest point is an edge, not the face centre
    const [n2] = depthRange([5, 5, 0], box);
    expect(n2).toBeCloseTo(Math.hypot(4, 4), 9);
    // the camera inside the box: near is 0
    expect(depthRange([0, 0, 0.5], box)[0]).toBe(0);
  });
});
