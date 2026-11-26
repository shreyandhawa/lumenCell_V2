import { describe, expect, it } from "vitest";
import { clampTarget, depthRange, FLAT, fitZoom, lerpCam, TILTED, wrapAngle } from "@/lib/camera";

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

  it("brings angles into [-180, 180)", () => {
    expect([190, -190, 720, -32, 180, -180].map(wrapAngle)).toEqual([-170, 170, 0, -32, -180, -180]);
  });

  it("turns the orbit the shorter way round, and lands exactly on the end angle", () => {
    const tilt = { target: [0, 0, 0] as [number, number, number], zoom: 8, ...TILTED };
    // A camera 600 degrees round from the tilt (dragged there, or spun before auto-rotate wrapped):
    // the short way back is 120 degrees on, not 600 back.
    const spun = { ...tilt, rotationOrbit: TILTED.rotationOrbit + 50 * 12 };
    const turned = (p: number) => lerpCam(spun, tilt, p).rotationOrbit! - spun.rotationOrbit;
    expect(turned(0.5)).toBeCloseTo(60, 9);
    expect(lerpCam(spun, tilt, 1).rotationOrbit).toBe(TILTED.rotationOrbit);
    // Never more than half a turn, whichever two angles, and always ending on the one asked for.
    const pairs: [number, number][] = [[170, -170], [-170, 170], [0, 180], [-1000, 1000], [359.5, 0]];
    for (const [from, to] of pairs) {
      const [a, b] = [{ ...tilt, rotationOrbit: from }, { ...tilt, rotationOrbit: to }];
      expect(Math.abs(lerpCam(a, b, 0.999).rotationOrbit! - from)).toBeLessThanOrEqual(180);
      expect(wrapAngle(lerpCam(a, b, 0.999).rotationOrbit! - to)).toBeCloseTo(0, 0);
      expect(lerpCam(a, b, 1).rotationOrbit).toBe(to);
    }
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

  it("keeps the target inside the cells' box on x and y, and leaves z alone", () => {
    const box = [
      [-1, -1, 0],
      [1, 1, 0],
    ] as [[number, number, number], [number, number, number]];
    expect(clampTarget([0.2, -0.4, 0], box)).toEqual([0.2, -0.4, 0]);
    expect(clampTarget([3, -5, 0.7], box)).toEqual([1, -1, 0.7]);
  });
});
