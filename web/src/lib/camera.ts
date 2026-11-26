import type { OrbitViewState } from "@deck.gl/core";
import type { Vec3 } from "@/lib/embedding";

export type Camera = Pick<OrbitViewState, "target" | "zoom" | "rotationX" | "rotationOrbit">;

export const FLAT = { rotationX: 0, rotationOrbit: 0 };
export const TILTED = { rotationX: 18, rotationOrbit: -32 };

export const lerp = (a: number, b: number, t: number) => a + (b - a) * t;

/** An angle in degrees, brought into [-180, 180). */
export const wrapAngle = (deg: number) => ((((deg + 180) % 360) + 360) % 360) - 180;

/**
 * The camera at t between a and b. The orbit turns the shorter way round, at most 180 degrees, so a
 * camera that has spun for a while never unwinds every turn; at t = 1 it is exactly b's angle.
 */
export const lerpCam = (a: Camera, b: Camera, t: number): Camera => {
  const from = a.rotationOrbit ?? 0;
  const to = b.rotationOrbit ?? 0;
  return {
    target: [0, 1, 2].map((k) => lerp(a.target[k] ?? 0, b.target[k] ?? 0, t)) as Vec3,
    zoom: lerp(a.zoom as number, b.zoom as number, t),
    rotationX: lerp(a.rotationX ?? 0, b.rotationX ?? 0, t),
    rotationOrbit: t === 1 ? to : from + wrapAngle(to - from) * t,
  };
};

/** Zoom at which the unit-span data fills 80% of the free stage (log2 of pixels per unit). */
export function fitZoom(width: number, height: number): number {
  return Math.log2(0.8 * Math.max(120, Math.min(width, height)));
}

/**
 * The camera's target, held inside the box that holds the cells on x and y, so a zoom or a pan cannot
 * carry the view off the plot into empty space. z is left as it is.
 */
export function clampTarget(target: ArrayLike<number>, [lo, hi]: readonly [Vec3, Vec3]): Vec3 {
  const at = (k: 0 | 1) => Math.min(Math.max(target[k]!, lo[k]!), hi[k]!);
  return [at(0), at(1), target[2] ?? 0];
}

/**
 * The nearest and farthest camera distances of the box that holds every cell: the depth fog runs
 * from clear at the first to faded at the second.
 */
export function depthRange(camera: ArrayLike<number>, [lo, hi]: readonly [Vec3, Vec3]): [number, number] {
  let near = 0;
  let far = 0;
  for (let k = 0; k < 3; k++) {
    const c = camera[k]!;
    const inside = Math.min(Math.max(c, lo[k]!), hi[k]!); // the box's closest point on this axis
    near += (c - inside) ** 2;
    far += Math.max((c - lo[k]!) ** 2, (c - hi[k]!) ** 2);
  }
  return [Math.sqrt(near), Math.sqrt(far)];
}
