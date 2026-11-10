/**
 * The motion system in one place, shared by CSS (index.css), Motion and deck.gl.
 *
 * UI moves on a strong ease-out and stays under 250 ms. Things that move across the
 * field (cells morphing between embeddings, the camera) use ease-in-out, a little
 * longer, because the eye has to follow them. Playback is linear: constant motion.
 */
import type { Transition } from "motion/react";

export const EASE_OUT = [0.23, 1, 0.32, 1] as const;
export const EASE_IN_OUT = [0.77, 0, 0.175, 1] as const;

export const DURATION = {
  press: 0.12,
  ui: 0.18,
  panel: 0.24,
  morph: 0.7,
  drawOn: 0.9,
} as const;

/** Panels and sheets: a critically damped spring, interruptible mid-flight. */
export const PANEL_SPRING: Transition = { type: "spring", duration: 0.38, bounce: 0 };
/** Releasing a drag: a little life, never a wobble. */
export const RELEASE_SPRING: Transition = { type: "spring", duration: 0.42, bounce: 0.12 };
export const FADE: Transition = { duration: DURATION.ui, ease: EASE_OUT };

/** A cubic-bezier easing function for deck.gl transitions (which take t -> t'). */
export function cubicBezier(x1: number, y1: number, x2: number, y2: number): (t: number) => number {
  const cx = 3 * x1;
  const bx = 3 * (x2 - x1) - cx;
  const ax = 1 - cx - bx;
  const cy = 3 * y1;
  const by = 3 * (y2 - y1) - cy;
  const ay = 1 - cy - by;
  const sampleX = (t: number) => ((ax * t + bx) * t + cx) * t;
  const sampleY = (t: number) => ((ay * t + by) * t + cy) * t;
  const slopeX = (t: number) => (3 * ax * t + 2 * bx) * t + cx;
  return (x: number) => {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    let t = x;
    for (let i = 0; i < 6; i++) {
      const err = sampleX(t) - x;
      const d = slopeX(t);
      if (Math.abs(err) < 1e-5 || Math.abs(d) < 1e-6) break;
      t -= err / d;
    }
    return sampleY(Math.min(1, Math.max(0, t)));
  };
}

export const easeInOut = cubicBezier(...EASE_IN_OUT);
export const easeOut = cubicBezier(...EASE_OUT);

export function prefersReducedMotion(): boolean {
  return typeof window !== "undefined" && window.matchMedia?.("(prefers-reduced-motion: reduce)").matches === true;
}
