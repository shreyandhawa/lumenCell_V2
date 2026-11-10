/** The time-lapse: a window of pseudotime swept from the root to the tips. */

export const PLAY_WINDOW = 0.14;
/** A full sweep, root to tips, in milliseconds. */
export const PLAY_DURATION = 7000;

/** The window at playback progress p in [0, 1]: it enters at the root and leaves at the tips. */
export function playWindow(p: number, width = PLAY_WINDOW): [number, number] {
  const t = Math.min(1, Math.max(0, p));
  const lo = t * (1 - width);
  return [round(lo), round(lo + width)];
}

/** Where to resume: the progress whose window starts at `lo` (0 when the window is the full range). */
export function progressFor(range: [number, number], width = PLAY_WINDOW): number {
  const [lo, hi] = range;
  if (lo <= 0 && hi >= 1) return 0;
  return Math.min(1, Math.max(0, lo / (1 - width)));
}

/** Soft edges inside the window, so cells brighten as it reaches them and dim as it leaves. */
export function softRange(range: [number, number]): [number, number] {
  const [lo, hi] = range;
  const edge = Math.min(0.03, (hi - lo) / 4);
  return [lo + edge, hi - edge];
}

function round(v: number): number {
  return Math.round(v * 1000) / 1000;
}
