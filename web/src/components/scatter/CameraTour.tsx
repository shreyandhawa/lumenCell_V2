import type { OrbitViewState } from "@deck.gl/core";
import { Clapperboard, Play, Square } from "lucide-react";
import { animate, type ValueAnimationTransition } from "motion";
import { useRef, useState } from "react";
import { CAMERA_SPRING } from "@/components/scatter/useStageMotion";
import { toolClass } from "@/components/ui";
import { type Camera, lerpCam } from "@/lib/camera";
import { useCamera } from "@/state/camera";

const MAX_STOPS = 8;
/** Seconds from one saved view to the next while the tour plays. */
const LEG = 1.6;

/**
 * Saved views of the 3-D stage, and a tour through them: the camera moves linearly from each to
 * the next, turning the shorter way round. Moving the view yourself stops the tour.
 */
export function CameraTour({ id, reduce, onStart }: { id: string; reduce: boolean; onStart: () => void }) {
  const [open, setOpen] = useState(false);
  const [stops, setStops] = useState<Camera[]>([]);
  const [playing, setPlaying] = useState(false);
  const stop = useRef<() => void>(() => {});

  // Moves the camera along a path (p from 0 to 1); anything else moving the camera stops it.
  const fly = (path: (p: number) => Camera, how: ValueAnimationTransition<number>, onDone?: () => void) => {
    stop.current();
    onStart();
    const { cameras, setCamera } = useCamera.getState();
    const base = cameras[id];
    if (!base) return;
    if (reduce) {
      setCamera(id, { ...base, ...path(1) });
      onDone?.();
      return;
    }
    let last: OrbitViewState = base;
    const c = animate(0, 1, {
      ...how,
      onUpdate: (p) => {
        if (useCamera.getState().cameras[id] !== last) {
          c.stop(); // the user moved it
          onDone?.();
          return;
        }
        last = { ...base, ...path(p) };
        setCamera(id, last);
      },
      onComplete: onDone,
    });
    stop.current = () => c.stop();
  };

  const goTo = (k: number) => {
    const from = useCamera.getState().cameras[id] as Camera | null;
    if (from) fly((p) => lerpCam(from, stops[k]!, p), CAMERA_SPRING);
  };
  const play = () => {
    const from = useCamera.getState().cameras[id] as Camera | null;
    if (!from || stops.length < 2) return;
    const legs = [from, ...stops];
    setPlaying(true);
    fly(
      (p) => {
        const x = p * (legs.length - 1);
        const k = Math.min(legs.length - 2, Math.floor(x));
        return lerpCam(legs[k]!, legs[k + 1]!, x - k);
      },
      { duration: LEG * (legs.length - 1), ease: "linear" },
      () => setPlaying(false),
    );
  };

  return (
    <div className="relative">
      {open && (
        <div
          id="camera-tour"
          role="group"
          aria-label="Camera tour"
          className="lift absolute right-0 bottom-full mb-2 w-60 rounded-xl bg-panel p-3 text-sm text-ink"
        >
          <p className="text-xs leading-snug text-ink-muted">Save views of the 3-D stage, then play through them.</p>
          <div className="mt-2 flex flex-wrap gap-1">
            {stops.map((_, k) => (
              <button
                key={k}
                type="button"
                aria-label={`Go to view ${k + 1}`}
                onClick={() => goTo(k)}
                className="press grid size-7 place-items-center rounded-md bg-panel-raised font-mono text-xs"
              >
                {k + 1}
              </button>
            ))}
            <button
              type="button"
              disabled={stops.length >= MAX_STOPS}
              onClick={() => {
                const c = useCamera.getState().cameras[id];
                if (c) setStops((s) => [...s, { target: c.target, zoom: c.zoom, rotationX: c.rotationX, rotationOrbit: c.rotationOrbit }]);
              }}
              className="press h-7 rounded-md bg-panel-raised px-2 text-xs disabled:text-ink-faint"
            >
              Save this view
            </button>
          </div>
          <div className="mt-2 flex justify-between">
            <button
              type="button"
              disabled={stops.length < 2 || reduce}
              title={reduce ? "Off while your system asks for reduced motion" : stops.length < 2 ? "Save two views first" : undefined}
              onClick={() => {
                if (!playing) return play();
                stop.current();
                setPlaying(false);
              }}
              className="press flex h-7 items-center gap-1.5 rounded-md px-2 text-xs disabled:text-ink-faint"
            >
              {playing ? <Square size={12} /> : <Play size={12} />} {playing ? "Stop tour" : "Play tour"}
            </button>
            {stops.length > 0 && (
              <button type="button" onClick={() => setStops([])} className="press h-7 px-2 text-xs text-ink-muted hover:text-ink">
                Clear
              </button>
            )}
          </div>
        </div>
      )}
      <button
        type="button"
        aria-label="Camera tour"
        aria-expanded={open}
        aria-controls="camera-tour"
        onClick={() => setOpen((o) => !o)}
        className={toolClass(open)}
      >
        <Clapperboard size={16} aria-hidden /> Tour
      </button>
    </div>
  );
}
