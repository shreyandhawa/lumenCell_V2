import { Pause, Play } from "lucide-react";
import { animate } from "motion";
import { useEffect } from "react";
import { rgbaCss, sequential } from "@/lib/colors";
import { useUi } from "@/state/ui";

/** Seconds for one sweep from the root to the end of pseudotime. */
const PLAY_SECONDS = 8;
const TICKS = [0, 0.25, 0.5, 0.75, 1];
const VIRIDIS = `linear-gradient(to right, ${TICKS.map((t) => rgbaCss(sequential(t, 255))).join(",")})`;

/** The bar and the gene trends share these columns, so their x axes line up. */
export const TIMELINE_COLS = "grid grid-cols-[52px_1fr_64px] gap-x-3 px-3";
/** Where the shared cursor sits inside the timeline, for a cursor in [0, 1]. */
export const cursorLeft = (t: number) => `calc(76px + ${t} * (100% - 164px))`;

/** Play and scrub pseudotime; one cursor shared with the gene trends. */
export function PseudotimeBar({ ready }: { ready: boolean }) {
  const cursor = useUi((s) => s.cursor);
  const setCursor = useUi((s) => s.setCursor);
  const playing = useUi((s) => s.playing);
  const setPlaying = useUi((s) => s.setPlaying);

  useEffect(() => {
    if (!playing) return;
    const now = useUi.getState().cursor ?? 0;
    const from = now >= 1 ? 0 : now;
    const controls = animate(from, 1, {
      duration: (1 - from) * PLAY_SECONDS,
      ease: "linear",
      onUpdate: setCursor,
      onComplete: () => setPlaying(false),
    });
    return () => controls.stop();
  }, [playing, setCursor, setPlaying]);

  return (
    <div aria-label="Pseudotime" role="group" className={`${TIMELINE_COLS} h-14 items-center rounded-xl bg-panel lift`}>
      <button
        type="button"
        aria-label={playing ? "Pause" : "Play"}
        disabled={!ready}
        onClick={() => setPlaying(!playing)}
        className="press grid size-8 place-items-center rounded-md bg-panel-raised text-ink disabled:text-ink-faint"
      >
        {playing ? <Pause size={16} fill="currentColor" /> : <Play size={16} />}
      </button>
      {ready ? (
        <div className="relative h-9">
          <div className="absolute top-2 right-0 left-0 h-1 rounded-full" style={{ background: VIRIDIS }} />
          <input
            type="range"
            className="scrub"
            aria-label="Pseudotime cursor"
            min={0}
            max={1}
            step={0.01}
            value={cursor ?? 0}
            onChange={(e) => {
              setPlaying(false);
              setCursor(Number(e.target.value));
            }}
          />
          {TICKS.map((t) => (
            <span
              key={t}
              aria-hidden
              className="pointer-events-none absolute top-5 font-mono text-2xs text-ink-muted"
              style={{ left: `${t * 100}%`, transform: `translateX(${t === 0 ? "0" : t === 1 ? "-100%" : "-50%"})` }}
            >
              {t}
            </span>
          ))}
        </div>
      ) : (
        <div>
          <div className="h-1 rounded-full bg-panel-raised" />
          <p className="mt-2 text-xs text-ink-muted">Pseudotime appears after a run</p>
        </div>
      )}
      <span className="text-right font-mono text-sm text-ink">pt {ready ? (cursor ?? 0).toFixed(2) : "–"}</span>
    </div>
  );
}
