import { Slider } from "@base-ui/react/slider";
import { Pause, Play, RotateCcw } from "lucide-react";
import { useEffect, useRef } from "react";
import type * as T from "@/api/generated";
import { cx, IconButton, Tip } from "@/components/ui";
import { PLAY_DURATION, playWindow, progressFor } from "@/lib/playback";
import { useUi } from "@/state/ui";

/** Pseudotime as the field's clock: scrub a window, or play it from the root to the tips. */
export function Timeline({ run }: { run: T.TrajectoryRun }) {
  const range = useUi((s) => s.ptRange);
  const setRange = useUi((s) => s.setPtRange);
  const playing = useUi((s) => s.playing);
  const setPlaying = useUi((s) => s.setPlaying);
  const [lo, hi] = range;
  const active = lo > 0 || hi < 1;
  usePlayback();

  const forks = [...new Set(run.divergence.map((d) => Math.round(d.pseudotime * 100) / 100))].filter(
    (p) => p > 0 && p < 1,
  );

  return (
    <section
      aria-label="Pseudotime window"
      className="flex h-[52px] items-center gap-3 rounded-2xl border border-line bg-panel pr-3 pl-2"
    >
      <Tip label={playing ? "Pause" : "Play the time-lapse"} shortcut="Space">
        <button
          type="button"
          aria-label={playing ? "Pause" : "Play the time-lapse"}
          onClick={() => setPlaying(!playing)}
          className="pressable grid size-9 shrink-0 place-items-center rounded-full bg-live text-live-ink [@media(hover:hover)]:hover:brightness-110"
        >
          {playing ? <Pause size={15} fill="currentColor" /> : <Play size={15} fill="currentColor" className="ml-0.5" />}
        </button>
      </Tip>

      <Slider.Root
        value={[lo, hi]}
        min={0}
        max={1}
        step={0.01}
        minStepsBetweenValues={1}
        onValueChange={(v) => {
          const [a = 0, b = 1] = v as number[];
          if (playing) setPlaying(false);
          setRange([a, b]);
        }}
        className="min-w-0 flex-1"
      >
        <Slider.Control className="relative flex h-9 touch-none items-center select-none">
          <Slider.Track className="relative h-1.5 w-full rounded-full bg-control">
            {forks.map((f) => (
              <span
                key={f}
                className="absolute top-1/2 h-3 w-px -translate-y-1/2 bg-ink-3"
                style={{ left: `${f * 100}%` }}
                title={`a lineage splits off at ${f.toFixed(2)}`}
                aria-hidden
              />
            ))}
            <Slider.Indicator className={cx("rounded-full", active || playing ? "bg-live" : "bg-line-strong")} />
            <Slider.Thumb
              index={0}
              aria-label="Pseudotime from"
              className="size-4 rounded-full border-2 border-panel bg-ink shadow-[0_1px_3px_rgb(0_0_0/0.4)] transition-transform duration-100 ease-out has-[:focus-visible]:outline-2 has-[:focus-visible]:outline-offset-2 has-[:focus-visible]:outline-focus data-dragging:scale-110"
            />
            <Slider.Thumb
              index={1}
              aria-label="Pseudotime to"
              className="size-4 rounded-full border-2 border-panel bg-ink shadow-[0_1px_3px_rgb(0_0_0/0.4)] transition-transform duration-100 ease-out has-[:focus-visible]:outline-2 has-[:focus-visible]:outline-offset-2 has-[:focus-visible]:outline-focus data-dragging:scale-110"
            />
          </Slider.Track>
        </Slider.Control>
      </Slider.Root>

      <span className="num w-[5.5rem] shrink-0 text-right text-xs text-ink-2" aria-live="off">
        {lo.toFixed(2)} – {hi.toFixed(2)}
      </span>
      <IconButton
        label="Reset"
        onClick={() => {
          setPlaying(false);
          setRange([0, 1]);
        }}
        className={cx(!active && "pointer-events-none opacity-0")}
        tabIndex={active ? 0 : -1}
      >
        <RotateCcw size={15} />
      </IconButton>
    </section>
  );
}

/** Drives the sweep while `playing` is on; it ends on the last window and stops. */
function usePlayback() {
  const playing = useUi((s) => s.playing);
  const frame = useRef(0);

  useEffect(() => {
    if (!playing) return;
    const { ptRange } = useUi.getState();
    const from = progressFor(ptRange) >= 0.999 ? 0 : progressFor(ptRange);
    const duration = PLAY_DURATION * (1 - from);
    const start = performance.now();
    const tick = (now: number) => {
      const p = from + (1 - from) * Math.min(1, (now - start) / duration);
      useUi.getState().setPtRange(playWindow(p));
      if (p < 1) frame.current = requestAnimationFrame(tick);
      else useUi.getState().setPlaying(false);
    };
    frame.current = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(frame.current);
  }, [playing]);
}
