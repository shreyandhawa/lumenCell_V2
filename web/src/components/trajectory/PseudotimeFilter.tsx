import { Button, Hint, Section } from "@/components/ui";
import { useUi } from "@/state/ui";

/** Show only cells within a pseudotime window (filtered on the GPU). */
export function PseudotimeFilter() {
  const [lo, hi] = useUi((s) => s.ptRange);
  const setRange = useUi((s) => s.setPtRange);
  const active = lo > 0 || hi < 1;
  return (
    <Section
      title="Pseudotime window"
      aside={active ? <Button onClick={() => setRange([0, 1])}>Reset</Button> : undefined}
    >
      <div className="relative h-5">
        <div className="absolute top-2 h-1 w-full rounded-full bg-line" />
        <div className="absolute top-2 h-1 rounded-full bg-accent" style={{ left: `${lo * 100}%`, width: `${(hi - lo) * 100}%` }} />
        <input
          type="range"
          className="range"
          aria-label="Pseudotime from"
          min={0}
          max={1}
          step={0.01}
          value={lo}
          onChange={(e) => setRange([Math.min(Number(e.target.value), hi - 0.01), hi])}
        />
        <input
          type="range"
          className="range"
          aria-label="Pseudotime to"
          min={0}
          max={1}
          step={0.01}
          value={hi}
          onChange={(e) => setRange([lo, Math.max(Number(e.target.value), lo + 0.01)])}
        />
      </div>
      <div className="mt-1 flex justify-between font-mono text-2xs text-ink-muted">
        <span>{lo.toFixed(2)}</span>
        <span>{hi.toFixed(2)}</span>
      </div>
      {active && <Hint>Cells without a pseudotime are hidden while a window is set.</Hint>}
    </Section>
  );
}
