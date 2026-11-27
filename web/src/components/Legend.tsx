import { rgbaCss, sequential } from "@/lib/colors";
import type { LegendSpec } from "@/components/useCellColors";

const MAX_ITEMS = 14;

export function Legend({ spec }: { spec: LegendSpec | null }) {
  if (!spec) return null;
  return (
    <div
      aria-label="Legend"
      className="lift pointer-events-auto flex max-w-[560px] flex-wrap items-center gap-x-3 gap-y-1 rounded-lg bg-panel px-2.5 py-1.5 text-xs"
    >
      <span className="text-ink-muted">{spec.title}</span>
      {spec.kind === "continuous" ? (
        <span className="flex items-center gap-2 font-mono text-ink">
          {spec.domain[0].toPrecision(3)}
          <span
            className="h-1.5 w-28 rounded-full"
            style={{
              background: `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map((t) => rgbaCss(sequential(t, 255))).join(",")})`,
            }}
          />
          {spec.domain[1].toPrecision(3)}
        </span>
      ) : (
        <>
          {spec.items.slice(0, MAX_ITEMS).map((item) => (
            <span key={item.label} className="flex items-center gap-1.5 text-ink">
              <span className="size-2 shrink-0 rounded-full" style={{ background: rgbaCss(item.color) }} />
              <span className="max-w-32 truncate font-mono">{item.label}</span>
            </span>
          ))}
          {spec.items.length > MAX_ITEMS && <span className="text-ink-muted">and {spec.items.length - MAX_ITEMS} more</span>}
        </>
      )}
    </div>
  );
}
