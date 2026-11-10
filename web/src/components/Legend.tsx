import { rgbaCss, sequential } from "@/lib/colors";
import type { LegendSpec } from "@/components/useCellColors";

const MAX_ITEMS = 14;

export function Legend({ spec }: { spec: LegendSpec | null }) {
  if (!spec) return null;
  return (
    <div className="pointer-events-auto max-w-56 rounded-md border border-line bg-panel/90 px-2.5 py-2 shadow-lg backdrop-blur">
      <div className="mb-1.5 truncate text-2xs font-semibold text-ink-muted" title={spec.title}>
        {spec.title}
      </div>
      {spec.kind === "continuous" ? (
        <div>
          <div
            className="h-2 w-40 rounded-sm"
            style={{
              background: `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map((t) => rgbaCss(sequential(t, 255))).join(",")})`,
            }}
          />
          <div className="mt-1 flex justify-between font-mono text-2xs text-ink-faint">
            <span>{spec.domain[0].toPrecision(3)}</span>
            <span>{spec.domain[1].toPrecision(3)}</span>
          </div>
        </div>
      ) : (
        <ul className="space-y-0.5">
          {spec.items.slice(0, MAX_ITEMS).map((item) => (
            <li key={item.label} className="flex items-center gap-1.5 text-2xs text-ink">
              <span className="size-2 shrink-0 rounded-full" style={{ background: rgbaCss(item.color) }} />
              <span className="truncate">{item.label}</span>
            </li>
          ))}
          {spec.items.length > MAX_ITEMS && (
            <li className="text-2xs text-ink-faint">and {spec.items.length - MAX_ITEMS} more</li>
          )}
        </ul>
      )}
    </div>
  );
}
