import { rgbaCss, sequential } from "@/lib/colors";
import type { LegendSpec } from "@/components/useCellColors";
import { cx } from "@/components/ui";
import { useUi } from "@/state/ui";

/** The colour key. Pointing at a category or lineage lights its cells on the field. */
export function Legend({ spec, max = 24, className }: { spec: LegendSpec | null; max?: number; className?: string }) {
  const setHighlight = useUi((s) => s.setHighlight);
  if (!spec) return null;
  return (
    <div aria-label="Legend" className={className}>
      <div className="mb-2 truncate text-xs font-medium text-ink-2" title={spec.title}>
        {spec.title}
      </div>
      {spec.kind === "continuous" ? (
        <div>
          <div
            className="h-2 rounded-full"
            style={{
              background: `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map((t) => rgbaCss(sequential(t, 255))).join(",")})`,
            }}
          />
          <div className="num mt-1.5 flex justify-between text-2xs text-ink-3">
            <span>{spec.domain[0].toPrecision(3)}</span>
            <span>{spec.domain[1].toPrecision(3)}</span>
          </div>
        </div>
      ) : (
        <ul className="-mx-1.5" onMouseLeave={() => setHighlight(null)}>
          {spec.items.slice(0, max).map((item) => (
            <li key={item.label}>
              <button
                type="button"
                disabled={!item.highlight}
                onMouseEnter={() => setHighlight(item.highlight)}
                onFocus={() => setHighlight(item.highlight)}
                onBlur={() => setHighlight(null)}
                className={cx(
                  "flex w-full items-center gap-2 rounded-md px-1.5 py-1 text-left text-xs text-ink transition-colors duration-150",
                  item.highlight && "[@media(hover:hover)]:hover:bg-control",
                  "disabled:cursor-default",
                )}
              >
                <span className="size-2.5 shrink-0 rounded-full" style={{ background: rgbaCss(item.color) }} />
                <span className="truncate">{item.label}</span>
              </button>
            </li>
          ))}
          {spec.items.length > max && (
            <li className="px-1.5 py-1 text-xs text-ink-3">and {spec.items.length - max} more</li>
          )}
        </ul>
      )}
    </div>
  );
}
