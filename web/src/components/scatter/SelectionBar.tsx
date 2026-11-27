import { fmt } from "@/components/ui";
import { composition, type SelectOp } from "@/lib/selection";
import { useSelection } from "@/state/selection";

/** What a selection is broken down by: the shown categorical colouring, lineage, or Leiden cluster. */
export interface GroupBy {
  title: string;
  codes: ArrayLike<number>;
  labels: string[];
}

const OPS: [SelectOp, string][] = [
  ["new", "New"],
  ["add", "Add"],
  ["subtract", "Subtract"],
  ["intersect", "Intersect"],
];
const TOP = 5;

/** The selection over the plot: how many cells, how the next lasso combines, and what they are. */
export function SelectionBar({ lasso, groupBy, total }: { lasso: boolean; groupBy: GroupBy | null; total: number }) {
  const { mask, count, op, setOp, setMask } = useSelection();
  if (!lasso && !mask) return null;

  const rows = mask && groupBy ? breakdown(mask, groupBy) : [];
  return (
    <div role="group" aria-label="Selection" className="lift w-72 rounded-xl bg-panel p-3 text-sm text-ink">
      <div className="flex items-baseline justify-between gap-2">
        <p aria-live="polite">
          {mask ? (
            <>
              <span className="font-mono">{fmt(count)}</span> cells selected
              <span className="text-ink-muted"> · {((100 * count) / total).toFixed(1)}%</span>
            </>
          ) : (
            "No cells selected"
          )}
        </p>
        {mask && (
          <button type="button" onClick={() => setMask(null)} className="press text-xs text-ink-muted hover:text-ink">
            Clear
          </button>
        )}
      </div>

      {lasso && (
        <>
          <div role="radiogroup" aria-label="Next lasso" className="mt-2 flex gap-0.5 rounded-lg bg-panel-raised p-0.5">
            {OPS.map(([v, label]) => (
              <button
                key={v}
                type="button"
                role="radio"
                aria-checked={op === v}
                onClick={() => setOp(v)}
                className={`press h-6 flex-1 rounded-md text-xs ${op === v ? "bg-panel text-ink" : "text-ink-muted hover:text-ink"}`}
              >
                {label}
              </button>
            ))}
          </div>
          <p className="mt-2 text-xs leading-snug text-ink-muted">
            Drag to draw around cells; in 3-D a lasso takes every cell behind it too, so turn the view (Shift +
            arrows) and Intersect to narrow it. Alt-click a cell for its whole {groupBy?.title ?? "group"}. Esc clears.
          </p>
        </>
      )}

      {rows.length > 0 && (
        <dl aria-label="Selection by group" className="mt-3 grid grid-cols-[1fr_auto] gap-x-3 gap-y-1 text-xs">
          {rows.map((r) => (
            <div key={r.label} className="contents">
              <dt className="truncate text-ink">{r.label}</dt>
              <dd className="text-right font-mono text-ink-muted">
                <span className="text-ink">{r.ofSelection}%</span> · {r.ofGroup}% of it
              </dd>
            </div>
          ))}
        </dl>
      )}
    </div>
  );
}

/** The largest groups in the selection: their share of the selection, and how much of each group it took. */
function breakdown(mask: Uint8Array, g: GroupBy) {
  const { selected, all } = composition(mask, g.codes, g.labels.length);
  const inside = selected.reduce((a, b) => a + b, 0) || 1;
  return selected
    .map((s, c) => ({ s, c }))
    .filter(({ s }) => s > 0)
    .sort((a, b) => b.s - a.s)
    .slice(0, TOP)
    .map(({ s, c }) => ({
      label: g.labels[c]!,
      ofSelection: percent(s / inside),
      ofGroup: percent(s / (all[c] || 1)),
    }));
}

/** A share as a whole percent; a few cells read "<1", never "0". */
const percent = (x: number) => (x > 0 && x < 0.005 ? "<1" : String(Math.round(100 * x)));
