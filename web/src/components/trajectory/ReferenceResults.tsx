import { Check, Minus, X } from "lucide-react";
import type * as T from "@/api/generated";

type Row = T.ReferenceRow;

const STATUS = {
  pass: { Icon: Check, text: "Pass" },
  fail: { Icon: X, text: "Fail" },
  unmeasured: { Icon: Minus, text: "Not measured" },
} as const;

const status = (r: Row) => (r.passed == null ? STATUS.unmeasured : r.passed ? STATUS.pass : STATUS.fail);

/** Rounded down to three places, so a score under a bound never reads as the bound. */
const formatScore = (v: number | null) =>
  v == null ? "–" : (Math.floor(Number((v * 1000).toFixed(6))) / 1000).toFixed(3).replace("-", "−");

/** The method's published benchmark results, one row per gate, grouped by dataset. */
export function ReferenceResults({ rows }: { rows: Row[] }) {
  if (rows.length === 0) return <p className="text-xs text-ink-muted">This method has no published benchmark results.</p>;
  const groups = new Map<string, Row[]>();
  for (const row of rows) groups.set(row.dataset, [...(groups.get(row.dataset) ?? []), row]);
  return (
    <>
      <p className="mb-3 text-xs leading-snug text-ink-muted">
        Reference results from this method's own benchmark. They are not a score for your dataset.
      </p>
      {[...groups].map(([dataset, group]) => (
        <div key={dataset} className="mb-3 last:mb-0">
          <h4 className="text-sm font-medium text-ink">{dataset}</h4>
          <p className="font-mono text-2xs break-words text-ink-faint">{group[0]?.source}</p>
          <ul aria-label={`${dataset}: reference results`} className="mt-1.5 space-y-2">
            {group.map((r) => {
              const { Icon, text } = status(r);
              return (
                <li key={r.gate} className="grid grid-cols-[16px_1fr_auto] gap-x-2 text-sm">
                  <Icon size={15} aria-hidden className="mt-0.5 text-ink" />
                  <span className="text-ink">
                    {r.metric} <span className="font-mono text-2xs text-ink-faint">{r.gate}</span>
                  </span>
                  <span className="font-mono text-xs text-ink">{formatScore(r.score)}</span>
                  <span className="col-start-2 col-end-4 text-xs leading-snug text-ink-muted">
                    <span className="font-medium text-ink">{text}</span> · needs {r.threshold}
                    {r.detail && ` · ${r.detail}`}
                  </span>
                </li>
              );
            })}
          </ul>
        </div>
      ))}
    </>
  );
}
