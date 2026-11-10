import { useMemo } from "react";
import { useDataset, useDoneRunId, useObs, useRunCells } from "@/api/queries";
import { Legend } from "@/components/Legend";
import { useCellColors } from "@/components/useCellColors";
import { Button, fmt, Hint, Kbd, Section } from "@/components/ui";
import { categoryLabel, LEIDEN } from "@/lib/colorBy";
import { useUi } from "@/state/ui";

/** The Cells panel: what the colours mean, and what is selected. */
export function Sidebar() {
  const dataset = useDataset().data;
  const { legend, loading } = useCellColors();
  if (!dataset) return null;

  return (
    <div className="flex h-full flex-col">
      <Section title="Cells">
        {loading ? (
          <div className="space-y-2" aria-hidden>
            <div className="skeleton h-3 w-24" />
            <div className="skeleton h-3 w-36" />
            <div className="skeleton h-3 w-28" />
          </div>
        ) : legend ? (
          <Legend spec={legend} />
        ) : (
          <Hint>Choose what to colour by in the toolbar, or search for a gene.</Hint>
        )}
      </Section>
      <SelectionSummary />
      <div className="mt-auto px-4 pt-3 pb-4">
        <Hint>
          <span className="num">{fmt(dataset.n_obs)}</span> cells · <span className="num">{fmt(dataset.n_vars)}</span>{" "}
          genes · X holds {dataset.x_kind === "counts" ? "raw counts" : "log-normalised values"}
        </Hint>
      </div>
    </div>
  );
}

function SelectionSummary() {
  const selection = useUi((s) => s.selection);
  const count = useUi((s) => s.selectionCount);
  const setSelection = useUi((s) => s.setSelection);
  const setLasso = useUi((s) => s.setLasso);
  const runId = useDoneRunId(useUi((s) => s.runId));
  const dataset = useDataset().data;
  const groupKey = dataset?.label_key ?? LEIDEN;
  const groups = useObs(selection ? groupKey : null).data;
  const cells = useRunCells(selection ? runId : null).data;

  const summary = useMemo(() => {
    if (!selection) return null;
    const tally = new Map<string, number>();
    let ptSum = 0;
    let ptN = 0;
    for (let i = 0; i < selection.length; i++) {
      if (!selection[i]) continue;
      if (groups?.codes) {
        const label = groups.categories[groups.codes[i]!] ?? "missing";
        tally.set(label, (tally.get(label) ?? 0) + 1);
      }
      const pt = cells?.pseudotime[i];
      if (pt !== undefined && Number.isFinite(pt)) {
        ptSum += pt;
        ptN++;
      }
    }
    const top = [...tally.entries()].sort((a, b) => b[1] - a[1]).slice(0, 5);
    return { top, meanPt: ptN ? ptSum / ptN : null };
  }, [selection, groups, cells]);

  if (!selection || !summary) {
    return (
      <Section title="Selection">
        <Hint>
          Draw around cells with the lasso <Kbd>L</Kbd> to see what they are.
        </Hint>
        <Button size="sm" tone="ghost" className="mt-2 -ml-2.5" onClick={() => setLasso(true)}>
          Start a lasso
        </Button>
      </Section>
    );
  }
  return (
    <Section
      title="Selection"
      aside={
        <Button size="sm" tone="ghost" onClick={() => setSelection(null)}>
          Clear
        </Button>
      }
    >
      <p className="mb-3 text-sm text-ink">
        <span className="num">{fmt(count)}</span> cells
        {summary.meanPt !== null && (
          <span className="text-ink-2">
            {" "}
            · mean pseudotime <span className="num">{summary.meanPt.toFixed(3)}</span>
          </span>
        )}
      </p>
      <ul className="space-y-2">
        {summary.top.map(([label, k]) => (
          <li key={label} className="text-xs text-ink">
            <div className="flex justify-between gap-2">
              <span className="truncate">{categoryLabel(groupKey, label)}</span>
              <span className="num text-ink-2">{((100 * k) / count).toFixed(0)}%</span>
            </div>
            <div className="mt-1 h-1 overflow-hidden rounded-full bg-control">
              <div
                className="h-1 origin-left rounded-full bg-focus transition-transform duration-300 ease-out"
                style={{ transform: `scaleX(${k / count})` }}
              />
            </div>
          </li>
        ))}
      </ul>
    </Section>
  );
}
