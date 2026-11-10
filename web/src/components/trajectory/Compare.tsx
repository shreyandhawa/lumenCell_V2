import { GitCompareArrows } from "lucide-react";
import { useState } from "react";
import { useComparison, useDoneRunId, useRun, useRuns } from "@/api/queries";
import { Button, fmt, Hint, Section, Select } from "@/components/ui";
import { useUi } from "@/state/ui";

/** How much the active run agrees with another one. */
export function Compare() {
  const active = useDoneRunId(useUi((s) => s.runId));
  const setColorBy = useUi((s) => s.setColorBy);
  const done = (useRuns().data?.runs ?? []).filter((r) => r.state === "done");
  const others = done.filter((r) => r.run_id !== active);
  const [other, setOther] = useState<string>("");
  const b = other && others.some((r) => r.run_id === other) ? other : (others[0]?.run_id ?? null);
  const cmp = useComparison(active, b);
  const runA = useRun(active).data;
  const runB = useRun(b).data;
  if (!active || others.length === 0) return null;

  return (
    <Section title="Compare runs">
      <Select
        label="Compare with"
        className="w-full"
        value={b}
        onChange={setOther}
        groups={[{ items: others.map((r) => ({ value: r.run_id, label: `${r.method} · ${r.run_id.slice(0, 6)}` })) }]}
      />
      {cmp.data && runA && runB && (
        <div className="mt-3 space-y-3">
          <p className="text-sm text-ink">
            Pseudotime rank agreement{" "}
            <span className="num font-semibold">{cmp.data.spearman === null ? "n/a" : cmp.data.spearman.toFixed(3)}</span>
            <span className="text-ink-3">
              {" "}
              over <span className="num">{fmt(cmp.data.n_common)}</span> cells in both
            </span>
          </p>
          <table className="w-full text-xs" aria-label="Matched lineages">
            <thead>
              <tr className="text-left text-ink-3">
                <th className="pb-1 font-medium">This run</th>
                <th className="pb-1 font-medium">{runB.provenance.method}</th>
                <th className="pb-1 text-right font-medium">Overlap</th>
              </tr>
            </thead>
            <tbody>
              {cmp.data.lineage_matches.map((m) => (
                <tr key={`${m.a}-${m.b}`} className="border-t border-line">
                  <td className="truncate py-1 text-ink">{runA.lineages[m.a]?.name}</td>
                  <td className="truncate py-1 text-ink">{runB.lineages[m.b]?.name}</td>
                  <td className="num py-1 text-right text-ink">{(m.jaccard * 100).toFixed(0)}%</td>
                </tr>
              ))}
            </tbody>
          </table>
          <Button size="sm" onClick={() => setColorBy({ kind: "disagreement", a: active, b: b! })}>
            <GitCompareArrows size={13} /> Colour cells by disagreement
          </Button>
          <Hint>Lineages are paired by cell overlap, so names and order may differ between methods.</Hint>
        </div>
      )}
    </Section>
  );
}
