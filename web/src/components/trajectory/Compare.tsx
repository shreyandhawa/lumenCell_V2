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
      <Select aria-label="Compare with" value={b ?? ""} onChange={(e) => setOther(e.target.value)}>
        {others.map((r) => (
          <option key={r.run_id} value={r.run_id}>
            {r.method} · {r.run_id.slice(0, 6)}
          </option>
        ))}
      </Select>
      {cmp.data && runA && runB && (
        <div className="mt-2 space-y-2">
          <p className="text-xs text-ink">
            Pseudotime rank agreement{" "}
            <span className="font-mono">{cmp.data.spearman === null ? "n/a" : cmp.data.spearman.toFixed(3)}</span>
            <span className="text-ink-faint"> over {fmt(cmp.data.n_common)} cells in both</span>
          </p>
          <table className="w-full text-2xs" aria-label="Matched lineages">
            <thead>
              <tr className="text-left text-ink-faint">
                <th className="font-normal">this run</th>
                <th className="font-normal">{runB.provenance.method}</th>
                <th className="text-right font-normal">overlap</th>
              </tr>
            </thead>
            <tbody>
              {cmp.data.lineage_matches.map((m) => (
                <tr key={`${m.a}-${m.b}`}>
                  <td className="truncate text-ink">{runA.lineages[m.a]?.name}</td>
                  <td className="truncate text-ink">{runB.lineages[m.b]?.name}</td>
                  <td className="text-right font-mono text-ink">{(m.jaccard * 100).toFixed(0)}%</td>
                </tr>
              ))}
            </tbody>
          </table>
          <Button onClick={() => setColorBy({ kind: "disagreement", a: active, b: b! })}>
            <GitCompareArrows size={12} /> Colour cells by disagreement
          </Button>
          <Hint>Lineages are paired by cell overlap, so names and order may differ between methods.</Hint>
        </div>
      )}
    </Section>
  );
}
