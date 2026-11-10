import { useMemo } from "react";
import { useDataset, useObs, useRunCells, useDoneRunId } from "@/api/queries";
import { GeneSearch } from "@/components/GeneSearch";
import { Badge, fmt, Hint, Section, Select } from "@/components/ui";
import { type ColorBy, useUi } from "@/state/ui";

const LEIDEN = "lumencell_leiden";

function encode(c: ColorBy | null): string {
  if (!c) return "";
  if (c.kind === "obs") return `obs:${c.column}`;
  if (c.kind === "gene") return `gene:${c.gene}`;
  if (c.kind === "disagreement") return `disagreement:${c.a}:${c.b}`;
  return c.kind;
}

function decode(v: string): ColorBy | null {
  if (v.startsWith("obs:")) return { kind: "obs", column: v.slice(4) };
  if (v.startsWith("gene:")) return { kind: "gene", gene: v.slice(5) };
  if (v === "pseudotime" || v === "lineage") return { kind: v };
  if (v.startsWith("disagreement:")) {
    const [, a = "", b = ""] = v.split(":");
    return { kind: "disagreement", a, b };
  }
  return null;
}

export function Sidebar() {
  const dataset = useDataset().data;
  const embedding = useUi((s) => s.embedding);
  const setEmbedding = useUi((s) => s.setEmbedding);
  const colorBy = useUi((s) => s.colorBy);
  const setColorBy = useUi((s) => s.setColorBy);
  const runId = useDoneRunId(useUi((s) => s.runId));
  if (!dataset) return null;

  const columns = dataset.obs.filter((c) => c.kind !== "string");
  const embeddingSource = dataset.embeddings.find((e) => e.name === embedding)?.source;

  return (
    <aside className="flex w-64 shrink-0 flex-col overflow-y-auto border-r border-line bg-panel" aria-label="View controls">
      <Section title="Embedding">
        <Select aria-label="Embedding" value={embedding ?? ""} onChange={(e) => setEmbedding(e.target.value)}>
          {dataset.embeddings.length === 0 && <option value="">waiting for preprocessing…</option>}
          {dataset.embeddings.map((e) => (
            <option key={e.name} value={e.name}>
              {e.name === "dc12" ? "DC1 vs DC2 (diffusion)" : e.name.toUpperCase()}
              {e.source === "computed" ? " · computed" : ""}
            </option>
          ))}
        </Select>
        <div className="mt-1.5">
          <Hint>
            {embeddingSource === "diffusion"
              ? "The first two diffusion components: the space the trajectory is inferred in."
              : "For display only. Trajectories are inferred in PCA / diffusion space and drawn here as a projection."}
          </Hint>
        </div>
      </Section>

      <Section title="Colour by">
        <Select aria-label="Colour by" value={encode(colorBy)} onChange={(e) => setColorBy(decode(e.target.value))}>
          <option value="">nothing</option>
          <optgroup label="Trajectory">
            <option value="pseudotime" disabled={!runId}>
              Pseudotime{runId ? "" : " (run a method first)"}
            </option>
            <option value="lineage" disabled={!runId}>
              Lineage{runId ? "" : " (run a method first)"}
            </option>
          </optgroup>
          <optgroup label="Annotations">
            {columns.map((c) => (
              <option key={c.name} value={`obs:${c.name}`}>
                {c.name === LEIDEN ? "Leiden clusters (computed)" : c.name}
              </option>
            ))}
          </optgroup>
          {colorBy?.kind === "gene" && (
            <optgroup label="Gene">
              <option value={encode(colorBy)}>{colorBy.gene}</option>
            </optgroup>
          )}
          {colorBy?.kind === "disagreement" && (
            <optgroup label="Comparison">
              <option value={encode(colorBy)}>Where the two runs disagree</option>
            </optgroup>
          )}
        </Select>
        <GeneSearch onPick={(gene) => setColorBy({ kind: "gene", gene })} />
      </Section>

      <SelectionSummary />

      <div className="mt-auto border-t border-line px-3 py-2">
        <Hint>
          {fmt(dataset.n_obs)} cells · {fmt(dataset.n_vars)} genes · X holds{" "}
          {dataset.x_kind === "counts" ? "raw counts" : "log-normalised values"}
        </Hint>
      </div>
    </aside>
  );
}

function SelectionSummary() {
  const selection = useUi((s) => s.selection);
  const count = useUi((s) => s.selectionCount);
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
        <Hint>Use the lasso (top left of the plot) to select cells.</Hint>
      </Section>
    );
  }
  return (
    <Section title="Selection" aside={<Badge tone="accent">{fmt(count)} cells</Badge>}>
      {summary.meanPt !== null && (
        <p className="mb-2 text-xs text-ink">
          mean pseudotime <span className="font-mono">{summary.meanPt.toFixed(3)}</span>
        </p>
      )}
      <ul className="space-y-1">
        {summary.top.map(([label, k]) => (
          <li key={label} className="text-2xs text-ink-muted">
            <div className="flex justify-between">
              <span className="truncate">{groupKey === LEIDEN ? `Cluster ${label}` : label}</span>
              <span className="font-mono">{((100 * k) / count).toFixed(0)}%</span>
            </div>
            <div className="mt-0.5 h-1 rounded-full bg-line">
              <div className="h-1 rounded-full bg-accent" style={{ width: `${(100 * k) / count}%` }} />
            </div>
          </li>
        ))}
      </ul>
    </Section>
  );
}
