import { Flag, Info, TriangleAlert } from "lucide-react";
import { useState } from "react";
import type { RootSpec } from "@/api/client";
import type * as T from "@/api/generated";
import { useDataset, useRootSuggestions, useSetRoot } from "@/api/queries";
import { Badge, Button, Hint, Section, Select } from "@/components/ui";

const LEIDEN = "lumencell_leiden";

const SOURCE: Record<T.RootSuggestion["source"], string> = {
  file: "stored in file",
  time_column: "earliest time point",
  gene_count: "heuristic",
};

function describe(spec: RootSpec): string {
  if (spec.kind === "cluster") return spec.key === LEIDEN ? `Cluster ${spec.value}` : `${spec.key} = ${spec.value}`;
  if (spec.kind === "cell") return `cell ${spec.obs_name}`;
  return "the root stored in the file";
}

/** Pick where the trajectory starts. Nothing is chosen for the user. */
export function RootPicker({ onDone }: { onDone?: () => void }) {
  const dataset = useDataset().data;
  const suggestions = useRootSuggestions().data?.suggestions ?? [];
  const setRoot = useSetRoot();
  const columns = (dataset?.obs ?? []).filter((c) => c.kind === "categorical" && c.categories?.length);
  const [column, setColumn] = useState<string>(dataset?.label_key ?? LEIDEN);
  const values = columns.find((c) => c.name === column)?.categories ?? [];
  const [value, setValue] = useState<string>("");

  const choose = (spec: RootSpec) =>
    setRoot.mutate(spec, {
      onSuccess: () => onDone?.(),
    });

  return (
    <Section title="Where does the trajectory start?">
      <p className="mb-3 text-xs leading-relaxed text-ink-muted">
        Pseudotime is measured from a starting population you choose. lumenCell never guesses it: a
        wrong guess silently reverses the trajectory.
      </p>

      {suggestions.length > 0 && (
        <div className="mb-3 space-y-2" aria-label="Suggested roots">
          {suggestions.map((s, i) => (
            <div key={i} className="rounded-md border border-line bg-panel-2 p-2">
              <div className="mb-1 flex items-center gap-1.5">
                <span className="text-xs font-medium text-ink">{describe(s.spec)}</span>
                <Badge tone={s.source === "gene_count" ? "warning" : "neutral"}>{SOURCE[s.source]}</Badge>
                {s.is_extremity && <Badge>end of graph</Badge>}
              </div>
              <Hint>{s.evidence}</Hint>
              <Button className="mt-1.5" onClick={() => choose(s.spec)} disabled={setRoot.isPending}>
                <Flag size={12} /> Start here
              </Button>
            </div>
          ))}
        </div>
      )}

      <div className="space-y-1.5" aria-label="Choose a cluster">
        <label className="block text-2xs text-ink-muted">Or choose a starting cluster</label>
        <div className="flex gap-1.5">
          <Select
            aria-label="Cluster column"
            value={column}
            onChange={(e) => {
              setColumn(e.target.value);
              setValue("");
            }}
          >
            {columns.map((c) => (
              <option key={c.name} value={c.name}>
                {c.name === LEIDEN ? "Leiden clusters" : c.name}
              </option>
            ))}
          </Select>
          <Select aria-label="Starting cluster" value={value} onChange={(e) => setValue(e.target.value)}>
            <option value="">pick…</option>
            {values.map((v) => (
              <option key={v} value={v}>
                {column === LEIDEN ? `Cluster ${v}` : v}
              </option>
            ))}
          </Select>
        </div>
        <Button
          tone="primary"
          disabled={!value || setRoot.isPending}
          onClick={() => choose({ kind: "cluster", key: column, value })}
        >
          Set root
        </Button>
      </div>
      {setRoot.error && (
        <p role="alert" className="mt-2 flex gap-1.5 text-2xs text-danger">
          <TriangleAlert size={12} className="mt-px shrink-0" />
          {setRoot.error.message}
        </p>
      )}
    </Section>
  );
}

export function RootSummary({ root, onChange }: { root: T.ResolvedRoot; onChange: () => void }) {
  return (
    <Section title="Root" aside={<Button onClick={onChange}>Change</Button>}>
      <p className="text-xs text-ink">
        <Flag size={12} className="mr-1 inline text-accent" />
        {describe(root.spec)}
      </p>
      <Hint>
        Cell <span className="font-mono">{root.obs_name}</span>: {root.resolved_by}.
      </Hint>
      {root.warnings.map((w) => (
        <p key={w.code} className="mt-1.5 flex gap-1.5 text-2xs text-warning">
          <Info size={12} className="mt-px shrink-0" />
          {w.message}
        </p>
      ))}
    </Section>
  );
}
