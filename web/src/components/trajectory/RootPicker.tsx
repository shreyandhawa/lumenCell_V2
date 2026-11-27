import { Check, Crosshair, TriangleAlert } from "lucide-react";
import { useState } from "react";
import type { RootSpec } from "@/api/client";
import type * as T from "@/api/generated";
import { useDataset, useRoot, useRootSuggestions, useSetRoot } from "@/api/queries";
import { Button, Select } from "@/components/ui";

const LEIDEN = "lumencell_leiden";

const SOURCE: Record<T.RootSuggestion["source"], string> = {
  file: "Stored in the file",
  time_column: "Earliest time point",
  gene_count: "Heuristic",
};

export function describeRoot(spec: RootSpec): string {
  if (spec.kind === "cluster") return spec.key === LEIDEN ? `Cluster ${spec.value}` : `${spec.key} = ${spec.value}`;
  if (spec.kind === "cell") return `cell ${spec.obs_name}`;
  return "the root stored in the file";
}

function same(a: RootSpec, b: RootSpec): boolean {
  if (a.kind === "cluster" && b.kind === "cluster") return a.key === b.key && a.value === b.value;
  return a.kind === b.kind;
}

/** Pick where the trajectory starts. Nothing is chosen for the user. */
export function RootPicker() {
  const root = useRoot().data?.root ?? null;
  const suggestions = useRootSuggestions().data?.suggestions ?? [];
  const setRoot = useSetRoot();
  // Any cluster of any categorical column, from the keyboard: suggestions miss declared roots.
  const dataset = useDataset().data;
  const columns = (dataset?.obs ?? []).filter((c) => c.kind === "categorical" && c.categories?.length);
  const [column, setColumn] = useState<string>(dataset?.label_key ?? LEIDEN);
  const values = columns.find((c) => c.name === column)?.categories ?? [];
  const [value, setValue] = useState("");
  // Once a root is chosen, only it stays in view, so the Run step does too.
  const [changing, setChanging] = useState(false);
  const all = root === null || changing;
  const choose = (spec: RootSpec) => setRoot.mutate(spec, { onSuccess: () => setChanging(false) });

  return (
    <>
      <p className="mb-3 text-sm leading-relaxed text-ink">
        <span>Where does the trajectory start?</span> It is never guessed.
      </p>
      {root && !(all && suggestions.some((s) => same(root.spec, s.spec))) && (
        <div className="mb-2 rounded-md border border-ink-muted bg-panel-raised p-2.5">
          <span className="flex items-center gap-2 text-sm text-ink">
            {describeRoot(root.spec)}
            <span className="ml-auto flex items-center gap-1 text-xs">
              <Check size={14} aria-hidden /> Chosen
            </span>
          </span>
          <span className="mt-1 block text-xs leading-snug text-ink-muted">Cell {root.obs_name}: {root.resolved_by}.</span>
        </div>
      )}
      {all && (
        <ul className="space-y-2" aria-label="Suggested roots">
          {suggestions.map((s, i) => {
            const chosen = root !== null && same(root.spec, s.spec);
            return (
              <li key={i}>
                <button
                  type="button"
                  aria-pressed={chosen}
                  disabled={setRoot.isPending}
                  onClick={() => choose(s.spec)}
                  className={`press w-full rounded-md border p-2.5 text-left ${
                    chosen ? "border-ink-muted bg-panel-raised" : "border-line hover:border-ink-muted"
                  }`}
                >
                  <span className="flex items-center gap-2 text-sm text-ink">
                    {describeRoot(s.spec)}
                    {chosen && (
                      <span className="ml-auto flex items-center gap-1 text-xs">
                        <Check size={14} aria-hidden /> Chosen
                      </span>
                    )}
                  </span>
                  <span className="mt-1 block text-xs leading-snug text-ink-muted">
                    {SOURCE[s.source]}: {s.evidence}
                  </span>
                </button>
              </li>
            );
          })}
        </ul>
      )}
      <p className="mt-3 flex items-center gap-2 text-xs text-ink-muted">
        <Crosshair size={14} aria-hidden /> Or click a cell to start from its cluster
      </p>
      {!all && (
        <Button tone="ghost" className="mt-1 -ml-2" onClick={() => setChanging(true)}>
          Other starts
        </Button>
      )}
      {all && (
        <div role="group" aria-label="Choose a cluster" className="mt-3">
          <p className="mb-1.5 text-xs text-ink-muted">Or choose a cluster</p>
          <div className="flex gap-2">
            <Select
              aria-label="Cluster column"
              className="min-w-0 flex-1"
              value={column}
              onChange={(e) => {
                setColumn(e.target.value);
                setValue("");
              }}
            >
              {columns.map((c) => (
                <option key={c.name} value={c.name}>
                  {c.name === LEIDEN ? "Leiden cluster" : c.name}
                </option>
              ))}
            </Select>
            <Select aria-label="Starting cluster" className="min-w-0 flex-1" value={value} onChange={(e) => setValue(e.target.value)}>
              <option value="">pick…</option>
              {values.map((v) => (
                <option key={v} value={v}>
                  {v}
                </option>
              ))}
            </Select>
          </div>
          <Button
            className="mt-2"
            disabled={!value || setRoot.isPending}
            onClick={() => choose({ kind: "cluster", key: column, value })}
          >
            Start here
          </Button>
        </div>
      )}
      {root?.warnings.map((w) => (
        <p key={w.code} className="mt-2 flex gap-1.5 text-xs text-ink">
          <TriangleAlert size={14} aria-hidden className="mt-px shrink-0 text-ink" />
          {w.message}
        </p>
      ))}
      {setRoot.error && (
        <p role="alert" className="mt-2 flex gap-1.5 text-xs text-ink">
          <TriangleAlert size={14} aria-hidden className="mt-px shrink-0 text-ink" />
          {setRoot.error.message}
        </p>
      )}
    </>
  );
}
