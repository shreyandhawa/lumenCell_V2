import { Flag } from "lucide-react";
import { useEffect, useState } from "react";
import type { RootSpec } from "@/api/client";
import type * as T from "@/api/generated";
import { useDataset, useRootSuggestions, useSetRoot } from "@/api/queries";
import { Button, Hint, Note, Section, Select, Tag } from "@/components/ui";
import { categoryLabel, columnLabel, LEIDEN } from "@/lib/colorBy";
import { useUi } from "@/state/ui";

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

/** Pick where the trajectory starts. Nothing is chosen for the user; pointing at a choice previews its cells. */
export function RootPicker({ onDone, onCancel }: { onDone?: () => void; onCancel?: () => void }) {
  const dataset = useDataset().data;
  const suggestions = useRootSuggestions().data?.suggestions ?? [];
  const setRoot = useSetRoot();
  const setHighlight = useUi((s) => s.setHighlight);
  const columns = (dataset?.obs ?? []).filter((c) => c.kind === "categorical" && c.categories?.length);
  const [column, setColumn] = useState<string>(dataset?.label_key ?? LEIDEN);
  const values = columns.find((c) => c.name === column)?.categories ?? [];
  const [value, setValue] = useState<string>("");

  // The chosen cluster stays lit on the field until it is set or the picker closes.
  useEffect(() => {
    setHighlight(value ? { kind: "category", column, value } : null);
  }, [column, value, setHighlight]);
  useEffect(() => () => setHighlight(null), [setHighlight]);

  const preview = (spec: RootSpec | null) => {
    if (spec?.kind === "cluster") setHighlight({ kind: "category", column: spec.key, value: spec.value });
    else setHighlight(value ? { kind: "category", column, value } : null);
  };
  const choose = (spec: RootSpec) => setRoot.mutate(spec, { onSuccess: () => onDone?.() });

  return (
    <Section title="Where does the trajectory start?">
      <p className="mb-4 text-sm text-ink-2">
        Pseudotime is measured from a starting population you choose. lumenCell never guesses it: a wrong guess
        silently reverses the trajectory.
      </p>

      {suggestions.length > 0 && (
        <ul className="-mx-4 mb-5 divide-y divide-line border-y border-line" aria-label="Suggested roots">
          {suggestions.map((s, i) => (
            <li
              key={i}
              className="px-4 py-3 transition-colors duration-150 [@media(hover:hover)]:hover:bg-control/60"
              onMouseEnter={() => preview(s.spec)}
              onMouseLeave={() => preview(null)}
            >
              <div className="mb-1 flex flex-wrap items-center gap-1.5">
                <span className="mr-1 text-sm font-semibold text-ink">{describe(s.spec)}</span>
                <Tag tone={s.source === "gene_count" ? "warn" : "neutral"}>{SOURCE[s.source]}</Tag>
                {s.is_extremity && <Tag>end of graph</Tag>}
              </div>
              <Hint>{s.evidence}</Hint>
              <Button
                size="sm"
                className="mt-2"
                onClick={() => choose(s.spec)}
                onFocus={() => preview(s.spec)}
                onBlur={() => preview(null)}
                disabled={setRoot.isPending}
              >
                <Flag size={12} /> Start here
              </Button>
            </li>
          ))}
        </ul>
      )}

      <div aria-label="Choose a cluster">
        <h3 className="mb-2 text-sm font-medium text-ink">Or choose a starting cluster</h3>
        <div className="flex gap-2">
          <Select
            label="Cluster column"
            className="w-1/2"
            value={column}
            onChange={(c) => {
              setColumn(c);
              setValue("");
            }}
            groups={[{ items: columns.map((c) => ({ value: c.name, label: columnLabel(c.name) })) }]}
          />
          <Select
            label="Starting cluster"
            className="w-1/2"
            value={value || null}
            placeholder="Pick one"
            onChange={setValue}
            groups={[{ items: values.map((v) => ({ value: v, label: categoryLabel(column, v) })) }]}
          />
        </div>
        <div className="mt-3 flex items-center gap-2">
          <Button
            tone="primary"
            disabled={!value || setRoot.isPending}
            onClick={() => choose({ kind: "cluster", key: column, value })}
          >
            Set as root
          </Button>
          {onCancel && (
            <Button tone="ghost" onClick={onCancel}>
              Keep the current root
            </Button>
          )}
        </div>
        {value && <Hint className="mt-2">Its cells are lit on the field.</Hint>}
      </div>
      {setRoot.error && (
        <div className="mt-3">
          <Note tone="error" role="alert">
            {setRoot.error.message}
          </Note>
        </div>
      )}
    </Section>
  );
}

export function RootSummary({ root, onChange }: { root: T.ResolvedRoot; onChange: () => void }) {
  return (
    <Section
      title="Root"
      aside={
        <Button size="sm" tone="ghost" onClick={onChange}>
          Change
        </Button>
      }
    >
      <p className="flex items-center gap-2 text-sm text-ink">
        <Flag size={13} className="shrink-0 text-live" />
        {describe(root.spec)}
      </p>
      <Hint className="mt-1">
        Cell <span className="num">{root.obs_name}</span>: {root.resolved_by}.
      </Hint>
      {root.warnings.length > 0 && (
        <div className="mt-2 space-y-1.5">
          {root.warnings.map((w) => (
            <Note key={w.code} tone="warning">
              {w.message}
            </Note>
          ))}
        </div>
      )}
    </Section>
  );
}
