import { Pencil } from "lucide-react";
import { useState } from "react";
import type * as T from "@/api/generated";
import { useRenameLineage } from "@/api/queries";
import { categorical, rgbaCss, SHARED } from "@/lib/colors";
import { fmt, Hint, Section } from "@/components/ui";
import { useUi } from "@/state/ui";

/** The pseudotime at which a lineage last splits from another one. */
function forkPoint(run: T.TrajectoryRun, id: number): number {
  let fork = 0;
  for (const d of run.divergence) if (d.a === id || d.b === id) fork = Math.max(fork, d.pseudotime);
  return fork;
}

export function LineageTree({ run }: { run: T.TrajectoryRun }) {
  const focus = useUi((s) => s.focusLineage);
  const setFocus = useUi((s) => s.setFocusLineage);
  return (
    <Section title="Lineages">
      <ul className="space-y-2" aria-label="Lineages">
        {run.lineages.map((lin) => {
          const fork = forkPoint(run, lin.id);
          const color = rgbaCss(categorical(lin.id, 255));
          const active = focus === lin.id;
          return (
            <li key={lin.id}>
              <div className="flex items-center gap-1.5">
                <button
                  type="button"
                  aria-pressed={active}
                  title={active ? "Show all lineages" : "Focus this lineage"}
                  onClick={() => setFocus(active ? null : lin.id)}
                  className={`size-3 shrink-0 rounded-full ring-offset-2 ring-offset-panel transition ${active ? "ring-2 ring-ink" : ""}`}
                  style={{ background: color }}
                />
                <LineageName runId={run.run_id} lineage={lin} />
                <span className="ml-auto font-mono text-2xs text-ink-faint">{fmt(lin.n_cells)}</span>
              </div>
              {/* shared stretch (grey) then this lineage's own stretch, on a 0-1 pseudotime axis */}
              <div className="mt-1 flex h-1.5 overflow-hidden rounded-full" aria-hidden>
                <div style={{ width: `${fork * 100}%`, background: rgbaCss(SHARED) }} />
                <div style={{ width: `${(1 - fork) * 100}%`, background: color }} />
              </div>
              <Hint>
                {fork > 0 ? `splits off at pseudotime ${fork.toFixed(2)}` : "no shared trunk"}
                {lin.terminal_purity != null && lin.terminal_purity < 0.5 ? " · mixed tip" : ""}
              </Hint>
            </li>
          );
        })}
      </ul>
    </Section>
  );
}

function LineageName({ runId, lineage }: { runId: string; lineage: T.Lineage }) {
  const rename = useRenameLineage(runId);
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(lineage.name);
  if (editing) {
    return (
      <input
        autoFocus
        aria-label={`Rename ${lineage.name}`}
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={() => setEditing(false)}
        onKeyDown={(e) => {
          if (e.key === "Escape") setEditing(false);
          if (e.key === "Enter" && draft.trim()) {
            rename.mutate({ lineage: lineage.id, name: draft.trim() }, { onSuccess: () => setEditing(false) });
          }
        }}
        className="h-6 min-w-0 flex-1 rounded border border-accent bg-panel-2 px-1.5 text-xs text-ink"
      />
    );
  }
  return (
    <span className="group flex min-w-0 items-center gap-1">
      <span className="truncate text-xs font-medium text-ink">{lineage.name}</span>
      <button
        type="button"
        aria-label={`Rename ${lineage.name}`}
        onClick={() => {
          setDraft(lineage.name);
          setEditing(true);
        }}
        className="text-ink-faint opacity-0 transition-opacity group-hover:opacity-100 focus-visible:opacity-100"
      >
        <Pencil size={11} />
      </button>
      {rename.error && <span className="text-2xs text-danger">{rename.error.message}</span>}
    </span>
  );
}
