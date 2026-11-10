import { Pencil } from "lucide-react";
import { useState } from "react";
import type * as T from "@/api/generated";
import { useRenameLineage } from "@/api/queries";
import { cx, fmt, Hint, Note, Section } from "@/components/ui";
import { categorical, rgbaCss, SHARED } from "@/lib/colors";
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
  const setHighlight = useUi((s) => s.setHighlight);
  return (
    <Section
      title="Lineages"
      aside={focus !== null ? <span className="text-xs text-ink-3">Showing one; click it again for all</span> : undefined}
    >
      <ul className="-mx-2 space-y-0.5" aria-label="Lineages" onMouseLeave={() => setHighlight(null)}>
        {run.lineages.map((lin) => {
          const fork = forkPoint(run, lin.id);
          const color = rgbaCss(categorical(lin.id, 255));
          const active = focus === lin.id;
          const dimmed = focus !== null && !active;
          return (
            <li
              key={lin.id}
              className={cx(
                "group relative rounded-lg px-2 py-2 transition-[background-color,opacity] duration-200",
                active ? "bg-control" : "[@media(hover:hover)]:hover:bg-control/60",
                dimmed && "opacity-55",
              )}
              onMouseEnter={() => setHighlight({ kind: "lineage", id: lin.id })}
            >
              {/* The whole row is the focus target; the rename control sits above it. */}
              <button
                type="button"
                aria-pressed={active}
                aria-label={active ? `Show all lineages` : `Focus ${lin.name}`}
                onClick={() => setFocus(active ? null : lin.id)}
                onFocus={() => setHighlight({ kind: "lineage", id: lin.id })}
                onBlur={() => setHighlight(null)}
                className="pressable absolute inset-0 rounded-lg"
              />
              <div className="pointer-events-none relative flex items-center gap-2">
                <span className="size-2.5 shrink-0 rounded-full" style={{ background: color }} />
                <LineageName runId={run.run_id} lineage={lin} />
                <span className="num ml-auto text-xs text-ink-3">{fmt(lin.n_cells)}</span>
              </div>
              {/* shared stretch then this lineage's own stretch, on the one 0-1 pseudotime axis */}
              <div className="pointer-events-none relative mt-2 flex h-1.5 overflow-hidden rounded-full" aria-hidden>
                <div style={{ width: `${fork * 100}%`, background: rgbaCss(SHARED) }} />
                <div style={{ width: `${(1 - fork) * 100}%`, background: color }} />
              </div>
              <div className="pointer-events-none relative mt-1.5">
                <Hint>
                  {fork > 0 ? (
                    <>
                      splits off at <span className="num">{fork.toFixed(2)}</span>
                    </>
                  ) : (
                    "no shared trunk"
                  )}
                  {lin.terminal_purity != null && lin.terminal_purity < 0.5 ? " · mixed tip" : ""}
                </Hint>
              </div>
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
      <span className="pointer-events-auto flex min-w-0 flex-1 flex-col gap-1">
        <input
          autoFocus
          autoComplete="off"
          spellCheck={false}
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
          className="h-7 min-w-0 rounded-md bg-control px-2 text-sm text-ink outline-2 outline-focus"
        />
        {rename.error && <Note tone="error">{rename.error.message}</Note>}
      </span>
    );
  }
  return (
    <span className="flex min-w-0 items-center gap-1">
      <span className="truncate text-sm font-medium text-ink">{lineage.name}</span>
      <button
        type="button"
        aria-label={`Rename ${lineage.name}`}
        onClick={() => {
          setDraft(lineage.name);
          setEditing(true);
        }}
        className="pressable pointer-events-auto grid size-6 place-items-center rounded-md text-ink-3 opacity-0 group-hover:opacity-100 focus-visible:opacity-100 [@media(hover:hover)]:hover:bg-control-hover [@media(hover:hover)]:hover:text-ink"
      >
        <Pencil size={12} />
      </button>
    </span>
  );
}
