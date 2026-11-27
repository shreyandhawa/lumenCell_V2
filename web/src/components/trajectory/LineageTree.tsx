import type * as T from "@/api/generated";
import { Label } from "@/components/ui";
import { useUi } from "@/state/ui";

/** The run's lineages; focusing one draws it in the accent on the stage and in the gene trends. */
export function LineageTree({ run }: { run: T.TrajectoryRun }) {
  const focus = useUi((s) => s.focusLineage);
  const setFocus = useUi((s) => s.setFocusLineage);
  return (
    <div>
      <Label>Lineages</Label>
      {/* Branches off one trunk, as lines leave an interchange; the focused one is the red station. */}
      <div className="relative">
        <span aria-hidden className="absolute top-4 bottom-4 left-[9.5px] w-[3px] rounded-full bg-ink-muted" />
      <ul className="space-y-1" aria-label="Lineages">
        {run.lineages.map((lin) => {
          const active = focus === lin.id;
          return (
            <li key={lin.id}>
              <button
                type="button"
                aria-pressed={active}
                title={`${lin.n_cells.toLocaleString("en-US")} cells`}
                onClick={() => setFocus(active ? null : lin.id)}
                className={`press flex h-8 w-full items-center gap-2 rounded-md pr-2 text-left text-sm text-ink ${
                  active ? "bg-panel-raised" : "hover:bg-panel-raised"
                }`}
              >
                <span aria-hidden className="relative h-full w-8 shrink-0">
                  <span className="absolute top-1/2 left-[11px] h-[3px] w-3.5 -translate-y-1/2 bg-ink-muted" />
                  <span
                    className={`absolute top-1/2 left-[22px] -translate-y-1/2 rounded-full ${
                      active ? "size-3 bg-accent ring-2 ring-ink" : "size-2 bg-ink-muted"
                    }`}
                  />
                </span>
                <span className="truncate">{lin.name}</span>
                {active && <span className="ml-auto text-xs text-ink-muted">focused</span>}
              </button>
            </li>
          );
        })}
      </ul>
      </div>
    </div>
  );
}
