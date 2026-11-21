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
                className={`press flex h-8 w-full items-center gap-3 rounded-md pr-2 pl-2 text-left text-[13px] text-ink ${
                  active ? "bg-panel-raised" : "hover:bg-panel-raised"
                }`}
              >
                <span aria-hidden className={`h-4 w-0.5 rounded-full ${active ? "bg-accent" : "bg-line"}`} />
                <span className="truncate">{lin.name}</span>
                {active && <span className="ml-auto text-xs text-ink-muted">focused</span>}
              </button>
            </li>
          );
        })}
      </ul>
    </div>
  );
}
