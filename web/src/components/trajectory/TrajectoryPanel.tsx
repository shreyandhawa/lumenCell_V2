import { useState } from "react";
import { useDataset, useRoot, useRun, useDoneRunId } from "@/api/queries";
import { LineageTree } from "@/components/trajectory/LineageTree";
import { PrepProgress } from "@/components/trajectory/PrepProgress";
import { PseudotimeFilter } from "@/components/trajectory/PseudotimeFilter";
import { RootPicker, RootSummary } from "@/components/trajectory/RootPicker";
import { MethodRunner, RunList } from "@/components/trajectory/Runs";
import { TrustPanel } from "@/components/trajectory/TrustPanel";
import { useUi } from "@/state/ui";

/** The trajectory workflow, in order: prepare -> choose a root -> run -> inspect. */
export function TrajectoryPanel() {
  const dataset = useDataset().data;
  const root = useRoot().data?.root ?? null;
  const runId = useDoneRunId(useUi((s) => s.runId));
  const run = useRun(runId).data;
  const [changingRoot, setChangingRoot] = useState(false);

  return (
    <aside className="flex w-80 shrink-0 flex-col overflow-y-auto border-l border-line bg-panel" aria-label="Trajectory">
      {dataset?.prep.state !== "ready" ? (
        <PrepProgress />
      ) : !root || changingRoot ? (
        <RootPicker onDone={() => setChangingRoot(false)} />
      ) : (
        <>
          <RootSummary root={root} onChange={() => setChangingRoot(true)} />
          <MethodRunner />
          <RunList />
          {run && (
            <>
              <LineageTree run={run} />
              <PseudotimeFilter />
              <TrustPanel run={run} />
            </>
          )}
        </>
      )}
    </aside>
  );
}
