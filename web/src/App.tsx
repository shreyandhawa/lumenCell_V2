import { LoaderCircle, PanelRightOpen, TriangleAlert } from "lucide-react";
import { AnimatePresence, MotionConfig } from "motion/react";
import { useEffect } from "react";
import { useDataset, useDoneRunId, useRun, useRuns } from "@/api/queries";
import { GenePanel } from "@/components/GenePanel";
import { Scatter } from "@/components/Scatter";
import { Sidebar } from "@/components/Sidebar";
import { cursorLeft, PseudotimeBar } from "@/components/trajectory/PseudotimeBar";
import { TrustPanel } from "@/components/trajectory/TrustPanel";
import { useSharedView } from "@/components/useSharedView";
import { useUi } from "@/state/ui";

const LEIDEN = "lumencell_leiden";

/** Pick sensible defaults once the data they depend on arrives; never overrides a user choice. */
function useDefaults() {
  const dataset = useDataset().data;
  const runs = useRuns().data?.runs;
  const { embedding, colorBy, runId, setEmbedding, setColorBy, setRun } = useUi();

  useEffect(() => {
    if (!dataset) return;
    const names = dataset.embeddings.filter((e) => e.dims === 2).map((e) => e.name); // the 3-D view has its own
    if ((embedding === null || !names.includes(embedding)) && dataset.default_embedding) {
      setEmbedding(dataset.default_embedding);
    }
  }, [dataset, embedding, setEmbedding]);

  useEffect(() => {
    if (!dataset || colorBy !== null || dataset.prep.state !== "ready") return;
    const label = dataset.obs.find((c) => c.name === dataset.label_key && c.kind === "categorical");
    setColorBy({ kind: "obs", column: label ? label.name : LEIDEN });
  }, [dataset, colorBy, setColorBy]);

  useEffect(() => {
    if (runId !== null || !runs) return;
    const latest = runs.find((r) => r.state === "done");
    if (latest) setRun(latest.run_id);
  }, [runs, runId, setRun]);
}

export default function App() {
  const dataset = useDataset();
  const run = useRun(useDoneRunId(useUi((s) => s.runId))).data;
  const cursor = useUi((s) => s.cursor);
  const sheetOpen = useUi((s) => s.sheetOpen);
  const setSheetOpen = useUi((s) => s.setSheetOpen);
  useDefaults();
  useSharedView(); // after the defaults, so a link's choices win

  if (dataset.isLoading) {
    return (
      <div className="grid h-full place-items-center text-[13px] text-ink-muted">
        <span className="flex items-center gap-2">
          <LoaderCircle size={16} className="animate-spin text-ink" /> Loading dataset…
        </span>
      </div>
    );
  }
  if (dataset.error) {
    return (
      <div className="grid h-full place-items-center p-8">
        <div role="alert" className="max-w-md rounded-xl bg-panel p-4 text-[13px]">
          <h1 className="mb-1 flex items-center gap-2 font-semibold text-ink">
            <TriangleAlert size={16} className="text-ink" /> Could not reach the lumenCell server
          </h1>
          <p className="text-ink-muted">{dataset.error.message}</p>
          <p className="mt-2 text-xs text-ink-muted">
            Start it with <code className="font-mono text-ink">lumencell launch your-data.h5ad</code>.
          </p>
        </div>
      </div>
    );
  }
  const trends = run !== undefined && cursor !== null;
  return (
    // Motion's own components drop transforms and keep fades under prefers-reduced-motion.
    <MotionConfig reducedMotion="user">
      <div className="relative h-full overflow-hidden bg-stage">
        <Scatter />
        <Sidebar />

        <div className="absolute right-[384px] bottom-3 left-[304px] z-10 flex flex-col gap-2">
          <AnimatePresence>{trends && <GenePanel key="trends" run={run} />}</AnimatePresence>
          <PseudotimeBar ready={run !== undefined} />
          {run && cursor !== null && (
            <div aria-hidden className="pointer-events-none absolute top-[30px] bottom-2 w-px bg-ink" style={{ left: cursorLeft(cursor) }}>
              <span className="absolute -top-1 -left-[3.5px] size-2 rounded-full bg-ink" />
            </div>
          )}
        </div>

        <AnimatePresence>
          {run && sheetOpen && <TrustPanel key="sheet" run={run} />}
        </AnimatePresence>
        {run && !sheetOpen && (
          <button
            type="button"
            aria-label="Show evidence"
            onClick={() => setSheetOpen(true)}
            className="press absolute top-3 right-3 z-10 flex h-[34px] items-center gap-2 rounded-lg bg-panel px-2.5 text-[13px] text-ink"
          >
            <PanelRightOpen size={16} /> Evidence
          </button>
        )}
      </div>
    </MotionConfig>
  );
}
