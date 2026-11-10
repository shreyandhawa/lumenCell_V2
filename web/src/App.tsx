import { useEffect } from "react";
import { Loader2 } from "lucide-react";
import { useDataset, useRuns } from "@/api/queries";
import { GenePanel } from "@/components/GenePanel";
import { Header } from "@/components/Header";
import { Scatter } from "@/components/Scatter";
import { Sidebar } from "@/components/Sidebar";
import { TrajectoryPanel } from "@/components/trajectory/TrajectoryPanel";
import { useUi } from "@/state/ui";

const LEIDEN = "lumencell_leiden";

/** Pick sensible defaults once the data they depend on arrives; never overrides a user choice. */
function useDefaults() {
  const dataset = useDataset().data;
  const runs = useRuns().data?.runs;
  const { embedding, colorBy, runId, setEmbedding, setColorBy, setRun } = useUi();

  useEffect(() => {
    if (!dataset) return;
    const names = dataset.embeddings.map((e) => e.name);
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
  useDefaults();

  if (dataset.isLoading) {
    return (
      <div className="grid h-full place-items-center text-sm text-ink-muted">
        <span className="flex items-center gap-2">
          <Loader2 size={16} className="animate-spin text-accent" /> Loading dataset…
        </span>
      </div>
    );
  }
  if (dataset.error) {
    return (
      <div className="grid h-full place-items-center p-8">
        <div role="alert" className="max-w-md rounded-lg border border-danger/40 bg-danger/10 p-4 text-sm">
          <h1 className="mb-1 font-semibold text-danger">Could not reach the lumenCell server</h1>
          <p className="text-ink-muted">{dataset.error.message}</p>
          <p className="mt-2 text-xs text-ink-faint">
            Start it with <code className="font-mono">lumencell launch your-data.h5ad</code>.
          </p>
        </div>
      </div>
    );
  }
  return (
    <div className="flex h-full flex-col">
      <Header />
      <div className="flex min-h-0 flex-1">
        <Sidebar />
        <main className="flex min-w-0 flex-1 flex-col">
          <div className="min-h-0 flex-1">
            <Scatter />
          </div>
          <GenePanel />
        </main>
        <TrajectoryPanel />
      </div>
    </div>
  );
}
