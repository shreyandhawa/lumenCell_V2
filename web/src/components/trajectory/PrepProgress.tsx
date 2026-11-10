import { Check, Circle, Loader2, TriangleAlert } from "lucide-react";
import { api } from "@/api/client";
import { useJobEvents } from "@/api/events";
import { useDataset, useStartPrepare } from "@/api/queries";
import { Button, Hint, Section } from "@/components/ui";

const LABELS: Record<string, string> = {
  normalize: "Normalise counts",
  hvg: "Select variable genes",
  pca: "PCA",
  neighbors: "Neighbour graph",
  leiden: "Leiden clusters",
  diffmap: "Diffusion map",
  paga: "PAGA cluster graph",
  umap: "UMAP (display only)",
};

export function PrepProgress() {
  const dataset = useDataset().data;
  const start = useStartPrepare();
  const prep = dataset?.prep;
  const live = useJobEvents(prep?.state === "running" ? prep.job_id : null);
  if (!prep) return null;

  return (
    <Section title="Preparing the dataset">
      <ol className="space-y-1" aria-label="Preprocessing stages">
        {prep.stages.map((s) => {
          const running = s.state === "running" || (live?.stage === s.name && live.state === "running");
          return (
            <li key={s.name} className="flex items-center gap-2 text-xs">
              {s.state === "done" || s.state === "skipped" ? (
                <Check size={13} className="text-success" />
              ) : running ? (
                <Loader2 size={13} className="animate-spin text-accent" />
              ) : s.state === "failed" ? (
                <TriangleAlert size={13} className="text-danger" />
              ) : (
                <Circle size={13} className="text-ink-faint" />
              )}
              <span className={s.state === "pending" ? "text-ink-faint" : "text-ink"}>
                {LABELS[s.name] ?? s.name}
              </span>
              {s.state === "skipped" && <span className="text-2xs text-ink-faint">(file has one)</span>}
              {s.seconds != null && s.state === "done" && (
                <span className="ml-auto font-mono text-2xs text-ink-faint">{s.seconds.toFixed(1)}s</span>
              )}
            </li>
          );
        })}
      </ol>
      <div className="mt-3 flex items-center gap-2">
        {prep.state === "missing" || prep.state === "failed" ? (
          <Button tone="primary" onClick={() => start.mutate()} disabled={start.isPending}>
            {prep.state === "failed" ? "Retry" : "Prepare now"}
          </Button>
        ) : (
          prep.job_id && (
            <Button tone="danger" onClick={() => void api.cancelJob(prep.job_id!)}>
              Cancel
            </Button>
          )
        )}
        <Hint>Runs once; later launches reuse the cache.</Hint>
      </div>
    </Section>
  );
}
