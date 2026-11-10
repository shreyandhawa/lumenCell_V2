import { Check, TriangleAlert } from "lucide-react";
import { api } from "@/api/client";
import { useJobEvents } from "@/api/events";
import { useDataset, useStartPrepare } from "@/api/queries";
import { Button, cx, Hint, Section } from "@/components/ui";

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
  const done = prep.stages.filter((s) => s.state === "done" || s.state === "skipped").length;

  return (
    <Section title="Preparing the dataset">
      <div className="mb-4 h-1 overflow-hidden rounded-full bg-control" aria-hidden>
        <div
          className="h-1 origin-left rounded-full bg-live transition-transform duration-500 ease-out"
          style={{ transform: `scaleX(${done / Math.max(prep.stages.length, 1)})` }}
        />
      </div>
      <ol className="space-y-0.5" aria-label="Preprocessing stages">
        {prep.stages.map((s) => {
          const running = s.state === "running" || (live?.stage === s.name && live.state === "running");
          const finished = s.state === "done" || s.state === "skipped";
          return (
            <li key={s.name} className="flex h-7 items-center gap-2.5 text-sm">
              <span className="grid size-4 shrink-0 place-items-center" aria-hidden>
                {finished ? (
                  <Check size={14} className="text-live" />
                ) : running ? (
                  <span className="size-2 animate-pulse rounded-full bg-live" />
                ) : s.state === "failed" ? (
                  <TriangleAlert size={14} className="text-danger" />
                ) : (
                  <span className="size-1.5 rounded-full bg-line-strong" />
                )}
              </span>
              <span className={cx(finished || running ? "text-ink" : "text-ink-3")}>{LABELS[s.name] ?? s.name}</span>
              {s.state === "skipped" && <span className="text-xs text-ink-3">the file has one</span>}
              {s.seconds != null && s.state === "done" && (
                <span className="num ml-auto text-xs text-ink-3">{s.seconds.toFixed(1)} s</span>
              )}
            </li>
          );
        })}
      </ol>
      <div className="mt-4 flex items-center gap-3">
        {prep.state === "missing" || prep.state === "failed" ? (
          <Button tone="primary" onClick={() => start.mutate()} disabled={start.isPending}>
            {prep.state === "failed" ? "Retry preparation" : "Prepare the dataset"}
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
