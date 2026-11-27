import { Check, CircleDashed, LoaderCircle, TriangleAlert } from "lucide-react";
import type * as T from "@/api/generated";
import { useJobEvents } from "@/api/events";
import { useStartPrepare } from "@/api/queries";
import { Button } from "@/components/ui";

const LABELS: Record<string, string> = {
  normalize: "Normalise",
  hvg: "Highly variable genes",
  pca: "PCA",
  neighbors: "Neighbours",
  leiden: "Leiden clusters",
  diffmap: "Diffusion map",
  paga: "PAGA",
  umap: "UMAP",
  cell_cycle: "Cell cycle",
};

/** "5 of 8": the stage running now (or the count done), for the step's aside. */
export function prepCount(prep: T.PrepStatus): string {
  const done = prep.stages.filter((s) => s.state === "done" || s.state === "skipped").length;
  const n = prep.stages.length;
  return `${prep.state === "running" ? Math.min(done + 1, n) : done} of ${n}`;
}

export function PrepProgress({ prep, hasFileEmbedding }: { prep: T.PrepStatus; hasFileEmbedding: boolean }) {
  const start = useStartPrepare();
  const live = useJobEvents(prep.state === "running" ? prep.job_id : null);

  return (
    <>
      <ol className="space-y-2.5" aria-label="Preprocessing stages">
        {prep.stages.map((s) => {
          const running = s.state === "running" || (live?.stage === s.name && live.state === "running");
          const done = s.state === "done" || s.state === "skipped";
          return (
            <li key={s.name} className="flex items-center gap-2.5 text-sm">
              {done ? (
                <Check size={15} aria-hidden className="text-ink-muted" />
              ) : running ? (
                <LoaderCircle size={15} aria-hidden className="animate-spin text-ink" />
              ) : s.state === "failed" ? (
                <TriangleAlert size={15} aria-hidden className="text-ink" />
              ) : (
                <CircleDashed size={15} aria-hidden className="text-ink-muted" />
              )}
              <span className={running ? "text-ink" : "text-ink-muted"}>{LABELS[s.name] ?? s.name}</span>
              <span className="ml-auto font-mono text-xs text-ink-muted">
                {s.state === "skipped" ? "in file" : running ? "running" : s.state === "pending" ? "" : s.state}
              </span>
            </li>
          );
        })}
      </ol>
      {prep.state === "missing" || prep.state === "failed" ? (
        <Button tone="primary" className="mt-3" onClick={() => start.mutate()} disabled={start.isPending}>
          {prep.state === "failed" ? "Retry" : "Prepare now"}
        </Button>
      ) : (
        <p className="mt-3 text-xs leading-relaxed text-ink-muted">
          Runs in the background.{" "}
          {hasFileEmbedding
            ? "The stage shows the embedding stored in the file meanwhile."
            : "The stage fills in when the UMAP is ready."}
        </p>
      )}
    </>
  );
}
