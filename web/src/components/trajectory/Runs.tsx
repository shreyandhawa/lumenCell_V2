import { LoaderCircle, TriangleAlert } from "lucide-react";
import { useState } from "react";
import { useJobEvents } from "@/api/events";
import { useCreateRun, useMethods, useRuns } from "@/api/queries";
import { Button, Mark, Select } from "@/components/ui";
import { useUi } from "@/state/ui";

export function MethodRunner() {
  const methods = useMethods().data?.methods ?? [];
  const create = useCreateRun();
  const runId = useUi((s) => s.runId);
  const setRun = useUi((s) => s.setRun);
  const setColorBy = useUi((s) => s.setColorBy);
  const [method, setMethod] = useState<string>("");
  const chosen = methods.find((m) => m.name === (method || methods.find((x) => x.available)?.name));
  const [jobId, setJobId] = useState<string | null>(null);
  const live = useJobEvents(jobId);
  const failed = useRuns().data?.runs.find((r) => r.run_id === runId && r.state === "failed");

  const run = () => {
    if (!chosen) return;
    create.mutate(
      { method: chosen.name },
      {
        onSuccess: (created) => {
          setJobId(created.job_id ?? null);
          setRun(created.run_id);
          setColorBy({ kind: "pseudotime" });
        },
      },
    );
  };

  return (
    <div className="space-y-3">
      <div>
        <label htmlFor="method" className="mb-1.5 block text-xs text-ink-muted">
          Method
        </label>
        <Select id="method" aria-label="Method" value={chosen?.name ?? ""} onChange={(e) => setMethod(e.target.value)}>
          {methods.map((m) => (
            <option key={m.name} value={m.name} disabled={!m.available}>
              {m.label}
              {m.available ? "" : " (not installed)"}
            </option>
          ))}
        </Select>
      </div>
      {chosen?.experimental && <Mark title={chosen.experimental_note ?? undefined}>Experimental</Mark>}
      {chosen && !chosen.available && chosen.install_hint && <p className="text-xs text-ink-muted">{chosen.install_hint}</p>}
      <div>
        <Button tone="primary" onClick={run} disabled={!chosen?.available || create.isPending}>
          Run {chosen?.label ?? ""}
        </Button>
      </div>
      {live && live.state === "running" && (
        <div aria-label="Run progress">
          <div className="mb-1 flex items-center gap-1.5 text-xs text-ink-muted">
            <LoaderCircle size={14} aria-hidden className="animate-spin text-ink" /> {live.stage ?? "starting"}
          </div>
          <div className="h-1 rounded-full bg-panel-raised">
            <div className="h-1 rounded-full bg-ink-muted" style={{ width: `${live.fraction * 100}%` }} />
          </div>
        </div>
      )}
      {(create.error || failed) && (
        <p role="alert" className="flex gap-1.5 text-xs text-ink">
          <TriangleAlert size={14} className="mt-px shrink-0 text-amber" />
          {create.error?.message ?? failed?.error ?? "The run failed."}
        </p>
      )}
    </div>
  );
}
