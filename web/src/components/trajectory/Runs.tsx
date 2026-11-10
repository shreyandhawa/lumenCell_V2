import { Loader2, Play, TriangleAlert } from "lucide-react";
import { useState } from "react";
import { useJobEvents } from "@/api/events";
import { useCreateRun, useMethods, useRuns } from "@/api/queries";
import { Badge, Button, Hint, Section, Select } from "@/components/ui";
import { useUi } from "@/state/ui";

export function MethodRunner() {
  const methods = useMethods().data?.methods ?? [];
  const create = useCreateRun();
  const setRun = useUi((s) => s.setRun);
  const setColorBy = useUi((s) => s.setColorBy);
  const [method, setMethod] = useState<string>("");
  const chosen = methods.find((m) => m.name === (method || methods.find((x) => x.available)?.name));
  const [jobId, setJobId] = useState<string | null>(null);
  const live = useJobEvents(jobId);

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
    <Section title="Method">
      <Select aria-label="Method" value={chosen?.name ?? ""} onChange={(e) => setMethod(e.target.value)}>
        {methods.map((m) => (
          <option key={m.name} value={m.name} disabled={!m.available}>
            {m.label}
            {m.available ? "" : " (not installed)"}
          </option>
        ))}
      </Select>
      {chosen && (
        <div className="mt-1.5 space-y-1">
          <Hint>{chosen.description}</Hint>
          {chosen.experimental && (
            <p className="flex gap-1.5 text-2xs text-warning">
              <TriangleAlert size={12} className="mt-px shrink-0" />
              Experimental: misses some of lumenCell's accuracy gates. See the trust panel after running.
            </p>
          )}
          {!chosen.available && chosen.install_hint && <Hint>{chosen.install_hint}</Hint>}
        </div>
      )}
      <Button tone="primary" className="mt-2" onClick={run} disabled={!chosen?.available || create.isPending}>
        <Play size={12} /> Run {chosen?.label ?? ""}
      </Button>
      {live && live.state === "running" && (
        <div className="mt-2" aria-label="Run progress">
          <div className="mb-1 flex items-center gap-1.5 text-2xs text-ink-muted">
            <Loader2 size={12} className="animate-spin text-accent" /> {live.stage ?? "starting"}
          </div>
          <div className="h-1 rounded-full bg-line">
            <div className="h-1 rounded-full bg-accent transition-[width]" style={{ width: `${live.fraction * 100}%` }} />
          </div>
        </div>
      )}
      {create.error && (
        <p role="alert" className="mt-2 text-2xs text-danger">
          {create.error.message}
        </p>
      )}
    </Section>
  );
}

export function RunList() {
  const runs = useRuns().data?.runs ?? [];
  const runId = useUi((s) => s.runId);
  const setRun = useUi((s) => s.setRun);
  if (runs.length === 0) return null;
  return (
    <Section title="Runs">
      <ul className="space-y-1">
        {runs.map((r) => (
          <li key={r.run_id}>
            <button
              type="button"
              onClick={() => setRun(r.run_id)}
              aria-current={r.run_id === runId}
              className={`flex w-full items-center gap-2 rounded-md border px-2 py-1.5 text-left text-xs transition-colors ${
                r.run_id === runId ? "border-accent/60 bg-accent/10" : "border-line hover:border-ink-faint"
              }`}
            >
              <span className="font-medium text-ink">{r.method}</span>
              <span className="font-mono text-2xs text-ink-faint">{r.run_id.slice(0, 6)}</span>
              <span className="ml-auto">
                {r.state === "done" ? (
                  <Badge>{r.n_lineages} lineages</Badge>
                ) : r.state === "failed" ? (
                  <Badge tone="warning">failed</Badge>
                ) : (
                  <Badge tone="accent">{r.state}</Badge>
                )}
              </span>
            </button>
            {r.state === "failed" && r.error && <Hint>{r.error}</Hint>}
          </li>
        ))}
      </ul>
    </Section>
  );
}
