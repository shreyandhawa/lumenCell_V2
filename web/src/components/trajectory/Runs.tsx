import { Play } from "lucide-react";
import { AnimatePresence, motion } from "motion/react";
import { useState } from "react";
import { useJobEvents } from "@/api/events";
import { useCreateRun, useDoneRunId, useMethods, useRun, useRuns } from "@/api/queries";
import { Button, cx, Hint, Note, Section, Select, Tag } from "@/components/ui";
import { EASE_OUT } from "@/lib/motion";
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
  const running = live?.state === "running";
  const activeRun = useRun(useDoneRunId(useUi((s) => s.runId))).data;
  // With a result on the field, running again is a secondary action: one row, quiet button.
  const compact = Boolean(activeRun) && !running;

  return (
    <Section title="Method">
      <div className={cx(compact && "flex items-center gap-2")}>
        <Select
          label="Method"
          className={compact ? "min-w-0 flex-1" : "w-full"}
          value={chosen?.name ?? null}
          onChange={setMethod}
          groups={[
            {
              items: methods.map((m) => ({
                value: m.name,
                label: m.label,
                disabled: !m.available,
                note: !m.available ? "not installed" : m.experimental ? "experimental" : undefined,
                noteTone: m.available && m.experimental ? ("warn" as const) : undefined,
              })),
            },
          ]}
        />
        {compact && (
          <Button onClick={run} disabled={!chosen?.available || create.isPending}>
            <Play size={13} fill="currentColor" />
            {chosen?.name === activeRun?.provenance.method ? "Run again" : `Run ${chosen?.label ?? ""}`}
          </Button>
        )}
      </div>
      {!compact && chosen && (
        <div className="mt-2 space-y-2">
          <Hint>{chosen.description}</Hint>
          {chosen.experimental && (
            <Note tone="warning">Experimental: misses some of lumenCell's accuracy gates. The trust panel says which.</Note>
          )}
          {!chosen.available && chosen.install_hint && <Hint>{chosen.install_hint}</Hint>}
        </div>
      )}
      {!compact && (
        <Button
          tone="primary"
          className="mt-3 w-full"
          onClick={run}
          disabled={!chosen?.available || create.isPending || running}
        >
          <Play size={13} fill="currentColor" /> Run {chosen?.label ?? ""}
        </Button>
      )}
      <AnimatePresence initial={false}>
        {running && live && (
          <motion.div
            key="progress"
            aria-label="Run progress"
            initial={{ opacity: 0, height: 0 }}
            animate={{ opacity: 1, height: "auto" }}
            exit={{ opacity: 0, height: 0 }}
            transition={{ duration: 0.2, ease: EASE_OUT }}
            className="overflow-hidden"
          >
            <div className="pt-3">
              <div className="mb-1.5 flex items-center justify-between text-xs text-ink-2" aria-live="polite">
                <span>{live.stage ?? "starting"}</span>
                <span className="num">{Math.round(live.fraction * 100)}%</span>
              </div>
              <div className="h-1 overflow-hidden rounded-full bg-control">
                <div
                  className="h-1 origin-left rounded-full bg-live transition-transform duration-300 ease-out"
                  style={{ transform: `scaleX(${live.fraction})` }}
                />
              </div>
            </div>
          </motion.div>
        )}
      </AnimatePresence>
      {create.error && (
        <div className="mt-2">
          <Note tone="error" role="alert">
            {create.error.message}
          </Note>
        </div>
      )}
    </Section>
  );
}

export function RunList() {
  const runs = useRuns().data?.runs ?? [];
  const runId = useUi((s) => s.runId);
  const setRun = useUi((s) => s.setRun);
  const seen = useUi((s) => s.seenRuns);
  if (runs.length === 0) return null;
  return (
    <Section title="Runs">
      <ul className="-mx-2 space-y-0.5">
        <AnimatePresence initial={false}>
          {runs.map((r) => {
            const active = r.run_id === runId;
            const fresh = r.state === "done" && !active && !seen.includes(r.run_id);
            return (
              <motion.li
                key={r.run_id}
                layout="position"
                initial={{ opacity: 0, y: -6 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.22, ease: EASE_OUT }}
              >
                <button
                  type="button"
                  onClick={() => setRun(r.run_id)}
                  aria-current={active}
                  className={cx(
                    "pressable flex h-9 w-full items-center gap-2 rounded-lg px-2 text-left text-sm",
                    active ? "bg-focus/12 text-ink" : "text-ink-2 [@media(hover:hover)]:hover:bg-control",
                  )}
                >
                  <span className={cx("size-1.5 shrink-0 rounded-full", fresh ? "bg-focus" : "bg-transparent")} aria-hidden />
                  <span className={cx("font-medium", active ? "text-ink" : "text-ink")}>{r.method}</span>
                  <span className="num text-xs text-ink-3">{r.run_id.slice(0, 6)}</span>
                  <span className="ml-auto">
                    {r.state === "done" ? (
                      <span className="num text-xs text-ink-3">
                        {r.n_lineages} lineage{r.n_lineages === 1 ? "" : "s"}
                      </span>
                    ) : r.state === "failed" ? (
                      <Tag tone="warn">failed</Tag>
                    ) : (
                      <Tag tone="live">{r.state}</Tag>
                    )}
                  </span>
                  {fresh && <span className="sr-only">new, not opened yet</span>}
                </button>
                {r.state === "failed" && r.error && <Hint className="px-2 pb-1">{r.error}</Hint>}
              </motion.li>
            );
          })}
        </AnimatePresence>
      </ul>
    </Section>
  );
}
