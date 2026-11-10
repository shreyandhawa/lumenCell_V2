import { Check } from "lucide-react";
import { AnimatePresence, motion } from "motion/react";
import { useState } from "react";
import { useDataset, useDoneRunId, useRoot, useRun } from "@/api/queries";
import { Compare } from "@/components/trajectory/Compare";
import { LineageTree } from "@/components/trajectory/LineageTree";
import { PrepProgress } from "@/components/trajectory/PrepProgress";
import { RootPicker, RootSummary } from "@/components/trajectory/RootPicker";
import { MethodRunner, RunList } from "@/components/trajectory/Runs";
import { TrustPanel } from "@/components/trajectory/TrustPanel";
import { cx } from "@/components/ui";
import { EASE_OUT } from "@/lib/motion";
import { useUi } from "@/state/ui";

const STEPS = ["Prepare", "Root", "Run", "Inspect"] as const;

/** The trajectory workflow, in order: prepare -> choose a root -> run -> inspect. */
export function TrajectoryPanel() {
  const dataset = useDataset().data;
  const rootQuery = useRoot();
  const root = rootQuery.data?.root ?? null;
  const runId = useDoneRunId(useUi((s) => s.runId));
  const run = useRun(runId).data;
  const [changingRoot, setChangingRoot] = useState(false);

  const prepared = dataset?.prep.state === "ready";
  const step = !prepared ? 0 : !root || changingRoot ? 1 : !run ? 2 : 3;
  // Until the saved root is known, show the panel's shape rather than flashing the picker.
  if (prepared && rootQuery.isPending) {
    return (
      <div aria-busy="true">
        <Stepper current={1} />
        <div className="space-y-3 px-4 pt-5">
          <div className="skeleton h-4 w-24" />
          <div className="skeleton h-3 w-48" />
          <div className="skeleton h-8 w-full" />
        </div>
      </div>
    );
  }

  return (
    <div className="pb-6">
      <Stepper current={step} />
      <AnimatePresence mode="popLayout" initial={false}>
        <motion.div
          key={step === 0 ? "prepare" : step === 1 ? "root" : "work"}
          initial={{ opacity: 0, y: 8 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -4, transition: { duration: 0.12 } }}
          transition={{ duration: 0.24, ease: EASE_OUT }}
        >
          {step === 0 ? (
            <PrepProgress />
          ) : step === 1 ? (
            <RootPicker onDone={() => setChangingRoot(false)} onCancel={root ? () => setChangingRoot(false) : undefined} />
          ) : (
            <>
              {root && <RootSummary root={root} onChange={() => setChangingRoot(true)} />}
              <MethodRunner />
              <RunList />
              {run && (
                <>
                  <LineageTree run={run} />
                  <TrustPanel run={run} />
                </>
              )}
              <Compare />
            </>
          )}
        </motion.div>
      </AnimatePresence>
    </div>
  );
}

function Stepper({ current }: { current: number }) {
  return (
    <ol className="flex items-center gap-1 border-b border-line px-4 py-3" aria-label="Workflow">
      {STEPS.map((label, i) => {
        const done = i < current;
        const now = i === current;
        return (
          <li key={label} className="flex min-w-0 items-center gap-1" aria-current={now ? "step" : undefined}>
            {i > 0 && <span className={cx("h-px w-3 shrink-0", done || now ? "bg-ink-3" : "bg-line-strong")} aria-hidden />}
            <span
              className={cx(
                "flex items-center gap-1.5 rounded-full py-0.5 pr-2 pl-1 text-xs font-medium whitespace-nowrap transition-colors duration-200",
                now ? "bg-control text-ink" : done ? "text-ink-2" : "text-ink-3",
              )}
            >
              <span
                className={cx(
                  "grid size-4 shrink-0 place-items-center rounded-full text-[10px] transition-colors duration-200",
                  done ? "bg-live text-live-ink" : now ? "bg-focus text-white" : "border border-line-strong",
                )}
                aria-hidden
              >
                {done ? <Check size={10} strokeWidth={3} /> : now ? i + 1 : null}
              </span>
              {label}
            </span>
          </li>
        );
      })}
    </ol>
  );
}
