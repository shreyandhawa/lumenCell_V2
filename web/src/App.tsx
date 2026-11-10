import { motion } from "motion/react";
import { useEffect, useLayoutEffect, useRef, useState, type ReactNode } from "react";
import { useDataset, useDoneRunId, useRun, useRuns } from "@/api/queries";
import { GenePanel } from "@/components/GenePanel";
import { CELLS_W, GAP, INSPECTOR_W, sheetHeight, stageInsets, TIMELINE_H } from "@/components/layout";
import { Scatter, isTyping } from "@/components/Scatter";
import { Sidebar } from "@/components/Sidebar";
import { Timeline } from "@/components/stage/Timeline";
import { Toolbar } from "@/components/Toolbar";
import { TrajectoryPanel } from "@/components/trajectory/TrajectoryPanel";
import { Button, cx, TooltipProvider } from "@/components/ui";
import { LEIDEN } from "@/lib/colorBy";
import { PANEL_SPRING } from "@/lib/motion";
import { useUi } from "@/state/ui";

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

/** Keyboard: Space plays, [ and ] toggle the panels, 1–9 pick an embedding. */
function useShortcuts() {
  const dataset = useDataset().data;
  const hasRun = useDoneRunId(useUi((s) => s.runId)) !== null;
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.metaKey || e.ctrlKey || e.altKey || isTyping(e.target)) return;
      const ui = useUi.getState();
      if (e.key === "[") ui.togglePanel("cells", false);
      else if (e.key === "]") ui.togglePanel("inspector", false);
      else if (e.key === " " && hasRun && !(e.target as HTMLElement)?.closest?.("button,[role=slider]")) {
        e.preventDefault();
        ui.setPlaying(!ui.playing);
      } else if (/^[1-9]$/.test(e.key) && dataset) {
        const target = dataset.embeddings[Number(e.key) - 1];
        if (target) ui.setEmbedding(target.name);
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [dataset, hasRun]);
}

function useTheme() {
  const theme = useUi((s) => s.theme);
  useLayoutEffect(() => {
    document.documentElement.dataset.theme = theme;
  }, [theme]);
}

function useStageSize(): [React.RefObject<HTMLDivElement | null>, { width: number; height: number }] {
  const ref = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState({ width: 1440, height: 800 });
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ro = new ResizeObserver(([e]) => e && setSize({ width: e.contentRect.width, height: e.contentRect.height }));
    ro.observe(el);
    return () => ro.disconnect();
  }, []);
  return [ref, size];
}

/** Below this stage width the panels float over the field as sheets and start closed. */
const NARROW = 900;
/** Below this stage width the Cells panel starts closed so the field keeps room. */
const COMPACT = 1200;

export default function App() {
  const dataset = useDataset();
  useDefaults();
  useShortcuts();
  useTheme();

  if (dataset.isLoading) return <Opening />;
  if (dataset.error) {
    return (
      <div className="grid h-full place-items-center bg-field p-8">
        <div role="alert" className="max-w-sm">
          <h1 className="text-lg font-semibold tracking-[-0.02em] text-ink">Can't reach the lumenCell server</h1>
          <p className="mt-2 text-sm text-ink-2">{dataset.error.message}</p>
          <p className="mt-3 text-sm text-ink-2">
            Start it from a terminal with <code className="num text-ink">lumencell launch your-data.h5ad</code>, then
            try again.
          </p>
          <Button tone="primary" className="mt-5" onClick={() => void dataset.refetch()}>
            Try again
          </Button>
        </div>
      </div>
    );
  }
  return (
    <TooltipProvider>
      <div className="flex h-full flex-col">
        <Toolbar />
        <Stage />
      </div>
    </TooltipProvider>
  );
}

function Stage() {
  const [ref, { width: stageWidth, height: stageHeight }] = useStageSize();
  const narrow = stageWidth < NARROW;
  const panels = useUi((s) => s.panels);
  const animate = useUi((s) => s.panelsAnimate);
  const detent = useUi((s) => s.sheet);
  const drag = useUi((s) => s.sheetDrag);
  const runId = useDoneRunId(useUi((s) => s.runId));
  const run = useRun(runId).data;
  const hasTimeline = Boolean(run);
  // On a narrow screen the panels cover the field rather than squeezing it.
  const insets = stageInsets(narrow ? { cells: false, inspector: false } : panels, detent, stageHeight, hasTimeline);
  const closedForNarrow = useRef(false);
  useEffect(() => {
    if (stageWidth >= COMPACT || closedForNarrow.current) return;
    closedForNarrow.current = true;
    useUi.setState({ panels: { cells: false, inspector: !narrow }, panelsAnimate: false });
  }, [stageWidth, narrow]);
  const transition = animate ? PANEL_SPRING : { duration: 0 };
  const sheetH = Math.max(44, Math.min(stageHeight * 0.8, sheetHeight(detent, stageHeight) - (drag ?? 0)));
  const sheetTransition = drag !== null ? { duration: 0 } : transition;
  const middle = { left: insets.left, right: insets.right };

  return (
    <main ref={ref} className="relative min-h-0 flex-1 overflow-hidden">
      <Scatter insets={insets} />

      <FloatingPanel
        side="left"
        open={panels.cells}
        width={Math.min(CELLS_W, stageWidth - 2 * GAP)}
        label="Cells panel"
        transition={transition}
      >
        <Sidebar />
      </FloatingPanel>
      <FloatingPanel
        side="right"
        open={panels.inspector}
        width={Math.min(INSPECTOR_W, stageWidth - 2 * GAP)}
        label="Trajectory"
        transition={transition}
      >
        <TrajectoryPanel />
      </FloatingPanel>

      {run && (
        <motion.div
          className="pointer-events-none absolute flex justify-center"
          initial={false}
          animate={{ ...middle, bottom: GAP + sheetH + GAP }}
          transition={sheetTransition}
          style={{ height: TIMELINE_H }}
        >
          <div className="pointer-events-auto w-full max-w-[44rem]">
            <Timeline run={run} />
          </div>
        </motion.div>
      )}

      <motion.div
        className="absolute"
        initial={false}
        animate={{ ...middle, bottom: GAP, height: sheetH }}
        transition={sheetTransition}
      >
        <GenePanel />
      </motion.div>
    </main>
  );
}

function FloatingPanel({
  side,
  open,
  width,
  label,
  transition,
  children,
}: {
  side: "left" | "right";
  open: boolean;
  width: number;
  label: string;
  transition: object;
  children: ReactNode;
}) {
  const away = (width + GAP * 2) * (side === "left" ? -1 : 1);
  return (
    <motion.aside
      aria-label={label}
      inert={!open}
      initial={false}
      animate={{ x: open ? 0 : away, opacity: open ? 1 : 0 }}
      transition={transition}
      className={cx(
        "absolute top-3 bottom-3 z-30 overflow-x-hidden overflow-y-auto overscroll-contain rounded-2xl border border-line bg-panel",
        side === "left" ? "left-3" : "right-3",
      )}
      style={{ width }}
    >
      {children}
    </motion.aside>
  );
}

/** The shell's shape while the dataset loads: no spinner in the middle of nothing. */
function Opening() {
  return (
    <div className="flex h-full flex-col bg-field" aria-busy="true" aria-label="Opening the dataset">
      <div className="flex h-12 items-center gap-3 border-b border-line bg-panel px-4">
        <h1 className="text-md font-semibold tracking-[-0.02em] text-ink">lumenCell</h1>
        <div className="skeleton h-4 w-40" />
        <div className="mx-auto skeleton h-8 w-80" />
      </div>
      <div className="relative flex-1">
        <div className="absolute top-3 bottom-3 left-3 w-[264px] rounded-2xl border border-line bg-panel p-4">
          <div className="skeleton h-4 w-20" />
          <div className="mt-4 skeleton h-3 w-36" />
          <div className="mt-2 skeleton h-3 w-28" />
        </div>
        <div className="absolute top-3 right-3 bottom-3 w-[336px] rounded-2xl border border-line bg-panel p-4">
          <div className="skeleton h-4 w-32" />
          <div className="mt-4 skeleton h-16 w-full" />
        </div>
      </div>
    </div>
  );
}
