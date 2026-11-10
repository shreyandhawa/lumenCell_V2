/**
 * Genes along pseudotime, in a bottom sheet. Every chart here shares the field's one
 * pseudotime axis: the playback window shows on it, and pointing at a time on any chart
 * lights the cells at that pseudotime on the field.
 */
import { Tabs } from "@base-ui/react/tabs";
import { useQueries } from "@tanstack/react-query";
import { scaleLinear } from "d3-scale";
import { line } from "d3-shape";
import { ChevronUp, X } from "lucide-react";
import { useEffect, useRef, useState, type PointerEvent as ReactPointerEvent } from "react";
import { api } from "@/api/client";
import type * as T from "@/api/generated";
import { useDoneRunId, useHeatmap, useRun } from "@/api/queries";
import { GeneSearch } from "@/components/GeneSearch";
import { cx, Hint, IconButton, Select } from "@/components/ui";
import { categorical, rgbaCss, sequential } from "@/lib/colors";
import { MAX_TREND_GENES, type SheetDetent, useUi } from "@/state/ui";

type Tab = "trends" | "heatmap";
const ORDER: SheetDetent[] = ["bar", "peek", "full"];

export function GenePanel() {
  const runId = useDoneRunId(useUi((s) => s.runId));
  const run = useRun(runId).data;
  const detent = useUi((s) => s.sheet);
  const setSheet = useUi((s) => s.setSheet);
  const [tab, setTab] = useState<Tab>("trends");
  const open = detent !== "bar";

  return (
    <section
      aria-label="Genes along pseudotime"
      className="@container flex h-full flex-col overflow-hidden rounded-2xl border border-line bg-panel"
    >
      <Tabs.Root
        value={tab}
        onValueChange={(t) => {
          setTab(t as Tab);
          if (!open) setSheet("peek");
        }}
        className="flex min-h-0 flex-1 flex-col"
      >
        <SheetHeader enabled={Boolean(run)}>
          <h2 className="shrink-0 text-sm font-semibold whitespace-nowrap text-ink @max-[36rem]:sr-only">Genes along pseudotime</h2>
          {run ? (
            <Tabs.List className="relative flex shrink-0 gap-0.5 @min-[36rem]:ml-3">
              <Tabs.Tab value="trends" className={TAB}>
                Trends
              </Tabs.Tab>
              <Tabs.Tab value="heatmap" className={TAB}>
                Top temporal genes
              </Tabs.Tab>
              <Tabs.Indicator className="absolute top-0 left-0 -z-0 h-full w-(--active-tab-width) translate-x-(--active-tab-left) rounded-md bg-control transition-[translate,width] duration-200 ease-out" />
            </Tabs.List>
          ) : (
            <span className="ml-3 text-xs text-ink-3">Available after a run</span>
          )}
          {run && (
            <IconButton
              label={open ? "Collapse gene panel" : "Expand gene panel"}
              onClick={() => setSheet(open ? "bar" : "peek")}
              className="ml-auto"
            >
              <ChevronUp size={16} className={cx("transition-transform duration-200 ease-out", open && "rotate-180")} />
            </IconButton>
          )}
        </SheetHeader>
        {run && (
          <div className={cx("min-h-0 flex-1 px-4 pb-4 transition-opacity duration-200", open ? "opacity-100" : "opacity-0")}>
            <Tabs.Panel value="trends" className="h-full rounded-lg outline-none focus-visible:outline-2 focus-visible:outline-focus">
              <Trends run={run} />
            </Tabs.Panel>
            <Tabs.Panel value="heatmap" className="h-full rounded-lg outline-none focus-visible:outline-2 focus-visible:outline-focus">
              <Heatmap run={run} />
            </Tabs.Panel>
          </div>
        )}
      </Tabs.Root>
    </section>
  );
}

const TAB =
  "pressable relative z-10 h-7 rounded-md px-2.5 text-xs font-medium whitespace-nowrap text-ink-2 outline-none data-active:text-ink [@media(hover:hover)]:hover:text-ink focus-visible:outline-2 focus-visible:outline-focus";

/** The sheet's top bar doubles as its drag handle: drag to resize, a quick flick is enough. */
function SheetHeader({ enabled, children }: { enabled: boolean; children: React.ReactNode }) {
  const setDrag = useUi((s) => s.setSheetDrag);
  const setSheet = useUi((s) => s.setSheet);
  const start = useRef<{ y: number; t: number } | null>(null);

  const onDown = (e: ReactPointerEvent<HTMLDivElement>) => {
    if (!enabled || (e.target as HTMLElement).closest("button,[role=tab]")) return;
    start.current = { y: e.clientY, t: performance.now() };
    e.currentTarget.setPointerCapture(e.pointerId);
  };
  const onMove = (e: ReactPointerEvent<HTMLDivElement>) => {
    if (start.current) setDrag(e.clientY - start.current.y);
  };
  const onUp = (e: ReactPointerEvent<HTMLDivElement>) => {
    const s = start.current;
    start.current = null;
    if (!s) return;
    const dy = e.clientY - s.y;
    const velocity = dy / Math.max(1, performance.now() - s.t);
    const current = ORDER.indexOf(useUi.getState().sheet);
    const step = Math.abs(dy) > 60 || Math.abs(velocity) > 0.11 ? (dy < 0 ? 1 : -1) : 0;
    setSheet(ORDER[Math.max(0, Math.min(ORDER.length - 1, current + step))]!);
  };

  return (
    <div
      onPointerDown={onDown}
      onPointerMove={onMove}
      onPointerUp={onUp}
      onPointerCancel={onUp}
      className={cx("relative flex h-11 shrink-0 touch-none items-center px-4 select-none", enabled && "cursor-row-resize")}
    >
      {enabled && <span className="absolute top-1.5 left-1/2 h-1 w-8 -translate-x-1/2 rounded-full bg-line-strong" aria-hidden />}
      {children}
    </div>
  );
}

/* ---- Trends -------------------------------------------------------------------------- */

function Trends({ run }: { run: T.TrajectoryRun }) {
  const genes = useUi((s) => s.trendGenes);
  const add = useUi((s) => s.addTrendGene);
  const remove = useUi((s) => s.removeTrendGene);
  const focus = useUi((s) => s.focusLineage);
  const lineages = run.lineages.filter((l) => focus === null || l.id === focus);
  const results = useQueries({
    queries: lineages.map((l) => ({
      queryKey: ["trends", run.run_id, l.id, genes.join(",")],
      queryFn: () => api.trends(run.run_id, l.id, genes),
      enabled: genes.length > 0,
      staleTime: Infinity,
    })),
  });

  return (
    <div className="flex h-full gap-5">
      <div className="w-52 shrink-0 space-y-3">
        <GeneSearch onPick={add} label="Add a gene" />
        {genes.length > 0 && (
          <ul className="flex flex-wrap gap-1" aria-label="Genes in the trends panel">
            {genes.map((g) => (
              <li key={g} translate="no" className="flex h-6 items-center gap-0.5 rounded-full bg-control pr-0.5 pl-2.5 text-xs text-ink">
                {g}
                <button
                  type="button"
                  aria-label={`Remove ${g}`}
                  onClick={() => remove(g)}
                  className="pressable grid size-5 place-items-center rounded-full text-ink-3 [@media(hover:hover)]:hover:bg-control-hover [@media(hover:hover)]:hover:text-ink"
                >
                  <X size={11} />
                </button>
              </li>
            ))}
          </ul>
        )}
        <Hint>
          Up to {MAX_TREND_GENES} genes. One line per lineage{focus !== null ? " (the focused one only)" : ""}: mean
          log-normalised expression, smoothed.
        </Hint>
      </div>
      {genes.length === 0 ? (
        <div className="grid flex-1 place-items-center text-sm text-ink-3">
          Add a gene, or pick one from Top temporal genes.
        </div>
      ) : (
        <div
          className="grid min-w-0 flex-1 auto-rows-[minmax(0,1fr)] grid-cols-[repeat(auto-fit,minmax(15rem,1fr))] gap-x-5 gap-y-2 overflow-y-auto"
          aria-label="Gene trends"
        >
          {genes.map((gene, k) => (
            <TrendChart
              key={gene}
              gene={gene}
              series={lineages.flatMap((l, i) => {
                const data = results[i]?.data;
                const trend = data?.genes[k];
                return data && trend ? [{ id: l.id, name: l.name, x: data.x, y: trend.y, r2: trend.r2 }] : [];
              })}
            />
          ))}
        </div>
      )}
    </div>
  );
}

interface Series {
  id: number;
  name: string;
  x: number[];
  y: (number | null)[];
  r2: number;
}

function useWidth(): [React.RefObject<HTMLDivElement | null>, number] {
  const ref = useRef<HTMLDivElement>(null);
  const [w, setW] = useState(240);
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ro = new ResizeObserver(([e]) => e && setW(Math.max(160, e.contentRect.width)));
    ro.observe(el);
    return () => ro.disconnect();
  }, []);
  return [ref, w];
}

function TrendChart({ gene, series }: { gene: string; series: Series[] }) {
  const [ref, w] = useWidth();
  const highlight = useUi((s) => s.highlight);
  const setHighlight = useUi((s) => s.setHighlight);
  const ptRange = useUi((s) => s.ptRange);
  const h = 150;
  const m = { top: 24, right: 8, bottom: 20, left: 32 };
  const yMax = Math.max(0.1, ...series.flatMap((s) => s.y.filter((v): v is number => v !== null)));
  const x = scaleLinear().domain([0, 1]).range([m.left, w - m.right]);
  const y = scaleLinear().domain([0, yMax]).range([h - m.bottom, m.top]).nice();
  const path = line<[number, number | null]>()
    .defined((d) => d[1] !== null)
    .x((d) => x(d[0]))
    .y((d) => y(d[1] ?? 0));
  const best = series.length ? Math.max(...series.map((s) => s.r2)) : 0;
  const at = highlight?.kind === "pseudotime" ? highlight.at : null;
  const windowed = ptRange[0] > 0 || ptRange[1] < 1;

  const onMove = (e: ReactPointerEvent<SVGSVGElement>) => {
    const r = e.currentTarget.getBoundingClientRect();
    const t = x.invert(e.clientX - r.left);
    setHighlight(t >= 0 && t <= 1 ? { kind: "pseudotime", at: t } : null);
  };

  return (
    <figure ref={ref} className="min-w-0 max-w-[36rem]" aria-label={`${gene} along pseudotime`}>
      <svg
        width={w}
        height={h}
        className="block overflow-visible"
        onPointerMove={onMove}
        onPointerLeave={() => setHighlight(null)}
      >
        <text x={m.left} y={12} className="fill-ink text-xs font-semibold">
          {gene}
        </text>
        <text x={w - m.right} y={12} textAnchor="end" className="num fill-ink-3 text-2xs">
          r² {best.toFixed(2)}
        </text>
        {windowed && (
          <rect
            x={x(ptRange[0])}
            y={m.top}
            width={Math.max(1, x(ptRange[1]) - x(ptRange[0]))}
            height={h - m.top - m.bottom}
            className="fill-live/10"
          />
        )}
        {y.ticks(3).map((t) => (
          <g key={t}>
            <line x1={m.left} x2={w - m.right} y1={y(t)} y2={y(t)} className="stroke-line" />
            <text x={m.left - 6} y={y(t) + 3} textAnchor="end" className="num fill-ink-3 text-2xs">
              {t}
            </text>
          </g>
        ))}
        {[0, 0.5, 1].map((t) => (
          <text key={t} x={x(t)} y={h - 4} textAnchor="middle" className="num fill-ink-3 text-2xs">
            {t}
          </text>
        ))}
        {series.map((s) => (
          <path
            key={s.id}
            d={path(s.x.map((xv, i) => [xv, s.y[i] ?? null])) ?? ""}
            pathLength={1}
            fill="none"
            stroke={rgbaCss(categorical(s.id, 255))}
            strokeWidth={2}
            strokeLinecap="round"
            className="draw-on"
          >
            <title>{`${s.name} (r² ${s.r2.toFixed(2)})`}</title>
          </path>
        ))}
        {at !== null && (
          <g className="pointer-events-none">
            <line x1={x(at)} x2={x(at)} y1={m.top} y2={h - m.bottom} className="stroke-ink-3" strokeDasharray="2 3" />
            {series.map((s) => {
              const v = valueAt(s, at);
              return v === null ? null : (
                <circle key={s.id} cx={x(at)} cy={y(v)} r={3.5} fill={rgbaCss(categorical(s.id, 255))} className="stroke-panel" strokeWidth={2} />
              );
            })}
          </g>
        )}
      </svg>
    </figure>
  );
}

function valueAt(s: Series, t: number): number | null {
  let best = -1;
  let dist = Infinity;
  for (let i = 0; i < s.x.length; i++) {
    const d = Math.abs(s.x[i]! - t);
    if (d < dist && s.y[i] !== null) {
      dist = d;
      best = i;
    }
  }
  return best >= 0 ? (s.y[best] ?? null) : null;
}

/* ---- Heatmap ------------------------------------------------------------------------- */

const ROW = 14;

function Heatmap({ run }: { run: T.TrajectoryRun }) {
  const focus = useUi((s) => s.focusLineage);
  const setFocus = useUi((s) => s.setFocusLineage);
  const add = useUi((s) => s.addTrendGene);
  const setColorBy = useUi((s) => s.setColorBy);
  const setHighlight = useUi((s) => s.setHighlight);
  const highlight = useUi((s) => s.highlight);
  // A focused lineage wins; otherwise the lineage picked here.
  const [picked, setPicked] = useState(0);
  const lineage = focus ?? picked;
  const setLineage = (id: number) => {
    setPicked(id);
    if (focus !== null) setFocus(id);
  };
  const heat = useHeatmap(run.run_id, lineage, 30);
  const canvas = useRef<HTMLCanvasElement>(null);
  const data = heat.data;
  const [row, setRow] = useState<number | null>(null);

  useEffect(() => {
    const el = canvas.current;
    if (!el || !data) return;
    const rows = data.genes.length;
    const cols = data.bins.length;
    el.width = cols;
    el.height = rows;
    const ctx = el.getContext("2d");
    if (!ctx) return;
    const img = ctx.createImageData(cols, rows);
    for (let r = 0; r < rows; r++) {
      for (let c = 0; c < cols; c++) {
        const [cr, cg, cb] = sequential(data.values[r * cols + c]!, 255);
        const o = (r * cols + c) * 4;
        img.data[o] = cr;
        img.data[o + 1] = cg;
        img.data[o + 2] = cb;
        img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }, [data]);

  const at = highlight?.kind === "pseudotime" ? highlight.at : null;

  return (
    <div className="flex h-full gap-5">
      <div className="w-52 shrink-0 space-y-3">
        <Select
          label="Lineage for the heatmap"
          className="w-full"
          value={String(lineage)}
          onChange={(v) => setLineage(Number(v))}
          groups={[{ items: run.lineages.map((l) => ({ value: String(l.id), label: l.name })) }]}
        />
        <Hint>
          The 30 genes whose expression changes most along this lineage, ordered by when they peak. Each row is scaled
          to its own range. Click a gene to plot it.
        </Hint>
      </div>
      {heat.isLoading ? (
        <div className="flex flex-1 flex-col gap-1.5 py-1" aria-label="Ranking genes">
          {Array.from({ length: 8 }, (_, i) => (
            <div key={i} className="skeleton h-3" style={{ width: `${92 - i * 4}%` }} />
          ))}
        </div>
      ) : data ? (
        <div className="flex min-w-0 flex-1 gap-3 overflow-y-auto" onMouseLeave={() => setRow(null)}>
          <ul className="shrink-0" aria-label="Top temporal genes" translate="no">
            {data.genes.map((g, i) => (
              <li key={g} style={{ height: ROW, lineHeight: `${ROW}px` }}>
                <button
                  type="button"
                  title={`r² ${data.r2[i]?.toFixed(2)}`}
                  onMouseEnter={() => setRow(i)}
                  onClick={() => {
                    add(g);
                    setColorBy({ kind: "gene", gene: g });
                  }}
                  className={cx(
                    "block text-2xs leading-[14px] transition-colors duration-100",
                    row === i ? "text-focus" : "text-ink-2 [@media(hover:hover)]:hover:text-focus",
                  )}
                >
                  {g}
                </button>
              </li>
            ))}
          </ul>
          <div className="relative min-w-0 flex-1">
            <canvas
              ref={canvas}
              aria-label="Heatmap of the top temporal genes"
              className="block w-full rounded-sm"
              style={{ height: data.genes.length * ROW, imageRendering: "pixelated" }}
              onPointerMove={(e) => {
                const r = e.currentTarget.getBoundingClientRect();
                setRow(Math.floor((e.clientY - r.top) / ROW));
                setHighlight({ kind: "pseudotime", at: Math.min(1, Math.max(0, (e.clientX - r.left) / r.width)) });
              }}
              onPointerLeave={() => setHighlight(null)}
            />
            {row !== null && (
              <div
                className="pointer-events-none absolute right-0 left-0 rounded-sm ring-1 ring-ink"
                style={{ top: row * ROW, height: ROW }}
                aria-hidden
              />
            )}
            {at !== null && (
              <div
                className="pointer-events-none absolute w-px bg-ink"
                style={{ left: `${at * 100}%`, top: 0, height: data.genes.length * ROW }}
                aria-hidden
              />
            )}
            <div className="num mt-1 flex justify-between text-2xs text-ink-3">
              <span>pseudotime 0</span>
              <span>1</span>
            </div>
          </div>
        </div>
      ) : null}
    </div>
  );
}
