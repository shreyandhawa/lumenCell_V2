import { useQueries } from "@tanstack/react-query";
import { scaleLinear } from "d3-scale";
import { line } from "d3-shape";
import { ChevronDown, ChevronUp, X } from "lucide-react";
import { useEffect, useRef, useState } from "react";
import { api } from "@/api/client";
import type * as T from "@/api/generated";
import { useDoneRunId, useHeatmap, useRun } from "@/api/queries";
import { GeneSearch } from "@/components/GeneSearch";
import { Hint, Select } from "@/components/ui";
import { categorical, rgbaCss, sequential } from "@/lib/colors";
import { MAX_TREND_GENES, useUi } from "@/state/ui";

type Tab = "trends" | "heatmap";

/** Gene expression along pseudotime, under the plot. Shown once a run has finished. */
export function GenePanel() {
  const runId = useDoneRunId(useUi((s) => s.runId));
  const run = useRun(runId).data;
  const [open, setOpen] = useState(true);
  const [tab, setTab] = useState<Tab>("trends");
  if (!run) return null;
  return (
    <section aria-label="Genes along pseudotime" className="shrink-0 border-t border-line bg-panel">
      <header className="flex h-9 items-center gap-1 px-3">
        <h2 className="mr-2 text-2xs font-semibold tracking-wide text-ink-muted uppercase">Genes along pseudotime</h2>
        <div role="tablist" className="flex gap-1">
          {(["trends", "heatmap"] as const).map((t) => (
            <button
              key={t}
              role="tab"
              type="button"
              aria-selected={tab === t}
              onClick={() => {
                setTab(t);
                setOpen(true);
              }}
              className={`rounded px-2 py-0.5 text-xs ${tab === t ? "bg-accent/15 text-accent" : "text-ink-muted hover:text-ink"}`}
            >
              {t === "trends" ? "Trends" : "Top temporal genes"}
            </button>
          ))}
        </div>
        <button
          type="button"
          aria-label={open ? "Collapse gene panel" : "Expand gene panel"}
          onClick={() => setOpen(!open)}
          className="ml-auto text-ink-muted hover:text-ink"
        >
          {open ? <ChevronDown size={16} /> : <ChevronUp size={16} />}
        </button>
      </header>
      {open && <div className="h-56 overflow-hidden px-3 pb-3">{tab === "trends" ? <Trends run={run} /> : <Heatmap run={run} />}</div>}
    </section>
  );
}

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
    <div className="flex h-full gap-3">
      <div className="w-52 shrink-0 space-y-2">
        <GeneSearch onPick={add} placeholder="Add a gene…" />
        <ul className="flex flex-wrap gap-1" aria-label="Genes in the trends panel">
          {genes.map((g) => (
            <li key={g} className="flex items-center gap-1 rounded border border-line bg-panel-2 py-0.5 pr-1 pl-1.5 font-mono text-2xs">
              {g}
              <button type="button" aria-label={`Remove ${g}`} onClick={() => remove(g)} className="text-ink-faint hover:text-ink">
                <X size={11} />
              </button>
            </li>
          ))}
        </ul>
        <Hint>
          Up to {MAX_TREND_GENES} genes. One line per lineage
          {focus !== null ? " (focused lineage only)" : ""}; mean log-normalised expression, smoothed.
        </Hint>
      </div>
      {genes.length === 0 ? (
        <div className="grid flex-1 place-items-center text-xs text-ink-faint">
          Add a gene, or pick one from “Top temporal genes”.
        </div>
      ) : (
        <div className="flex flex-1 gap-3 overflow-x-auto" aria-label="Gene trends">
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

function TrendChart({ gene, series }: { gene: string; series: Series[] }) {
  const w = 230;
  const h = 170;
  const m = { top: 18, right: 8, bottom: 20, left: 30 };
  const yMax = Math.max(0.1, ...series.flatMap((s) => s.y.filter((v): v is number => v !== null)));
  const x = scaleLinear().domain([0, 1]).range([m.left, w - m.right]);
  const y = scaleLinear().domain([0, yMax]).range([h - m.bottom, m.top]).nice();
  const path = line<[number, number | null]>()
    .defined((d) => d[1] !== null)
    .x((d) => x(d[0]))
    .y((d) => y(d[1] ?? 0));
  const best = series.length ? Math.max(...series.map((s) => s.r2)) : 0;
  return (
    <figure className="shrink-0" aria-label={`${gene} along pseudotime`}>
      <svg width={w} height={h} className="overflow-visible">
        <text x={m.left} y={11} className="fill-ink font-mono text-[11px] font-semibold">
          {gene}
        </text>
        <text x={w - m.right} y={11} textAnchor="end" className="fill-ink-faint text-[10px]">
          r² {best.toFixed(2)}
        </text>
        <line x1={m.left} x2={w - m.right} y1={h - m.bottom} y2={h - m.bottom} className="stroke-line" />
        <line x1={m.left} x2={m.left} y1={m.top} y2={h - m.bottom} className="stroke-line" />
        {[0, 0.5, 1].map((t) => (
          <text key={t} x={x(t)} y={h - 6} textAnchor="middle" className="fill-ink-faint text-[10px]">
            {t}
          </text>
        ))}
        {y.ticks(3).map((t) => (
          <text key={t} x={m.left - 4} y={y(t) + 3} textAnchor="end" className="fill-ink-faint text-[10px]">
            {t}
          </text>
        ))}
        {series.map((s) => (
          <path
            key={s.id}
            d={path(s.x.map((xv, i) => [xv, s.y[i] ?? null])) ?? ""}
            fill="none"
            stroke={rgbaCss(categorical(s.id, 255))}
            strokeWidth={2}
          >
            <title>{`${s.name} (r² ${s.r2.toFixed(2)})`}</title>
          </path>
        ))}
      </svg>
    </figure>
  );
}

function Heatmap({ run }: { run: T.TrajectoryRun }) {
  const focus = useUi((s) => s.focusLineage);
  const setFocus = useUi((s) => s.setFocusLineage);
  const add = useUi((s) => s.addTrendGene);
  const setColorBy = useUi((s) => s.setColorBy);
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

  return (
    <div className="flex h-full gap-3">
      <div className="w-52 shrink-0 space-y-2">
        <Select aria-label="Lineage for the heatmap" value={lineage} onChange={(e) => setLineage(Number(e.target.value))}>
          {run.lineages.map((l) => (
            <option key={l.id} value={l.id}>
              {l.name}
            </option>
          ))}
        </Select>
        <Hint>
          The 30 genes whose expression changes most along this lineage, ordered by when they peak. Each row is
          scaled to its own range. Click a gene to plot it.
        </Hint>
      </div>
      {heat.isLoading ? (
        <div className="grid flex-1 place-items-center text-xs text-ink-faint">Ranking genes…</div>
      ) : data ? (
        <div className="flex min-w-0 flex-1 gap-2 overflow-y-auto">
          <ul className="shrink-0 py-0" aria-label="Top temporal genes">
            {data.genes.map((g, i) => (
              <li key={g} className="h-[14px] leading-[14px]">
                <button
                  type="button"
                  title={`r² ${data.r2[i]?.toFixed(2)}`}
                  onClick={() => {
                    add(g);
                    setColorBy({ kind: "gene", gene: g });
                  }}
                  className="font-mono text-[10px] text-ink-muted hover:text-accent"
                >
                  {g}
                </button>
              </li>
            ))}
          </ul>
          <div className="min-w-0 flex-1">
            <canvas
              ref={canvas}
              aria-label="Heatmap of the top temporal genes"
              className="w-full"
              style={{ height: data.genes.length * 14, imageRendering: "pixelated" }}
            />
            <div className="flex justify-between font-mono text-2xs text-ink-faint">
              <span>pseudotime 0</span>
              <span>1</span>
            </div>
          </div>
        </div>
      ) : null}
    </div>
  );
}
