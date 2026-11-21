import { useQueries } from "@tanstack/react-query";
import { line } from "d3-shape";
import { motion } from "motion/react";
import { api } from "@/api/client";
import type * as T from "@/api/generated";
import { useHeatmap } from "@/api/queries";
import { TIMELINE_COLS } from "@/components/trajectory/PseudotimeBar";
import { Label } from "@/components/ui";
import { MAX_TREND_GENES, useUi } from "@/state/ui";

const H = 28;

/** Linear interpolation of a binned trend at pseudotime t; null where the bins are empty. */
function at(x: number[], y: (number | null)[], t: number): number | null {
  const k = x.findIndex((v) => v >= t);
  if (k === -1) return y[y.length - 1] ?? null;
  if (k === 0) return y[0] ?? null;
  const a = y[k - 1] ?? null;
  const b = y[k] ?? null;
  if (a === null || b === null) return a ?? b;
  return a + ((b - a) * (t - x[k - 1]!)) / (x[k]! - x[k - 1]!);
}

/**
 * Gene expression along pseudotime, on the bar's x axis. One line per lineage (the focused one in
 * the accent). Genes: the ones the user picked, then the run's top temporal genes.
 */
export function GenePanel({ run }: { run: T.TrajectoryRun }) {
  const picked = useUi((s) => s.trendGenes);
  const focus = useUi((s) => s.focusLineage);
  const cursor = useUi((s) => s.cursor) ?? 0;
  const top = useHeatmap(run.run_id, focus ?? 0, MAX_TREND_GENES);
  const genes = [...new Set([...picked, ...(top.data?.genes ?? [])])].slice(0, MAX_TREND_GENES);
  const results = useQueries({
    queries: run.lineages.map((l) => ({
      queryKey: ["trends", run.run_id, l.id, genes.join(",")],
      queryFn: () => api.trends(run.run_id, l.id, genes),
      enabled: genes.length > 0,
      staleTime: Infinity,
    })),
  });
  const focusName = focus === null ? null : run.lineages[focus]?.name;

  return (
    <motion.section
      aria-label="Gene trends"
      initial={{ opacity: 0, y: 8 }}
      animate={{ opacity: 1, y: 0 }}
      exit={{ opacity: 0, y: 8, transition: { duration: 0.15 } }}
      transition={{ duration: 0.2, ease: [0.23, 1, 0.32, 1] }}
      className="rounded-xl bg-panel pt-3 pb-2"
    >
      <div className="px-3">
        <Label>Gene trends · {focusName ?? "all lineages"} · same axis as the bar</Label>
      </div>
      {genes.length === 0 ? (
        <p className="px-3 pb-2 text-xs text-ink-muted">{top.isError ? "Could not rank genes." : "Ranking genes…"}</p>
      ) : (
        <ul aria-label="Genes">
          {genes.map((gene, k) => {
            const series = run.lineages.flatMap((l, i) => {
              const data = results[i]?.data;
              const y = data?.genes[k]?.y;
              return data && y ? [{ id: l.id, x: data.x, y }] : [];
            });
            const yMax = Math.max(1e-6, ...series.flatMap((s) => s.y.filter((v): v is number => v !== null)));
            const path = line<[number, number | null]>()
              .defined((d) => d[1] !== null)
              .x((d) => d[0] * 100)
              .y((d) => H - 2 - ((d[1] ?? 0) / yMax) * (H - 4));
            const shown = focus === null ? series : series.filter((s) => s.id === focus);
            const values = shown.map((s) => at(s.x, s.y, cursor)).filter((v): v is number => v !== null);
            const value = values.length ? values.reduce((a, b) => a + b, 0) / values.length : null;
            return (
              <li key={gene} className={`${TIMELINE_COLS} h-10 items-center`}>
                <span className="truncate text-right font-mono text-2xs text-ink" title={gene}>
                  {gene}
                </span>
                <svg viewBox={`0 0 100 ${H}`} preserveAspectRatio="none" className="h-7 w-full" aria-hidden>
                  <line x1={0} x2={100} y1={H - 2} y2={H - 2} className="stroke-line" vectorEffect="non-scaling-stroke" />
                  {/* the focused lineage last, so it draws on top */}
                  {[...series.filter((s) => s.id !== focus), ...series.filter((s) => s.id === focus)].map((s) => (
                      <path
                        key={s.id}
                        d={path(s.x.map((xv, i) => [xv, s.y[i] ?? null])) ?? ""}
                        fill="none"
                        vectorEffect="non-scaling-stroke"
                        strokeWidth={s.id === focus ? 2 : 1.25}
                        className={s.id === focus ? "stroke-accent" : "stroke-ink-muted"}
                      />
                    ))}
                </svg>
                <span
                  className="font-mono text-xs text-ink"
                  title={focus === null ? "Mean of the lineage trends at the cursor" : "At the cursor"}
                >
                  {value === null ? "–" : value.toFixed(2)}
                </span>
              </li>
            );
          })}
        </ul>
      )}
    </motion.section>
  );
}
