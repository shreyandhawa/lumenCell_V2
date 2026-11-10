import { useMemo } from "react";
import { useDataset, useDisagreement, useGeneExpression, useObs, useRun, useRunCells, useDoneRunId } from "@/api/queries";
import {
  categoryPalette,
  categorical,
  colorsForCodes,
  colorsForLineages,
  colorsForValues,
  MISSING,
  robustDomain,
  type RGBA,
  SHARED,
  uniform,
} from "@/lib/colors";
import { useUi } from "@/state/ui";

export type LegendSpec =
  | { kind: "categorical"; title: string; items: { label: string; color: RGBA }[] }
  | { kind: "continuous"; title: string; domain: [number, number] };

export interface CellColors {
  colors: Uint8Array | null;
  legend: LegendSpec | null;
  /** Human-readable value of cell i for the tooltip. */
  describe: (i: number) => string | null;
  loading: boolean;
}

const DEFAULT: RGBA = [148, 163, 184, 200];
const LEIDEN = "lumencell_leiden";

export function useCellColors(): CellColors {
  const colorBy = useUi((s) => s.colorBy);
  const runId = useDoneRunId(useUi((s) => s.runId));
  const dataset = useDataset().data;
  const n = dataset?.n_obs ?? 0;

  const obsColumn = colorBy?.kind === "obs" ? colorBy.column : null;
  const gene = colorBy?.kind === "gene" ? colorBy.gene : null;
  const needsRun = colorBy?.kind === "pseudotime" || colorBy?.kind === "lineage";
  const obs = useObs(obsColumn);
  const expr = useGeneExpression(gene);
  const cells = useRunCells(needsRun ? runId : null);
  const run = useRun(colorBy?.kind === "lineage" ? runId : null);
  const disagreement = useDisagreement(
    colorBy?.kind === "disagreement" ? colorBy.a : null,
    colorBy?.kind === "disagreement" ? colorBy.b : null,
  );

  return useMemo<CellColors>(() => {
    const none: CellColors = { colors: n ? uniform(n, DEFAULT) : null, legend: null, describe: () => null, loading: false };
    if (!colorBy || !n) return none;

    if (colorBy.kind === "obs") {
      const v = obs.data;
      if (!v) return { ...none, loading: obs.isLoading };
      if (v.kind === "continuous" && v.values) {
        const domain = robustDomain(v.values);
        const values = v.values;
        return {
          colors: colorsForValues(values, domain),
          legend: { kind: "continuous", title: colorBy.column, domain },
          describe: (i) => (Number.isFinite(values[i]!) ? values[i]!.toPrecision(3) : "missing"),
          loading: false,
        };
      }
      const info = dataset?.obs.find((c) => c.name === colorBy.column);
      const palette = categoryPalette(v.categories.length, info?.colors);
      const codes = v.codes!;
      const leiden = colorBy.column === LEIDEN;
      const label = (c: string) => (leiden ? `Cluster ${c}` : c);
      return {
        colors: colorsForCodes(codes, palette),
        legend: {
          kind: "categorical",
          title: leiden ? "Leiden clusters" : colorBy.column,
          items: v.categories.map((c, i) => ({ label: label(c), color: palette[i]! })),
        },
        describe: (i) => (codes[i]! >= 0 ? label(v.categories[codes[i]!]!) : "missing"),
        loading: false,
      };
    }

    if (colorBy.kind === "gene") {
      const values = expr.data;
      if (!values) return { ...none, loading: expr.isLoading };
      const domain = robustDomain(values);
      return {
        colors: colorsForValues(values, domain),
        legend: { kind: "continuous", title: `${colorBy.gene} (log-normalised)`, domain },
        describe: (i) => values[i]!.toFixed(2),
        loading: false,
      };
    }

    if (colorBy.kind === "disagreement") {
      const diff = disagreement.data;
      if (!diff) return { ...none, loading: disagreement.isLoading };
      const domain = robustDomain(diff);
      return {
        colors: colorsForValues(diff, domain),
        legend: { kind: "continuous", title: "rank difference between runs", domain },
        describe: (i) => (Number.isFinite(diff[i]!) ? `ranks differ by ${(diff[i]! * 100).toFixed(1)}%` : "not in both runs"),
        loading: false,
      };
    }

    const c = cells.data;
    if (!c) return { ...none, loading: cells.isLoading };
    if (colorBy.kind === "pseudotime") {
      const pt = c.pseudotime;
      return {
        colors: colorsForValues(pt, [0, 1]),
        legend: { kind: "continuous", title: "pseudotime", domain: [0, 1] },
        describe: (i) => (Number.isFinite(pt[i]!) ? `pseudotime ${pt[i]!.toFixed(3)}` : "no pseudotime"),
        loading: false,
      };
    }
    const lineages = run.data?.lineages ?? [];
    const primary = c.primaryLineage;
    return {
      colors: colorsForLineages(primary),
      legend: {
        kind: "categorical",
        title: "lineage",
        items: [
          ...lineages.map((l) => ({ label: l.name, color: categorical(l.id) })),
          { label: "shared trunk", color: SHARED },
          { label: "unassigned", color: MISSING },
        ],
      },
      describe: (i) => {
        const p = primary[i]!;
        return p >= 0 ? (lineages[p]?.name ?? `lineage ${p + 1}`) : p === -1 ? "shared trunk" : "unassigned";
      },
      loading: false,
    };
  }, [
    colorBy,
    n,
    obs.data,
    obs.isLoading,
    expr.data,
    expr.isLoading,
    cells.data,
    cells.isLoading,
    run.data,
    dataset,
    disagreement.data,
    disagreement.isLoading,
  ]);
}

/** Dim cells outside the selection, so the selection reads at a glance. */
export function withSelection(colors: Uint8Array, selection: Uint8Array | null): Uint8Array {
  if (!selection) return colors;
  const out = colors.slice();
  for (let i = 0; i < selection.length; i++) {
    if (!selection[i]) out[i * 4 + 3] = Math.round(out[i * 4 + 3]! * 0.18);
  }
  return out;
}
