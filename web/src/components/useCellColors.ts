import { useMemo } from "react";
import { useDataset, useDisagreement, useGeneExpression, useObs, useRun, useRunCells, useDoneRunId } from "@/api/queries";
import { categoryLabel, columnLabel } from "@/lib/colorBy";
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

export interface LegendItem {
  label: string;
  color: RGBA;
  /** What hovering this item highlights on the field. */
  highlight: { kind: "category"; column: string; value: string } | { kind: "lineage"; id: number } | null;
}

export type LegendSpec =
  | { kind: "categorical"; title: string; items: LegendItem[] }
  | { kind: "continuous"; title: string; domain: [number, number] };

export interface CellColors {
  colors: Uint8Array | null;
  legend: LegendSpec | null;
  /** Human-readable value of cell i for the hover card. */
  describe: (i: number) => string | null;
  loading: boolean;
}

const UNCOLOURED: Record<"dark" | "light", RGBA> = {
  dark: [150, 162, 182, 190],
  light: [92, 102, 120, 190],
};

export function useCellColors(): CellColors {
  const colorBy = useUi((s) => s.colorBy);
  const theme = useUi((s) => s.theme);
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
    const none: CellColors = {
      colors: n ? uniform(n, UNCOLOURED[theme]) : null,
      legend: null,
      describe: () => null,
      loading: false,
    };
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
      const column = colorBy.column;
      return {
        colors: colorsForCodes(codes, palette),
        legend: {
          kind: "categorical",
          title: columnLabel(column),
          items: v.categories.map((c, i) => ({
            label: categoryLabel(column, c),
            color: palette[i]!,
            highlight: { kind: "category", column, value: c },
          })),
        },
        describe: (i) => (codes[i]! >= 0 ? categoryLabel(column, v.categories[codes[i]!]!) : "missing"),
        loading: false,
      };
    }

    if (colorBy.kind === "gene") {
      const values = expr.data;
      if (!values) return { ...none, loading: expr.isLoading };
      const domain = robustDomain(values);
      return {
        colors: colorsForValues(values, domain),
        legend: { kind: "continuous", title: `${colorBy.gene}, log-normalised`, domain },
        describe: (i) => `${colorBy.gene} ${values[i]!.toFixed(2)}`,
        loading: false,
      };
    }

    if (colorBy.kind === "disagreement") {
      const diff = disagreement.data;
      if (!diff) return { ...none, loading: disagreement.isLoading };
      const domain = robustDomain(diff);
      return {
        colors: colorsForValues(diff, domain),
        legend: { kind: "continuous", title: "Rank difference between runs", domain },
        describe: (i) =>
          Number.isFinite(diff[i]!) ? `ranks differ by ${(diff[i]! * 100).toFixed(1)}%` : "not in both runs",
        loading: false,
      };
    }

    const c = cells.data;
    if (!c) return { ...none, loading: cells.isLoading };
    if (colorBy.kind === "pseudotime") {
      const pt = c.pseudotime;
      return {
        colors: colorsForValues(pt, [0, 1]),
        legend: { kind: "continuous", title: "Pseudotime", domain: [0, 1] },
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
        title: "Lineage",
        items: [
          ...lineages.map((l) => ({
            label: l.name,
            color: categorical(l.id),
            highlight: { kind: "lineage" as const, id: l.id },
          })),
          { label: "Shared trunk", color: SHARED, highlight: null },
          { label: "Unassigned", color: MISSING, highlight: null },
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
    theme,
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

/** Which cells the current hover highlight names (1 = named), or null when nothing is highlighted. */
export function useHighlightMask(): Uint8Array | null {
  const highlight = useUi((s) => s.highlight);
  const runId = useDoneRunId(useUi((s) => s.runId));
  const n = useDataset().data?.n_obs ?? 0;
  const obs = useObs(highlight?.kind === "category" ? highlight.column : null);
  const cells = useRunCells(highlight?.kind === "lineage" || highlight?.kind === "pseudotime" ? runId : null);

  return useMemo(() => {
    if (!highlight || !n) return null;
    const mask = new Uint8Array(n);
    if (highlight.kind === "category") {
      const v = obs.data;
      const code = v?.categories.indexOf(highlight.value) ?? -1;
      if (!v?.codes || code < 0) return null;
      for (let i = 0; i < n; i++) mask[i] = v.codes[i] === code ? 1 : 0;
      return mask;
    }
    const c = cells.data;
    if (!c) return null;
    if (highlight.kind === "lineage") {
      return c.membership.slice(highlight.id * n, highlight.id * n + n);
    }
    const half = 0.035;
    for (let i = 0; i < n; i++) {
      const pt = c.pseudotime[i]!;
      mask[i] = Math.abs(pt - highlight.at) <= half ? 1 : 0;
    }
    return mask;
  }, [highlight, n, obs.data, cells.data]);
}

/**
 * Dim what is not in focus: cells outside the selection, and cells the hover highlight
 * does not name. The colour attribute transitions on the GPU, so the dimming eases in.
 */
export function withEmphasis(
  colors: Uint8Array,
  selection: Uint8Array | null,
  highlight: Uint8Array | null,
): Uint8Array {
  if (!selection && !highlight) return colors;
  const out = colors.slice();
  const n = out.length / 4;
  for (let i = 0; i < n; i++) {
    let a = out[i * 4 + 3]!;
    if (selection && !selection[i]) a *= 0.16;
    if (highlight && !highlight[i]) a *= 0.12;
    out[i * 4 + 3] = Math.round(a);
  }
  return out;
}

