import { OrthographicView, type OrthographicViewState } from "@deck.gl/core";
import { DataFilterExtension } from "@deck.gl/extensions";
import { PathLayer, ScatterplotLayer, TextLayer } from "@deck.gl/layers";
import DeckGL, { type DeckGLRef } from "@deck.gl/react";
import { Lasso, Maximize2, X } from "lucide-react";
import { useEffect, useMemo, useRef, useState, type PointerEvent } from "react";
import { useCurves, useDataset, useEmbedding, useRun, useRunCells, useDoneRunId } from "@/api/queries";
import { Legend } from "@/components/Legend";
import { useCellColors, withSelection } from "@/components/useCellColors";
import { categorical } from "@/lib/colors";
import { bounds, fitView, selectInPolygon, type Point } from "@/lib/geometry";
import { useUi } from "@/state/ui";

const VIEW = new OrthographicView({ id: "ortho", flipY: false });
const FILTER = new DataFilterExtension({ filterSize: 2 });

function useSize(): [React.RefObject<HTMLDivElement | null>, { width: number; height: number }] {
  const ref = useRef<HTMLDivElement>(null);
  const [size, setSize] = useState({ width: 0, height: 0 });
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const ro = new ResizeObserver(([entry]) => {
      if (entry) setSize({ width: entry.contentRect.width, height: entry.contentRect.height });
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);
  return [ref, size];
}

function pointRadius(n: number): number {
  if (n < 5_000) return 3;
  if (n < 50_000) return 2;
  return 1.4;
}

export function Scatter() {
  const embedding = useUi((s) => s.embedding);
  const runId = useDoneRunId(useUi((s) => s.runId));
  const ptRange = useUi((s) => s.ptRange);
  const focusLineage = useUi((s) => s.focusLineage);
  const selection = useUi((s) => s.selection);
  const lassoActive = useUi((s) => s.lassoActive);
  const setLasso = useUi((s) => s.setLasso);
  const setSelection = useUi((s) => s.setSelection);

  const n = useDataset().data?.n_obs ?? 0;
  const positions = useEmbedding(embedding).data;
  const { colors, legend, describe } = useCellColors();
  const cells = useRunCells(runId).data;
  const run = useRun(runId).data;
  const curves = useCurves(runId, embedding).data;

  const [containerRef, size] = useSize();
  const deckRef = useRef<DeckGLRef>(null);
  const [viewState, setViewState] = useState<OrthographicViewState | null>(null);

  // Fit the view whenever the embedding (or the canvas size, before the first fit) changes.
  const fitted = useRef<string | null>(null);
  useEffect(() => {
    if (!positions || size.width === 0 || fitted.current === embedding) return;
    fitted.current = embedding;
    setViewState({ ...fitView(bounds(positions), size.width, size.height, 0.12), minZoom: -10, maxZoom: 20 });
  }, [positions, embedding, size.width, size.height]);

  const resetView = () => {
    if (positions && size.width) setViewState({ ...fitView(bounds(positions), size.width, size.height, 0.12) });
  };

  // [pseudotime, on-focused-lineage] per cell, read on the GPU by DataFilterExtension.
  const filterValues = useMemo(() => {
    const out = new Float32Array(n * 2);
    for (let i = 0; i < n; i++) {
      const pt = cells?.pseudotime[i];
      out[i * 2] = pt === undefined ? 0.5 : Number.isFinite(pt) ? pt : -0.5;
      out[i * 2 + 1] =
        focusLineage === null || !cells ? 1 : cells.membership[focusLineage * n + i]! ? 1 : 0;
    }
    return out;
  }, [n, cells, focusLineage]);

  const ptFiltered = cells !== undefined && (ptRange[0] > 0 || ptRange[1] < 1);
  const filterRange: [number, number][] = [ptFiltered ? ptRange : [-1, 2], [0.5, 1.5]];

  const fill = useMemo(() => (colors ? withSelection(colors, selection) : null), [colors, selection]);

  const layers = [];
  if (positions && fill && positions.length === n * 2) {
    layers.push(
      new ScatterplotLayer({
        id: "cells",
        data: {
          length: n,
          attributes: {
            getPosition: { value: positions, size: 2 },
            getFillColor: { value: fill, size: 4, normalized: true },
            getFilterValue: { value: filterValues, size: 2 },
          },
        },
        radiusUnits: "pixels",
        getRadius: pointRadius(n),
        stroked: false,
        pickable: true,
        autoHighlight: true,
        highlightColor: [255, 255, 255, 220],
        extensions: [FILTER],
        filterRange,
      }),
    );
  }
  if (curves && run) {
    const paths = curves.curves
      .filter((c) => c.points.length > 1)
      .map((c) => ({ lineage: c.lineage, path: c.points as [number, number][] }));
    const dim = (id: number) => (focusLineage === null || focusLineage === id ? 255 : 60);
    layers.push(
      new PathLayer({
        id: "curve-casing",
        data: paths,
        getPath: (d) => d.path,
        getColor: (d) => [10, 12, 18, Math.round(dim(d.lineage) * 0.8)],
        getWidth: 6,
        widthUnits: "pixels",
        jointRounded: true,
        capRounded: true,
        updateTriggers: { getColor: focusLineage },
      }),
      new PathLayer({
        id: "curves",
        data: paths,
        getPath: (d) => d.path,
        getColor: (d) => {
          const c = categorical(d.lineage);
          return [c[0], c[1], c[2], dim(d.lineage)];
        },
        getWidth: 3,
        widthUnits: "pixels",
        jointRounded: true,
        capRounded: true,
        updateTriggers: { getColor: focusLineage },
      }),
    );
    const rootId = run.milestones.find((m) => m.kind === "root")?.id;
    layers.push(
      new ScatterplotLayer({
        id: "milestones",
        data: curves.milestones,
        getPosition: (m) => [m.x, m.y],
        getRadius: (m) => (m.id === rootId ? 7 : 4),
        radiusUnits: "pixels",
        getFillColor: (m) => (m.id === rootId ? [255, 255, 255, 255] : [235, 238, 245, 230]),
        getLineColor: [10, 12, 18, 255],
        lineWidthUnits: "pixels",
        getLineWidth: (m) => (m.id === rootId ? 3 : 1.5),
        stroked: true,
      }),
      new TextLayer({
        id: "lineage-labels",
        data: paths.map((p) => ({ ...p, name: run.lineages[p.lineage]?.name ?? "" })),
        getPosition: (d) => d.path[d.path.length - 1]!,
        getText: (d) => d.name,
        getSize: 12,
        getColor: [240, 242, 248, 255],
        getTextAnchor: "middle",
        getAlignmentBaseline: "bottom",
        getPixelOffset: [0, -10],
        background: true,
        getBackgroundColor: [16, 18, 26, 210],
        backgroundPadding: [4, 2],
        fontFamily: "Inter, system-ui, sans-serif",
        fontWeight: 600,
      }),
    );
  }

  return (
    <div ref={containerRef} className="relative h-full w-full overflow-hidden bg-surface" data-testid="scatter">
      {viewState && (
        <DeckGL
          ref={deckRef}
          views={VIEW}
          viewState={viewState}
          onViewStateChange={({ viewState: vs }) => setViewState(vs as OrthographicViewState)}
          controller={{ dragPan: !lassoActive, scrollZoom: { smooth: true } }}
          layers={layers}
          getTooltip={({ index, layer }) =>
            layer?.id === "cells" && index >= 0 ? { text: describe(index) ?? `cell ${index}` } : null
          }
          getCursor={({ isDragging }) => (lassoActive ? "crosshair" : isDragging ? "grabbing" : "default")}
        />
      )}
      {!positions && (
        <div className="absolute inset-0 grid place-items-center text-sm text-ink-faint">
          {embedding ? "Loading embedding…" : "No embedding yet: preprocessing will compute one."}
        </div>
      )}
      {lassoActive && positions && (
        <LassoOverlay
          onDone={(screen) => {
            setLasso(false);
            const viewport = deckRef.current?.deck?.getViewports()[0];
            if (!viewport || screen.length < 3) return;
            const world = screen.map((p) => viewport.unproject(p as [number, number]).slice(0, 2) as unknown as Point);
            setSelection(selectInPolygon(positions, world));
          }}
        />
      )}
      <div className="absolute top-3 left-3 flex gap-1.5">
        <ToolButton label="Lasso select" active={lassoActive} onClick={() => setLasso(!lassoActive)}>
          <Lasso size={14} />
        </ToolButton>
        <ToolButton label="Fit view" onClick={resetView}>
          <Maximize2 size={14} />
        </ToolButton>
        {selection && (
          <ToolButton label="Clear selection" onClick={() => setSelection(null)}>
            <X size={14} />
          </ToolButton>
        )}
      </div>
      <div className="pointer-events-none absolute right-3 bottom-3">
        <Legend spec={legend} />
      </div>
    </div>
  );
}

function ToolButton({
  label,
  active = false,
  onClick,
  children,
}: {
  label: string;
  active?: boolean;
  onClick: () => void;
  children: React.ReactNode;
}) {
  return (
    <button
      type="button"
      aria-label={label}
      aria-pressed={active}
      title={label}
      onClick={onClick}
      className={`grid size-8 place-items-center rounded-md border shadow-sm backdrop-blur transition-colors ${
        active ? "border-accent bg-accent/20 text-accent" : "border-line bg-panel/90 text-ink-muted hover:text-ink"
      }`}
    >
      {children}
    </button>
  );
}

function LassoOverlay({ onDone }: { onDone: (points: [number, number][]) => void }) {
  const [points, setPoints] = useState<[number, number][]>([]);
  const drawing = useRef(false);
  const local = (e: PointerEvent<HTMLDivElement>): [number, number] => {
    const r = e.currentTarget.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  };
  return (
    <div
      className="absolute inset-0 cursor-crosshair"
      onPointerDown={(e) => {
        drawing.current = true;
        e.currentTarget.setPointerCapture(e.pointerId);
        setPoints([local(e)]);
      }}
      onPointerMove={(e) => {
        if (!drawing.current) return;
        // Read the event now: React clears currentTarget before a deferred updater runs.
        const point = local(e);
        setPoints((p) => [...p, point]);
      }}
      onPointerUp={() => {
        drawing.current = false;
        onDone(points);
        setPoints([]);
      }}
    >
      <svg className="pointer-events-none absolute inset-0 size-full">
        <polygon
          points={points.map((p) => p.join(",")).join(" ")}
          className="fill-accent/10 stroke-accent"
          strokeWidth={1.5}
          strokeDasharray="4 3"
        />
      </svg>
    </div>
  );
}
