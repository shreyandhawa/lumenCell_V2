/**
 * The field: every cell on one shared stage, drawn on the GPU.
 *
 * One orthographic orbit camera serves every embedding, so switching embeddings is a
 * morph of the same cells (their positions transition on the GPU) and DC 3D lifts the
 * diffusion map into its third component. Panels float over the field; the view's
 * padding keeps the cells centred in the area they leave free.
 */
import { LinearInterpolator, OrbitView, type OrbitViewState } from "@deck.gl/core";
import { DataFilterExtension } from "@deck.gl/extensions";
import { PathLayer, ScatterplotLayer, TextLayer } from "@deck.gl/layers";
import DeckGL from "@deck.gl/react";
import { Lasso, Maximize2, X } from "lucide-react";
import { animate } from "motion";
import { useEffect, useEffectEvent, useMemo, useRef, useState, type PointerEvent } from "react";
import { useCurves, useDataset, useDoneRunId, useEmbedding, useRun, useRunCells } from "@/api/queries";
import { useCellColors, useHighlightMask, withEmphasis } from "@/components/useCellColors";
import { cx, IconButton } from "@/components/ui";
import { embeddingLabel } from "@/lib/colorBy";
import { categorical, type RGBA } from "@/lib/colors";
import { DepthCueExtension } from "@/lib/depthCue";
import { fitStage, type Insets, onStage, selectInPolygon, type Stage } from "@/lib/geometry";
import { DURATION, easeInOut, easeOut, prefersReducedMotion } from "@/lib/motion";
import { softRange } from "@/lib/playback";
import { useUi } from "@/state/ui";

const FILTER = new DataFilterExtension({ filterSize: 2 });
const DEPTH = new DepthCueExtension();
const CELL_EXTENSIONS = [FILTER, DEPTH];
const CAMERA = new LinearInterpolator(["target", "zoom", "rotationX", "rotationOrbit"]);
const NO_PT = -2;
const TILT_3D = { rotationX: 22, rotationOrbit: -34 };

type ViewState = OrbitViewState & {
  transitionDuration?: number;
  transitionInterpolator?: LinearInterpolator;
  transitionEasing?: (t: number) => number;
};

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
  if (n < 5_000) return 2.8;
  if (n < 50_000) return 1.9;
  return 1.3;
}

/** A number that eases to its target (Motion drives it; reduced motion jumps). */
function useTween(target: number, duration: number, ease: (t: number) => number = easeOut): number {
  const [value, setValue] = useState(target);
  const current = useRef(target);
  const instant = duration <= 0 || prefersReducedMotion();
  useEffect(() => {
    if (instant) {
      current.current = target;
      return;
    }
    const controls = animate(current.current, target, {
      duration: duration / 1000,
      ease,
      onUpdate: (v) => {
        current.current = v;
        setValue(v);
      },
    });
    return () => controls.stop();
  }, [target, duration, ease, instant]);
  return instant ? target : value;
}

function useTweenedInsets(insets: Insets, animateIt: boolean): Insets {
  const ms = animateIt ? 320 : 0;
  return {
    top: useTween(insets.top, ms),
    right: useTween(insets.right, ms),
    bottom: useTween(insets.bottom, ms),
    left: useTween(insets.left, ms),
  };
}

export function Scatter({ insets }: { insets: Insets }) {
  const embedding = useUi((s) => s.embedding);
  const runId = useDoneRunId(useUi((s) => s.runId));
  const ptRange = useUi((s) => s.ptRange);
  const playing = useUi((s) => s.playing);
  const focusLineage = useUi((s) => s.focusLineage);
  const selection = useUi((s) => s.selection);
  const lassoActive = useUi((s) => s.lassoActive);
  const setLasso = useUi((s) => s.setLasso);
  const setSelection = useUi((s) => s.setSelection);
  const setHovered = useUi((s) => s.setHovered);
  const hovered = useUi((s) => s.hovered);
  const panelsAnimate = useUi((s) => s.panelsAnimate);
  const theme = useUi((s) => s.theme);

  const n = useDataset().data?.n_obs ?? 0;
  const stage = useEmbedding(embedding).data;
  const { colors, describe } = useCellColors();
  const highlightMask = useHighlightMask();
  const cells = useRunCells(runId).data;
  const run = useRun(runId).data;
  const curves = useCurves(runId, embedding).data;

  const [containerRef, size] = useSize();
  const [viewState, setViewState] = useState<ViewState | null>(null);
  const padding = useTweenedInsets(insets, panelsAnimate);
  const reduced = prefersReducedMotion();

  const view = useMemo(
    () => new OrbitView({ id: "stage", orthographic: true, orbitAxis: "Y", padding }),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [padding.top, padding.right, padding.bottom, padding.left],
  );

  // ---- camera -------------------------------------------------------------------------
  const fitFor = (s: Stage): ViewState => {
    const { zoom } = fitStage(s.radius, size.width, size.height, insets, s.dims);
    return {
      target: [0, 0, 0],
      zoom,
      minZoom: zoom - 3,
      maxZoom: zoom + 9,
      ...(s.dims === 3 ? TILT_3D : { rotationX: 0, rotationOrbit: 0 }),
    };
  };
  const fly = (to: ViewState, ms: number, easing = easeInOut) =>
    setViewState({
      ...to,
      transitionDuration: reduced ? 0 : ms,
      transitionInterpolator: CAMERA,
      transitionEasing: easing,
    });

  const placed = useRef<{ embedding: string | null; dims: number } | null>(null);
  useEffect(() => {
    if (!stage || size.width === 0 || stage.positions.length !== n * 3) return;
    const prev = placed.current;
    if (prev && prev.embedding === embedding) return;
    placed.current = { embedding, dims: stage.dims };
    // The first placement is instant; later ones fly with the cells' morph.
    if (!prev) setViewState(fitFor(stage));
    else fly(fitFor(stage), DURATION.morph * 1000);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [stage, embedding, size.width, n]);

  const resetView = () => stage && fly(fitFor(stage), 520, easeOut);

  // When the window is resized, keep the whole stage in view.
  const lastSize = useRef(size);
  useEffect(() => {
    const prev = lastSize.current;
    lastSize.current = size;
    if (!stage || !placed.current || prev.width === 0) return;
    if (Math.abs(prev.width - size.width) < 2 && Math.abs(prev.height - size.height) < 2) return;
    const id = window.setTimeout(() => setViewState(fitFor(stage)), 120);
    return () => window.clearTimeout(id);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [size.width, size.height]);

  // Frame the focused lineage's cells; release back to the whole stage. This answers a
  // store change, so it subscribes to the store rather than running as an effect.
  const frameLineage = (focus: number | null) => {
    if (!stage || !cells) return;
    if (focus === null) {
      resetView();
      return;
    }
    const focusLineage = focus;
    let r = 0;
    let cx = 0;
    let cy = 0;
    let cz = 0;
    let k = 0;
    const p = stage.positions;
    for (let i = 0; i < n; i++) {
      if (!cells.membership[focusLineage * n + i]) continue;
      cx += p[i * 3]!;
      cy += p[i * 3 + 1]!;
      cz += p[i * 3 + 2]!;
      k++;
    }
    if (!k) return;
    cx /= k;
    cy /= k;
    cz /= k;
    for (let i = 0; i < n; i++) {
      if (!cells.membership[focusLineage * n + i]) continue;
      r = Math.max(r, Math.hypot(p[i * 3]! - cx, p[i * 3 + 1]! - cy, p[i * 3 + 2]! - cz));
    }
    const base = fitFor(stage);
    const { zoom } = fitStage(Math.max(r, 0.15), size.width, size.height, insets, 3);
    fly({ ...base, target: [cx, cy, cz], zoom: Math.min(zoom, (base.zoom ?? 0) + 2.5) }, 620);
  };
  const onFocusChange = useEffectEvent(frameLineage);
  useEffect(
    () =>
      useUi.subscribe((s, prev) => {
        if (s.focusLineage !== prev.focusLineage) onFocusChange(s.focusLineage);
      }),
    [],
  );

  // Keyboard: F fits, L toggles the lasso, Esc cancels the lasso or clears the selection.
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.metaKey || e.ctrlKey || e.altKey || isTyping(e.target)) return;
      if (e.key === "f" || e.key === "F") resetView();
      else if (e.key === "l" || e.key === "L") setLasso(!useUi.getState().lassoActive);
      else if (e.key === "Escape") {
        if (useUi.getState().lassoActive) setLasso(false);
        else setSelection(null);
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  });

  // ---- what each cell shows -------------------------------------------------------------
  const fill = useMemo(
    () => (colors ? withEmphasis(colors, selection, highlightMask) : null),
    [colors, selection, highlightMask],
  );

  // [pseudotime, on-the-focused-lineage] per cell, read on the GPU by DataFilterExtension.
  const filterValues = useMemo(() => {
    const out = new Float32Array(n * 2);
    for (let i = 0; i < n; i++) {
      const pt = cells?.pseudotime[i];
      out[i * 2] = pt === undefined ? 0.5 : Number.isFinite(pt) ? pt : NO_PT;
      out[i * 2 + 1] = focusLineage === null || !cells ? 1 : cells.membership[focusLineage * n + i]! ? 1 : 0;
    }
    return out;
  }, [n, cells, focusLineage]);

  const windowed = cells !== undefined && (ptRange[0] > 0 || ptRange[1] < 1);
  const ready = stage && fill && stage.positions.length === n * 3 && filterValues.length === n * 2;
  const radius = pointRadius(n);
  const is3D = stage?.dims === 3;

  // Curves draw on from the root the first time a run is shown; they fade back after a morph.
  const drawn = useRef(new Set<string>());
  const [drawKey, setDrawKey] = useState<string | null>(null);
  useEffect(() => {
    if (runId && curves && !drawn.current.has(runId)) {
      drawn.current.add(runId);
      setDrawKey(runId);
    }
  }, [runId, curves]);
  const drawProgress = useTween(drawKey ? 1 : 0, drawKey ? DURATION.drawOn * 1000 : 0);
  const curveShow = useCurveVisibility(embedding);
  const progress = drawKey === runId ? drawProgress : 1;

  const fieldColor: RGBA = theme === "dark" ? [6, 8, 11, 255] : [230, 234, 240, 255];
  const inkColor: RGBA = theme === "dark" ? [230, 234, 240, 255] : [16, 20, 26, 255];
  // A faint rim under every cell keeps the embedding's shape when a colour runs dark
  // (viridis' low end on darkfield, its high end on brightfield). Opaque, so it never stacks.
  const rimColor: RGBA = theme === "dark" ? [30, 36, 46, 255] : [203, 209, 219, 255];

  // Cells the time-lapse window has passed: their own colour mixed most of the way to the
  // field, drawn opaque so overlapping cells cannot add up to full brightness.
  // One data object per dataset: deck.gl compares `data` by identity, and a new object
  // on every render would recompute every per-cell attribute (80 ms at 100k cells).
  const cellData = useMemo(() => ({ length: n }), [n]);

  const glowFill = useMemo(() => {
    if (!fill) return null;
    const out = new Uint8Array(fill.length);
    for (let i = 0; i < fill.length; i += 4) {
      // A light field swallows tints faster than a dark one, so brightfield keeps more colour.
      const k = (theme === "dark" ? 0.2 : 0.4) * (fill[i + 3]! / 255);
      out[i] = Math.round(fieldColor[0] + (fill[i]! - fieldColor[0]) * k);
      out[i + 1] = Math.round(fieldColor[1] + (fill[i + 1]! - fieldColor[1]) * k);
      out[i + 2] = Math.round(fieldColor[2] + (fill[i + 2]! - fieldColor[2]) * k);
      out[i + 3] = 255;
    }
    return out;
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [fill, theme]);

  const layers = [];
  if (ready) {
    const data = cellData;
    const positionOf = (_: unknown, { index, target }: { index: number; target: number[] }) => {
      target[0] = stage.positions[index * 3]!;
      target[1] = stage.positions[index * 3 + 1]!;
      target[2] = stage.positions[index * 3 + 2]!;
      return target as [number, number, number];
    };
    const colourFrom =
      (buf: Uint8Array) =>
      (_: unknown, { index, target }: { index: number; target: number[] }) => {
        target[0] = buf[index * 4]!;
        target[1] = buf[index * 4 + 1]!;
        target[2] = buf[index * 4 + 2]!;
        target[3] = buf[index * 4 + 3]!;
        return target as unknown as RGBA;
      };
    const colourOf = colourFrom(fill);
    const filterOf = (_: unknown, { index, target }: { index: number; target: number[] }) => {
      target[0] = filterValues[index * 2]!;
      target[1] = filterValues[index * 2 + 1]!;
      return target as [number, number];
    };
    const shared = {
      data,
      getPosition: positionOf,
      getFillColor: colourOf,
      getFilterValue: filterOf,
      radiusUnits: "pixels" as const,
      getRadius: radius,
      billboard: true,
      stroked: false,
      extensions: CELL_EXTENSIONS,
      depthCue: is3D ? 1 : 0,
      depthCueRadius: Math.max(stage.radius, 0.5),
      parameters: { depthCompare: "always" as const },
      updateTriggers: { getPosition: stage, getFillColor: fill, getFilterValue: filterValues },
      transitions: {
        getPosition: { duration: reduced ? 0 : DURATION.morph * 1000, easing: easeInOut },
        getFillColor: { duration: 220, easing: easeOut },
      },
    };
    // What is on the field at all: every cell, or while the time-lapse plays, the cells
    // the window has reached so far (cells ahead of it stay dark).
    const shown: [number, number] = windowed ? [0, playing ? ptRange[1] : 1] : [-3, 3];
    layers.push(
      new ScatterplotLayer({
        ...shared,
        id: "rim",
        getRadius: radius + 1.1,
        getFillColor: rimColor,
        pickable: false,
        filterRange: [shown, [0.5, 1.5]],
        updateTriggers: { ...shared.updateTriggers, getFillColor: theme },
        transitions: { getPosition: shared.transitions.getPosition },
      }),
    );
    if (windowed && glowFill) {
      // The time-lapse afterglow: cells the window has passed stay dimly lit; while paused,
      // every cell outside the window stays dimly visible for context.
      layers.push(
        new ScatterplotLayer({
          ...shared,
          id: "afterglow",
          getFillColor: colourFrom(glowFill),
          pickable: false,
          filterRange: [
            [0, playing ? ptRange[0] : 1],
            [0.5, 1.5],
          ],
          updateTriggers: { ...shared.updateTriggers, getFillColor: glowFill },
        }),
      );
    }
    layers.push(
      new ScatterplotLayer({
        ...shared,
        id: "cells",
        pickable: !lassoActive,
        filterRange: [windowed ? ptRange : [-3, 3], [0.5, 1.5]],
        filterSoftRange: windowed ? [softRange(ptRange), [0.5, 1.5]] : null,
        filterTransformColor: true,
      }),
    );
    if (playing && theme === "dark") {
      // The exposure front: cells brighten as the window reaches them, then settle.
      // Darkfield only: brightfield has no fluorescence to excite, and ink there would darken.
      const front = Math.max(ptRange[0], ptRange[1] - 0.035);
      layers.push(
        new ScatterplotLayer({
          ...shared,
          id: "exposure",
          getFillColor: [inkColor[0], inkColor[1], inkColor[2], 110],
          getRadius: radius + 0.6,
          pickable: false,
          filterRange: [
            [front, ptRange[1]],
            [0.5, 1.5],
          ],
          filterSoftRange: [
            [front, ptRange[1] - 0.012],
            [0.5, 1.5],
          ],
          filterTransformColor: true,
          updateTriggers: { ...shared.updateTriggers, getFillColor: theme },
          transitions: { getPosition: shared.transitions.getPosition },
        }),
      );
    }
    if (hovered && hovered.index < n) {
      const i = hovered.index;
      layers.push(
        new ScatterplotLayer({
          id: "hover-ring",
          data: [[stage.positions[i * 3]!, stage.positions[i * 3 + 1]!, stage.positions[i * 3 + 2]!]],
          getPosition: (d: number[]) => d as [number, number, number],
          radiusUnits: "pixels",
          getRadius: radius + 4,
          filled: false,
          stroked: true,
          lineWidthUnits: "pixels",
          getLineWidth: 1.5,
          getLineColor: inkColor,
          billboard: true,
          parameters: { depthCompare: "always" },
        }),
      );
    }
  }

  if (ready && curves && run && curveShow > 0) {
    const t = stage.transform;
    const paths = curves.curves
      .filter((c) => c.points.length > 1)
      .map((c) => ({ lineage: c.lineage, path: trimPath(c.points.map((p) => onStage(p, t)), progress) }));
    const dim = (id: number) => (focusLineage === null || focusLineage === id ? 1 : 0.22);
    const alpha = (id: number, a = 255) => Math.round(a * dim(id) * curveShow);
    layers.push(
      new PathLayer({
        id: "curve-casing",
        data: paths,
        getPath: (d) => d.path,
        getColor: (d) => [fieldColor[0], fieldColor[1], fieldColor[2], alpha(d.lineage, 210)],
        getWidth: 7,
        widthUnits: "pixels",
        jointRounded: true,
        capRounded: true,
        billboard: true,
        parameters: { depthCompare: "always" },
        updateTriggers: { getColor: [focusLineage, curveShow, theme], getPath: progress },
      }),
      new PathLayer({
        id: "curves",
        data: paths,
        getPath: (d) => d.path,
        getColor: (d) => {
          const c = categorical(d.lineage, 255);
          return [c[0], c[1], c[2], alpha(d.lineage)];
        },
        getWidth: 3,
        widthUnits: "pixels",
        jointRounded: true,
        capRounded: true,
        billboard: true,
        parameters: { depthCompare: "always" },
        updateTriggers: { getColor: [focusLineage, curveShow], getPath: progress },
      }),
    );
    const rootId = run.milestones.find((m) => m.kind === "root")?.id;
    const marks = curves.milestones.map((m) => ({
      id: m.id,
      root: m.id === rootId,
      position: onStage([m.x, m.y, ...(m.z == null ? [] : [m.z])], t),
    }));
    const root = marks.find((m) => m.root);
    if (root && progress < 1) {
      // The root announces itself once, as the curves start to draw.
      layers.push(
        new ScatterplotLayer({
          id: "root-pulse",
          data: [root],
          getPosition: (m: (typeof marks)[number]) => m.position,
          radiusUnits: "pixels",
          getRadius: 7 + progress * 22,
          filled: false,
          stroked: true,
          lineWidthUnits: "pixels",
          getLineWidth: 2,
          getLineColor: [inkColor[0], inkColor[1], inkColor[2], Math.round(255 * (1 - progress))],
          billboard: true,
          parameters: { depthCompare: "always" },
          updateTriggers: { getRadius: progress, getLineColor: progress },
        }),
      );
    }
    layers.push(
      new ScatterplotLayer({
        id: "milestones",
        data: marks.filter((m) => m.root || progress >= 1),
        getPosition: (m: (typeof marks)[number]) => m.position,
        getRadius: (m: (typeof marks)[number]) => (m.root ? 6 : 3.5),
        radiusUnits: "pixels",
        getFillColor: [inkColor[0], inkColor[1], inkColor[2], Math.round(255 * curveShow)],
        getLineColor: [fieldColor[0], fieldColor[1], fieldColor[2], Math.round(255 * curveShow)],
        lineWidthUnits: "pixels",
        getLineWidth: (m: (typeof marks)[number]) => (m.root ? 3 : 1.5),
        stroked: true,
        billboard: true,
        parameters: { depthCompare: "always" },
        updateTriggers: { getFillColor: [curveShow, theme], getLineColor: [curveShow, theme] },
      }),
      new TextLayer({
        id: "lineage-labels",
        data: paths.filter((p) => p.path.length > 0),
        getPosition: (d) => d.path[d.path.length - 1]!,
        getText: (d) => run.lineages[d.lineage]?.name ?? "",
        getSize: 12,
        getColor: (d) => [inkColor[0], inkColor[1], inkColor[2], Math.round(alpha(d.lineage) * labelFade(progress))],
        getTextAnchor: "middle",
        getAlignmentBaseline: "bottom",
        getPixelOffset: [0, -12],
        fontSettings: { sdf: true, fontSize: 64, buffer: 10, radius: 12 },
        outlineWidth: 5,
        outlineColor: [fieldColor[0], fieldColor[1], fieldColor[2], 235],
        fontFamily: "Atkinson Hyperlegible Next Variable, system-ui, sans-serif",
        fontWeight: 600,
        billboard: true,
        parameters: { depthCompare: "always" },
        updateTriggers: { getColor: [focusLineage, curveShow, progress, theme] },
      }),
    );
  }

  return (
    <div ref={containerRef} className="absolute inset-0 overflow-hidden bg-field" data-testid="scatter">
      <Fade show={Boolean(viewState && ready)}>
        {viewState && (
          <DeckGL
            views={view}
            viewState={viewState}
            onViewStateChange={({ viewState: vs }) => setViewState(vs as ViewState)}
            controller={{
              dragPan: !lassoActive,
              dragRotate: is3D && !lassoActive,
              dragMode: is3D ? "rotate" : "pan",
              scrollZoom: { smooth: true, speed: 0.012 },
              inertia: 320,
              keyboard: false,
            }}
            layers={layers}
            onHover={({ index, layer, x, y }) => {
              const next = layer?.id === "cells" && index >= 0 ? { index, x, y } : null;
              const prev = useUi.getState().hovered;
              if (next?.index !== prev?.index || (next && (next.x !== prev?.x || next.y !== prev?.y))) setHovered(next);
            }}
            getCursor={({ isDragging }) => (lassoActive ? "crosshair" : isDragging ? "grabbing" : is3D ? "grab" : "default")}
          />
        )}
      </Fade>

      {!stage && (
        <div className="absolute inset-0 grid place-items-center" style={{ paddingLeft: insets.left, paddingRight: insets.right, paddingBottom: insets.bottom }}>
          <p className="text-sm text-ink-3">
            {embedding ? "Loading the embedding…" : "The field fills once preprocessing has made an embedding."}
          </p>
        </div>
      )}

      {lassoActive && stage && (
        <LassoOverlay
          onDone={(poly) => {
            setLasso(false);
            // Project through the view itself: the DeckGL ref can hold a finalized instance in StrictMode.
            if (!viewState || poly.length < 3) return;
            const viewport = view.makeViewport({ width: size.width, height: size.height, viewState });
            if (!viewport) return;
            const screen = new Float32Array(n * 2);
            const p = stage.positions;
            for (let i = 0; i < n; i++) {
              const [sx, sy] = viewport.project([p[i * 3]!, p[i * 3 + 1]!, p[i * 3 + 2]!]);
              screen[i * 2] = sx!;
              screen[i * 2 + 1] = sy!;
            }
            setSelection(selectInPolygon(screen, poly));
          }}
        />
      )}

      {hovered && !lassoActive && <HoverCard x={hovered.x} y={hovered.y} index={hovered.index} describe={describe} pt={cells?.pseudotime[hovered.index]} />}

      <div className="pointer-events-none absolute space-y-1.5" style={{ top: insets.top + 6, left: insets.left + 6 }}>
        {stage && embedding && (
          <p
            className="truncate text-2xs text-ink-3"
            style={{ maxWidth: Math.max(120, size.width - insets.left - insets.right - (selection ? 112 : 78) - 18) }}
          >
            {embeddingLabel(embedding)}
            {stage.dims === 3
              ? " · the inference space; drag to orbit, shift-drag to pan"
              : embedding.startsWith("dc")
                ? " · the inference space"
                : " · display only; trajectories are inferred in diffusion space"}
          </p>
        )}
        {cells && <PseudotimeStamp range={ptRange} playing={playing} active={windowed || playing} />}
      </div>

      {is3D && viewState && size.width > 0 && (
        <AxisTriad
          viewport={view.makeViewport({ width: size.width, height: size.height, viewState })}
          style={{ left: insets.left + 14, bottom: insets.bottom + 18 }}
        />
      )}

      <div className="absolute flex gap-1" style={{ top: insets.top, right: insets.right }}>
        <div className="flex gap-0.5 rounded-xl border border-line bg-panel p-0.5">
          <IconButton label="Lasso select" shortcut="L" pressed={lassoActive} onClick={() => setLasso(!lassoActive)}>
            <Lasso size={16} />
          </IconButton>
          <IconButton label="Fit view" shortcut="F" onClick={resetView}>
            <Maximize2 size={15} />
          </IconButton>
          {selection && (
            <IconButton label="Clear selection" shortcut="Esc" onClick={() => setSelection(null)}>
              <X size={16} />
            </IconButton>
          )}
        </div>
      </div>

    </div>
  );
}

function labelFade(progress: number): number {
  return Math.max(0, Math.min(1, (progress - 0.75) / 0.25));
}

/** The leading part of a path, `progress` of the way along it (by point count). */
function trimPath(path: [number, number, number][], progress: number): [number, number, number][] {
  if (progress >= 1 || path.length < 2) return path;
  const exact = progress * (path.length - 1);
  const whole = Math.floor(exact);
  const out = path.slice(0, whole + 1);
  const f = exact - whole;
  const a = path[whole]!;
  const b = path[whole + 1];
  if (b && f > 0) out.push([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]);
  return out;
}

/** Hide curves while cells morph to a new embedding, then fade them back in. */
function useCurveVisibility(embedding: string | null): number {
  const [target, setTarget] = useState(1);
  const last = useRef(embedding);
  useEffect(() => {
    if (last.current === embedding) return;
    const first = last.current === null;
    last.current = embedding;
    if (first || prefersReducedMotion()) return;
    setTarget(0);
    const id = window.setTimeout(() => setTarget(1), DURATION.morph * 1000 - 80);
    return () => window.clearTimeout(id);
  }, [embedding]);
  return useTween(target, target === 0 ? 0 : 260);
}

function Fade({ show, children }: { show: boolean; children: React.ReactNode }) {
  return (
    <div
      className={cx("absolute inset-0 transition-opacity duration-500 ease-out", show ? "opacity-100" : "opacity-0")}
    >
      {children}
    </div>
  );
}

/** The acquisition stamp every time-lapse frame carries: quiet at rest, lit while the window moves. */
function PseudotimeStamp({ range, playing, active }: { range: [number, number]; playing: boolean; active: boolean }) {
  const [lo, hi] = range;
  return (
    <div
      aria-hidden
      className={cx(
        "num w-[11ch] text-[1.75rem] leading-none font-semibold tracking-[-0.03em] whitespace-nowrap transition-colors duration-300",
        active ? "text-ink" : "text-ink-3",
      )}
    >
      pt {playing ? ((lo + hi) / 2).toFixed(2) : `${lo.toFixed(2)}–${hi.toFixed(2)}`}
    </div>
  );
}

/** DC1, DC2 and DC3 drawn from the camera's point of view, so a still frame reads as 3-D. */
function AxisTriad({
  viewport,
  style,
}: {
  viewport: { project: (p: number[]) => number[] } | null;
  style: React.CSSProperties;
}) {
  if (!viewport) return null;
  const o = viewport.project([0, 0, 0]);
  const axes = [
    { label: "DC1", p: viewport.project([1, 0, 0]) },
    { label: "DC2", p: viewport.project([0, 1, 0]) },
    { label: "DC3", p: viewport.project([0, 0, 1]) },
  ].map(({ label, p }) => ({ label, dx: p[0]! - o[0]!, dy: p[1]! - o[1]!, depth: p[2]! - o[2]! }));
  const longest = Math.max(...axes.map((a) => Math.hypot(a.dx, a.dy)), 1e-6);
  const L = 30;
  return (
    <svg
      width={104}
      height={104}
      viewBox="-52 -52 104 104"
      className="pointer-events-none absolute"
      style={style}
      aria-hidden
    >
      {[...axes]
        .sort((a, b) => b.depth - a.depth)
        .map((a) => {
          const x = (a.dx / longest) * L;
          const y = (a.dy / longest) * L;
          return (
            <g key={a.label} className={a.depth < 0 ? "text-ink" : "text-ink-2"}>
              <line x1={0} y1={0} x2={x} y2={y} stroke="currentColor" strokeWidth={1.75} strokeLinecap="round" />
              <text x={x * 1.36} y={y * 1.36 + 4} textAnchor="middle" className="num fill-current text-[11px] font-semibold">
                {a.label}
              </text>
            </g>
          );
        })}
    </svg>
  );
}

function HoverCard({
  x,
  y,
  index,
  describe,
  pt,
}: {
  x: number;
  y: number;
  index: number;
  describe: (i: number) => string | null;
  pt: number | undefined;
}) {
  const value = describe(index);
  const showPt = pt !== undefined && !(value ?? "").startsWith("pseudotime");
  return (
    <div
      role="tooltip"
      className="pointer-events-none absolute z-10 rounded-lg bg-raised px-2.5 py-1.5 text-xs shadow-ply"
      style={{ left: x + 14, top: y + 14 }}
    >
      <div className="font-medium text-ink">{value ?? `cell ${index}`}</div>
      <div className="num text-2xs text-ink-3">
        cell {index}
        {showPt && (Number.isFinite(pt) ? ` · pt ${pt!.toFixed(3)}` : " · no pseudotime")}
      </div>
    </div>
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
      className="absolute inset-0 z-[5] cursor-crosshair"
      onPointerDown={(e) => {
        drawing.current = true;
        e.currentTarget.setPointerCapture(e.pointerId);
        setPoints([local(e)]);
      }}
      onPointerMove={(e) => {
        if (!drawing.current) return;
        // Read the event now: React clears currentTarget before a deferred updater runs.
        const point = local(e);
        setPoints((p) => {
          const last = p[p.length - 1];
          return last && Math.hypot(point[0] - last[0], point[1] - last[1]) < 3 ? p : [...p, point];
        });
      }}
      onPointerUp={() => {
        drawing.current = false;
        onDone(points);
        setPoints([]);
      }}
    >
      <svg className="pointer-events-none absolute inset-0 size-full">
        <path d={smoothPath(points)} className="fill-focus/10 stroke-focus" strokeWidth={1.5} strokeLinejoin="round" />
      </svg>
    </div>
  );
}

/** A closed path through the pointer's samples, smoothed with quadratic curves through midpoints. */
function smoothPath(p: [number, number][]): string {
  if (p.length < 2) return "";
  let d = `M${p[0]![0]},${p[0]![1]}`;
  for (let i = 1; i < p.length - 1; i++) {
    const [x, y] = p[i]!;
    const [nx, ny] = p[i + 1]!;
    d += ` Q${x},${y} ${(x + nx) / 2},${(y + ny) / 2}`;
  }
  const last = p[p.length - 1]!;
  return `${d} L${last[0]},${last[1]} Z`;
}

function isTyping(target: EventTarget | null): boolean {
  const el = target as HTMLElement | null;
  if (!el) return false;
  return el.isContentEditable || el.tagName === "INPUT" || el.tagName === "TEXTAREA" || el.tagName === "SELECT";
}

export { isTyping };
