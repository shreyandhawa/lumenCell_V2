import { OrbitView, type OrbitViewState } from "@deck.gl/core";
import { DataFilterExtension } from "@deck.gl/extensions";
import { PathLayer, ScatterplotLayer } from "@deck.gl/layers";
import DeckGL from "@deck.gl/react";
import { Rotate3d, RotateCcw, SlidersHorizontal } from "lucide-react";
import { animate, cubicBezier } from "motion";
import { motion, useReducedMotion } from "motion/react";
import { useEffect, useMemo, useRef, useState } from "react";
import {
  useCurves,
  useDataset,
  useDoneRunId,
  useEmbedding,
  useObs,
  useRoot,
  useRun,
  useRunCells,
  useRuns,
  useSetRoot,
} from "@/api/queries";
import { Legend } from "@/components/Legend";
import { useCellColors } from "@/components/useCellColors";
import { Gimbal } from "@/components/scatter/Gimbal";
import { useStageSize } from "@/components/scatter/useStageSize";
import { describeRoot } from "@/components/trajectory/RootPicker";
import { fmt } from "@/components/ui";
import { type Camera, FLAT, fitZoom as fitZoomFor, lerp, lerpCam, TILTED } from "@/lib/camera";
import { normalise, type Vec3 } from "@/lib/embedding";
import { DepthFadeExtension } from "@/lib/gl/depthFade";
import { useCamera } from "@/state/camera";
import { useUi } from "@/state/ui";

const LEIDEN = "lumencell_leiden";
const FILTER = new DataFilterExtension({ filterSize: 1 });
const FADE = new DepthFadeExtension();
// The stage fills the window; the data is centred in the space between the rail, the sheet and the bar.
const RAIL = 304;
const SHEET = 384;
const BAR = 80;
const EASE_OUT = [0.23, 1, 0.32, 1] as const;
const EASE_IN_OUT = [0.77, 0, 0.175, 1] as const;
const MORPH_SECONDS = 0.6;
// The GPU moves the cells on the same curve as the camera, so the root pin and curves stay on them.
const MORPH = { duration: MORPH_SECONDS * 1000, easing: cubicBezier(...EASE_IN_OUT) };
/** Auto-rotate speed, degrees per second. */
const SPIN = 12;
/** How far the farthest cells fade toward the stage in 3-D (0 = not at all, 1 = into the stage). */
const FOG = 0.6;
/** The stage colour per theme (index.css --color-stage): what fog and faint cells fade toward. */
const STAGE = { dark: [0, 0, 0], light: [244, 246, 243] } as const;
const VIEW_ID = "stage";

/** On-screen point radius at the fitted zoom, in pixels. */
function pointRadius(n: number): number {
  if (n < 5_000) return 3;
  if (n < 50_000) return 2;
  return 1.4;
}

export function Scatter() {
  const embedding = useUi((s) => s.embedding);
  const uiRunId = useUi((s) => s.runId);
  const runId = useDoneRunId(uiRunId);
  const view = useUi((s) => s.view);
  const setView = useUi((s) => s.setView);
  const theme = useUi((s) => s.theme);
  const cursor = useUi((s) => s.cursor);
  const focus = useUi((s) => s.focusLineage);
  const sheetOpen = useUi((s) => s.sheetOpen);
  const colorBy = useUi((s) => s.colorBy);
  const pointSize = useUi((s) => s.pointSize);
  const setPointSize = useUi((s) => s.setPointSize);
  const fog = useUi((s) => s.fog);
  const setFog = useUi((s) => s.setFog);
  const [optionsOpen, setOptionsOpen] = useState(false);

  const dataset = useDataset().data;
  const n = dataset?.n_obs ?? 0;
  const root = useRoot().data?.root ?? null;
  const run = useRun(runId).data;
  const runState = useRuns().data?.runs.find((r) => r.run_id === uiRunId)?.state;
  const flatRaw = useEmbedding(embedding).data;
  // DC1-DC3 come from preparation, so 3-D (and picking a root in it) opens before any run.
  const deepRaw = useEmbedding(dataset?.prep.state === "ready" ? "dc123" : null).data;
  const curves2 = useCurves(runId, embedding).data;
  const curves3 = useCurves(runId, runId ? "dc123" : null).data;
  const cells = useRunCells(runId).data;
  const { colors, legend, describe } = useCellColors();
  const setRoot = useSetRoot();
  const reduce = useReducedMotion() ?? false;

  // Before a run, a click on a cell starts the trajectory from that cell's cluster (the one shown).
  const picking = dataset?.prep.state === "ready" && !runId;
  const shownColumn = colorBy?.kind === "obs" ? dataset?.obs.find((c) => c.name === colorBy.column) : undefined;
  const clusterColumn = shownColumn && shownColumn.kind !== "continuous" ? shownColumn.name : LEIDEN;
  const clusters = useObs(picking ? clusterColumn : null).data;

  const flat = useMemo(() => (flatRaw && flatRaw.length === n * 2 ? normalise(flatRaw, 2) : null), [flatRaw, n]);
  const deep = useMemo(() => (deepRaw && deepRaw.length === n * 3 ? normalise(deepRaw, 3) : null), [deepRaw, n]);

  const [containerRef, size] = useStageSize();
  const sheetVisible = run !== undefined && sheetOpen;
  const [padTween, setPadTween] = useState(0);
  const padRight = reduce ? (sheetVisible ? SHEET : 0) : padTween;
  const padding = { left: RAIL, right: padRight, top: 0, bottom: BAR };
  const stageView = useMemo(
    () => new OrbitView({ id: VIEW_ID, orbitAxis: "Y", padding: { left: RAIL, right: padRight, top: 0, bottom: BAR } }),
    [padRight],
  );
  const fitZoom = fitZoomFor(size.width - padding.left - padding.right, size.height - padding.bottom);
  const userViewState = useCamera((s) => s.cameras[VIEW_ID] ?? null);
  const setCamera = useCamera((s) => s.setCamera);
  const setViewState = (vs: OrbitViewState) => setCamera(VIEW_ID, vs);
  // The first fit, until the user or an animation moves the camera.
  const viewState = useMemo<OrbitViewState | null>(
    () =>
      userViewState ??
      (flat && size.width ? { target: [0, 0, 0], ...FLAT, zoom: fitZoom, minZoom: fitZoom - 3, maxZoom: fitZoom + 6 } : null),
    [userViewState, flat, size.width, fitZoom],
  );
  const [t, setT] = useState(0); // 0 = 2-D embedding, 1 = DC1-DC3; in between while morphing

  // Live values for the animations below, which must not restart when these change.
  const live = useRef({ viewState, t, fitZoom, padRight });
  useEffect(() => {
    live.current = { viewState, t, fitZoom, padRight };
  });

  // The sheet opening or closing moves the stage centre with it.
  useEffect(() => {
    if (reduce) return;
    const c = animate(live.current.padRight, sheetVisible ? SHEET : 0, { duration: 0.24, ease: EASE_OUT, onUpdate: setPadTween });
    return () => c.stop();
  }, [sheetVisible, reduce]);

  const tween = (to: Camera, toT: number, seconds: number, ease: readonly number[]) => {
    const from = live.current.viewState;
    const t0 = live.current.t;
    if (!from) return () => {};
    if (reduce) {
      setT(toT);
      setViewState({ ...from, ...to });
      return () => {};
    }
    const c = animate(0, 1, {
      duration: seconds,
      ease: [...ease] as [number, number, number, number],
      onUpdate: (p) => {
        setT(lerp(t0, toT, p));
        setViewState({ ...from, ...lerpCam(from as Camera, to, p) });
      },
    });
    return () => c.stop();
  };

  // The 2-D <-> 3-D morph: positions and camera move together.
  useEffect(() => {
    const toT = view === "3d" ? 1 : 0;
    if (live.current.t === toT || (toT === 1 && !deep)) return;
    const cam = { target: [0, 0, 0] as Vec3, zoom: live.current.fitZoom, ...(toT ? TILTED : FLAT) };
    return tween(cam, toT, MORPH_SECONDS, EASE_IN_OUT);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- tween reads live values; only these start a morph
  }, [view, deep]);

  // Auto-rotate, in 3-D only; any drag, key or reset stops it.
  const [spinning, setSpinning] = useState(false);
  const spin = spinning && t === 1 && !reduce;
  useEffect(() => {
    if (!spin) return;
    let last = performance.now();
    let id = requestAnimationFrame(function step(now) {
      const vs = live.current.viewState;
      if (vs) setCamera(VIEW_ID, { ...vs, rotationOrbit: (vs.rotationOrbit ?? 0) + ((now - last) / 1000) * SPIN });
      last = now;
      id = requestAnimationFrame(step);
    });
    return () => cancelAnimationFrame(id);
  }, [spin, setCamera]);

  // deck.gl makes its root focusable for the keyboard controls once it has loaded; name it and say
  // which keys work.
  const labelPlot = () => {
    const root = containerRef.current?.querySelector<HTMLElement>('[tabindex="0"]');
    root?.setAttribute("role", "application");
    root?.setAttribute("aria-label", "Plot. Arrow keys pan, minus and equals zoom, Shift with arrows orbits in 3-D, R resets.");
  };

  const reset = () => {
    setSpinning(false);
    tween({ target: [0, 0, 0], zoom: fitZoom, ...(view === "3d" ? TILTED : FLAT) }, t, 0.3, EASE_OUT);
  };

  // The set piece: when a run finishes, pseudotime washes out from the root, once. It holds at 0
  // from the moment the run is done, so the finished colours never flash before the wash starts.
  const [wash, setWash] = useState(1);
  const [washing, setWashing] = useState(false);
  const [seenState, setSeenState] = useState(runState);
  if (runState !== seenState) {
    setSeenState(runState);
    const finished = runState === "done" && (seenState === "queued" || seenState === "running") && !reduce;
    setWashing(finished);
    setWash(finished ? 0 : 1);
  }
  useEffect(() => {
    if (!cells || !washing) return;
    const c = animate(0, 1, { duration: 0.9, ease: EASE_OUT, onUpdate: setWash, onComplete: () => setWashing(false) });
    return () => c.stop();
  }, [cells, washing]);

  // One attribute object per embedding: deck.gl starts the GPU morph when the object changes, so a
  // recolour (which rebuilds the layer data) never restarts it.
  const flatAttr = useMemo(() => (flat ? { value: flat.xyz, size: 3 } : null), [flat]);
  const deepAttr = useMemo(() => (deep ? { value: deep.xyz, size: 3 } : null), [deep]);
  const positions = view === "3d" && deepAttr ? deepAttr : flatAttr;
  /** Where cell i is now, mid-morph included: the GPU's interpolation, repeated for one cell (the root pin). */
  const at = (i: number): Vec3 | null => {
    if (!flat || i < 0 || i >= n) return null;
    const p = [0, 1, 2].map((k) => (deep ? lerp(flat.xyz[i * 3 + k]!, deep.xyz[i * 3 + k]!, t) : flat.xyz[i * 3 + k]!));
    return p.every(Number.isFinite) ? (p as Vec3) : null;
  };

  const viewport = useMemo(
    () => (viewState && size.width ? stageView.makeViewport({ width: size.width, height: size.height, viewState }) : null),
    [stageView, viewState, size.width, size.height],
  );

  const ptValues = useMemo(() => {
    const out = new Float32Array(n).fill(-1);
    if (cells) for (let i = 0; i < n; i++) if (Number.isFinite(cells.pseudotime[i]!)) out[i] = cells.pseudotime[i]!;
    return out;
  }, [n, cells]);
  // Cells past the cursor (or the wash) keep only their halo.
  const limit = cells ? (cursor ?? (washing ? wash : 2)) : 2;

  const paths = useMemo(() => {
    const shape = (cs: typeof curves2, norm: typeof flat) =>
      cs && norm ? cs.curves.filter((c) => c.points.length > 1).map((c) => ({ lineage: c.lineage, path: c.points.map(norm.map) })) : null;
    const a = shape(curves2, flat);
    const b = shape(curves3, deep);
    if (t === 0 || !b) return a ?? [];
    if (t === 1 || !a) return b;
    return a.map((c, i) => {
      const d = b[i];
      if (!d || d.path.length !== c.path.length) return t < 0.5 ? c : (d ?? c);
      return { lineage: c.lineage, path: c.path.map((p, k) => [0, 1, 2].map((j) => lerp(p[j]!, d.path[k]![j]!, t)) as Vec3) };
    });
  }, [curves2, curves3, flat, deep, t]);

  const haloData = useMemo(() => ({ length: n, attributes: { getPosition: positions! } }), [n, positions]);
  const cellData = useMemo(
    () => ({
      length: n,
      attributes: {
        getPosition: positions!,
        getFillColor: { value: colors!, size: 4, normalized: true },
        getFilterValue: { value: ptValues, size: 1 },
      },
    }),
    [n, positions, colors, ptValues],
  );
  const ink: [number, number, number, number] = theme === "dark" ? [243, 245, 242, 190] : [20, 23, 22, 170];
  const layers = [];
  if (positions && colors) {
    // Cells are balls of one size in embedding units (near ones look bigger in 3-D), sized so they
    // are r pixels at the fitted zoom, and never smaller than 0.6 r or larger than 3 r before perspective.
    const r = pointRadius(n) * pointSize;
    const radius = r / 2 ** fitZoom;
    // The fog spans the box that holds the cells, which moves with them while morphing.
    const box = (k: 0 | 1) =>
      [0, 1, 2].map((j) => (flat && deep ? lerp(flat.bounds[k][j]!, deep.bounds[k][j]!, t) : (flat?.bounds[k][j] ?? 0))) as Vec3;
    const transitions = reduce ? undefined : { getPosition: MORPH };
    layers.push(
      // A faint ghost of every cell, so cells past the cursor still show where they are. Drawn first,
      // without depth, so cells always cover it.
      new ScatterplotLayer({
        id: "halo",
        data: haloData,
        getFillColor: theme === "dark" ? [255, 255, 255, 10] : [0, 0, 0, 14],
        radiusUnits: "pixels",
        getRadius: r + 1.5,
        billboard: true,
        parameters: { depthCompare: "always", depthWriteEnabled: false },
        transitions,
        updateTriggers: { getFillColor: theme, getRadius: r },
      }),
      new ScatterplotLayer({
        id: "cells",
        data: cellData,
        radiusUnits: "common",
        getRadius: radius,
        radiusMinPixels: 0.6 * r,
        radiusMaxPixels: 3 * r,
        billboard: true,
        pickable: true,
        autoHighlight: true,
        highlightColor: theme === "dark" ? [255, 255, 255, 220] : [0, 0, 0, 200],
        extensions: [FILTER, FADE],
        filterRange: [-2, limit],
        filterSoftRange: [-2, limit - 0.06],
        depthFadeStage: STAGE[theme],
        depthFadeStrength: fog ? FOG * t : 0,
        depthFadeBounds: [box(0), box(1)],
        transitions,
        updateTriggers: { getRadius: radius },
      }),
    );
  }
  if (paths.length) {
    const ordered = [...paths.filter((p) => p.lineage !== focus), ...paths.filter((p) => p.lineage === focus)];
    layers.push(
      new PathLayer({
        id: "curves",
        data: ordered,
        getPath: (d) => d.path,
        getColor: (d) => (d.lineage === focus ? [224, 51, 143, 255] : ink),
        getWidth: (d) => (d.lineage === focus ? 3 : 1.5),
        widthUnits: "pixels",
        billboard: true,
        jointRounded: true,
        capRounded: true,
        parameters: { depthCompare: "always" },
        updateTriggers: { getColor: [focus, theme], getWidth: focus },
      }),
    );
  }

  // The root pin follows the root cell through pans, zooms and the morph.
  const rootXyz = root ? at(root.cell_index) : null;
  const pin = viewport && rootXyz ? viewport.project(rootXyz) : null;

  const embeddingInfo = dataset?.embeddings.find((e) => e.name === embedding);
  const p = run?.provenance;
  const caption =
    view === "3d"
      ? `DC1–DC3 · ${p?.inference_rep === "X_diffmap" ? `3 of ${p.inference_dims}` : "the first 3"} diffusion components${fog ? " · far cells fade" : ""}`
      : !embedding
        ? ""
        : paths.length && p
          ? `${embedding} · curves are projections of the tree inferred in ${p.inference_rep}`
          : legend
            ? `${embedding} · coloured by ${legend.title === "Leiden" ? "Leiden cluster" : legend.title}`
            : `${embedding}, ${embeddingInfo?.source === "file" ? "stored in the file" : "computed by lumenCell"} · ${fmt(n)} cells`;
  const right = Math.max(12, padRight);
  const raised = cursor !== null && run !== undefined; // the gene trends sit above the bar

  return (
    <div
      ref={containerRef}
      className="absolute inset-0 overflow-hidden bg-stage"
      data-testid="scatter"
      onKeyDown={(e) => {
        // R resets the view, from the plot or its controls (never while typing in a field).
        const field = e.target instanceof HTMLInputElement || e.target instanceof HTMLSelectElement;
        if (e.key.toLowerCase() === "r" && !field && !e.metaKey && !e.ctrlKey && !e.altKey && viewState) reset();
      }}
    >
      {viewState && (
        <DeckGL
          views={stageView}
          viewState={viewState}
          onLoad={labelPlot}
          onViewStateChange={({ viewState: vs }) => {
            if (t !== 0 && t !== 1) return;
            setSpinning(false);
            setViewState(vs as OrbitViewState);
          }}
          // Arrows pan, -/= zoom; with Shift they orbit, in 3-D only.
          controller={
            t === 0
              ? { dragMode: "pan", dragRotate: false, keyboard: { rotateSpeedX: 0, rotateSpeedY: 0 } }
              : t === 1
                ? { dragMode: "rotate" }
                : false
          }
          layers={layers}
          onClick={(info) => {
            if (!picking || info.layer?.id !== "cells" || info.index < 0 || !clusters?.codes) return;
            const code = clusters.codes[info.index]!;
            if (code >= 0) setRoot.mutate({ kind: "cluster", key: clusterColumn, value: clusters.categories[code]! });
          }}
          getTooltip={({ index, layer }) =>
            layer?.id === "cells" && index >= 0
              ? {
                  text: describe(index) ?? `cell ${index}`,
                  style: {
                    background: "var(--color-panel-raised)",
                    color: "var(--color-ink)",
                    fontFamily: "var(--font-sans)",
                    fontSize: "12px",
                    borderRadius: "4px",
                    padding: "4px 6px",
                  },
                }
              : null
          }
          getCursor={({ isDragging, isHovering }) =>
            isDragging ? "grabbing" : isHovering && picking ? "crosshair" : t === 1 ? "grab" : "default"
          }
        />
      )}

      {!flat && (
        <div className="absolute inset-0 grid place-items-center pl-[304px] text-[13px] text-ink-muted">
          {embedding ? "Loading the embedding…" : "No embedding yet: preparation will compute one."}
        </div>
      )}

      {root && pin && (
        <div className="pointer-events-none absolute" style={{ left: pin[0], top: pin[1] }}>
          <motion.span
            key={root.obs_name}
            initial={{ opacity: 0, scale: 0.9 }}
            animate={{ opacity: 1, scale: 1 }}
            transition={{ duration: 0.2, ease: EASE_OUT }}
            className="absolute -top-2 -left-2 size-4 rounded-full border-2 border-stage bg-accent"
          />
          <span className="absolute -top-[11px] left-4 rounded bg-panel-raised px-1.5 py-0.5 text-xs whitespace-nowrap text-ink">
            Root · {describeRoot(root.spec)}
          </span>
        </div>
      )}

      <p className="pointer-events-none absolute top-[22px] left-[308px] text-xs text-ink-muted">{caption}</p>

      <div
        role="group"
        aria-label="View"
        className="absolute top-3 z-10 flex gap-0.5 rounded-lg bg-panel p-1"
        // clear of the "Show evidence" button when the sheet is hidden
        style={{ right: Math.max(padRight + 6, run ? 126 : 18) }}
      >
        {(["2d", "3d"] as const).map((v) => (
          <button
            key={v}
            type="button"
            aria-pressed={view === v}
            disabled={v === "3d" && !deep}
            title={v === "3d" && !deep ? "3-D opens once preparation has finished" : undefined}
            onClick={() => {
              setSpinning(false);
              setView(v);
            }}
            className={`press h-[26px] rounded-md px-2.5 text-[13px] disabled:text-ink-faint ${
              view === v ? "bg-panel-raised text-ink" : "text-ink-muted hover:text-ink"
            }`}
          >
            {v === "2d" ? "2-D" : "3-D"}
          </button>
        ))}
      </div>

      {legend && colorBy?.kind !== "pseudotime" && (
        <div className="absolute left-[304px] z-10" style={{ bottom: raised ? 316 : 94 }}>
          <Legend spec={legend} />
        </div>
      )}

      {viewState && (
        <div className="absolute z-10 flex items-end gap-3" style={{ right, bottom: raised ? 316 : 120 }}>
          <div className="relative flex flex-col gap-1.5">
            {optionsOpen && (
              <div
                id="display-options"
                role="group"
                aria-label="Display options"
                className="absolute right-0 bottom-full mb-2 w-56 rounded-xl bg-panel p-3 text-[13px] text-ink"
              >
                <label className="flex flex-col gap-1.5">
                  <span className="flex justify-between">
                    Point size <span className="font-mono text-xs text-ink-muted">{pointSize.toFixed(2)}×</span>
                  </span>
                  <input
                    type="range"
                    min={0.5}
                    max={3}
                    step={0.25}
                    value={pointSize}
                    onChange={(e) => setPointSize(Number(e.target.value))}
                    className="accent-ink"
                  />
                </label>
                <label className="mt-3 flex items-center gap-2">
                  <input type="checkbox" checked={fog} onChange={(e) => setFog(e.target.checked)} className="accent-ink" />
                  Fade far cells in depth
                </label>
              </div>
            )}
            <button
              type="button"
              aria-label="Display options"
              aria-expanded={optionsOpen}
              aria-controls="display-options"
              onClick={() => setOptionsOpen((o) => !o)}
              className={`press grid size-8 place-items-center rounded-md text-ink ${optionsOpen ? "bg-panel-raised" : "bg-panel"}`}
            >
              <SlidersHorizontal size={16} />
            </button>
            {t === 1 && (
              <button
                type="button"
                aria-label="Auto-rotate"
                aria-pressed={spin}
                disabled={reduce}
                title={reduce ? "Off while your system asks for reduced motion" : undefined}
                onClick={() => setSpinning((s) => !s)}
                className={`press grid size-8 place-items-center rounded-md text-ink disabled:text-ink-faint ${spin ? "bg-panel-raised" : "bg-panel"}`}
              >
                <Rotate3d size={16} />
              </button>
            )}
            <button
              type="button"
              aria-label="Reset view"
              title="Reset view (R)"
              onClick={reset}
              className="press grid size-8 place-items-center rounded-md bg-panel text-ink"
            >
              <RotateCcw size={16} />
            </button>
          </div>
          {t > 0 && viewport && <Gimbal m={viewport.viewMatrix as unknown as number[]} />}
        </div>
      )}
    </div>
  );
}
