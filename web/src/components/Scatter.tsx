import { type Layer, OrbitView, type OrbitViewState, type Viewport } from "@deck.gl/core";
import { DataFilterExtension } from "@deck.gl/extensions";
import { ScatterplotLayer } from "@deck.gl/layers";
import DeckGL, { type DeckGLRef } from "@deck.gl/react";
import { ImageDown, Lasso, Rotate3d, RotateCcw, SlidersHorizontal } from "lucide-react";
import { animate, cubicBezier } from "motion";
import { motion, useReducedMotion } from "motion/react";
import { useEffect, useMemo, useRef, useState } from "react";
import {
  useCurves,
  useDataset,
  useDoneRunId,
  useEmbedding,
  useKnnGraph,
  useMethods,
  useObs,
  usePagaGraph,
  useRoot,
  useRun,
  useRunCells,
  useRuns,
  useSetRoot,
} from "@/api/queries";
import { Legend } from "@/components/Legend";
import { useCellColors } from "@/components/useCellColors";
import { CameraTour } from "@/components/scatter/CameraTour";
import { Gimbal } from "@/components/scatter/Gimbal";
import { type Ends, knnLayer, pagaLayers, type PlacedNode } from "@/components/scatter/graphLayers";
import { type GroupBy, SelectionBar } from "@/components/scatter/SelectionBar";
import { trajectoryLayers } from "@/components/scatter/trajectoryLayers";
import { useStageSize } from "@/components/scatter/useStageSize";
import { describeRoot } from "@/components/trajectory/RootPicker";
import { fmt, Select } from "@/components/ui";
import { axesLabel, axesName, setAxis } from "@/lib/axes";
import { type Camera, clampTarget, FLAT, fitZoom as fitZoomFor, lerp, lerpCam, TILTED } from "@/lib/camera";
import { normalise, type Vec3 } from "@/lib/embedding";
import { type Point, selectInPolygon } from "@/lib/geometry";
import { exportPng } from "@/lib/exportPng";
import { DepthFadeExtension } from "@/lib/gl/depthFade";
import { edgeEnds } from "@/lib/graph";
import { combine, projectAll } from "@/lib/selection";
import { centroid, cursorPoint, type PlacedMilestone, placeMilestones, prefixTo } from "@/lib/trajectory";
import { useCamera } from "@/state/camera";
import { useSelection } from "@/state/selection";
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
/** The split view: the 2-D embedding on the left, DC1-DC3 on the right, each with its own camera. */
const PANES = ["flat", "deep"] as const;
type Pane = (typeof PANES)[number];
type Positions = { value: Float32Array; size: number };
const ON_TOP = { depthCompare: "always", depthWriteEnabled: false } as const;
/** A layer id without its pane: "deep:cells" -> "cells". */
const kind = (id: string | undefined) => id?.slice(id.indexOf(":") + 1);
/** 2-D: drag and arrows pan, never rotate. */
const flatControl = (keys: boolean) => ({
  dragMode: "pan" as const,
  dragRotate: false,
  keyboard: keys && { rotateSpeedX: 0, rotateSpeedY: 0 },
});
/** 3-D: drag orbits. With the lasso on, the left drag draws, so the right drag (deck's alternate mode) orbits. */
const deepControl = (keys: boolean, lasso: boolean) => ({ dragMode: lasso ? ("pan" as const) : ("rotate" as const), keyboard: keys });

/** ``fn``, kept per positions object: deck.gl uploads every attribute again when the data object changes. */
function perPositions<T>(fn: (p: Positions) => T): (p: Positions) => T {
  const seen = new WeakMap<Positions, T>();
  return (p) => {
    let v = seen.get(p);
    if (v === undefined) seen.set(p, (v = fn(p)));
    return v;
  };
}

/**
 * How far past the fitted view the zoom may go. Deeper, the view sits in the empty space between
 * cells and draws nothing; the target is also kept over the cells (see onViewStateChange).
 */
const MAX_ZOOM_PAST_FIT = 4;

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
  const axes = useUi((s) => s.axes);
  const setAxes = useUi((s) => s.setAxes);
  const overlays = useUi((s) => s.overlays);
  const setOverlay = useUi((s) => s.setOverlay);
  const [optionsOpen, setOptionsOpen] = useState(false);

  const dataset = useDataset().data;
  const n = dataset?.n_obs ?? 0;
  const root = useRoot().data?.root ?? null;
  const run = useRun(runId).data;
  const methods = useMethods().data?.methods;
  const runState = useRuns().data?.runs.find((r) => r.run_id === uiRunId)?.state;
  const flatRaw = useEmbedding(embedding).data;
  // DC1-DC3 come from preparation, so 3-D (and picking a root in it) opens before any run.
  // The 3-D view shows the chosen three diffusion components; while new ones load, the old stay up.
  const deepName = axesName(axes);
  const deepQuery = useEmbedding(dataset?.prep.state === "ready" ? deepName : null);
  const deepRaw = deepQuery.data;
  const deepStale = deepQuery.isPlaceholderData; // drawn, but not the components asked for yet
  const components = dataset?.embeddings.find((e) => e.name === "dc123")?.components ?? null;
  const curves2 = useCurves(runId, embedding).data;
  const curves3 = useCurves(runId, runId ? deepName : null).data;
  const pagaFlat = usePagaGraph(overlays.paga ? embedding : null).data;
  const pagaDeep = usePagaGraph(overlays.paga && deepRaw ? deepName : null).data;
  const knnData = useKnnGraph(overlays.knn).data;
  const knn = overlays.knn ? knnData : undefined; // a disabled query keeps its cached edges
  const cells = useRunCells(runId).data;
  const { colors, legend, describe } = useCellColors();
  const setRoot = useSetRoot();
  const reduce = useReducedMotion() ?? false;

  // Before a run, a click on a cell starts the trajectory from that cell's cluster (the one shown).
  const picking = dataset?.prep.state === "ready" && !runId;
  const shownColumn = colorBy?.kind === "obs" ? dataset?.obs.find((c) => c.name === colorBy.column) : undefined;
  const clusterColumn = shownColumn && shownColumn.kind !== "continuous" ? shownColumn.name : LEIDEN;
  // The lasso tool: a left drag draws, other buttons and keys still move the view.
  const [lasso, setLasso] = useState(false);
  const [drawing, setDrawing] = useState<Point[] | null>(null);
  // Side by side, the hovered cell is ringed in both panes.
  const [hovered, setHovered] = useState<number | null>(null);
  const { mask, op, setMask } = useSelection();
  const deckRef = useRef<DeckGLRef>(null);
  const clusters = useObs(picking || lasso || mask ? clusterColumn : null).data;

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
  // The first fit, until the user or an animation moves the camera.
  const viewState = useMemo<OrbitViewState | null>(
    () =>
      userViewState ??
      (flat && size.width ? { target: [0, 0, 0], ...FLAT, zoom: fitZoom, minZoom: fitZoom - 3, maxZoom: fitZoom + MAX_ZOOM_PAST_FIT } : null),
    [userViewState, flat, size.width, fitZoom],
  );
  const [t, setT] = useState(0); // 0 = 2-D embedding, 1 = DC1-DC3; in between while morphing
  // The box that holds the cells at this point of the morph, the same blend the fog uses.
  const morphBox = (flat && deep
    ? [0, 1].map((k) => [0, 1, 2].map((j) => lerp(flat.bounds[k as 0 | 1][j]!, deep.bounds[k as 0 | 1][j]!, t))) as [Vec3, Vec3]
    : flat?.bounds) as [Vec3, Vec3] | undefined;

  // Side by side: two panes share the free stage, each fitted and moved on its own.
  const split = view === "split" && flat !== null && deep !== null;
  const half = Math.max(0, Math.floor((size.width - RAIL - padRight) / 2));
  const paneFit = fitZoomFor(half, size.height - BAR);
  const cameras = useCamera((s) => s.cameras);
  const paneStates = useMemo<Record<Pane, OrbitViewState> | null>(() => {
    if (!split) return null;
    const fit = (cam: typeof FLAT): OrbitViewState => ({
      target: [0, 0, 0],
      ...cam,
      zoom: paneFit,
      minZoom: paneFit - 3,
      maxZoom: paneFit + MAX_ZOOM_PAST_FIT,
    });
    return { flat: cameras.flat ?? fit(FLAT), deep: cameras.deep ?? fit(TILTED) };
  }, [split, cameras.flat, cameras.deep, paneFit]);
  // deck.gl gives every key to every view, so only the pane the pointer was over last has keys.
  const [active, setActive] = useState<Pane>("deep");
  const paneAt = (x: number): Pane => (x < RAIL + half ? "flat" : "deep");
  const paneViews = useMemo(
    () =>
      split
        ? PANES.map(
            (id, k) =>
              new OrbitView({
                id,
                orbitAxis: "Y",
                x: RAIL + k * half,
                width: half,
                padding: { bottom: BAR },
                controller: id === "flat" ? flatControl(active === id) : deepControl(active === id, lasso),
              }),
          )
        : null,
    [split, half, active, lasso],
  );
  const paneViewports = useMemo(
    () =>
      paneViews && paneStates && size.width
        ? paneViews.map((v) => v.makeViewport({ width: size.width, height: size.height, viewState: paneStates[v.id as Pane] }))
        : null,
    [paneViews, paneStates, size.width, size.height],
  );

  // Live values for the animations below, which must not restart when these change.
  const live = useRef({ viewState, t, fitZoom, padRight, paneStates });
  useEffect(() => {
    live.current = { viewState, t, fitZoom, padRight, paneStates };
  });

  // The sheet opening or closing moves the stage centre with it.
  useEffect(() => {
    if (reduce) return;
    const c = animate(live.current.padRight, sheetVisible ? SHEET : 0, { duration: 0.24, ease: EASE_OUT, onUpdate: setPadTween });
    return () => c.stop();
  }, [sheetVisible, reduce]);

  /** Move a camera (the stage's unless a pane is named) and the morph toward toT together. */
  const tween = (
    to: Camera,
    toT: number,
    seconds: number,
    ease: readonly number[],
    id: string = VIEW_ID,
    from: OrbitViewState | null = live.current.viewState,
  ) => {
    const t0 = live.current.t;
    if (!from) return () => {};
    if (reduce) {
      setT(toT);
      setCamera(id, { ...from, ...to });
      return () => {};
    }
    const c = animate(0, 1, {
      duration: seconds,
      ease: [...ease] as [number, number, number, number],
      onUpdate: (p) => {
        setT(lerp(t0, toT, p));
        setCamera(id, { ...from, ...lerpCam(from as Camera, to, p) });
      },
    });
    return () => c.stop();
  };

  // The 2-D <-> 3-D morph: positions and camera move together. It waits for the stage's first
  // camera and, into 3-D, for the diffusion components; it runs again when they arrive.
  const staged = viewState !== null;
  useEffect(() => {
    if (view === "split") return; // the panes do not morph; the stage morphs when it is shown again
    const toT = view === "3d" ? 1 : 0;
    if (live.current.t === toT || (toT === 1 && !deep) || !staged) return;
    // A shared link's 3-D camera, once; otherwise the default tilt.
    const { arrival, setArrival } = useCamera.getState();
    if (arrival && toT === 1) setArrival(null);
    const cam = (toT === 1 && arrival) || { target: [0, 0, 0] as Vec3, zoom: live.current.fitZoom, ...(toT ? TILTED : FLAT) };
    return tween(cam, toT, MORPH_SECONDS, EASE_IN_OUT);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- tween reads live values; only these start a morph
  }, [view, deep, staged]);

  // Auto-rotate, in 3-D only (the right pane side by side); any drag, key or reset stops it.
  const [spinning, setSpinning] = useState(false);
  const spin = spinning && (split || t === 1) && !reduce;
  useEffect(() => {
    if (!spin) return;
    let last = performance.now();
    let id = requestAnimationFrame(function step(now) {
      const panes = live.current.paneStates;
      const vs = panes ? panes.deep : live.current.viewState;
      if (vs) setCamera(panes ? "deep" : VIEW_ID, { ...vs, rotationOrbit: (vs.rotationOrbit ?? 0) + ((now - last) / 1000) * SPIN });
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
    root?.setAttribute(
      "aria-label",
      "Plot. Arrow keys pan, minus and equals zoom, Shift with arrows orbits in 3-D, R resets. Side by side, keys move the side the pointer was over last.",
    );
  };

  const reset = () => {
    setSpinning(false);
    if (paneStates) {
      for (const id of PANES)
        tween({ target: [0, 0, 0], zoom: paneFit, ...(id === "deep" ? TILTED : FLAT) }, t, 0.3, EASE_OUT, id, paneStates[id]);
      return;
    }
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
  const at = (i: number, tt = t): Vec3 | null => {
    if (!flat || i < 0 || i >= n) return null;
    const p = [0, 1, 2].map((k) => (deep ? lerp(flat.xyz[i * 3 + k]!, deep.xyz[i * 3 + k]!, tt) : flat.xyz[i * 3 + k]!));
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

  const shapes = useMemo(() => {
    const shape = (cs: typeof curves2, norm: typeof flat) =>
      cs && norm ? cs.curves.filter((c) => c.points.length > 1).map((c) => ({ lineage: c.lineage, path: c.points.map(norm.map) })) : null;
    return { a: shape(curves2, flat), b: deepStale ? null : shape(curves3, deep) };
  }, [curves2, curves3, flat, deep, deepStale]);
  /** The curves at tt: 0 on the 2-D embedding, 1 on DC1-DC3, in between while morphing. */
  const pathsAt = (tt: number) => {
    const { a, b } = shapes;
    if (tt === 0 || !b) return tt === 1 ? [] : (a ?? []); // never the 2-D curves on the 3-D view
    if (tt === 1 || !a) return b;
    return a.map((c, i) => {
      const d = b[i];
      if (!d || d.path.length !== c.path.length) return tt < 0.5 ? c : (d ?? c);
      return { lineage: c.lineage, path: c.path.map((p, k) => [0, 1, 2].map((j) => lerp(p[j]!, d.path[k]![j]!, tt)) as Vec3) };
    });
  };
  const paths = pathsAt(t);

  const mixAt = (tt: number) => (a: Vec3, b: Vec3) => a.map((v, k) => lerp(v, b[k]!, tt)) as Vec3;
  /** Where each item (a milestone, a PAGA cluster) is in one view. */
  const where = (items: readonly { id: number; x: number; y: number; z?: number | null }[] | undefined, norm: typeof flat) =>
    items && norm ? new Map(items.map((m) => [m.id, norm.map([m.x, m.y, m.z ?? 0])] as const)) : null;
  const milestonesAt = (tt: number) =>
    run
      ? placeMilestones(
          run.milestones,
          tt < 1 ? where(curves2?.milestones, flat) : null,
          tt > 0 && !deepStale ? where(curves3?.milestones, deep) : null,
          mixAt(tt),
        )
      : [];
  /** PAGA's clusters and edges at tt, like the milestones. */
  const pagaAt = (tt: number) => {
    const g = pagaFlat ?? pagaDeep;
    if (!g) return null;
    const a = tt < 1 ? where(pagaFlat?.nodes, flat) : null;
    const b = tt > 0 && !deepStale ? where(pagaDeep?.nodes, deep) : null;
    return { nodes: placeMilestones(g.nodes, a, b, mixAt(tt)), edges: g.edges };
  };

  // The playhead: where the focused lineage's cells are at the cursor, in each view.
  const heads = useMemo(() => {
    if (cursor === null || focus === null || !cells || !flat || focus >= cells.nLineages) return null;
    const member = cells.membership.subarray(focus * n, (focus + 1) * n);
    const a = cursorPoint(cells.pseudotime, member, flat.xyz, cursor);
    const b = deep ? cursorPoint(cells.pseudotime, member, deep.xyz, cursor) : null;
    return { a, b };
  }, [cursor, focus, cells, flat, deep, n]);
  /** The focused lineage's curve up to where its cells are at the cursor, at tt. */
  const travelledAt = (tt: number) => {
    const head = heads && (heads.a && heads.b ? mixAt(tt)(heads.a, heads.b) : tt < 0.5 ? heads.a : heads.b);
    const path = pathsAt(tt).find((p) => p.lineage === focus)?.path;
    return head && path ? prefixTo(path, head) : null;
  };
  const deepFocus = shapes.b?.find((p) => p.lineage === focus)?.path;

  // Focusing a lineage in 3-D (the right pane side by side) turns the camera around it.
  useEffect(() => {
    const centre = deepFocus && centroid(deepFocus);
    const pane = paneStates?.deep;
    const vs = pane ?? live.current.viewState;
    if ((!pane && live.current.t !== 1) || !vs || !centre) return;
    const to = { target: centre, zoom: vs.zoom, rotationX: vs.rotationX, rotationOrbit: vs.rotationOrbit };
    return pane ? tween(to, live.current.t, 0.45, EASE_OUT, "deep", pane) : tween(to, 1, 0.45, EASE_OUT);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- only a focus change moves the camera
  }, [focus]);

  const haloOf = useMemo(() => perPositions((pos) => ({ length: n, attributes: { getPosition: pos } })), [n]);
  const selected = useMemo(() => (mask ? Float32Array.from(mask) : null), [mask]);
  const cellsOf = useMemo(
    () =>
      perPositions((pos) => ({
        length: n,
        attributes: {
          getPosition: pos,
          getFillColor: { value: colors!, size: 4, normalized: true },
          getFilterValue: { value: ptValues, size: 1 },
          ...(selected ? { getSelected: { value: selected, size: 1 } } : {}),
        },
      })),
    [n, colors, ptValues, selected],
  );
  // The neighbour graph's line ends, per view, from the same positions as the cells.
  const endsOf = useMemo(
    () =>
      perPositions((pos): Ends => {
        const { source, target } = edgeEnds(knn!.src, knn!.dst, pos.value);
        return {
          length: knn!.src.length,
          attributes: { getSourcePosition: { value: source, size: 3 }, getTargetPosition: { value: target, size: 3 } },
        };
      }),
    [knn],
  );
  const ink: [number, number, number, number] = theme === "dark" ? [243, 245, 242, 190] : [20, 23, 22, 170];
  // Cells are balls of one size in embedding units (near ones look bigger in 3-D), sized so they
  // are r pixels at the fitted zoom, and never smaller than 0.6 r or larger than 3 r before perspective.
  const r = pointRadius(n) * pointSize;
  /** Everything drawn in one view, at tt between the 2-D embedding (0) and DC1-DC3 (1). */
  const scene = (prefix: string, tt: number, pos: Positions | null, zoom: number): Layer[] => {
    if (!pos || !colors) return [];
    const radius = r / 2 ** zoom;
    // The fog spans the box that holds the cells, which moves with them while morphing.
    const box = (k: 0 | 1) =>
      [0, 1, 2].map((j) => (flat && deep ? lerp(flat.bounds[k][j]!, deep.bounds[k][j]!, tt) : (flat?.bounds[k][j] ?? 0))) as Vec3;
    const transitions = reduce ? undefined : { getPosition: MORPH };
    const paga = pagaAt(tt);
    const out: Layer[] = [
      // A faint ghost of every cell, so cells past the cursor still show where they are. Drawn first,
      // without depth, so cells always cover it. Off while a selection shrinks the unselected cells,
      // where it would ring each of them.
      new ScatterplotLayer({
        id: `${prefix}halo`,
        visible: mask === null,
        data: haloOf(pos),
        getFillColor: theme === "dark" ? [255, 255, 255, 10] : [0, 0, 0, 14],
        radiusUnits: "pixels",
        getRadius: r + 1.5,
        billboard: true,
        parameters: { depthCompare: "always", depthWriteEnabled: false },
        transitions,
        updateTriggers: { getFillColor: theme, getRadius: r },
      }),
      ...(knn ? [knnLayer(prefix, endsOf(pos), ink, reduce ? undefined : MORPH)] : []),
      new ScatterplotLayer({
        id: `${prefix}cells`,
        data: cellsOf(pos),
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
        depthFadeStrength: fog ? FOG * tt : 0,
        depthFadeBounds: [box(0), box(1)],
        depthFadeDim: mask !== null,
        transitions,
        updateTriggers: { getRadius: radius, getSelected: mask },
      }),
      ...(paga ? pagaLayers(prefix, paga.nodes, paga.edges, ink, STAGE[theme]) : []),
      ...trajectoryLayers({
        prefix,
        paths: pathsAt(tt),
        milestones: milestonesAt(tt),
        focus,
        travelled: travelledAt(tt),
        ink,
        stage: STAGE[theme],
      }),
    ];
    if (paneStates && hovered !== null)
      out.push(
        new ScatterplotLayer({
          id: `${prefix}hover`,
          data: [hovered],
          getPosition: (i: number) => [pos.value[i * 3]!, pos.value[i * 3 + 1]!, pos.value[i * 3 + 2]!],
          filled: false,
          stroked: true,
          getLineColor: [ink[0], ink[1], ink[2], 255],
          lineWidthUnits: "pixels",
          getLineWidth: 1.5,
          radiusUnits: "pixels",
          getRadius: r + 5,
          billboard: true,
          parameters: ON_TOP,
          updateTriggers: { getPosition: hovered },
        }),
      );
    return out;
  };
  const layers = paneStates
    ? [...scene("flat:", 0, flatAttr, paneFit), ...scene("deep:", 1, deepAttr, paneFit)]
    : scene("", t, positions, fitZoom);
  const lineageName = (id: number) => run?.lineages.find((l) => l.id === id)?.name ?? `lineage ${id + 1}`;
  // The shown colour's value, then pseudotime and lineage when they are not already that value.
  const describeCell = (i: number) => {
    const parts = [describe(i) ?? `cell ${i}`];
    const pt = cells?.pseudotime[i];
    const primary = cells?.primaryLineage[i];
    if (colorBy?.kind !== "pseudotime" && pt !== undefined && Number.isFinite(pt)) parts.push(`pseudotime ${pt.toFixed(2)}`);
    if (colorBy?.kind !== "lineage" && primary !== undefined && primary > -2)
      parts.push(primary === -1 ? "shared trunk" : lineageName(primary));
    return parts.join(" · ");
  };
  const describeMilestone = (m: PlacedMilestone) => {
    const end = run?.lineages.find((l) => l.terminal_milestone === m.id);
    return m.kind === "tip" && end
      ? `End of ${end.name} · pseudotime ${m.pseudotime.toFixed(2)}`
      : `Fork · pseudotime ${m.pseudotime.toFixed(2)}`;
  };

  // What a selection is broken down by, and what Alt-click selects: the shown grouping.
  const leiden = clusterColumn === LEIDEN;
  const groupBy: GroupBy | null =
    colorBy?.kind === "lineage" && cells && run
      ? { title: "lineage", codes: cells.primaryLineage, labels: run.lineages.map((l) => l.name) }
      : clusters?.codes
        ? {
            title: leiden ? "Leiden cluster" : clusterColumn,
            codes: clusters.codes,
            labels: clusters.categories.map((c) => (leiden ? `Cluster ${c}` : c)),
          }
        : null;

  // A lasso takes every drawn cell inside it on screen, at any depth; cells past the cursor are not drawn.
  // Side by side, it takes cells in the pane it starts in; each pane projects into its own pixels.
  const finishLasso = (poly: Point[]) => {
    const pane = paneViewports ? paneAt(poly[0]![0]) : null;
    const vp = paneViewports ? paneViewports[pane === "flat" ? 0 : 1]! : viewport;
    if (!vp || !flat) return;
    const xyz = (pane ? pane === "deep" : t === 1) && deep ? deep.xyz : flat.xyz;
    const onPane = poly.map(([x, y]): Point => [x - vp.x, y - vp.y]);
    const next = selectInPolygon(projectAll(xyz, vp.pixelProjectionMatrix as unknown as number[]), onPane);
    if (limit < 2) for (let i = 0; i < n; i++) if (ptValues[i]! > limit) next[i] = 0;
    setMask(combine(mask, next, op));
  };
  const selectGroupAt = (x: number, y: number) => {
    const layerIds = paneStates ? ["flat:cells", "deep:cells"] : ["cells"];
    const i = deckRef.current?.pickObject({ x, y, radius: 4, layerIds })?.index ?? -1;
    const code = groupBy && i >= 0 ? groupBy.codes[i]! : -1;
    if (!groupBy || code < 0) return;
    const next = new Uint8Array(n);
    for (let k = 0; k < n; k++) if (groupBy.codes[k] === code) next[k] = 1;
    setMask(combine(mask, next, op));
  };
  const local = (e: React.PointerEvent<HTMLDivElement>): Point => {
    const r = e.currentTarget.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  };

  // The root pin follows the root cell through pans, zooms and the morph; side by side, in each pane.
  const pinOn = (vp: Viewport | null | undefined, xyz: Vec3 | null): Point | null => {
    if (!vp || !xyz) return null;
    const [x, y] = vp.project(xyz) as [number, number];
    return x >= 0 && x <= vp.width && y >= 0 && y <= vp.height ? [vp.x + x, vp.y + y] : null;
  };
  const pins = !root
    ? []
    : paneViewports
      ? paneViewports.map((vp, k) => pinOn(vp, at(root.cell_index, k)))
      : [pinOn(viewport, at(root.cell_index))];

  const embeddingInfo = dataset?.embeddings.find((e) => e.name === embedding);
  const p = run?.provenance;
  const dcCaption = `${axesLabel(axes)} · ${components ? `3 of ${components}` : "3"} diffusion components${fog ? " · far cells fade" : ""}`;
  const flatCaption = `${embedding} · a 2-D display projection`;
  const caption =
    view === "3d"
      ? dcCaption
      : !embedding
        ? ""
        : paths.length && p
          ? `${embedding} · curves are projections of the tree inferred in ${p.inference_rep}`
          : legend
            ? `${embedding} · coloured by ${legend.title === "Leiden" ? "Leiden cluster" : legend.title}`
            : `${embedding}, ${embeddingInfo?.source === "file" ? "stored in the file" : "computed by lumenCell"} · ${fmt(n)} cells`;
  const right = Math.max(12, padRight);
  // The view buttons sit clear of the "Show evidence" button when the sheet is hidden.
  const toggleRight = Math.max(padRight + 6, run ? 126 : 18);
  // Side by side, the 3-D caption sits beside the view buttons when there is room, else below them.
  const captionRoom = size.width - toggleRight - 172 - (RAIL + half + 12);
  const gimbalAt = paneViewports ? paneViewports[1] : t > 0 ? viewport : null;
  const axisNames = [`DC${axes[0]}`, `DC${axes[1]}`, `DC${axes[2]}`] as const;

  // The plot as a PNG, with its caption, root, legend, gimbal and provenance.
  const savePng = async () => {
    const deck = deckRef.current?.deck;
    const gl = containerRef.current?.querySelector("canvas");
    if (!deck || !gl || !dataset) return;
    deck.redraw("export"); // draw now: the canvas is read in this same task
    const css = getComputedStyle(document.documentElement);
    const token = (name: string) => css.getPropertyValue(`--color-${name}`).trim();
    const method = run && methods?.find((m) => m.name === run.provenance.method);
    const footer = [
      `lumenCell ${dataset.versions.lumencell ?? ""}`,
      `${dataset.name} · ${dataset.fingerprint}`,
      run ? `${method?.label ?? run.provenance.method}${method?.experimental ? " (experimental)" : ""} · run ${run.run_id}` : "no run",
      new Date().toLocaleDateString("en-CA"), // YYYY-MM-DD, local
    ].join(" · ");
    const blob = await exportPng({
      gl,
      crop: { x: RAIL, y: 0, width: size.width - RAIL - padRight, height: size.height - BAR },
      colors: { stage: token("stage"), ink: token("ink"), muted: token("ink-muted"), accent: token("accent") },
      font: getComputedStyle(document.body).fontFamily,
      captions: paneStates ? [flatCaption, dcCaption] : [caption],
      pins: pins.flatMap((pin) => (pin ? [{ x: pin[0], y: pin[1] }] : [])),
      pinLabel: root ? `Root · ${describeRoot(root.spec)}` : "",
      legend,
      gimbal: gimbalAt ? { m: gimbalAt.viewMatrix, names: axisNames } : null,
      footer,
    });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `lumencell-${dataset.name}-${view}.png`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  const raised = cursor !== null && run !== undefined; // the gene trends sit above the bar

  return (
    <div
      ref={containerRef}
      className="absolute inset-0 overflow-hidden bg-stage"
      data-testid="scatter"
      onPointerMove={(e) => {
        const pane = paneStates ? paneAt(local(e)[0]) : active;
        if (pane !== active) setActive(pane);
      }}
      onKeyDown={(e) => {
        // R resets the view, from the plot or its controls (never while typing in a field). Esc clears the selection.
        const field = e.target instanceof HTMLInputElement || e.target instanceof HTMLSelectElement;
        if (e.key.toLowerCase() === "r" && !field && !e.metaKey && !e.ctrlKey && !e.altKey && viewState) reset();
        if (e.key === "Escape") {
          setDrawing(null);
          setMask(null);
        }
      }}
      // In lasso mode a left drag on the plot draws, and never reaches deck.gl (capture phase).
      onPointerDownCapture={(e) => {
        if (!lasso || e.button !== 0 || !(e.target instanceof HTMLCanvasElement) || (!paneStates && t !== 0 && t !== 1)) return;
        e.stopPropagation();
        e.currentTarget.setPointerCapture(e.pointerId);
        setDrawing([local(e)]);
      }}
      onPointerMoveCapture={(e) => {
        if (!drawing) return;
        e.stopPropagation();
        const p = local(e);
        const last = drawing[drawing.length - 1]!;
        if (Math.hypot(p[0] - last[0], p[1] - last[1]) >= 3) setDrawing([...drawing, p]);
      }}
      onPointerUpCapture={(e) => {
        if (!drawing) return;
        e.stopPropagation();
        if (drawing.length >= 3) finishLasso(drawing);
        else if (e.altKey) selectGroupAt(...drawing[0]!);
        setDrawing(null);
      }}
    >
      {viewState && (
        <DeckGL
          ref={deckRef}
          views={paneViews ?? stageView}
          viewState={paneStates ?? viewState}
          onLoad={labelPlot}
          onViewStateChange={({ viewId, viewState: vs }) => {
            if (!paneStates && t !== 0 && t !== 1) return;
            setSpinning(false);
            // The target stays over the cells, so a zoom or pan cannot carry the view off the plot.
            const next = vs as OrbitViewState;
            const box = paneStates ? (viewId === "deep" ? deep?.bounds : flat?.bounds) : morphBox;
            const target = box ? clampTarget(next.target as Vec3, box) : (next.target as Vec3);
            setCamera(paneStates ? viewId : VIEW_ID, { ...next, target });
          }}
          // Arrows pan, -/= zoom; with Shift they orbit, in 3-D only. Side by side, each pane has its own.
          controller={paneStates ? undefined : t === 0 ? flatControl(true) : t === 1 ? deepControl(true, lasso) : false}
          layerFilter={paneStates ? ({ layer, viewport: vp }) => layer.id.startsWith(`${vp.id}:`) : undefined}
          onHover={paneStates ? (info) => setHovered(kind(info.layer?.id) === "cells" && info.index >= 0 ? info.index : null) : undefined}
          layers={layers}
          onClick={(info) => {
            if (!picking || kind(info.layer?.id) !== "cells" || info.index < 0 || !clusters?.codes) return;
            const code = clusters.codes[info.index]!;
            if (code >= 0) setRoot.mutate({ kind: "cluster", key: clusterColumn, value: clusters.categories[code]! });
          }}
          getTooltip={({ index, layer, object }) => {
            const text =
              kind(layer?.id) === "milestones"
                ? describeMilestone(object as PlacedMilestone)
                : kind(layer?.id) === "paga-nodes"
                  ? `Cluster ${(object as PlacedNode).id} · ${fmt((object as PlacedNode).n_cells)} cells · PAGA`
                  : kind(layer?.id) === "cells" && index >= 0
                    ? describeCell(index)
                    : null;
            return text
              ? {
                  text,
                  style: {
                    background: "var(--color-panel-raised)",
                    color: "var(--color-ink)",
                    fontFamily: "var(--font-sans)",
                    fontSize: "12px",
                    borderRadius: "4px",
                    padding: "4px 6px",
                  },
                }
              : null;
          }}
          getCursor={({ isDragging, isHovering }) =>
            lasso
              ? "crosshair"
              : isDragging
                ? "grabbing"
                : isHovering && picking
                  ? "crosshair"
                  : (paneStates ? active === "deep" : t === 1)
                    ? "grab"
                    : "default"
          }
        />
      )}

      {!flat && (
        <div className="absolute inset-0 grid place-items-center pl-[304px] text-[13px] text-ink-muted">
          {embedding ? "Loading the embedding…" : "No embedding yet: preparation will compute one."}
        </div>
      )}

      {root &&
        pins.map(
          (pin, k) =>
            pin && (
              <div key={k} className="pointer-events-none absolute" style={{ left: pin[0], top: pin[1] }}>
                <motion.span
                  key={root.obs_name}
                  initial={{ opacity: 0, scale: 0.9 }}
                  animate={{ opacity: 1, scale: 1 }}
                  transition={{ duration: 0.2, ease: EASE_OUT }}
                  className="absolute -top-2 -left-2 size-4 rounded-full border-2 border-stage bg-accent"
                />
                {k === 0 && (
                  <span className="absolute -top-[11px] left-4 rounded bg-panel-raised px-1.5 py-0.5 text-xs whitespace-nowrap text-ink">
                    Root · {describeRoot(root.spec)}
                  </span>
                )}
              </div>
            ),
        )}

      {paneStates ? (
        <>
          <div aria-hidden className="pointer-events-none absolute inset-y-0 w-px bg-line" style={{ left: RAIL + half }} />
          <p className="pointer-events-none absolute top-[22px] text-xs text-ink-muted" style={{ left: RAIL + 4, maxWidth: half - 16 }}>
            {flatCaption}
          </p>
          <p
            className="pointer-events-none absolute text-xs text-ink-muted"
            style={
              captionRoom >= 200
                ? { top: 22, left: RAIL + half + 12, maxWidth: captionRoom }
                : { top: 56, left: RAIL + half + 12, maxWidth: half - 24 }
            }
          >
            {dcCaption}
          </p>
        </>
      ) : (
        <p className="pointer-events-none absolute top-[22px] left-[308px] text-xs text-ink-muted">{caption}</p>
      )}

      {drawing && drawing.length > 1 && (
        <svg aria-hidden className="pointer-events-none absolute inset-0 size-full">
          <polygon
            points={drawing.map((p) => p.join(",")).join(" ")}
            fill="var(--color-ink)"
            fillOpacity={0.08}
            stroke="var(--color-ink)"
            strokeWidth={1.25}
            strokeDasharray="4 3"
          />
        </svg>
      )}

      <div className="absolute top-[48px] left-[304px] z-10">
        <SelectionBar lasso={lasso} groupBy={groupBy} total={n} />
      </div>

      <div
        role="group"
        aria-label="View"
        className="absolute top-3 z-10 flex gap-0.5 rounded-lg bg-panel p-1"
        style={{ right: toggleRight }}
      >
        {(["2d", "3d", "split"] as const).map((v) => (
          <button
            key={v}
            type="button"
            aria-pressed={view === v}
            disabled={v !== "2d" && !deep}
            title={
              v !== "2d" && !deep
                ? "3-D opens once preparation has finished"
                : v === "split"
                  ? "The 2-D embedding and DC1–DC3 side by side"
                  : undefined
            }
            onClick={() => {
              setSpinning(false);
              setView(v);
            }}
            className={`press h-[26px] rounded-md px-2.5 text-[13px] disabled:text-ink-faint ${
              view === v ? "bg-panel-raised text-ink" : "text-ink-muted hover:text-ink"
            }`}
          >
            {{ "2d": "2-D", "3d": "3-D", split: "Split" }[v]}
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
                className="absolute right-0 bottom-full mb-2 w-72 rounded-xl bg-panel p-3 text-[13px] text-ink"
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
                {components && (
                  <fieldset className="mt-3">
                    <legend>3-D axes</legend>
                    <div className="mt-1.5 grid grid-cols-3 gap-1.5">
                      {(["x", "y", "z"] as const).map((name, k) => (
                        <label key={name} className="flex flex-col gap-0.5 text-xs text-ink-muted">
                          {name}
                          <Select
                            aria-label={`${name.toUpperCase()} axis`}
                            value={axes[k]}
                            onChange={(e) => setAxes(setAxis(axes, k, Number(e.target.value)))}
                          >
                            {Array.from({ length: components }, (_, i) => (
                              <option key={i} value={i + 1}>
                                DC{i + 1}
                              </option>
                            ))}
                          </Select>
                        </label>
                      ))}
                    </div>
                  </fieldset>
                )}
                <fieldset className="mt-3">
                  <legend>Graphs drawn over this view</legend>
                  <label className="mt-1.5 flex items-center gap-2">
                    <input
                      type="checkbox"
                      checked={overlays.paga}
                      onChange={(e) => setOverlay("paga", e.target.checked)}
                      className="accent-ink"
                    />
                    PAGA cluster graph
                  </label>
                  <label className="mt-1 flex items-center gap-2">
                    <input
                      type="checkbox"
                      checked={overlays.knn}
                      onChange={(e) => setOverlay("knn", e.target.checked)}
                      className="accent-ink"
                    />
                    Neighbour graph
                  </label>
                  <p className="mt-1.5 text-xs leading-snug text-ink-muted">
                    Both are built from neighbours in PCA space. A line joins where its ends sit in this view, so its length means nothing.
                    {knn?.sampled && ` Showing ${fmt(knn.src.length)} of ${fmt(knn.total)} neighbour edges, a random sample.`}
                  </p>
                </fieldset>
              </div>
            )}
            <button
              type="button"
              aria-label="Lasso select"
              aria-pressed={lasso}
              title="Lasso select: drag to draw around cells"
              onClick={() => {
                setLasso((on) => !on);
                setDrawing(null);
              }}
              className={`press grid size-8 place-items-center rounded-md text-ink ${lasso ? "bg-panel-raised" : "bg-panel"}`}
            >
              <Lasso size={16} />
            </button>
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
            <button
              type="button"
              aria-label="Export PNG"
              title="Save the plot as a PNG, with its legend and where it came from"
              onClick={() => void savePng()}
              className="press grid size-8 place-items-center rounded-md bg-panel text-ink"
            >
              <ImageDown size={16} />
            </button>
            {t === 1 && !paneStates && <CameraTour id={VIEW_ID} reduce={reduce} onStart={() => setSpinning(false)} />}
            {(paneStates || t === 1) && (
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
          {gimbalAt && <Gimbal m={gimbalAt.viewMatrix as unknown as number[]} names={axisNames} />}
        </div>
      )}
    </div>
  );
}
