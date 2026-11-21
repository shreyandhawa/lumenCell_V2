import { OrbitView, type OrbitViewState } from "@deck.gl/core";
import { DataFilterExtension } from "@deck.gl/extensions";
import { PathLayer, ScatterplotLayer } from "@deck.gl/layers";
import DeckGL from "@deck.gl/react";
import { RotateCcw } from "lucide-react";
import { animate } from "motion";
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
import { describeRoot } from "@/components/trajectory/RootPicker";
import { fmt } from "@/components/ui";
import { useUi } from "@/state/ui";

const LEIDEN = "lumencell_leiden";
const FILTER = new DataFilterExtension({ filterSize: 1 });
// The stage fills the window; the data is centred in the space between the rail, the sheet and the bar.
const RAIL = 304;
const SHEET = 384;
const BAR = 80;
const FLAT = { rotationX: 0, rotationOrbit: 0 };
const TILTED = { rotationX: 18, rotationOrbit: -32 };
const EASE_OUT = [0.23, 1, 0.32, 1] as const;
const EASE_IN_OUT = [0.77, 0, 0.175, 1] as const;
/** How much the farthest cells fade in 3-D (share of their alpha taken away). */
const FAR_FADE = 0.75;

type Camera = Pick<OrbitViewState, "target" | "zoom" | "rotationX" | "rotationOrbit">;
type Vec3 = [number, number, number];

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

/** Centre an (n x dims) embedding on 0 and scale its longest side to 1, as xyz (z = 0 in 2-D). */
function normalise(src: Float32Array, dims: 2 | 3) {
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < src.length; i++) {
    const v = src[i]!;
    const k = i % dims;
    if (Number.isFinite(v)) {
      lo[k] = Math.min(lo[k]!, v);
      hi[k] = Math.max(hi[k]!, v);
    }
  }
  const mid = [0, 1, 2].map((k) => (k < dims ? (lo[k]! + hi[k]!) / 2 : 0));
  const span = Math.max(...[0, 1, 2].map((k) => (k < dims ? hi[k]! - lo[k]! : 0))) || 1;
  const map = (p: ArrayLike<number>): Vec3 => [0, 1, 2].map((k) => (k < dims ? (p[k]! - mid[k]!) / span : 0)) as Vec3;
  const n = src.length / dims;
  const xyz = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) xyz.set(map(src.subarray(i * dims, i * dims + dims)), i * 3);
  return { xyz, map };
}

const lerp = (a: number, b: number, t: number) => a + (b - a) * t;
const lerpCam = (a: Camera, b: Camera, t: number): Camera => ({
  target: [0, 1, 2].map((k) => lerp(a.target[k] ?? 0, b.target[k] ?? 0, t)) as Vec3,
  zoom: lerp(a.zoom as number, b.zoom as number, t),
  rotationX: lerp(a.rotationX ?? 0, b.rotationX ?? 0, t),
  rotationOrbit: lerp(a.rotationOrbit ?? 0, b.rotationOrbit ?? 0, t),
});

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

  const dataset = useDataset().data;
  const n = dataset?.n_obs ?? 0;
  const root = useRoot().data?.root ?? null;
  const run = useRun(runId).data;
  const runState = useRuns().data?.runs.find((r) => r.run_id === uiRunId)?.state;
  const flatRaw = useEmbedding(embedding).data;
  const deepRaw = useEmbedding(runId ? "dc123" : null).data;
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

  const [containerRef, size] = useSize();
  const sheetVisible = run !== undefined && sheetOpen;
  const [padTween, setPadTween] = useState(0);
  const padRight = reduce ? (sheetVisible ? SHEET : 0) : padTween;
  const padding = { left: RAIL, right: padRight, top: 0, bottom: BAR };
  const stageView = useMemo(
    () => new OrbitView({ id: "stage", orbitAxis: "Y", padding: { left: RAIL, right: padRight, top: 0, bottom: BAR } }),
    [padRight],
  );
  const fitZoom = Math.log2(
    0.8 * Math.max(120, Math.min(size.width - padding.left - padding.right, size.height - padding.bottom)),
  );
  const [userViewState, setViewState] = useState<OrbitViewState | null>(null);
  // The first fit, until the user or an animation moves the camera.
  const viewState = useMemo<OrbitViewState | null>(
    () =>
      userViewState ??
      (flat && size.width ? { target: [0, 0, 0], ...FLAT, zoom: fitZoom, minZoom: fitZoom - 3, maxZoom: fitZoom + 6 } : null),
    [userViewState, flat, size.width, fitZoom],
  );
  const [t, setT] = useState(0); // 0 = 2-D embedding, 1 = DC1-DC3; in between while morphing
  const [wash, setWash] = useState(1);

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

  // 3-D needs a finished run.
  useEffect(() => {
    if (!runId && view === "3d") setView("2d");
  }, [runId, view, setView]);

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
    return tween(cam, toT, 0.6, EASE_IN_OUT);
    // eslint-disable-next-line react-hooks/exhaustive-deps -- tween reads live values; only these start a morph
  }, [view, deep]);

  const reset = () => tween({ target: [0, 0, 0], zoom: fitZoom, ...(view === "3d" ? TILTED : FLAT) }, t, 0.3, EASE_OUT);

  // The set piece: when a run finishes, pseudotime washes out from the root, once.
  const washPending = useRef(false);
  const lastState = useRef(runState);
  useEffect(() => {
    if (runState === "done" && (lastState.current === "queued" || lastState.current === "running")) washPending.current = true;
    lastState.current = runState;
  }, [runState]);
  useEffect(() => {
    if (!cells || !washPending.current) return;
    washPending.current = false;
    if (reduce) return;
    const c = animate(0, 1, { duration: 0.9, ease: EASE_OUT, onUpdate: setWash });
    return () => {
      c.stop();
      setWash(1);
    };
  }, [cells, reduce]);

  const positions = useMemo(() => {
    if (!flat) return null;
    if (!deep || t === 0) return flat.xyz;
    if (t === 1) return deep.xyz;
    const out = new Float32Array(flat.xyz.length);
    for (let i = 0; i < out.length; i++) out[i] = lerp(flat.xyz[i]!, deep.xyz[i]!, t);
    return out;
  }, [flat, deep, t]);

  const viewport = useMemo(
    () => (viewState && size.width ? stageView.makeViewport({ width: size.width, height: size.height, viewState }) : null),
    [stageView, viewState, size.width, size.height],
  );

  // In 3-D, far cells fade: alpha falls with view-space depth.
  const fill = useMemo(() => {
    if (!colors || !positions || t === 0 || !viewport) return colors;
    const m = viewport.viewMatrix as unknown as number[];
    const depth = new Float32Array(n);
    let near = -Infinity;
    let far = Infinity;
    for (let i = 0; i < n; i++) {
      const d = m[2]! * positions[i * 3]! + m[6]! * positions[i * 3 + 1]! + m[10]! * positions[i * 3 + 2]! + m[14]!;
      depth[i] = d;
      near = Math.max(near, d);
      far = Math.min(far, d);
    }
    const span = near - far || 1;
    const out = colors.slice();
    for (let i = 0; i < n; i++) out[i * 4 + 3] = Math.round(colors[i * 4 + 3]! * (1 - FAR_FADE * t * ((near - depth[i]!) / span)));
    return out;
  }, [colors, positions, t, viewport, n]);

  const ptValues = useMemo(() => {
    const out = new Float32Array(n).fill(-1);
    if (cells) for (let i = 0; i < n; i++) if (Number.isFinite(cells.pseudotime[i]!)) out[i] = cells.pseudotime[i]!;
    return out;
  }, [n, cells]);
  // Cells past the cursor (or the wash) keep only their halo.
  const limit = cells ? (cursor ?? (wash < 1 ? wash : 2)) : 2;

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

  const ink: [number, number, number, number] = theme === "dark" ? [243, 245, 242, 190] : [20, 23, 22, 170];
  const layers = [];
  if (positions && fill) {
    const r = pointRadius(n);
    const overdraw = { depthCompare: "always", depthWriteEnabled: false } as const;
    layers.push(
      new ScatterplotLayer({
        id: "halo",
        data: { length: n, attributes: { getPosition: { value: positions, size: 3 } } },
        getFillColor: theme === "dark" ? [255, 255, 255, 18] : [0, 0, 0, 22],
        radiusUnits: "pixels",
        getRadius: r + 1.5,
        parameters: overdraw,
        updateTriggers: { getFillColor: theme },
      }),
      new ScatterplotLayer({
        id: "cells",
        data: {
          length: n,
          attributes: {
            getPosition: { value: positions, size: 3 },
            getFillColor: { value: fill, size: 4, normalized: true },
            getFilterValue: { value: ptValues, size: 1 },
          },
        },
        radiusUnits: "pixels",
        getRadius: r,
        pickable: true,
        autoHighlight: true,
        highlightColor: theme === "dark" ? [255, 255, 255, 220] : [0, 0, 0, 200],
        parameters: overdraw,
        extensions: [FILTER],
        filterRange: [-2, limit],
        filterSoftRange: [-2, limit - 0.06],
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
        jointRounded: true,
        capRounded: true,
        parameters: { depthCompare: "always" },
        updateTriggers: { getColor: [focus, theme], getWidth: focus },
      }),
    );
  }

  // The root pin follows the root cell through pans, zooms and the morph.
  const rootXyz = root && positions ? Array.from(positions.subarray(root.cell_index * 3, root.cell_index * 3 + 3)) : null;
  const pin = viewport && rootXyz?.every(Number.isFinite) ? viewport.project(rootXyz) : null;

  const embeddingInfo = dataset?.embeddings.find((e) => e.name === embedding);
  const p = run?.provenance;
  const caption =
    view === "3d"
      ? `DC1–DC3 · ${p?.inference_rep === "X_diffmap" ? `3 of ${p.inference_dims}` : "the first 3"} diffusion components · far cells fade`
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
    <div ref={containerRef} className="absolute inset-0 overflow-hidden bg-stage" data-testid="scatter">
      {viewState && (
        <DeckGL
          views={stageView}
          viewState={viewState}
          onViewStateChange={({ viewState: vs }) => {
            if (t === 0 || t === 1) setViewState(vs as OrbitViewState);
          }}
          controller={t === 0 ? { dragMode: "pan", dragRotate: false } : t === 1 ? { dragMode: "rotate" } : false}
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
            disabled={v === "3d" && !runId}
            title={v === "3d" && !runId ? "3-D opens once a run has finished" : undefined}
            onClick={() => setView(v)}
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
          <button
            type="button"
            aria-label="Reset view"
            onClick={reset}
            className="press grid size-8 place-items-center rounded-md bg-panel text-ink"
          >
            <RotateCcw size={16} />
          </button>
          {t > 0 && viewport && <Gimbal m={viewport.viewMatrix as unknown as number[]} />}
        </div>
      )}
    </div>
  );
}

/** Which way DC1, DC2 and DC3 point on screen. */
function Gimbal({ m }: { m: number[] }) {
  const axes = (["DC1", "DC2", "DC3"] as const).map((name, k) => {
    const [x, y, z] = [m[k * 4]!, m[k * 4 + 1]!, m[k * 4 + 2]!];
    const len = Math.hypot(x, y, z) || 1;
    return { name, x: (x / len) * 30, y: (-y / len) * 30 };
  });
  return (
    <svg
      role="img"
      aria-label="Orientation of DC1, DC2 and DC3"
      viewBox="-48 -48 96 96"
      className="size-24 rounded-full bg-stage/80 text-ink-muted"
    >
      <circle r={47} fill="none" stroke="currentColor" strokeOpacity={0.5} />
      {axes.map((a) => (
        <g key={a.name}>
          <line x1={0} y1={0} x2={a.x} y2={a.y} stroke="var(--color-ink)" strokeWidth={1.25} />
          <text x={a.x * 1.3} y={a.y * 1.3 + 3} textAnchor="middle" className="fill-ink-muted font-mono text-[9px]">
            {a.name}
          </text>
        </g>
      ))}
      <circle r={2.5} fill="var(--color-ink)" />
    </svg>
  );
}
