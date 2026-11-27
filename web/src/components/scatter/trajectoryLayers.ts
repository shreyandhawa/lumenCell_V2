/**
 * The trajectory over the cells: curves, milestones and the pseudotime playhead.
 *
 * Curves run through the middle of their cells, so a plain depth test would hide most of them in
 * 3-D. Each curve is drawn twice: faintly on top of everything (the stretches behind cells), then
 * fully where it is in front. The focused lineage is drawn fully on top: it is what is being studied.
 */
import type { Layer } from "@deck.gl/core";
import { PathLayer, ScatterplotLayer } from "@deck.gl/layers";
import type { Vec3 } from "@/lib/embedding";
import type { PlacedMilestone } from "@/lib/trajectory";

type RGBA = [number, number, number, number];

export interface TrajectoryLayersInput {
  paths: { lineage: number; path: Vec3[] }[];
  milestones: PlacedMilestone[];
  focus: number | null;
  /** The focused lineage's curve up to the cursor, ending where its cells are at the cursor. */
  travelled: Vec3[] | null;
  ink: RGBA;
  /** The focused lineage and the playhead. */
  accent: RGBA;
  stage: readonly [number, number, number];
  /** Prepended to every layer id: the pane the layers belong to in the split view. */
  prefix?: string;
}

const ON_TOP = { depthCompare: "always", depthWriteEnabled: false } as const;

export function trajectoryLayers({ paths, milestones, focus, travelled, ink, accent, stage, prefix = "" }: TrajectoryLayersInput): Layer[] {
  const layers: Layer[] = [];
  const ordered = [...paths.filter((p) => p.lineage !== focus), ...paths.filter((p) => p.lineage === focus)];
  const color = (d: { lineage: number }, alpha: number): RGBA =>
    d.lineage === focus ? [accent[0], accent[1], accent[2], alpha] : [ink[0], ink[1], ink[2], Math.round((ink[3] * alpha) / 255)];
  const width = (d: { lineage: number }) => (d.lineage === focus ? 3 : 1.5);
  const common = {
    data: ordered,
    getPath: (d: { path: Vec3[] }) => d.path,
    getWidth: width,
    widthUnits: "pixels" as const,
    billboard: true,
    jointRounded: true,
    capRounded: true,
  };
  if (ordered.length) {
    layers.push(
      new PathLayer({
        ...common,
        id: `${prefix}curves-behind`,
        getColor: (d) => color(d, d.lineage === focus ? 255 : 70),
        parameters: ON_TOP,
        updateTriggers: { getColor: [focus, ink, accent], getWidth: focus },
      }),
      new PathLayer({
        ...common,
        id: `${prefix}curves`,
        getColor: (d) => color(d, 255),
        updateTriggers: { getColor: [focus, ink, accent], getWidth: focus },
      }),
    );
  }
  if (travelled && travelled.length > 1) {
    layers.push(
      new PathLayer({
        id: `${prefix}travelled`,
        data: [travelled],
        getPath: (d: Vec3[]) => d,
        getColor: accent,
        getWidth: 5,
        widthUnits: "pixels",
        billboard: true,
        jointRounded: true,
        capRounded: true,
        parameters: ON_TOP,
      }),
      new ScatterplotLayer({
        id: `${prefix}playhead`,
        data: [travelled[travelled.length - 1]!],
        getPosition: (d: Vec3) => d,
        getFillColor: accent,
        getLineColor: [...stage, 255],
        stroked: true,
        lineWidthUnits: "pixels",
        getLineWidth: 2,
        radiusUnits: "pixels",
        getRadius: 6,
        billboard: true,
        parameters: ON_TOP,
        updateTriggers: { getLineColor: stage },
      }),
    );
  }
  // Forks are rings, tips are solid dots; the root has its own pin.
  const marks = milestones.filter((m) => m.kind === "fork" || m.kind === "tip");
  if (marks.length) {
    layers.push(
      new ScatterplotLayer({
        id: `${prefix}milestones`,
        data: marks,
        getPosition: (d: PlacedMilestone) => d.position,
        getFillColor: (d: PlacedMilestone) => (d.kind === "tip" ? ink : [...stage, 255]),
        getLineColor: ink,
        stroked: true,
        lineWidthUnits: "pixels",
        getLineWidth: 1.5,
        radiusUnits: "pixels",
        getRadius: (d: PlacedMilestone) => (d.kind === "tip" ? 4 : 5),
        billboard: true,
        pickable: true,
        parameters: ON_TOP,
        updateTriggers: { getFillColor: [ink, stage], getLineColor: ink },
      }),
    );
  }
  return layers;
}
