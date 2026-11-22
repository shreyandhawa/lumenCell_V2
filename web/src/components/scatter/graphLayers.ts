/**
 * Graphs drawn over the plot. Both come from the cells' neighbours in PCA space; here a line only
 * joins where its two ends sit in this view, so a long line says nothing about distance.
 */
import type { Layer } from "@deck.gl/core";
import { LineLayer, ScatterplotLayer } from "@deck.gl/layers";
import type * as T from "@/api/generated";
import type { Vec3 } from "@/lib/embedding";

type RGBA = [number, number, number, number];
export type PlacedNode = T.GraphNode & { position: Vec3 };
type Attr = { value: Float32Array; size: number };
export type Ends = { length: number; attributes: { getSourcePosition: Attr; getTargetPosition: Attr } };

const ON_TOP = { depthCompare: "always", depthWriteEnabled: false } as const;

/** PAGA: each Leiden cluster at the mean of its cells, sized by its cells; lines by connectivity. */
export function pagaLayers(
  prefix: string,
  nodes: PlacedNode[],
  edges: T.GraphEdge[],
  ink: RGBA,
  stage: readonly [number, number, number],
): Layer[] {
  const at = new Map(nodes.map((n) => [n.id, n.position]));
  const lines = edges.filter((e) => at.has(e.source) && at.has(e.target));
  const biggest = Math.max(1, ...nodes.map((n) => n.n_cells));
  return [
    new LineLayer<T.GraphEdge>({
      id: `${prefix}paga-edges`,
      data: lines,
      getSourcePosition: (e) => at.get(e.source)!,
      getTargetPosition: (e) => at.get(e.target)!,
      getColor: (e) => [ink[0], ink[1], ink[2], Math.round(40 + 215 * e.weight)],
      getWidth: (e) => 1 + 5 * e.weight,
      widthUnits: "pixels",
      parameters: ON_TOP,
      updateTriggers: { getColor: ink },
    }),
    new ScatterplotLayer<PlacedNode>({
      id: `${prefix}paga-nodes`,
      data: nodes,
      getPosition: (n) => n.position,
      getRadius: (n) => 3 + 9 * Math.sqrt(n.n_cells / biggest),
      radiusUnits: "pixels",
      getFillColor: [...stage, 255],
      getLineColor: [ink[0], ink[1], ink[2], 255],
      stroked: true,
      lineWidthUnits: "pixels",
      getLineWidth: 1.5,
      billboard: true,
      pickable: true,
      parameters: ON_TOP,
      updateTriggers: { getFillColor: stage, getLineColor: ink },
    }),
  ];
}

/** The neighbour graph, faint, under the cells: it shows in the gaps, where its long-range lines are. */
export function knnLayer(prefix: string, ends: Ends, ink: RGBA, transitions: object | undefined): Layer {
  return new LineLayer({
    id: `${prefix}knn`,
    data: ends,
    getColor: [ink[0], ink[1], ink[2], 40],
    getWidth: 1,
    widthUnits: "pixels",
    parameters: ON_TOP,
    transitions: transitions && { getSourcePosition: transitions, getTargetPosition: transitions },
  });
}
