import type { ColorBy } from "@/state/ui";

export const LEIDEN = "lumencell_leiden";
export const NONE = "none";

/** A ColorBy as one string, for selects and for tests. */
export function encodeColorBy(c: ColorBy | null): string {
  if (!c) return NONE;
  if (c.kind === "obs") return `obs:${c.column}`;
  if (c.kind === "gene") return `gene:${c.gene}`;
  if (c.kind === "disagreement") return `disagreement:${c.a}:${c.b}`;
  return c.kind;
}

export function decodeColorBy(v: string): ColorBy | null {
  if (v.startsWith("obs:")) return { kind: "obs", column: v.slice(4) };
  if (v.startsWith("gene:")) return { kind: "gene", gene: v.slice(5) };
  if (v === "pseudotime" || v === "lineage") return { kind: v };
  if (v.startsWith("disagreement:")) {
    const [, a = "", b = ""] = v.split(":");
    return { kind: "disagreement", a, b };
  }
  return null;
}

export function columnLabel(name: string): string {
  return name === LEIDEN ? "Leiden clusters" : name;
}

export function categoryLabel(column: string, value: string): string {
  return column === LEIDEN ? `Cluster ${value}` : value;
}

const EMBEDDING_LABELS: Record<string, string> = {
  umap: "UMAP",
  tsne: "t-SNE",
  draw_graph_fa: "FA",
  dc12: "DC 2D",
  dc123: "DC 3D",
};

export function embeddingLabel(name: string): string {
  return EMBEDDING_LABELS[name] ?? name.toUpperCase();
}
