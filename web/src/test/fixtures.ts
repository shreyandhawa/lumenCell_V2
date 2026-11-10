import type * as T from "@/api/generated";

const stages = ["normalize", "hvg", "pca", "neighbors", "leiden", "diffmap", "paga", "umap"];

export function dataset(prep: T.PrepStatus["state"], root: T.ResolvedRoot | null = null): T.DatasetInfo {
  const ready = prep === "ready";
  return {
    schema_version: 1,
    name: "synthetic",
    n_obs: 4,
    n_vars: 3,
    x_kind: "counts",
    gene_symbol_column: null,
    obs: [
      { name: "cell_type", kind: "categorical", categories: ["stem", "A", "B"], colors: null, min: null, max: null },
      ...(ready
        ? [{ name: "lumencell_leiden", kind: "categorical" as const, categories: ["0", "1"], colors: null, min: null, max: null }]
        : []),
    ],
    embeddings: ready ? [{ name: "umap", dims: 2, source: "computed" }] : [],
    default_embedding: ready ? "umap" : null,
    cluster_key: ready ? "lumencell_leiden" : null,
    label_key: "cell_type",
    prep: {
      state: prep,
      stages: stages.map((name, i) => ({
        name,
        state: ready || i < 2 ? "done" : i === 2 && prep === "running" ? "running" : "pending",
        seconds: ready || i < 2 ? 0.5 : null,
        message: null,
      })),
      job_id: prep === "running" ? "job-prep" : null,
    },
    root,
    versions: { lumencell: "0.1.0" },
  };
}

export const resolvedRoot: T.ResolvedRoot = {
  spec: { kind: "cluster", key: "cell_type", value: "stem" },
  cell_index: 0,
  obs_name: "cell0",
  cluster_key: "cell_type",
  cluster_label: "stem",
  resolved_by: "the cell of cell_type = stem farthest from the rest",
  warnings: [],
};

export const suggestions: T.RootSuggestions = {
  suggestions: [
    {
      spec: { kind: "cluster", key: "cell_type", value: "stem" },
      source: "time_column",
      score: 0.9,
      evidence: "90% of its cells come from the earliest day (day 0)",
      is_extremity: true,
    },
  ],
};

export const methods: T.MethodList = {
  methods: [
    {
      name: "paga_dpt",
      label: "PAGA + DPT",
      family: "cluster graph + diffusion pseudotime",
      description: "Diffusion pseudotime with lineages along the PAGA tree.",
      available: true,
      install_hint: null,
      experimental: true,
      params_schema: { properties: {} },
      reference_scores: {},
    },
  ],
};

export const run: T.TrajectoryRun = {
  run_id: "run1",
  state: "done",
  topology_kind: "tree",
  milestones: [
    { id: 0, kind: "root", label: "0", pseudotime: 0 },
    { id: 1, kind: "tip", label: "1", pseudotime: 1 },
  ],
  segments: [{ id: 0, source: 0, target: 1, length: 1 }],
  lineages: [
    {
      id: 0,
      name: "Erythroid",
      auto_name: "Erythroid",
      terminal_milestone: 1,
      terminal_label: "A",
      terminal_purity: 0.9,
      n_cells: 4,
      segments: [0],
    },
  ],
  divergence: [],
  warnings: [{ code: "disconnected_cells", severity: "info", message: "2 cells are not connected." }],
  n_unassigned: 0,
  runtime_seconds: 0.1,
  cell_arrays: ["pseudotime", "membership", "primary_lineage", "segment", "progress"],
  provenance: {
    method: "paga_dpt",
    method_version: "1.12.5",
    lumencell_version: "0.1.0",
    params: { n_dcs: 10 },
    param_hash: "run1",
    dataset_fingerprint: "fp",
    inference_rep: "X_diffmap",
    inference_dims: 10,
    root: resolvedRoot,
    terminals: ["A"],
    terminals_source: "auto",
    created_at: "2026-11-10T00:00:00Z",
  },
};
