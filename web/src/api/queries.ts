/**
 * Server state, via TanStack Query. Query keys are the resource identity, so a
 * 100k-cell embedding or gene column is fetched once and shared by every view
 * that asks for it (in-flight requests are deduplicated too).
 */
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api, type RootSpec } from "@/api/client";
import type * as T from "@/api/generated";

const forever = { staleTime: Infinity, gcTime: 30 * 60_000 } as const;

export const keys = {
  dataset: ["dataset"] as const,
  embedding: (name: string) => ["embedding", name] as const,
  obs: (column: string) => ["obs", column] as const,
  expression: (gene: string) => ["expression", gene.toLowerCase()] as const,
  methods: ["methods"] as const,
  root: ["root"] as const,
  suggestions: ["root", "suggestions"] as const,
  runs: ["runs"] as const,
  run: (id: string) => ["run", id] as const,
  cells: (id: string) => ["cells", id] as const,
  curves: (id: string, embedding: string) => ["curves", id, embedding] as const,
  genes: (q: string) => ["genes", q] as const,
};

export function useDataset() {
  return useQuery({
    queryKey: keys.dataset,
    queryFn: api.dataset,
    // Poll while preprocessing runs: embeddings and clusters appear when it finishes.
    refetchInterval: (q) => (q.state.data && q.state.data.prep.state !== "ready" ? 1500 : false),
  });
}

export function usePrepReady(): boolean {
  return useDataset().data?.prep.state === "ready";
}

export function useEmbedding(name: string | null) {
  return useQuery({
    queryKey: keys.embedding(name ?? ""),
    queryFn: async () => (await api.embedding(name!)).get("coords", "f32").data,
    enabled: name !== null,
    ...forever,
  });
}

export interface ObsVector {
  kind: "categorical" | "boolean" | "continuous";
  codes?: Int16Array | Int32Array;
  values?: Float32Array;
  categories: string[];
}

export function useObs(column: string | null) {
  return useQuery({
    queryKey: keys.obs(column ?? ""),
    queryFn: async (): Promise<ObsVector> => {
      const c = await api.obs(column!);
      const kind = c.meta.kind as ObsVector["kind"];
      if (kind === "continuous") return { kind, values: c.get("values", "f32").data, categories: [] };
      const arr = c.arrays.get("codes")!;
      return { kind, codes: arr.data as Int16Array | Int32Array, categories: c.meta.categories as string[] };
    },
    enabled: column !== null,
    ...forever,
  });
}

export function useGeneExpression(gene: string | null) {
  return useQuery({
    queryKey: keys.expression(gene ?? ""),
    queryFn: async () => (await api.expression([gene!])).get("expression", "f32").data,
    enabled: gene !== null,
    ...forever,
  });
}

export function useGeneSearch(q: string) {
  return useQuery({ queryKey: keys.genes(q), queryFn: () => api.genes(q, 12), enabled: q.length > 0 });
}

export function useMethods() {
  return useQuery({ queryKey: keys.methods, queryFn: api.methods, ...forever });
}

export function useRoot() {
  const ready = usePrepReady();
  return useQuery({ queryKey: keys.root, queryFn: api.root, enabled: ready });
}

export function useRootSuggestions() {
  const ready = usePrepReady();
  return useQuery({ queryKey: keys.suggestions, queryFn: api.rootSuggestions, enabled: ready });
}

export function useSetRoot() {
  const qc = useQueryClient();
  return useMutation({
    mutationFn: (spec: RootSpec) => api.setRoot(spec),
    onSuccess: (resolved) => {
      qc.setQueryData<T.RootState>(keys.root, { root: resolved });
      void qc.invalidateQueries({ queryKey: keys.dataset });
    },
  });
}

export function useResolveRoot() {
  return useMutation({ mutationFn: (spec: RootSpec) => api.resolveRoot(spec) });
}

export function useRuns() {
  const ready = usePrepReady();
  return useQuery({
    queryKey: keys.runs,
    queryFn: api.runs,
    enabled: ready,
    refetchInterval: (q) =>
      q.state.data?.runs.some((r) => r.state === "queued" || r.state === "running") ? 1000 : false,
  });
}

/**
 * The active run's id once it has finished, else null. Run details, cells and
 * curves only exist for finished runs, so every per-run query keys off this.
 */
export function useDoneRunId(runId: string | null): string | null {
  const runs = useRuns().data?.runs;
  return runId !== null && runs?.some((r) => r.run_id === runId && r.state === "done") ? runId : null;
}

export function useCreateRun() {
  const qc = useQueryClient();
  return useMutation({
    mutationFn: (req: T.TrajectoryRequest) => api.createRun(req),
    onSuccess: () => void qc.invalidateQueries({ queryKey: keys.runs }),
  });
}

export function useRun(id: string | null) {
  return useQuery({ queryKey: keys.run(id ?? ""), queryFn: () => api.run(id!), enabled: id !== null });
}

export interface RunCells {
  pseudotime: Float32Array;
  membership: Uint8Array;
  primaryLineage: Int16Array;
  nLineages: number;
}

export function useRunCells(id: string | null) {
  return useQuery({
    queryKey: keys.cells(id ?? ""),
    queryFn: async (): Promise<RunCells> => {
      const c = await api.cells(id!);
      const membership = c.get("membership", "u8");
      return {
        pseudotime: c.get("pseudotime", "f32").data,
        membership: membership.data,
        primaryLineage: c.get("primary_lineage", "i16").data,
        nLineages: membership.shape[0] ?? 0,
      };
    },
    enabled: id !== null,
    ...forever,
  });
}

export function useCurves(id: string | null, embedding: string | null) {
  return useQuery({
    queryKey: keys.curves(id ?? "", embedding ?? ""),
    queryFn: () => api.curves(id!, embedding!),
    enabled: id !== null && embedding !== null,
    ...forever,
  });
}

export function useRenameLineage(runId: string) {
  const qc = useQueryClient();
  return useMutation({
    mutationFn: ({ lineage, name }: { lineage: number; name: string }) => api.renameLineage(runId, lineage, name),
    onSuccess: () => void qc.invalidateQueries({ queryKey: keys.run(runId) }),
  });
}

export function useStartPrepare() {
  const qc = useQueryClient();
  return useMutation({
    mutationFn: api.startPrepare,
    onSuccess: () => void qc.invalidateQueries({ queryKey: keys.dataset }),
  });
}
