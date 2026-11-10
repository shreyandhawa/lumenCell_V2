/**
 * Thin typed wrappers over /api/v1. JSON shapes come from the generated
 * OpenAPI types; per-cell data comes back as LCB1 containers.
 */
import type * as T from "@/api/generated";
import { decodeLcb, type LcbContainer } from "@/api/lcb1";

export const API_BASE = "/api/v1";

export class ApiError extends Error {
  override name = "ApiError";
  constructor(
    readonly status: number,
    readonly code: string,
    message: string,
    readonly detail?: Record<string, unknown> | null,
  ) {
    super(message);
  }
}

async function request(path: string, init?: RequestInit): Promise<Response> {
  const res = await fetch(`${API_BASE}${path}`, init);
  if (!res.ok) {
    let body: T.ErrorBody | undefined;
    try {
      body = (await res.json()) as T.ErrorBody;
    } catch {
      // not a JSON error body
    }
    throw new ApiError(
      res.status,
      body?.error.code ?? "http_error",
      body?.error.message ?? `${res.status} ${res.statusText}`,
      body?.error.detail,
    );
  }
  return res;
}

async function getJson<R>(path: string): Promise<R> {
  return (await (await request(path)).json()) as R;
}

async function sendJson<R>(method: string, path: string, body?: unknown): Promise<R> {
  const init: RequestInit = { method, headers: { "content-type": "application/json" } };
  if (body !== undefined) init.body = JSON.stringify(body);
  return (await (await request(path, init)).json()) as R;
}

async function getLcb(path: string, init?: RequestInit): Promise<LcbContainer> {
  return decodeLcb(await (await request(path, init)).arrayBuffer());
}

const enc = encodeURIComponent;

export type RootSpec = T.RootRequest["root"];

export const api = {
  dataset: () => getJson<T.DatasetInfo>("/dataset"),
  embedding: (name: string) => getLcb(`/embeddings/${enc(name)}`),
  obs: (column: string) => getLcb(`/obs/${enc(column)}`),
  expression: (genes: string[]) => getLcb(`/expression?genes=${enc(genes.join(","))}`),
  genes: (q: string, limit = 20) => getJson<T.GeneSearchResult>(`/genes?q=${enc(q)}&limit=${limit}`),
  prepStatus: () => getJson<T.PrepStatus>("/prepare"),
  startPrepare: () => sendJson<T.JobInfo>("POST", "/prepare"),
  methods: () => getJson<T.MethodList>("/methods"),
  root: () => getJson<T.RootState>("/root"),
  setRoot: (root: RootSpec) => sendJson<T.ResolvedRoot>("PUT", "/root", { root }),
  resolveRoot: (root: RootSpec) => sendJson<T.ResolvedRoot>("POST", "/root/resolve", { root }),
  rootSuggestions: () => getJson<T.RootSuggestions>("/root/suggestions"),
  createRun: (req: T.TrajectoryRequest) => sendJson<T.RunCreated>("POST", "/trajectories", req),
  runs: () => getJson<T.RunList>("/trajectories"),
  run: (id: string) => getJson<T.TrajectoryRun>(`/trajectories/${enc(id)}`),
  cells: (id: string) => getLcb(`/trajectories/${enc(id)}/cells`),
  curves: (id: string, embedding: string) =>
    getJson<T.Curves>(`/trajectories/${enc(id)}/curves?embedding=${enc(embedding)}`),
  renameLineage: (id: string, lineage: number, name: string) =>
    sendJson<T.Lineage>("PATCH", `/trajectories/${enc(id)}/lineages/${lineage}`, { name }),
  job: (id: string) => getJson<T.JobInfo>(`/jobs/${enc(id)}`),
  cancelJob: (id: string) => sendJson<T.JobInfo>("DELETE", `/jobs/${enc(id)}`),
};

export function jobEventsUrl(jobId: string): string {
  return `${API_BASE}/jobs/${enc(jobId)}/events`;
}
