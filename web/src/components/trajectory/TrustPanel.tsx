import type * as T from "@/api/generated";
import { useCurves, useMethods } from "@/api/queries";
import { cx, Disclosure, Hint, Note, Section } from "@/components/ui";
import { embeddingLabel, LEIDEN } from "@/lib/colorBy";
import { useUi } from "@/state/ui";

const FIDELITY_WARN = 0.8;

/** Everything needed to judge a result: what may be wrong with it first, then where it came from. */
export function TrustPanel({ run }: { run: T.TrajectoryRun }) {
  const embedding = useUi((s) => s.embedding);
  const curves = useCurves(run.run_id, embedding).data;
  const method = useMethods().data?.methods.find((m) => m.name === run.provenance.method);
  const p = run.provenance;
  const notices = [...run.warnings, ...(curves?.warnings ?? [])];

  return (
    <Section title="How far to trust this">
      <div className="space-y-2">
        {method?.experimental && (
          <Note tone="warning">
            {method.label} is experimental. {method.experimental_note}
          </Note>
        )}
        {notices.length > 0 ? (
          <ul className="space-y-2" aria-label="Warnings">
            {notices.map((w, i) => (
              <li key={`${w.code}-${i}`}>
                <Note tone={w.severity}>{w.message}</Note>
              </li>
            ))}
          </ul>
        ) : (
          <Hint>No warnings for this run.</Hint>
        )}
      </div>

      {curves && embedding && (
        <div className="mt-4">
          <h3 className="mb-1.5 text-xs font-medium text-ink-2">Curve fidelity on {embeddingLabel(embedding)}</h3>
          <ul className="space-y-1.5">
            {curves.curves.map((c) => {
              const lin = run.lineages[c.lineage];
              const f = c.fidelity;
              const low = f != null && f < FIDELITY_WARN;
              return (
                <li key={c.lineage} className="grid grid-cols-[1fr_4.5rem_2.5rem] items-center gap-2 text-xs">
                  <span className="truncate text-ink-2">{lin?.name}</span>
                  <span className="h-1 overflow-hidden rounded-full bg-control" aria-hidden>
                    <span
                      className={cx("block h-1 origin-left rounded-full", low ? "bg-warn" : "bg-ink-3")}
                      style={{ transform: `scaleX(${Math.max(0, f ?? 0)})` }}
                    />
                  </span>
                  <span className={cx("num text-right", low ? "text-warn" : "text-ink")}>
                    {f == null ? "n/a" : f.toFixed(2)}
                  </span>
                </li>
              );
            })}
          </ul>
        </div>
      )}

      <dl className="mt-4 grid grid-cols-[auto_1fr] gap-x-4 gap-y-1.5 text-xs">
        <dt className="text-ink-3">Method</dt>
        <dd className="text-ink">
          {method?.label ?? p.method} <span className="num text-ink-3">{p.method_version}</span>
        </dd>
        <dt className="text-ink-3">Root</dt>
        <dd className="text-ink">
          {p.root.cluster_key === LEIDEN ? `Cluster ${p.root.cluster_label}` : (p.root.cluster_label ?? p.root.obs_name)}
          <span className="text-ink-3"> · cell </span>
          <span className="num text-ink-3">{p.root.obs_name}</span>
        </dd>
        <dt className="text-ink-3">Inferred in</dt>
        <dd className="num text-ink">
          {p.inference_rep} · {p.inference_dims} dims
        </dd>
        <dt className="text-ink-3">Drawn on</dt>
        <dd className="text-ink">
          {embedding ? embeddingLabel(embedding) : "…"}{" "}
          <span className="text-ink-3">(projection, {curves?.curve_source ?? "…"} curves)</span>
        </dd>
        <dt className="text-ink-3">Unassigned</dt>
        <dd className="text-ink">
          <span className="num">{run.n_unassigned.toLocaleString()}</span> cells
        </dd>
      </dl>
      <div className="mt-2">
        <Disclosure label="Parameters">
          <p className="num text-xs break-all text-ink-2">
            {Object.entries(p.params)
              .map(([k, v]) => `${k}=${String(v)}`)
              .join("  ") || "defaults"}
          </p>
        </Disclosure>
      </div>
    </Section>
  );
}
