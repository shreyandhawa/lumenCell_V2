import { Info, ShieldAlert, TriangleAlert } from "lucide-react";
import type * as T from "@/api/generated";
import { useCurves, useMethods } from "@/api/queries";
import { Hint, Section } from "@/components/ui";
import { useUi } from "@/state/ui";

const ICON = { info: Info, warning: TriangleAlert, error: ShieldAlert } as const;
const TONE = { info: "text-ink-muted", warning: "text-warning", error: "text-danger" } as const;

/** Everything needed to judge a result: where it came from, and what may be wrong with it. */
export function TrustPanel({ run }: { run: T.TrajectoryRun }) {
  const embedding = useUi((s) => s.embedding);
  const curves = useCurves(run.run_id, embedding).data;
  const method = useMethods().data?.methods.find((m) => m.name === run.provenance.method);
  const p = run.provenance;
  const notices = [...run.warnings, ...(curves?.warnings ?? [])];

  return (
    <Section title="How far to trust this">
      <dl className="grid grid-cols-[auto_1fr] gap-x-3 gap-y-1 text-2xs">
        <dt className="text-ink-faint">method</dt>
        <dd className="text-ink">
          {method?.label ?? p.method} <span className="text-ink-faint">({p.method_version})</span>
        </dd>
        <dt className="text-ink-faint">root</dt>
        <dd className="text-ink">
          {p.root.cluster_key === "lumencell_leiden" ? `Cluster ${p.root.cluster_label}` : (p.root.cluster_label ?? p.root.obs_name)}
          <span className="text-ink-faint"> · cell {p.root.obs_name}</span>
        </dd>
        <dt className="text-ink-faint">inferred in</dt>
        <dd className="font-mono text-ink">
          {p.inference_rep} · {p.inference_dims} dims
        </dd>
        <dt className="text-ink-faint">drawn on</dt>
        <dd className="text-ink">
          {embedding} <span className="text-ink-faint">(projection, {curves?.curve_source ?? "…"})</span>
        </dd>
        <dt className="text-ink-faint">params</dt>
        <dd className="font-mono break-all text-ink">
          {Object.entries(p.params)
            .map(([k, v]) => `${k}=${String(v)}`)
            .join(" ")}
        </dd>
        <dt className="text-ink-faint">unassigned</dt>
        <dd className="text-ink">{run.n_unassigned.toLocaleString()} cells</dd>
      </dl>

      {curves && (
        <div className="mt-2">
          <div className="mb-1 text-2xs text-ink-faint">curve fidelity on this embedding</div>
          {curves.curves.map((c) => {
            const lin = run.lineages[c.lineage];
            const f = c.fidelity;
            return (
              <div key={c.lineage} className="flex items-center justify-between text-2xs">
                <span className="truncate text-ink-muted">{lin?.name}</span>
                <span className={`font-mono ${f != null && f < 0.8 ? "text-warning" : "text-ink"}`}>
                  {f == null ? "n/a" : f.toFixed(2)}
                </span>
              </div>
            );
          })}
        </div>
      )}

      {method?.experimental && (
        <p className="mt-2 flex gap-1.5 text-2xs text-warning">
          <TriangleAlert size={12} className="mt-px shrink-0" />
          {method.label} is experimental: on lumenCell's synthetic benchmark it misses some accuracy gates
          (branch assignment near forks, cross-lineage timing). Pseudotime within each lineage scored 0.95+.
        </p>
      )}

      {notices.length > 0 ? (
        <ul className="mt-2 space-y-1.5" aria-label="Warnings">
          {notices.map((w, i) => {
            const Icon = ICON[w.severity];
            return (
              <li key={`${w.code}-${i}`} className={`flex gap-1.5 text-2xs leading-relaxed ${TONE[w.severity]}`}>
                <Icon size={12} className="mt-0.5 shrink-0" />
                <span>{w.message}</span>
              </li>
            );
          })}
        </ul>
      ) : (
        <div className="mt-2">
          <Hint>No warnings for this run.</Hint>
        </div>
      )}
    </Section>
  );
}
