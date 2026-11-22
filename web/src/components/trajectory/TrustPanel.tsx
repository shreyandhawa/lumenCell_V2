import { Check, LoaderCircle, Minus, PanelRightClose, TriangleAlert, X } from "lucide-react";
import { motion } from "motion/react";
import type * as T from "@/api/generated";
import { useCurves, useDataset, useMethods } from "@/api/queries";
import { describeRoot } from "@/components/trajectory/RootPicker";
import { ReferenceResults } from "@/components/trajectory/ReferenceResults";
import { KeyValues, Label, Mark } from "@/components/ui";
import { axesLabel, axesName } from "@/lib/axes";
import { useUi } from "@/state/ui";

// These labels describe C1 (connectivity) and C5 (clonal support) only, and hold only while
// TAU_TREE, TAU_LIN and TAU_PT are None (C2-C4 excluded); revisit if a tau is ever set.
// The payload has no verified flag, so "supported" reads as provisional (the rule is unverified).
const LEVEL_LABEL: Record<T.TreeSupport["level"], string> = {
  not_supported: "Not supported: cells left without a pseudotime",
  weak: "Clones do not back these lineages",
  provisional: "Not shown to be a tree",
  supported: "Not shown to be a tree",
};

const STATE_ICON = { pass: Check, fail: X, pending: LoaderCircle } as const;

const minus = (s: string) => s.replace("-", "−");
const num = (v: number | null | undefined, digits = 3) => (v == null ? "–" : minus(v.toFixed(digits)));

function list(items: string[]): string {
  return items.length < 2 ? (items[0] ?? "") : `${items.slice(0, -1).join(", ")} and ${items[items.length - 1]}`;
}

/** The verdict in plain words: is it a real failure or a cap, and which checks it reads. */
function explain(ts: T.TreeSupport): string {
  const lead =
    ts.level === "weak" || ts.level === "not_supported"
      ? "A real failure, not a cap."
      : ts.checks.some((c) => c.state === "pending")
        ? "Some checks are still running."
        : ts.checks.some((c) => c.state === "unavailable")
          ? "Some checks could not be measured."
          : "Capped: the verdict rule failed its own check on fresh seeds.";
  const read = ts.checks.filter((c) => c.state !== "excluded" && c.state !== "not_applicable");
  const scope =
    read.length < ts.checks.length
      ? ` The verdict reads ${list(read.map((c) => `${c.name} (${c.id})`))} only, and is not evidence of a tree.`
      : "";
  return lead + scope;
}

function checkLine(c: T.SupportCheck, afterExcluded: boolean, clonal: T.ClonalAgreement | null): string {
  if (c.state === "excluded")
    return afterExcluded ? "Left out for the same reason" : "Left out: did not separate trees from the controls";
  if (c.state === "pending") return "Still being computed";
  if (c.state === "unavailable") return "Could not be measured";
  if (c.id === "C1") return "Share of cells with a pseudotime";
  if (c.id === "C5") {
    const p = clonal?.primary;
    if (c.state === "not_applicable") return c.reason ?? "No clone labels are set";
    if (clonal?.status !== "ok" || p?.p == null || c.threshold == null) return c.reason ?? "Clones could not be scored";
    const meets = p.p <= c.threshold;
    const gain = (p.delta_eta ?? 0) > 0;
    return `p ${p.p.toPrecision(3)} ${meets ? "meets" : "misses"} ${c.threshold}, ${meets === gain ? "and" : "but"} ${
      gain ? "a gain" : "no gain"
    } over the clustering baseline`;
  }
  return c.reason ?? "";
}

function checkValue(c: T.SupportCheck, clonal: T.ClonalAgreement | null): string {
  if (c.id === "C1" && c.value != null && c.threshold != null)
    return `${c.value.toFixed(2)} ${c.state === "pass" ? "≥" : "<"} ${c.threshold.toFixed(2)}`;
  if (c.id === "C5") return clonal?.primary?.delta_eta == null ? "–" : `Δη ${num(clonal.primary.delta_eta)}`;
  return num(c.value);
}

function Checks({ ts, clonal }: { ts: T.TreeSupport; clonal: T.ClonalAgreement | null }) {
  return (
    <ul className="space-y-2" aria-label="Checks">
      {ts.checks.map((c, i) => {
        const Icon = c.state in STATE_ICON ? STATE_ICON[c.state as keyof typeof STATE_ICON] : Minus;
        const afterExcluded = ts.checks.slice(0, i).some((x) => x.state === "excluded");
        return (
          <li key={c.id} className="grid grid-cols-[16px_24px_1fr_auto] gap-x-2 text-[13px]" title={c.reason ?? undefined}>
            <Icon size={15} aria-label={c.state.replace("_", " ")} className={`mt-0.5 text-ink ${c.state === "pending" ? "animate-spin" : ""}`} />
            <span className="mt-px font-mono text-xs text-ink-muted">{c.id}</span>
            <span className="text-ink">{c.name}</span>
            <span className="font-mono text-xs text-ink">{checkValue(c, clonal)}</span>
            <span className="col-start-3 col-end-5 mt-0.5 text-xs leading-snug text-ink-muted">
              {checkLine(c, afterExcluded, clonal)}
            </span>
          </li>
        );
      })}
    </ul>
  );
}

/** η against each null and the k-means baseline, on one line from 0. */
function NullLine({ score }: { score: T.ClonalScore }) {
  const marks = [
    ...score.nulls.flatMap((n) => (n.eta == null ? [] : [{ label: n.null, v: n.eta, base: false }])),
    ...(score.eta_base == null ? [] : [{ label: "k-means", v: score.eta_base, base: true }]),
  ].sort((a, b) => a.v - b.v);
  const max = Math.max(0.5, ...marks.map((m) => m.v));
  const near = (i: number, j: number) => marks[j] !== undefined && Math.abs(marks[i]!.v - marks[j]!.v) / max < 0.08;
  // Close pairs sit side by side; a third close mark drops to a second row; ends stay inside.
  const anchor = (i: number, f: number) =>
    f > 0.92 ? "-100%" : f < 0.08 ? "0%" : near(i, i + 1) ? "-100%" : near(i, i - 1) ? "0%" : "-50%";
  return (
    <div className="relative mb-4 h-9 font-mono text-2xs text-ink-muted" aria-hidden>
      <span className="absolute top-0 left-0">0</span>
      <span className="absolute top-0 right-0">{max === 0.5 ? "0.5" : max.toFixed(2)}</span>
      <div className="absolute top-[13px] right-0 left-0 h-px bg-line" />
      {marks.map((m, i) => {
        const x = `${(m.v / max) * 100}%`;
        const row = near(i, i - 1) && near(i - 1, i - 2) ? 11 : 0;
        return (
          <div key={m.label}>
            {m.base ? (
              <span className="absolute top-[7px] h-3 w-0.5 -translate-x-1/2 bg-ink" style={{ left: x }} />
            ) : (
              <span
                className={`absolute top-[10px] size-[7px] -translate-x-1/2 rounded-full ${m.label === "N1" ? "bg-ink" : "bg-ink-muted"}`}
                style={{ left: x }}
              />
            )}
            <span
              className="absolute px-0.5 whitespace-nowrap"
              style={{ left: x, top: 20 + row, transform: `translateX(${anchor(i, m.v / max)})` }}
            >
              {m.label}
            </span>
          </div>
        );
      })}
    </div>
  );
}

function Clones({ clonal }: { clonal: T.ClonalAgreement | null }) {
  const p = clonal?.primary;
  if (!clonal || clonal.status === "not_configured")
    return <p className="text-xs text-ink-muted">No clone labels are set for this dataset, so clonal support does not apply.</p>;
  if (clonal.status === "pending") return <p className="text-xs text-ink-muted">Clones are still being scored.</p>;
  if (clonal.status !== "ok" || !p) return <p className="text-xs text-ink-muted">Not scored: {clonal.notes[0] ?? clonal.status}</p>;
  const s = clonal.settings;
  const nullWord: Record<string, string> = {
    N1: s?.single_donor ? "one donor" : `within ${s?.donor_key ?? "donor"}`,
    N2: "+ pseudotime decile",
    N3: `+ ${s?.strata_keys.join(", ") ?? "strata"}`,
  };
  return (
    <>
      <NullLine score={p} />
      <KeyValues
        label="Clones against nulls"
        rows={[
          ...p.nulls.map((n): [string, string] => [`η against ${n.null}`, `${num(n.eta)} · ${nullWord[n.null]}`]),
          ["η, k-means clusters", num(p.eta_base)],
          ["Δη, lineages − clusters", `${num(p.delta_eta)} · p ${p.p_delta == null ? "–" : p.p_delta.toPrecision(3)}`],
          ["Clones", `${p.n_clones.toLocaleString("en-US")} (${num(p.effective_clones, 1)} effective)`],
        ]}
      />
    </>
  );
}

/** The Evidence sheet: verdict, checks, clones against nulls, provenance. */
export function TrustPanel({ run }: { run: T.TrajectoryRun }) {
  const embedding = useUi((s) => s.embedding);
  const setSheetOpen = useUi((s) => s.setSheetOpen);
  const curves = useCurves(run.run_id, embedding).data;
  const axes = useUi((s) => s.axes);
  const curves3 = useCurves(run.run_id, axesName(axes)).data;
  const dataset = useDataset().data;
  const fidelity = (cs: T.Curves | undefined, lineage: number) => {
    const f = cs?.curves.find((c) => c.lineage === lineage)?.fidelity;
    return f === null || f === undefined ? "–" : f.toFixed(2);
  };
  const method = useMethods().data?.methods.find((m) => m.name === run.provenance.method);
  const p = run.provenance;
  const ts = run.trust?.tree_support ?? null;
  const clonal = run.trust?.clonal ?? null;
  const notices = [...run.warnings, ...(curves?.warnings ?? [])];
  const rule = "border-t border-line pt-3 mt-3";

  return (
    <motion.aside
      aria-label="Evidence"
      initial={{ opacity: 0, x: 16 }}
      animate={{ opacity: 1, x: 0 }}
      exit={{ opacity: 0, x: 16, transition: { duration: 0.16, ease: [0.23, 1, 0.32, 1] } }}
      transition={{ duration: 0.24, ease: [0.23, 1, 0.32, 1] }}
      className="absolute top-3 right-3 bottom-3 z-10 w-[360px] overflow-y-auto rounded-xl bg-panel px-4 pt-3 pb-4"
    >
      <header className="mb-3 flex items-center justify-between">
        <h2 className="text-base font-semibold">Evidence</h2>
        <button
          type="button"
          aria-label="Hide evidence"
          onClick={() => setSheetOpen(false)}
          className="press grid size-8 place-items-center rounded-md bg-panel-raised text-ink"
        >
          <PanelRightClose size={16} />
        </button>
      </header>

      <section aria-label="Verdict">
        <Label>Verdict</Label>
        {ts ? (
          <>
            <Mark>{LEVEL_LABEL[ts.level]}</Mark>
            <p className="mt-2 text-xs leading-snug text-ink-muted" title={ts.reasons.join("\n")}>
              {explain(ts)}
            </p>
          </>
        ) : (
          <p className="text-xs text-ink-muted">The checks are still running.</p>
        )}
      </section>

      {ts && (
        <section aria-label="Checks" className={rule}>
          <Label>Checks</Label>
          <Checks ts={ts} clonal={clonal} />
        </section>
      )}

      <section aria-label="Clones against nulls" className={rule}>
        <Label>Clones against nulls</Label>
        <Clones clonal={clonal} />
      </section>

      <section aria-label="Reference results" className={rule}>
        <Label>Reference results</Label>
        <ReferenceResults rows={method?.reference_scores ?? []} />
      </section>

      {run.lineages.length > 0 && (
        <section aria-label="Curve fidelity" className={rule}>
          <Label>Curve fidelity</Label>
          <p className="mb-2 text-xs leading-snug text-ink-muted">
            How well each drawn curve keeps its cells' pseudotime order (Spearman, 1 is perfect). A low value means
            the view bends the lineage, not that the lineage is wrong.
          </p>
          <KeyValues
            label="Curve fidelity"
            rows={run.lineages.map((l) => [
              l.name,
              `${embedding ?? "2-D"} ${fidelity(curves, l.id)} · ${axesLabel(axes)} ${fidelity(curves3, l.id)}`,
            ])}
          />
        </section>
      )}

      <section aria-label="Provenance" className={rule}>
        <Label>Provenance</Label>
        <KeyValues
          label="Provenance"
          rows={[
            ["Method", `${method?.label ?? p.method}${method?.experimental ? " (experimental)" : ""}`],
            ["Inferred in", `${p.inference_rep}, ${p.inference_dims} components`],
            ["Root", describeRoot(p.root.spec)],
            ["Data", `${dataset?.name ?? ""} · ${p.dataset_fingerprint.slice(0, 8)}`],
            ["lumenCell", p.lumencell_version],
          ]}
        />
      </section>

      {notices.length > 0 && (
        <section aria-label="Warnings" className={rule}>
          <Label>Warnings</Label>
          <ul className="space-y-1.5">
            {notices.map((w, i) => (
              <li key={`${w.code}-${i}`} className="flex gap-1.5 text-xs leading-relaxed text-ink">
                <TriangleAlert size={14} aria-hidden className="mt-0.5 shrink-0 text-ink" />
                {w.message}
              </li>
            ))}
          </ul>
        </section>
      )}
    </motion.aside>
  );
}
