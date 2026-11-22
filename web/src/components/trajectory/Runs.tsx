import { LoaderCircle, TriangleAlert } from "lucide-react";
import { useState } from "react";
import { ApiError } from "@/api/client";
import { useJobEvents } from "@/api/events";
import { useCreateRun, useDataset, useMethods, useRuns } from "@/api/queries";
import { Button, Mark, Select } from "@/components/ui";
import { checkParams, initialValues, ParamForm, type ParamSchema, type Values } from "@/components/trajectory/ParamForm";
import { TerminalPicker } from "@/components/trajectory/TerminalPicker";
import { useUi } from "@/state/ui";

/** What the user has entered for one method. Kept per method, so switching back restores it. */
interface Form {
  values: Values;
  terminals: string[];
}

/** The cluster column the terminal names refer to when a method has no cluster_key of its own. */
const TERMINAL_COLUMN = "lumencell_leiden";

export function MethodRunner() {
  const methods = useMethods().data?.methods ?? [];
  const dataset = useDataset().data;
  const create = useCreateRun();
  const runId = useUi((s) => s.runId);
  const setRun = useUi((s) => s.setRun);
  const setColorBy = useUi((s) => s.setColorBy);
  const [method, setMethod] = useState<string>("");
  const chosen = methods.find((m) => m.name === (method || methods.find((x) => x.available)?.name));
  const [forms, setForms] = useState<Record<string, Form>>({});
  const [jobId, setJobId] = useState<string | null>(null);
  const live = useJobEvents(jobId);
  const failed = useRuns().data?.runs.find((r) => r.run_id === runId && r.state === "failed");

  const schema = (chosen?.params_schema ?? {}) as ParamSchema;
  const props = schema.properties ?? {};
  const form: Form = forms[chosen?.name ?? ""] ?? { values: initialValues(schema), terminals: [] };
  const setForm = (patch: Partial<Form>) =>
    chosen && setForms((all) => ({ ...all, [chosen.name]: { ...(all[chosen.name] ?? form), ...patch } }));
  const { params, errors } = checkParams(schema, form.values);
  const invalid = Object.keys(errors).length > 0;

  // Terminal names are labels of the method's cluster column; a label the column lacks is dropped.
  const column = form.values.cluster_key || TERMINAL_COLUMN;
  const labels = dataset?.obs.find((c) => c.name === column && c.kind === "categorical")?.categories ?? null;
  const terminals = labels ? form.terminals.filter((t) => labels.includes(t)) : [];

  // Validation errors from the server, under their field; any other error is listed below the form.
  const serverErrors = create.error instanceof ApiError && create.error.code === "invalid_params" ? create.error.detail : null;
  const detailErrors = ((serverErrors?.errors ?? []) as { loc: string[]; msg: string }[]).filter((e) => e.loc.length > 0);
  const byField = Object.fromEntries(detailErrors.filter((e) => props[e.loc[0]!]).map((e) => [e.loc[0], e.msg]));
  const otherErrors = detailErrors.filter((e) => !props[e.loc[0]!]).map((e) => e.msg);

  const run = () => {
    if (!chosen || invalid) return;
    create.mutate(
      { method: chosen.name, params, terminals: terminals.length > 0 ? terminals : null },
      {
        onSuccess: (created) => {
          setJobId(created.job_id ?? null);
          setRun(created.run_id);
          setColorBy({ kind: "pseudotime" });
        },
      },
    );
  };

  return (
    <div className="space-y-3">
      <div>
        <label htmlFor="method" className="mb-1.5 block text-xs text-ink-muted">
          Method
        </label>
        <Select id="method" aria-label="Method" value={chosen?.name ?? ""} onChange={(e) => setMethod(e.target.value)}>
          {methods.map((m) => (
            <option key={m.name} value={m.name} disabled={!m.available}>
              {m.label}
              {m.available ? "" : " (not installed)"}
            </option>
          ))}
        </Select>
      </div>
      {chosen?.experimental && <Mark title={chosen.experimental_note ?? undefined}>Experimental</Mark>}
      {chosen && !chosen.available && chosen.install_hint && <p className="text-xs text-ink-muted">{chosen.install_hint}</p>}
      {chosen?.available && "cluster_key" in props && (
        <TerminalPicker column={column} labels={labels} value={terminals} onChange={(next) => setForm({ terminals: next })} />
      )}
      {chosen?.available && Object.keys(props).length > 0 && (
        <ParamForm
          schema={schema}
          values={form.values}
          errors={{ ...errors, ...byField }}
          onChange={(name, raw) => setForm({ values: { ...form.values, [name]: raw } })}
        />
      )}
      {invalid && <p className="text-xs text-ink-muted">Fix the values above to run.</p>}
      <div>
        <Button tone="primary" onClick={run} disabled={!chosen?.available || create.isPending || invalid}>
          Run {chosen?.label ?? ""}
        </Button>
      </div>
      {live && live.state === "running" && (
        <div aria-label="Run progress">
          <div className="mb-1 flex items-center gap-1.5 text-xs text-ink-muted">
            <LoaderCircle size={14} aria-hidden className="animate-spin text-ink" /> {live.stage ?? "starting"}
          </div>
          <div className="h-1 rounded-full bg-panel-raised">
            <div className="h-1 rounded-full bg-ink-muted" style={{ width: `${live.fraction * 100}%` }} />
          </div>
        </div>
      )}
      {(create.error || failed) && (
        <p role="alert" className="flex gap-1.5 text-xs text-ink">
          <TriangleAlert size={14} aria-hidden className="mt-px shrink-0 text-ink" />
          {[create.error?.message ?? failed?.error ?? "The run failed.", ...otherErrors].join(" ")}
        </p>
      )}
    </div>
  );
}
