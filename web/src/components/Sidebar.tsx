import { Database, Moon, Sun } from "lucide-react";
import { useState } from "react";
import { useDataset, useDoneRunId, useMethods, useRoot, useRun } from "@/api/queries";
import { GeneSearch } from "@/components/GeneSearch";
import { LineageTree } from "@/components/trajectory/LineageTree";
import { PrepProgress, prepCount } from "@/components/trajectory/PrepProgress";
import { describeRoot, RootPicker } from "@/components/trajectory/RootPicker";
import { MethodRunner } from "@/components/trajectory/Runs";
import { fmt, Mark, Select, Step } from "@/components/ui";
import { colorByKey, parseColorBy, useUi } from "@/state/ui";

const LEIDEN = "lumencell_leiden";

/** The steps rail: Prepare, Pick a root, Run, Inspect. */
export function Sidebar() {
  const dataset = useDataset().data;
  const root = useRoot().data?.root ?? null;
  const run = useRun(useDoneRunId(useUi((s) => s.runId))).data;
  const method = useMethods().data?.methods.find((m) => m.name === run?.provenance.method);
  const theme = useUi((s) => s.theme);
  const setTheme = useUi((s) => s.setTheme);
  const colorBy = useUi((s) => s.colorBy);
  const setColorBy = useUi((s) => s.setColorBy);
  const addTrendGene = useUi((s) => s.addTrendGene);
  // A done step stays closed until its header is clicked.
  const [reopen, setReopen] = useState<number | null>(null);
  const toggle = (n: number) => () => setReopen(reopen === n ? null : n);
  if (!dataset) return null;

  const prepDone = dataset.prep.state === "ready";
  const runDone = run !== undefined;
  // The station you are at: the first step not yet done. Earlier ones are done, later ones ahead.
  const current = !prepDone ? 1 : runDone ? 4 : root === null ? 2 : 3;
  const stateOf = (n: number) => (n < current ? "done" : n === current ? "active" : "pending");
  const columns = dataset.obs.filter((c) => c.kind !== "string");

  return (
    <aside aria-label="Steps" className="shell absolute top-3 bottom-3 left-3 z-10 flex w-[280px] flex-col overflow-y-auto rounded-xl bg-panel p-4">
      <header className="flex items-center justify-between">
        <h1 className="text-xl font-semibold tracking-tight">lumenCell</h1>
        <button
          type="button"
          aria-label={theme === "dark" ? "Switch to light theme" : "Switch to dark theme"}
          onClick={() => setTheme(theme === "dark" ? "light" : "dark")}
          className="press grid size-8 place-items-center rounded-md bg-panel-raised text-ink"
        >
          {theme === "dark" ? <Moon size={16} /> : <Sun size={16} />}
        </button>
      </header>
      <div className="mt-5 border-b border-line pb-3">
        <p className="flex items-center gap-2 text-sm font-medium text-ink">
          <Database size={15} aria-hidden className="shrink-0 text-ink-muted" />
          <span className="truncate" title={dataset.name}>
            {dataset.name}
          </span>
        </p>
        <p className="mt-1 font-mono text-xs text-ink-muted">
          {fmt(dataset.n_obs)} cells · {fmt(dataset.n_vars)} genes
        </p>
      </div>

      <div className="mt-2">
        <Step n={1} title="Prepare" state={stateOf(1)} aside={prepDone ? "done" : prepCount(dataset.prep)} open={!prepDone}>
          <PrepProgress prep={dataset.prep} hasFileEmbedding={dataset.embeddings.some((e) => e.source === "file")} />
        </Step>

        <Step
          n={2}
          title="Pick a root"
          state={stateOf(2)}
          aside={!prepDone ? "after prep" : root ? describeRoot(root.spec) : null}
          open={prepDone && (current === 2 || reopen === 2)}
          onToggle={current > 2 ? toggle(2) : undefined}
        >
          <RootPicker />
        </Step>

        <Step
          n={3}
          title="Run"
          state={stateOf(3)}
          aside={runDone ? (method?.label ?? run.provenance.method) : prepDone && root ? "ready" : null}
          open={prepDone && root !== null && (!runDone || reopen === 3 || method?.experimental === true)}
          onToggle={runDone ? toggle(3) : undefined}
        >
          {runDone && reopen !== 3 ? (
            <Mark title={method?.experimental_note ?? undefined}>Experimental</Mark>
          ) : (
            <MethodRunner />
          )}
        </Step>

        <Step
          n={4}
          last
          title="Inspect"
          state={stateOf(4)}
          aside={runDone ? `${run.lineages.length} lineages` : null}
          open={runDone}
        >
          {run && (
            <div className="space-y-5">
              <p className="text-xs leading-snug text-ink-muted">
                Focus a lineage, then play pseudotime below. Evidence says how far to trust it.
              </p>
              <div>
                <label htmlFor="colour-by" className="mb-1.5 block text-xs text-ink-muted">
                  Colour by
                </label>
                <Select id="colour-by" aria-label="Colour by" value={colorByKey(colorBy)} onChange={(e) => setColorBy(parseColorBy(e.target.value))}>
                  <option value="">Nothing</option>
                  <optgroup label="Trajectory">
                    <option value="pseudotime">Pseudotime</option>
                    <option value="lineage">Lineage</option>
                    {run.cell_arrays.includes("extras.entropy") && <option value="entropy">Fate entropy</option>}
                  </optgroup>
                  <optgroup label="Annotations">
                    {columns.map((c) => (
                      <option key={c.name} value={`obs:${c.name}`}>
                        {c.name === LEIDEN ? "Leiden cluster" : c.name}
                      </option>
                    ))}
                  </optgroup>
                  {colorBy?.kind === "gene" && (
                    <optgroup label="Gene">
                      <option value={colorByKey(colorBy)}>{colorBy.gene}</option>
                    </optgroup>
                  )}
                </Select>
                <GeneSearch
                  onPick={(gene) => {
                    setColorBy({ kind: "gene", gene });
                    addTrendGene(gene);
                  }}
                />
              </div>
              <LineageTree run={run} />
            </div>
          )}
        </Step>
      </div>
    </aside>
  );
}
