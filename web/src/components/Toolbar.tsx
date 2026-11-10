import { Moon, PanelLeft, PanelRight, Sun } from "lucide-react";
import { useDataset, useDoneRunId } from "@/api/queries";
import { GeneSearch } from "@/components/GeneSearch";
import { cx, fmt, IconButton, Segmented, Select, type SelectGroup } from "@/components/ui";
import { columnLabel, decodeColorBy, embeddingLabel, encodeColorBy, NONE } from "@/lib/colorBy";
import { useUi } from "@/state/ui";

const STAGE_LABELS: Record<string, string> = {
  normalize: "normalising counts",
  hvg: "selecting variable genes",
  pca: "PCA",
  neighbors: "neighbour graph",
  leiden: "Leiden clusters",
  diffmap: "diffusion map",
  paga: "PAGA graph",
  umap: "UMAP",
};

const EMBEDDING_HINT: Record<string, string> = {
  file: "From the file, for display only. Trajectories are inferred in diffusion space and drawn here as a projection.",
  computed: "Computed by lumenCell, for display only. Trajectories are inferred in diffusion space and drawn here as a projection.",
  diffusion: "Diffusion components: the space the trajectory is inferred in.",
};

export function Toolbar() {
  const dataset = useDataset().data;
  const panels = useUi((s) => s.panels);
  const togglePanel = useUi((s) => s.togglePanel);
  const theme = useUi((s) => s.theme);
  const setTheme = useUi((s) => s.setTheme);

  return (
    <header className="relative z-20 flex min-h-12 shrink-0 flex-wrap items-center gap-x-3 gap-y-2 border-b border-line bg-panel px-3 py-2 md:flex-nowrap md:py-0">
      <h1 translate="no" className="pl-1 text-md font-semibold tracking-[-0.02em] text-ink">
        lumenCell
      </h1>
      {dataset && (
        <div className="flex min-w-0 items-baseline gap-2 max-lg:hidden">
          <span className="truncate text-sm text-ink-2" title={dataset.name}>
            {dataset.name}
          </span>
          <span className="num shrink-0 text-xs text-ink-3">
            {fmt(dataset.n_obs)} cells × {fmt(dataset.n_vars)} genes
          </span>
        </div>
      )}

      <div className="order-last flex w-full min-w-0 items-center gap-2 md:order-none md:mx-auto md:w-auto">
        <EmbeddingSwitch />
        <ColourBy />
        <GeneSearch
          className="w-44 max-md:hidden"
          onPick={(gene) => {
            useUi.getState().setColorBy({ kind: "gene", gene });
            useUi.getState().addTrendGene(gene);
          }}
        />
      </div>

      <PrepStatus />
      <div className="ml-auto flex items-center gap-0.5 md:ml-0">
        <IconButton label="Cells panel" shortcut="[" pressed={panels.cells} onClick={() => togglePanel("cells")}>
          <PanelLeft size={16} />
        </IconButton>
        <IconButton
          label="Workflow panel"
          shortcut="]"
          pressed={panels.inspector}
          onClick={() => togglePanel("inspector")}
        >
          <PanelRight size={16} />
        </IconButton>
        <IconButton
          label={theme === "dark" ? "Switch to brightfield" : "Switch to darkfield"}
          onClick={() => setTheme(theme === "dark" ? "light" : "dark")}
        >
          {theme === "dark" ? <Sun size={16} /> : <Moon size={16} />}
        </IconButton>
      </div>
    </header>
  );
}

function EmbeddingSwitch() {
  const dataset = useDataset().data;
  const embedding = useUi((s) => s.embedding);
  const setEmbedding = useUi((s) => s.setEmbedding);
  const embeddings = dataset?.embeddings ?? [];
  if (embeddings.length === 0) {
    return <div className="skeleton h-8 w-44" aria-label="Embeddings appear once preprocessing has run" />;
  }
  if (embeddings.length > 5) {
    return (
      <Select
        label="Embedding"
        className="w-36"
        value={embedding}
        onChange={setEmbedding}
        groups={[{ items: embeddings.map((e) => ({ value: e.name, label: embeddingLabel(e.name) })) }]}
      />
    );
  }
  return (
    <Segmented
      label="Embedding"
      value={embedding}
      onChange={setEmbedding}
      options={embeddings.map((e) => ({ value: e.name, label: embeddingLabel(e.name), title: EMBEDDING_HINT[e.source] }))}
    />
  );
}

function ColourBy() {
  const dataset = useDataset().data;
  const colorBy = useUi((s) => s.colorBy);
  const setColorBy = useUi((s) => s.setColorBy);
  const runId = useDoneRunId(useUi((s) => s.runId));
  if (!dataset) return null;

  const groups: SelectGroup[] = [
    { items: [{ value: NONE, label: "No colour" }] },
    {
      label: "Trajectory",
      items: [
        { value: "pseudotime", label: "Pseudotime", disabled: !runId, note: runId ? undefined : "after a run" },
        { value: "lineage", label: "Lineage", disabled: !runId, note: runId ? undefined : "after a run" },
      ],
    },
    {
      label: "Annotations",
      items: dataset.obs
        .filter((c) => c.kind !== "string")
        .map((c) => ({ value: `obs:${c.name}`, label: columnLabel(c.name) })),
    },
  ];
  if (colorBy?.kind === "gene") {
    groups.push({ label: "Gene", items: [{ value: encodeColorBy(colorBy), label: colorBy.gene }] });
  }
  if (colorBy?.kind === "disagreement") {
    groups.push({ label: "Comparison", items: [{ value: encodeColorBy(colorBy), label: "Disagreement between runs" }] });
  }
  return (
    <Select
      label="Colour by"
      className="w-48 max-md:min-w-0 max-md:flex-1"
      value={encodeColorBy(colorBy)}
      groups={groups}
      onChange={(v) => setColorBy(decodeColorBy(v))}
    />
  );
}

function PrepStatus() {
  const prep = useDataset().data?.prep;
  if (!prep || prep.state === "ready") return null;
  const running = prep.stages.find((s) => s.state === "running");
  const failed = prep.state === "failed";
  return (
    <span
      role="status"
      className={cx("flex shrink-0 items-center gap-2 text-xs", failed ? "text-warn" : "text-ink-2")}
    >
      <span className={cx("size-1.5 rounded-full", failed ? "bg-warn" : "animate-pulse bg-live")} aria-hidden />
      {failed ? "Preparation failed" : `Preparing${running ? `: ${STAGE_LABELS[running.name] ?? running.name}` : "…"}`}
    </span>
  );
}
