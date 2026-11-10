import { useDataset } from "@/api/queries";
import { Badge, fmt } from "@/components/ui";

export function Header() {
  const dataset = useDataset().data;
  const prep = dataset?.prep.state;
  return (
    <header className="flex h-11 shrink-0 items-center gap-3 border-b border-line bg-panel px-4">
      <h1 className="text-sm font-semibold tracking-tight">
        lumen<span className="text-accent">Cell</span>
      </h1>
      {dataset && (
        <>
          <span className="h-4 w-px bg-line" aria-hidden />
          <span className="truncate text-sm text-ink" title={dataset.name}>
            {dataset.name}
          </span>
          <span className="font-mono text-2xs text-ink-faint">
            {fmt(dataset.n_obs)} cells × {fmt(dataset.n_vars)} genes
          </span>
          <span className="ml-auto">
            {prep === "ready" ? (
              <Badge>prepared</Badge>
            ) : prep === "failed" ? (
              <Badge tone="warning">prep failed</Badge>
            ) : (
              <Badge tone="accent">preparing…</Badge>
            )}
          </span>
        </>
      )}
    </header>
  );
}
