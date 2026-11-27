import { Search } from "lucide-react";
import { useState } from "react";
import { useGeneSearch } from "@/api/queries";

export function GeneSearch({
  onPick,
  placeholder = "Colour by gene…",
}: {
  onPick: (gene: string) => void;
  placeholder?: string;
}) {
  const [q, setQ] = useState("");
  const results = useGeneSearch(q.trim());
  return (
    <div className="relative mt-2">
      <Search size={14} aria-hidden className="pointer-events-none absolute top-2.5 left-2.5 text-ink-muted" />
      <input
        type="search"
        aria-label={placeholder}
        placeholder={placeholder}
        value={q}
        onChange={(e) => setQ(e.target.value)}
        onKeyDown={(e) => {
          const first = results.data?.genes[0];
          if (e.key === "Enter" && first) {
            onPick(first);
            setQ("");
          }
        }}
        className="h-9 w-full rounded-md bg-panel-raised pr-2 pl-8 text-sm text-ink placeholder:text-ink-muted focus-visible:outline-2 focus-visible:outline-ink"
      />
      {q && results.data && (
        <ul role="listbox" aria-label="Genes" className="absolute z-10 mt-1 max-h-56 w-full overflow-y-auto rounded-md border border-line bg-panel-raised py-1 shadow-xl">
          {results.data.genes.length === 0 && <li className="px-2 py-1 text-xs text-ink-muted">no match</li>}
          {results.data.genes.map((g) => (
            <li key={g}>
              <button
                type="button"
                role="option"
                aria-selected={false}
                className="w-full px-2.5 py-1 text-left font-mono text-xs text-ink hover:bg-line"
                onClick={() => {
                  onPick(g);
                  setQ("");
                }}
              >
                {g}
              </button>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
