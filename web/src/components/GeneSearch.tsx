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
      <Search size={13} className="pointer-events-none absolute top-2 left-2 text-ink-faint" />
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
        className="h-7 w-full rounded-md border border-line bg-panel-2 pr-2 pl-7 text-xs text-ink placeholder:text-ink-faint focus-visible:outline-2 focus-visible:outline-accent"
      />
      {q && results.data && (
        <ul role="listbox" aria-label="Genes" className="absolute z-10 mt-1 max-h-56 w-full overflow-y-auto rounded-md border border-line bg-panel-2 py-1 shadow-xl">
          {results.data.genes.length === 0 && <li className="px-2 py-1 text-2xs text-ink-faint">no match</li>}
          {results.data.genes.map((g) => (
            <li key={g}>
              <button
                type="button"
                role="option"
                aria-selected={false}
                className="w-full px-2 py-1 text-left font-mono text-xs text-ink hover:bg-accent/15"
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
