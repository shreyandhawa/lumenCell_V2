import { Combobox } from "@base-ui/react/combobox";
import { Search } from "lucide-react";
import { useState } from "react";
import { useGeneSearch } from "@/api/queries";
import { cx } from "@/components/ui";

/** Type a gene, pick it with the arrow keys or the pointer. Results come from the server. */
export function GeneSearch({
  onPick,
  label = "Colour by gene",
  className,
}: {
  onPick: (gene: string) => void;
  /** The accessible name; the placeholder adds an ellipsis. */
  label?: string;
  className?: string;
}) {
  const [q, setQ] = useState("");
  const query = q.trim();
  const results = useGeneSearch(query);
  const genes = results.data?.genes ?? [];

  return (
    <Combobox.Root
      items={genes}
      filter={null}
      value={null}
      inputValue={q}
      onInputValueChange={(next) => setQ(next)}
      onValueChange={(gene) => {
        if (typeof gene === "string") {
          onPick(gene);
          setQ("");
        }
      }}
    >
      <div className={cx("relative", className)}>
        <Search size={14} className="pointer-events-none absolute top-1/2 left-2.5 -translate-y-1/2 text-ink-3" aria-hidden />
        <Combobox.Input
          aria-label={label}
          placeholder={`${label}…`}
          spellCheck={false}
          autoComplete="off"
          translate="no"
          className="h-8 w-full rounded-lg bg-control pr-2 pl-8 text-sm text-ink transition-colors outline-none placeholder:text-ink-3 [@media(hover:hover)]:hover:bg-control-hover focus-visible:outline-2 focus-visible:outline-offset-0 focus-visible:outline-focus"
        />
      </div>
      <Combobox.Portal>
        <Combobox.Positioner className="z-50 outline-none" sideOffset={6}>
          <Combobox.Popup
            aria-busy={results.isFetching || undefined}
            className="popup w-[var(--anchor-width)] min-w-48 rounded-xl bg-raised p-1 shadow-ply outline-none"
          >
            <Combobox.Empty>
              {query && !results.isFetching ? (
                <div className="px-2.5 py-1.5 text-xs text-ink-3">No gene matches “{query}”.</div>
              ) : null}
            </Combobox.Empty>
            <Combobox.List className="max-h-[min(20rem,var(--available-height))] overflow-y-auto">
              {(gene: string) => (
                <Combobox.Item
                  key={gene}
                  value={gene}
                  translate="no"
                  className="cursor-default rounded-md px-2.5 py-1.5 text-sm text-ink outline-none select-none data-highlighted:bg-focus/15"
                >
                  {gene}
                </Combobox.Item>
              )}
            </Combobox.List>
          </Combobox.Popup>
        </Combobox.Positioner>
      </Combobox.Portal>
    </Combobox.Root>
  );
}
