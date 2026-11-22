/**
 * Terminal states, picked from the labels of the cluster column the method reads
 * (the method's cluster_key). Only labels the dataset has are offered; none picked
 * means the method finds the terminals itself.
 */
export function TerminalPicker({
  column,
  labels,
  value,
  onChange,
}: {
  column: string;
  labels: string[] | null;
  value: string[];
  onChange: (next: string[]) => void;
}) {
  if (labels === null) {
    return (
      <p className="text-xs leading-snug text-ink-muted">
        {column} has no labels in this dataset, so the method finds the terminals itself.
      </p>
    );
  }
  const toggle = (label: string) => {
    const next = new Set(value);
    if (next.has(label)) next.delete(label);
    else next.add(label);
    onChange(labels.filter((l) => next.has(l)));
  };
  return (
    <fieldset>
      <legend className="mb-1 text-xs text-ink-muted">Terminal states, labels of {column}</legend>
      <p className="mb-1.5 text-xs leading-snug text-ink-muted">Leave all unticked and the method finds them.</p>
      <div className="grid max-h-40 grid-cols-2 gap-x-3 gap-y-1 overflow-y-auto">
        {labels.map((label) => (
          <label key={label} className="flex items-center gap-2 text-[13px] text-ink">
            <input type="checkbox" checked={value.includes(label)} onChange={() => toggle(label)} />
            {label}
          </label>
        ))}
      </div>
    </fieldset>
  );
}
