import { Check, ChevronDown, TriangleAlert } from "lucide-react";
import type { ButtonHTMLAttributes, ReactNode, SelectHTMLAttributes } from "react";

type StepState = "active" | "done" | "pending";

/** The line's stretch already travelled, and the stretch ahead. */
const travelled = "bg-ink-muted";
const ahead = "bg-ink-faint";

/**
 * One step of the rail: Prepare, Pick a root, Run, Inspect, drawn as stations on one line. Done
 * stations are filled, the current one is the red "you are here", the rest are hollow. A done
 * step's header reopens it.
 */
export function Step({
  n,
  title,
  state,
  aside,
  open,
  onToggle,
  last = false,
  children,
}: {
  n: number;
  title: string;
  state: StepState;
  aside?: ReactNode;
  open: boolean;
  onToggle?: () => void;
  /** The last station: the line ends here. */
  last?: boolean;
  children?: ReactNode;
}) {
  const active = state === "active";
  const head = (
    <>
      <span
        aria-hidden
        className={`relative grid size-[22px] shrink-0 place-items-center rounded-full text-2xs font-semibold ${
          state === "done"
            ? "bg-ink text-panel"
            : active
              ? "bg-accent text-on-accent ring-2 ring-ink"
              : "border-2 border-ink-faint bg-panel text-ink-muted"
        }`}
      >
        {state === "done" ? <Check size={13} strokeWidth={3} /> : n}
      </span>
      <span className={`shrink-0 ${active ? "text-base font-semibold text-ink" : "text-sm font-medium text-ink"}`}>{title}</span>
      {aside != null && <span className="ml-auto min-w-0 pl-2 text-right text-xs leading-tight text-balance text-ink-muted">{aside}</span>}
    </>
  );
  return (
    <section aria-label={title} aria-current={active ? "step" : undefined} className="relative py-2.5">
      {/* The line through the stations, drawn under them: in from the last station, out to the next. */}
      {n > 1 && <span aria-hidden className={`absolute top-0 left-[9.5px] h-[22px] w-[3px] ${state === "pending" ? ahead : travelled}`} />}
      {!last && <span aria-hidden className={`absolute top-[22px] bottom-0 left-[9.5px] w-[3px] ${state === "done" ? travelled : ahead}`} />}
      {onToggle ? (
        <button type="button" aria-expanded={open} onClick={onToggle} className="flex w-full items-center gap-2.5 text-left">
          {head}
        </button>
      ) : (
        <div className="flex items-center gap-2.5">{head}</div>
      )}
      {open && children && <div className="mt-3 pl-8">{children}</div>}
    </section>
  );
}

type Tone = "primary" | "secondary" | "ghost";

export function Button({
  tone = "secondary",
  className = "",
  ...props
}: ButtonHTMLAttributes<HTMLButtonElement> & { tone?: Tone }) {
  const tones: Record<Tone, string> = {
    // Ink, as on a station sign (white on the blue line): the red is kept for "you are here".
    primary: "bg-ink px-3 text-panel hover:opacity-90 disabled:bg-panel-raised",
    secondary: "bg-panel-raised px-3 text-ink hover:bg-line",
    ghost: "px-2 text-ink-muted hover:text-ink",
  };
  return (
    <button
      type="button"
      className={`press inline-flex h-8 items-center justify-center gap-1.5 rounded text-sm font-medium focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-ink disabled:cursor-not-allowed disabled:text-ink-faint ${tones[tone]} ${className}`}
      {...props}
    />
  );
}

export function Select({ className = "", ...props }: SelectHTMLAttributes<HTMLSelectElement>) {
  return (
    <div className={`relative ${className}`}>
      <select
        className="h-9 w-full appearance-none rounded-md bg-panel-raised pr-8 pl-2.5 text-sm text-ink focus-visible:outline-2 focus-visible:outline-ink"
        {...props}
      />
      <ChevronDown size={16} aria-hidden className="pointer-events-none absolute top-2.5 right-2.5 text-ink-muted" />
    </div>
  );
}

/** The amber mark: experimental, provisional and weak, always with its icon and a label. */
export function Mark({ children, title }: { children: ReactNode; title?: string }) {
  return (
    <span
      title={title}
      className="inline-flex items-center gap-1.5 rounded border border-amber bg-amber-tint px-1.5 py-0.5 text-xs text-ink"
    >
      <TriangleAlert size={14} aria-hidden className="shrink-0 text-amber" />
      {children}
    </span>
  );
}

export function Label({ children }: { children: ReactNode }) {
  return <h3 className="mb-2 text-2xs font-semibold tracking-[0.06em] text-ink-muted uppercase">{children}</h3>;
}

/** Key-value rows: a muted key, a mono value. */
export function KeyValues({ rows, label }: { rows: [string, ReactNode][]; label: string }) {
  return (
    <dl aria-label={label} className="grid grid-cols-[120px_1fr] gap-x-3 gap-y-1.5 text-xs">
      {rows.map(([k, v]) => (
        <div key={k} className="contents">
          <dt className="text-ink-muted">{k}</dt>
          <dd className="font-mono break-words text-ink">{v}</dd>
        </div>
      ))}
    </dl>
  );
}

export function fmt(n: number): string {
  return n.toLocaleString("en-US");
}
