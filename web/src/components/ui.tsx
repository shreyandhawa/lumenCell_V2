import type { ButtonHTMLAttributes, ReactNode, SelectHTMLAttributes } from "react";

export function Section({
  title,
  aside,
  children,
}: {
  title: string;
  aside?: ReactNode;
  children: ReactNode;
}) {
  return (
    <section className="border-b border-line px-3 py-3" aria-label={title}>
      <header className="mb-2 flex items-center justify-between gap-2">
        <h2 className="text-2xs font-semibold tracking-wide text-ink-muted uppercase">{title}</h2>
        {aside}
      </header>
      {children}
    </section>
  );
}

type Tone = "primary" | "quiet" | "danger";

export function Button({
  tone = "quiet",
  className = "",
  ...props
}: ButtonHTMLAttributes<HTMLButtonElement> & { tone?: Tone }) {
  const tones: Record<Tone, string> = {
    primary: "bg-accent text-accent-ink hover:brightness-110 disabled:opacity-40",
    quiet: "border border-line bg-panel-2 text-ink hover:border-ink-faint disabled:opacity-40",
    danger: "border border-danger/50 text-danger hover:bg-danger/10",
  };
  return (
    <button
      type="button"
      className={`inline-flex h-7 items-center justify-center gap-1.5 rounded-md px-2.5 text-xs font-medium transition-[filter,border-color,background-color] focus-visible:outline-2 focus-visible:outline-offset-1 focus-visible:outline-accent disabled:cursor-not-allowed ${tones[tone]} ${className}`}
      {...props}
    />
  );
}

export function Select({ className = "", ...props }: SelectHTMLAttributes<HTMLSelectElement>) {
  return (
    <select
      className={`h-7 w-full rounded-md border border-line bg-panel-2 px-2 text-xs text-ink focus-visible:outline-2 focus-visible:outline-accent ${className}`}
      {...props}
    />
  );
}

export function Badge({ tone = "neutral", children }: { tone?: "neutral" | "warning" | "accent"; children: ReactNode }) {
  const tones = {
    neutral: "border-line text-ink-muted",
    warning: "border-warning/40 text-warning",
    accent: "border-accent/40 text-accent",
  };
  return (
    <span className={`inline-flex items-center rounded border px-1.5 py-px text-2xs font-medium whitespace-nowrap ${tones[tone]}`}>
      {children}
    </span>
  );
}

export function Hint({ children }: { children: ReactNode }) {
  return <p className="text-2xs leading-relaxed text-ink-faint">{children}</p>;
}

export function fmt(n: number): string {
  return n.toLocaleString("en-US");
}
