/**
 * The control vocabulary. Every control is a standard one (Base UI underneath, for
 * keyboard and screen-reader behaviour); the world only sets its type, colour and motion.
 */
import { Collapsible } from "@base-ui/react/collapsible";
import { Select as BaseSelect } from "@base-ui/react/select";
import { Toggle } from "@base-ui/react/toggle";
import { ToggleGroup } from "@base-ui/react/toggle-group";
import { Tooltip } from "@base-ui/react/tooltip";
import { Check, ChevronDown, ChevronRight, CircleAlert, Info, TriangleAlert } from "lucide-react";
import { motion } from "motion/react";
import { type ButtonHTMLAttributes, type ReactElement, type ReactNode, useId } from "react";
import { PANEL_SPRING } from "@/lib/motion";

export function cx(...parts: (string | false | null | undefined)[]): string {
  return parts.filter(Boolean).join(" ");
}

export function fmt(n: number): string {
  return n.toLocaleString("en-US");
}

/* ---- Buttons -------------------------------------------------------------------------- */

type Tone = "primary" | "secondary" | "ghost" | "danger";
type Size = "sm" | "md";

const TONES: Record<Tone, string> = {
  primary: "bg-live text-live-ink font-semibold [@media(hover:hover)]:hover:brightness-110 disabled:bg-control disabled:text-ink-3",
  secondary: "bg-control text-ink [@media(hover:hover)]:hover:bg-control-hover disabled:text-ink-3",
  ghost: "text-ink-2 [@media(hover:hover)]:hover:bg-control [@media(hover:hover)]:hover:text-ink disabled:text-ink-3",
  danger: "bg-control text-danger [@media(hover:hover)]:hover:bg-control-hover",
};
const SIZES: Record<Size, string> = {
  sm: "h-7 gap-1.5 rounded-md px-2.5 text-xs",
  md: "h-8 gap-2 rounded-lg px-3 text-sm",
};

export function Button({
  tone = "secondary",
  size = "md",
  className,
  ...props
}: ButtonHTMLAttributes<HTMLButtonElement> & { tone?: Tone; size?: Size }) {
  return (
    <button
      type="button"
      className={cx(
        "pressable inline-flex shrink-0 items-center justify-center font-medium whitespace-nowrap select-none disabled:cursor-not-allowed",
        SIZES[size],
        TONES[tone],
        className,
      )}
      {...props}
    />
  );
}

/** A square icon button with a tooltip; `label` is also its accessible name. */
export function IconButton({
  label,
  shortcut,
  pressed,
  className,
  children,
  ...props
}: ButtonHTMLAttributes<HTMLButtonElement> & { label: string; shortcut?: string; pressed?: boolean }) {
  return (
    <Tip label={label} shortcut={shortcut}>
      <button
        type="button"
        aria-label={label}
        aria-pressed={pressed}
        className={cx(
          "pressable grid size-8 shrink-0 place-items-center rounded-lg",
          pressed ? "bg-focus/15 text-focus" : "text-ink-2 [@media(hover:hover)]:hover:bg-control [@media(hover:hover)]:hover:text-ink",
          className,
        )}
        {...props}
      >
        {children}
      </button>
    </Tip>
  );
}

/* ---- Tooltip -------------------------------------------------------------------------- */

export function TooltipProvider({ children }: { children: ReactNode }) {
  return (
    <Tooltip.Provider delay={500} timeout={500}>
      {children}
    </Tooltip.Provider>
  );
}

export function Tip({ label, shortcut, children }: { label: string; shortcut?: string; children: ReactElement }) {
  return (
    <Tooltip.Root>
      <Tooltip.Trigger render={children} />
      <Tooltip.Portal>
        <Tooltip.Positioner sideOffset={8} className="z-50">
          <Tooltip.Popup className="popup flex items-center gap-2 rounded-md bg-raised px-2 py-1 text-xs text-ink shadow-ply">
            {label}
            {shortcut && <Kbd>{shortcut}</Kbd>}
          </Tooltip.Popup>
        </Tooltip.Positioner>
      </Tooltip.Portal>
    </Tooltip.Root>
  );
}

export function Kbd({ children }: { children: ReactNode }) {
  return (
    <kbd className="rounded border border-line-strong px-1 font-sans text-2xs leading-4 text-ink-2">{children}</kbd>
  );
}

/* ---- Segmented control ---------------------------------------------------------------- */

export interface SegmentOption<T extends string> {
  value: T;
  label: string;
  disabled?: boolean;
  title?: string;
}

/** One-of-n choice with a sliding indicator. Never empty: re-pressing the active segment does nothing. */
export function Segmented<T extends string>({
  label,
  value,
  options,
  onChange,
  className,
}: {
  label: string;
  value: T | null;
  options: SegmentOption<T>[];
  onChange: (v: T) => void;
  className?: string;
}) {
  const id = useId();
  return (
    <ToggleGroup
      aria-label={label}
      value={value ? [value] : []}
      onValueChange={(next) => {
        const v = next[0] as T | undefined;
        if (v && v !== value) onChange(v);
      }}
      className={cx("flex rounded-lg bg-control p-0.5", className)}
    >
      {options.map((o) => {
        const on = o.value === value;
        return (
          <Toggle
            key={o.value}
            value={o.value}
            disabled={o.disabled}
            title={o.title}
            className={cx(
              "pressable relative h-7 rounded-md px-2.5 text-xs font-medium whitespace-nowrap select-none data-disabled:cursor-not-allowed data-disabled:text-ink-3",
              on ? "text-ink" : "text-ink-2 [@media(hover:hover)]:hover:text-ink",
            )}
          >
            {on && (
              <motion.span
                layoutId={`seg-${id}`}
                transition={PANEL_SPRING}
                className="absolute inset-0 rounded-md bg-raised shadow-[0_1px_2px_rgb(0_0_0/0.35)]"
                aria-hidden
              />
            )}
            <span className="relative">{o.label}</span>
          </Toggle>
        );
      })}
    </ToggleGroup>
  );
}

/* ---- Select --------------------------------------------------------------------------- */

export interface SelectItem {
  value: string;
  label: string;
  disabled?: boolean;
  /** Shown faintly after the label, e.g. why an item is disabled. */
  note?: string;
  /** "warn" sets the note in the amber channel (warnings and the experimental flag). */
  noteTone?: "quiet" | "warn";
}
export interface SelectGroup {
  label?: string;
  items: SelectItem[];
}

export function Select({
  label,
  value,
  groups,
  onChange,
  placeholder = "Choose…",
  className,
}: {
  label: string;
  value: string | null;
  groups: SelectGroup[];
  onChange: (value: string) => void;
  placeholder?: string;
  className?: string;
}) {
  const items = Object.fromEntries(groups.flatMap((g) => g.items.map((i) => [i.value, i.label])));
  return (
    <BaseSelect.Root items={items} value={value} onValueChange={(v) => v !== null && onChange(v as string)}>
      <BaseSelect.Trigger
        aria-label={label}
        className={cx(
          "pressable flex h-8 min-w-0 items-center justify-between gap-2 rounded-lg bg-control pr-2 pl-3 text-left text-sm text-ink select-none [@media(hover:hover)]:hover:bg-control-hover data-popup-open:bg-control-hover",
          className,
        )}
      >
        <BaseSelect.Value className="truncate data-placeholder:text-ink-3" placeholder={placeholder} />
        <BaseSelect.Icon className="shrink-0 text-ink-3">
          <ChevronDown size={14} />
        </BaseSelect.Icon>
      </BaseSelect.Trigger>
      <BaseSelect.Portal>
        <BaseSelect.Positioner className="z-50 outline-none" sideOffset={6} alignItemWithTrigger={false}>
          <BaseSelect.Popup className="popup min-w-[var(--anchor-width)] rounded-xl bg-raised p-1 text-sm shadow-ply outline-none">
            <BaseSelect.List className="max-h-[min(24rem,var(--available-height))] overflow-y-auto">
              {groups.map((g, gi) => (
                <BaseSelect.Group key={g.label ?? gi} className={gi > 0 ? "mt-1 border-t border-line pt-1" : ""}>
                  {g.label && (
                    <BaseSelect.GroupLabel className="px-2.5 pt-1.5 pb-1 text-2xs font-semibold text-ink-3">
                      {g.label}
                    </BaseSelect.GroupLabel>
                  )}
                  {g.items.map((it) => (
                    <BaseSelect.Item
                      key={it.value}
                      value={it.value}
                      disabled={it.disabled}
                      className="grid cursor-default grid-cols-[1rem_1fr] items-center gap-2 rounded-md py-1.5 pr-3 pl-2 text-ink outline-none select-none data-disabled:text-ink-3 data-highlighted:bg-focus/15"
                    >
                      <BaseSelect.ItemIndicator className="col-start-1 text-focus">
                        <Check size={14} />
                      </BaseSelect.ItemIndicator>
                      <span className="col-start-2 flex min-w-0 items-baseline gap-2">
                        <BaseSelect.ItemText className="truncate">{it.label}</BaseSelect.ItemText>
                        {it.note && (
                          <span className={cx("truncate text-xs", it.noteTone === "warn" ? "text-warn" : "text-ink-3")}>
                            {it.note}
                          </span>
                        )}
                      </span>
                    </BaseSelect.Item>
                  ))}
                </BaseSelect.Group>
              ))}
            </BaseSelect.List>
          </BaseSelect.Popup>
        </BaseSelect.Positioner>
      </BaseSelect.Portal>
    </BaseSelect.Root>
  );
}

/* ---- Grouping without boxes ------------------------------------------------------------ */

/** A titled group inside a panel. Hierarchy comes from the type scale, not from a border. */
export function Section({
  title,
  aside,
  children,
  className,
}: {
  title: string;
  aside?: ReactNode;
  children: ReactNode;
  className?: string;
}) {
  return (
    <section className={cx("px-4 pt-5 pb-1", className)} aria-label={title}>
      <header className="mb-2.5 flex min-h-7 items-center justify-between gap-2">
        <h2 className="text-md font-semibold tracking-[-0.01em] text-balance text-ink">{title}</h2>
        {aside}
      </header>
      {children}
    </section>
  );
}

/** A disclosure for secondary detail (provenance, parameters). */
export function Disclosure({ label, children }: { label: string; children: ReactNode }) {
  return (
    <Collapsible.Root>
      <Collapsible.Trigger className="pressable group flex h-7 items-center gap-1 rounded-md text-xs font-medium text-ink-2 [@media(hover:hover)]:hover:text-ink">
        <ChevronRight size={14} className="transition-transform duration-150 ease-out group-data-panel-open:rotate-90" />
        {label}
      </Collapsible.Trigger>
      <Collapsible.Panel className="collapse">
        <div className="pt-1 pb-2">{children}</div>
      </Collapsible.Panel>
    </Collapsible.Root>
  );
}

export function Hint({ children, className }: { children: ReactNode; className?: string }) {
  return <p className={cx("text-xs text-ink-3", className)}>{children}</p>;
}

const NOTE = {
  info: { icon: Info, className: "text-ink-2" },
  warning: { icon: TriangleAlert, className: "text-warn" },
  error: { icon: CircleAlert, className: "text-danger" },
} as const;

/** A notice in its state colour: icon plus text, no box, no coloured side border. */
export function Note({
  tone,
  children,
  role,
}: {
  tone: keyof typeof NOTE;
  children: ReactNode;
  role?: "alert" | "status";
}) {
  const { icon: Icon, className } = NOTE[tone];
  return (
    <p role={role} className={cx("flex gap-2 text-xs leading-[1.125rem]", className)}>
      <Icon size={14} className="mt-0.5 shrink-0" aria-hidden />
      <span>{children}</span>
    </p>
  );
}

/** A small state marker; colour only where it carries meaning. */
export function Tag({ tone = "neutral", children }: { tone?: "neutral" | "live" | "warn" | "focus"; children: ReactNode }) {
  const tones = {
    neutral: "bg-control text-ink-2",
    live: "bg-live/15 text-live",
    warn: "bg-warn/15 text-warn",
    focus: "bg-focus/15 text-focus",
  };
  return (
    <span className={cx("inline-flex h-5 items-center rounded-full px-2 text-2xs font-semibold whitespace-nowrap", tones[tone])}>
      {children}
    </span>
  );
}
