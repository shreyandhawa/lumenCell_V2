import { Select } from "@/components/ui";

/**
 * The part of JSON Schema that method parameters use (see engine/methods/*.py):
 * integers, numbers, strings, enums, and an optional value written as anyOf [T, null].
 */
export interface Prop {
  title?: string;
  description?: string;
  type?: string;
  enum?: string[];
  default?: unknown;
  minimum?: number;
  maximum?: number;
  exclusiveMinimum?: number;
  anyOf?: Prop[];
}

export interface ParamSchema {
  properties?: Record<string, Prop>;
}

export type Values = Record<string, string>;

/** The field's own type, and whether null is allowed (an empty field sends null: the method's rule). */
function shape(prop: Prop): { base: Prop; nullable: boolean } {
  const options = prop.anyOf ?? [];
  const option = options.find((o) => o.type !== "null");
  return { base: option ? { ...prop, ...option } : prop, nullable: options.some((o) => o.type === "null") };
}

function bounds(p: Prop): string | undefined {
  if (p.minimum !== undefined && p.maximum !== undefined) return `${p.minimum} to ${p.maximum}`;
  if (p.exclusiveMinimum !== undefined) return `more than ${p.exclusiveMinimum}`;
  if (p.minimum !== undefined) return `at least ${p.minimum}`;
  return undefined;
}

/** Each field starts at its schema default, as text, so the inputs can hold what the user is typing. */
export function initialValues(schema: ParamSchema): Values {
  return Object.fromEntries(Object.entries(schema.properties ?? {}).map(([name, p]) => [name, p.default == null ? "" : String(p.default)]));
}

/** Why this text can't be sent, or null when it can. */
export function fieldError(prop: Prop, raw: string): string | null {
  const { base, nullable } = shape(prop);
  if (raw.trim() === "") return nullable ? null : "Enter a value.";
  if (base.enum) return base.enum.includes(raw) ? null : "Choose one of the listed values.";
  if (base.type === "string") return null;
  const n = Number(raw);
  if (!Number.isFinite(n)) return "Enter a number.";
  if (base.type === "integer" && !Number.isInteger(n)) return "Enter a whole number.";
  if (base.minimum !== undefined && n < base.minimum) return `Must be at least ${base.minimum}.`;
  if (base.maximum !== undefined && n > base.maximum) return `Must be at most ${base.maximum}.`;
  if (base.exclusiveMinimum !== undefined && n <= base.exclusiveMinimum) return `Must be more than ${base.exclusiveMinimum}.`;
  return null;
}

/** The parameters to send, or the error of each field that cannot be sent. */
export function checkParams(
  schema: ParamSchema,
  values: Values,
): { params: Record<string, string | number | null>; errors: Record<string, string> } {
  const params: Record<string, string | number | null> = {};
  const errors: Record<string, string> = {};
  for (const [name, prop] of Object.entries(schema.properties ?? {})) {
    const raw = values[name] ?? "";
    const error = fieldError(prop, raw);
    if (error) {
      errors[name] = error;
    } else if (raw.trim() === "") {
      params[name] = null;
    } else {
      const { base } = shape(prop);
      params[name] = base.enum || base.type === "string" ? raw : Number(raw);
    }
  }
  return { params, errors };
}

const INPUT = "h-9 w-full rounded-md bg-panel-raised px-2.5 text-sm text-ink focus-visible:outline-2 focus-visible:outline-ink";

/** One input per parameter, typed from its schema. Errors are shown under the field they belong to. */
export function ParamForm({
  schema,
  values,
  errors,
  onChange,
}: {
  schema: ParamSchema;
  values: Values;
  errors: Record<string, string>;
  onChange: (name: string, raw: string) => void;
}) {
  return (
    <div className="space-y-3">
      {Object.entries(schema.properties ?? {}).map(([name, prop]) => {
        const { base } = shape(prop);
        const id = `param-${name}`;
        const error = errors[name];
        const help = [prop.description, bounds(base)].filter(Boolean).join(" · ");
        return (
          <div key={name}>
            <label htmlFor={id} className="mb-1 block text-xs text-ink-muted">
              {prop.title ?? name}
            </label>
            {base.enum ? (
              <Select
                id={id}
                value={values[name] ?? ""}
                aria-invalid={error !== undefined}
                onChange={(e) => onChange(name, e.target.value)}
              >
                {base.enum.map((v) => (
                  <option key={v} value={v}>
                    {v}
                  </option>
                ))}
              </Select>
            ) : (
              <input
                id={id}
                type={base.type === "string" ? "text" : "number"}
                step={base.type === "integer" ? 1 : "any"}
                min={base.minimum}
                max={base.maximum}
                className={INPUT}
                value={values[name] ?? ""}
                aria-invalid={error !== undefined}
                aria-describedby={help ? `${id}-help` : undefined}
                onChange={(e) => onChange(name, e.target.value)}
              />
            )}
            {help && (
              <p id={`${id}-help`} className="mt-1 text-xs leading-snug text-ink-muted">
                {help}
              </p>
            )}
            {error && <p className="mt-1 text-xs text-ink">{error}</p>}
          </div>
        );
      })}
    </div>
  );
}
