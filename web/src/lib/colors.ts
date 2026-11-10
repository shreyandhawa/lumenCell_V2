/**
 * Per-cell colour buffers (RGBA, Uint8, 4 bytes per cell) handed to deck.gl as
 * binary attributes, so recolouring 100k cells never allocates JS objects.
 */

export type RGBA = readonly [number, number, number, number];

/** Distinct, colour-blind-aware categorical palette (Okabe-Ito first, then Tableau). */
export const CATEGORICAL = [
  "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#7F3C8D",
  "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2", "#EDC948", "#B07AA1", "#FF9DA7",
  "#9C755F", "#BAB0AC", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#A6761D",
] as const;

/** Cells with no value (NaN, missing category, no pseudotime). */
export const MISSING: RGBA = [128, 128, 140, 70];
/** Cells on the shared trunk (more than one lineage). */
export const SHARED: RGBA = [170, 170, 182, 200];

// Viridis, sampled at 10 stops: perceptually uniform and readable in grey.
const VIRIDIS = [
  "#440154", "#482878", "#3e4989", "#31688e", "#26828e",
  "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725",
].map(hexToRgb);

export function hexToRgb(hex: string): [number, number, number] {
  const h = hex.replace("#", "");
  const v = Number.parseInt(h.length === 3 ? [...h].map((c) => c + c).join("") : h.slice(0, 6), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

export function sequential(t: number, alpha = 230): RGBA {
  const x = Math.min(1, Math.max(0, t)) * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x));
  const f = x - i;
  const a = VIRIDIS[i]!;
  const b = VIRIDIS[i + 1]!;
  return [
    Math.round(a[0] + (b[0] - a[0]) * f),
    Math.round(a[1] + (b[1] - a[1]) * f),
    Math.round(a[2] + (b[2] - a[2]) * f),
    alpha,
  ];
}

export function categorical(i: number, alpha = 220): RGBA {
  const [r, g, b] = hexToRgb(CATEGORICAL[i % CATEGORICAL.length]!);
  return [r, g, b, alpha];
}

function write(out: Uint8Array, i: number, c: RGBA): void {
  out[i * 4] = c[0];
  out[i * 4 + 1] = c[1];
  out[i * 4 + 2] = c[2];
  out[i * 4 + 3] = c[3];
}

export function uniform(n: number, c: RGBA): Uint8Array {
  const out = new Uint8Array(n * 4);
  for (let i = 0; i < n; i++) write(out, i, c);
  return out;
}

/** Colours per category: the file's own colours when it has one per category, else the palette. */
export function categoryPalette(n: number, fileColors?: readonly string[] | null): RGBA[] {
  if (fileColors && fileColors.length === n) {
    return fileColors.map((hex) => {
      const [r, g, b] = hexToRgb(hex);
      return [r, g, b, 220] as const;
    });
  }
  return Array.from({ length: n }, (_, i) => categorical(i));
}

export function colorsForCodes(codes: ArrayLike<number>, palette: readonly RGBA[]): Uint8Array {
  const out = new Uint8Array(codes.length * 4);
  for (let i = 0; i < codes.length; i++) {
    const c = codes[i]!;
    write(out, i, c >= 0 && c < palette.length ? palette[c]! : MISSING);
  }
  return out;
}

/** 1st-99th percentile of the finite values: one outlier cannot wash out the scale. */
export function robustDomain(values: ArrayLike<number>): [number, number] {
  const finite: number[] = [];
  for (let i = 0; i < values.length; i++) {
    const v = values[i]!;
    if (Number.isFinite(v)) finite.push(v);
  }
  if (finite.length === 0) return [0, 1];
  finite.sort((a, b) => a - b);
  const at = (q: number) => finite[Math.min(finite.length - 1, Math.floor(q * (finite.length - 1)))]!;
  const lo = at(0.01);
  const hi = at(0.99);
  return hi > lo ? [lo, hi] : [finite[0]!, finite[finite.length - 1]! + 1e-9];
}

export function colorsForValues(values: ArrayLike<number>, domain: [number, number]): Uint8Array {
  const [lo, hi] = domain;
  const span = hi - lo || 1;
  const out = new Uint8Array(values.length * 4);
  for (let i = 0; i < values.length; i++) {
    const v = values[i]!;
    write(out, i, Number.isFinite(v) ? sequential((v - lo) / span) : MISSING);
  }
  return out;
}

/** primary_lineage: index, -1 shared trunk, -2 unassigned. */
export function colorsForLineages(primary: ArrayLike<number>): Uint8Array {
  const out = new Uint8Array(primary.length * 4);
  for (let i = 0; i < primary.length; i++) {
    const p = primary[i]!;
    write(out, i, p >= 0 ? categorical(p) : p === -1 ? SHARED : MISSING);
  }
  return out;
}

export function rgbaCss(c: RGBA): string {
  return `rgba(${c[0]}, ${c[1]}, ${c[2]}, ${(c[3] / 255).toFixed(2)})`;
}
