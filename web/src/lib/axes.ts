/** Which three diffusion components the 3-D view shows, in x, y, z order. */
export type Axes = [number, number, number];

export const DEFAULT_AXES: Axes = [1, 2, 3];

const isDefault = (a: Axes) => a[0] === 1 && a[1] === 2 && a[2] === 3;

/** The server's name for these components: dc123 for DC1-DC3, else dc:i,j,k. */
export const axesName = (a: Axes) => (isDefault(a) ? "dc123" : `dc:${a.join(",")}`);

/** DC1–DC3 for the first three, else each one: DC1, DC2, DC4. */
export const axesLabel = (a: Axes) => (isDefault(a) ? "DC1–DC3" : a.map((i) => `DC${i}`).join(", "));

/** Put component `dc` on axis k. The axis that had it takes k's old component, so the three stay different. */
export function setAxis(a: Axes, k: number, dc: number): Axes {
  const next: Axes = [...a];
  const j = a.indexOf(dc);
  if (j >= 0) next[j] = a[k]!;
  next[k] = dc;
  return next;
}
