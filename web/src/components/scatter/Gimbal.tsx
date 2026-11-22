/** Which way the three shown diffusion components (DC1, DC2 and DC3 by default) point on screen. */
export function Gimbal({ m, names }: { m: number[]; names: readonly [string, string, string] }) {
  const axes = names.map((name, k) => {
    const [x, y, z] = [m[k * 4]!, m[k * 4 + 1]!, m[k * 4 + 2]!];
    const len = Math.hypot(x, y, z) || 1;
    return { name, x: (x / len) * 30, y: (-y / len) * 30 };
  });
  return (
    <svg
      role="img"
      aria-label={`Orientation of ${names[0]}, ${names[1]} and ${names[2]}`}
      viewBox="-48 -48 96 96"
      className="size-24 rounded-full bg-stage/80 text-ink-muted"
    >
      <circle r={47} fill="none" stroke="currentColor" strokeOpacity={0.5} />
      {axes.map((a) => (
        <g key={a.name}>
          <line x1={0} y1={0} x2={a.x} y2={a.y} stroke="var(--color-ink)" strokeWidth={1.25} />
          <text x={a.x * 1.3} y={a.y * 1.3 + 3} textAnchor="middle" className="fill-ink-muted font-mono text-[9px]">
            {a.name}
          </text>
        </g>
      ))}
      <circle r={2.5} fill="var(--color-ink)" />
    </svg>
  );
}
