/**
 * The plot as a PNG: the drawn cells cut from deck.gl's canvas, on the theme's stage colour, with
 * what the screen shows around them (caption, root pin, legend, gimbal) and where it came from.
 */
import type { LegendSpec } from "@/components/useCellColors";
import { rgbaCss, sequential } from "@/lib/colors";

export interface PngInput {
  /** deck.gl's canvas, drawn in this same task (its buffer is not kept between frames). */
  gl: HTMLCanvasElement;
  /** The plot area to keep, in CSS pixels of the canvas. */
  crop: { x: number; y: number; width: number; height: number };
  colors: { stage: string; ink: string; muted: string; accent: string };
  font: string;
  captions: string[];
  /** Root pins in CSS pixels of the canvas; the first is labelled. */
  pins: { x: number; y: number }[];
  pinLabel: string;
  legend: LegendSpec | null;
  gimbal: { m: ArrayLike<number>; names: readonly string[] } | null;
  /** The provenance line along the bottom. */
  footer: string;
}

const MAX_ITEMS = 14;

export function exportPng(p: PngInput): Promise<Blob> {
  const scale = p.gl.width / p.gl.clientWidth || 1;
  const { x, y, width: w, height: h } = p.crop;
  const out = document.createElement("canvas");
  out.width = Math.round(w * scale);
  out.height = Math.round(h * scale);
  const g = out.getContext("2d")!;
  g.fillStyle = p.colors.stage;
  g.fillRect(0, 0, out.width, out.height);
  g.drawImage(p.gl, x * scale, y * scale, w * scale, h * scale, 0, 0, out.width, out.height);
  g.scale(scale, scale);
  g.textBaseline = "middle";
  const text = (s: string, tx: number, ty: number, color: string, size = 12) => {
    g.font = `${size}px ${p.font}`;
    g.fillStyle = color;
    g.fillText(s, tx, ty);
    return g.measureText(s).width;
  };

  p.captions.forEach((c, i) => text(c, 16, 22 + i * 18, p.colors.muted));

  p.pins.forEach((pin, i) => {
    const [px, py] = [pin.x - x, pin.y - y];
    g.beginPath();
    g.arc(px, py, 7, 0, 2 * Math.PI);
    g.fillStyle = p.colors.stage;
    g.fill();
    g.beginPath();
    g.arc(px, py, 5, 0, 2 * Math.PI);
    g.fillStyle = p.colors.accent;
    g.fill();
    if (i === 0) text(p.pinLabel, px + 12, py, p.colors.ink);
  });

  if (p.legend) drawLegend(g, p.legend, 16, h - 52, w - 160, text, p.colors.muted);
  if (p.gimbal) drawGimbal(g, p.gimbal, w - 56, h - 76, p.colors, text);
  text(p.footer, 16, h - 18, p.colors.muted, 11);

  return new Promise((resolve, reject) => out.toBlob((b) => (b ? resolve(b) : reject(new Error("the browser made no PNG"))), "image/png"));
}

type Text = (s: string, x: number, y: number, color: string, size?: number) => number;

function drawLegend(g: CanvasRenderingContext2D, spec: LegendSpec, x0: number, y: number, maxWidth: number, text: Text, muted: string) {
  let x = x0 + text(spec.title, x0, y, muted) + 12;
  if (spec.kind === "continuous") {
    x += text(spec.domain[0].toPrecision(3), x, y, muted) + 8;
    const ramp = g.createLinearGradient(x, 0, x + 112, 0);
    for (const t of [0, 0.25, 0.5, 0.75, 1]) ramp.addColorStop(t, rgbaCss(sequential(t, 255)));
    g.fillStyle = ramp;
    g.fillRect(x, y - 3, 112, 6);
    text(spec.domain[1].toPrecision(3), x + 120, y, muted);
    return;
  }
  for (const item of spec.items.slice(0, MAX_ITEMS)) {
    // text() left the font at 12px, so this measures the label as it will be drawn
    if (x + 20 + g.measureText(item.label).width > x0 + maxWidth) {
      text(`and ${spec.items.length - spec.items.indexOf(item)} more`, x, y, muted);
      return;
    }
    g.beginPath();
    g.arc(x + 4, y, 4, 0, 2 * Math.PI);
    g.fillStyle = rgbaCss(item.color);
    g.fill();
    x += 12 + text(item.label, x + 12, y, muted) + 12;
  }
  if (spec.items.length > MAX_ITEMS) text(`and ${spec.items.length - MAX_ITEMS} more`, x, y, muted);
}

/** The same orientation the on-screen gimbal shows: each axis's direction in the view matrix. */
function drawGimbal(
  g: CanvasRenderingContext2D,
  { m, names }: { m: ArrayLike<number>; names: readonly string[] },
  cx: number,
  cy: number,
  colors: { ink: string; muted: string },
  text: Text,
) {
  g.beginPath();
  g.arc(cx, cy, 40, 0, 2 * Math.PI);
  g.strokeStyle = colors.muted;
  g.lineWidth = 1;
  g.stroke();
  names.forEach((name, k) => {
    const [ax, ay, az] = [m[k * 4]!, m[k * 4 + 1]!, m[k * 4 + 2]!];
    const len = Math.hypot(ax, ay, az) || 1;
    const [dx, dy] = [(ax / len) * 26, (-ay / len) * 26];
    g.beginPath();
    g.moveTo(cx, cy);
    g.lineTo(cx + dx, cy + dy);
    g.strokeStyle = colors.ink;
    g.lineWidth = 1.25;
    g.stroke();
    g.textAlign = "center";
    text(name, cx + dx * 1.3, cy + dy * 1.3, colors.muted, 9);
    g.textAlign = "start";
  });
}
