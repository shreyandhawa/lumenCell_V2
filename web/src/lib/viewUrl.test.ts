import { describe, expect, it } from "vitest";
import { decodeView, encodeView, type SharedView } from "@/lib/viewUrl";

const FP = "9d273491";
const view: SharedView = {
  run: "palantir-1a2b",
  colorBy: { kind: "obs", column: "lumencell_leiden" },
  view: "3d",
  cameras: { stage: { target: [0.1, -0.2, 0.05], zoom: 8.25, rotationX: 18, rotationOrbit: -32.5 } },
  focus: 1,
  cursor: 0.4,
  axes: [1, 2, 4],
  pointSize: 1.25,
};

describe("view in the URL hash", () => {
  it("round-trips everything it stores", () => {
    const hash = encodeView(view, FP);
    expect(decodeView(`#${hash}`, FP)).toEqual(view);
    expect(decodeView(encodeView({ ...view, colorBy: { kind: "gene", gene: "CD8A" }, view: "split", cameras: {} }, FP), FP)).toMatchObject({
      colorBy: { kind: "gene", gene: "CD8A" },
      view: "split",
    });
  });

  it("leaves out defaults, so a plain view makes a short link", () => {
    const hash = encodeView(
      { ...view, run: null, colorBy: null, view: "2d", cameras: {}, focus: null, cursor: null, axes: [1, 2, 3], pointSize: 1 },
      FP,
    );
    expect(hash).toBe(`v=1&data=${FP}&view=2d`);
  });

  it("is ignored on another dataset, another version, or no hash", () => {
    const hash = encodeView(view, FP);
    expect(decodeView(hash, "other")).toBeNull();
    expect(decodeView(hash.replace("v=1", "v=2"), FP)).toBeNull();
    expect(decodeView("", FP)).toBeNull();
  });

  it("drops each bad field instead of guessing", () => {
    const bad = `v=1&data=${FP}&view=4d&run=../etc&color=obs:&cam.stage=1,2,x,4,5,6&cam.flat=1,2&focus=-1&t=7&axes=1,1,2&size=40`;
    expect(decodeView(bad, FP)).toEqual({
      run: null,
      colorBy: null,
      view: "2d",
      cameras: {},
      focus: null,
      cursor: null,
      axes: [1, 2, 3],
      pointSize: 1,
    });
  });
});
