import { describe, expect, it, vi } from "vitest";
import { DepthFadeExtension, depthFadeModule, UNIFORM_BLOCK } from "@/lib/gl/depthFade";

const GLSL_TYPE: Record<string, string> = { f32: "float", "vec3<f32>": "vec3" };

describe("depthFade shader module", () => {
  it("declares its uniforms in the same order and types as the GLSL block (std140 offsets follow it)", () => {
    const fields = [...UNIFORM_BLOCK.matchAll(/^\s+(\w+)\s+(\w+);$/gm)].map((m) => [m[2], m[1]]);
    const declared = Object.entries(depthFadeModule.uniformTypes).map(([name, type]) => [name, GLSL_TYPE[type]]);
    expect(fields).toEqual(declared);
  });

  it("feeds the fog range from the camera and the stage colour as 0-1 RGB", () => {
    const setShaderModuleProps = vi.fn();
    const layer = {
      props: {
        depthFadeStage: [255, 0, 51] as const,
        depthFadeStrength: 0.7,
        depthFadeBounds: [
          [-1, -1, -1],
          [1, 1, 1],
        ],
      },
      setShaderModuleProps,
    };
    const draw = DepthFadeExtension.prototype.draw as (this: unknown, p: unknown) => void;
    draw.call(layer, { context: { viewport: { cameraPosition: [0, 0, 5] } } });
    expect(setShaderModuleProps).toHaveBeenCalledWith({
      depthFade: { stage: [1, 0, 0.2], near: 4, far: expect.closeTo(Math.hypot(1, 1, 6), 9), strength: 0.7, dim: 0 },
    });
  });

  it("reads the selection once per cell, not per vertex of each cell's quad", () => {
    // deck.gl's default step mode is per vertex: every cell would read cells 0-3's selection.
    const add = vi.fn();
    const init = DepthFadeExtension.prototype.initializeState as unknown as (this: unknown) => void;
    init.call({ getAttributeManager: () => ({ add }) });
    expect(add).toHaveBeenCalledWith({
      depthFadeSelected: expect.objectContaining({ stepMode: "dynamic", accessor: "getSelected", size: 1 }),
    });
  });
});
