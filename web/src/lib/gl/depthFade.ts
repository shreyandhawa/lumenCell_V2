/**
 * Real depth for dense point clouds. Cells are drawn opaque with the depth buffer on, so the nearest
 * cell wins instead of the last one drawn, and far cells fade toward the stage colour on the GPU
 * (fog), so orbiting never touches a per-cell JS buffer.
 *
 * Colour alpha is kept as meaning, not transparency: a cell's colour is blended toward the stage by
 * its alpha (missing values stay faint) before it is made opaque.
 */
import { type Layer, LayerExtension, project } from "@deck.gl/core";
import { depthRange } from "@/lib/camera";
import type { Vec3 } from "@/lib/embedding";

// Field order must match uniformTypes below (std140 offsets follow it).
export const UNIFORM_BLOCK = /* glsl */ `\
layout(std140) uniform depthFadeUniforms {
  vec3 stage;
  float near;
  float far;
  float strength;
} depthFade;
`;

export type DepthFadeUniforms = {
  stage: [number, number, number];
  near: number;
  far: number;
  strength: number;
};

export const depthFadeModule = {
  name: "depthFade",
  dependencies: [project],
  vs: UNIFORM_BLOCK,
  fs: UNIFORM_BLOCK,
  inject: {
    // Fade in the vertex stage, so the picking highlight (fragment stage) is never fogged.
    "vs:DECKGL_FILTER_COLOR": /* glsl */ `
    color.rgb = mix(depthFade.stage, color.rgb, color.a);
    color.a = 1.0;
    float depthFade_t = clamp(
      (distance(geometry.position.xyz, project.cameraPosition) - depthFade.near) / max(depthFade.far - depthFade.near, 1e-6),
      0.0, 1.0
    );
    color.rgb = mix(color.rgb, depthFade.stage, depthFade.strength * depthFade_t);
    `,
    // Antialiased rims would write depth while translucent and punch holes in cells behind them.
    "fs:DECKGL_FILTER_COLOR": /* glsl */ `
    if (color.a < 0.5) discard;
    color.a = 1.0;
    `,
  },
  getUniforms: (props?: Partial<DepthFadeUniforms>) => props ?? {},
  uniformTypes: {
    stage: "vec3<f32>",
    near: "f32",
    far: "f32",
    strength: "f32",
  },
} as const;

export interface DepthFadeProps {
  /** Stage colour, 0-255 RGB: what far and faint cells fade toward. */
  depthFadeStage?: readonly [number, number, number];
  /** 0 = no fog (2-D), up to 1 = the farthest cells reach the stage colour. */
  depthFadeStrength?: number;
  /** The box that holds every cell (min corner, max corner): the fog spans it. */
  depthFadeBounds?: readonly [Vec3, Vec3];
}

const UNIT_BOX: [Vec3, Vec3] = [
  [-0.5, -0.5, -0.5],
  [0.5, 0.5, 0.5],
];

export class DepthFadeExtension extends LayerExtension {
  static defaultProps = { depthFadeStage: [0, 0, 0], depthFadeStrength: 0, depthFadeBounds: UNIT_BOX };
  static extensionName = "DepthFadeExtension";

  getShaders() {
    return { modules: [depthFadeModule] };
  }

  draw(this: Layer<DepthFadeProps>, params: { context: { viewport: { cameraPosition: number[] } } }) {
    const { depthFadeStage = [0, 0, 0], depthFadeStrength = 0, depthFadeBounds = UNIT_BOX } = this.props;
    const [near, far] = depthRange(params.context.viewport.cameraPosition, depthFadeBounds);
    const uniforms: DepthFadeUniforms = {
      stage: depthFadeStage.map((c) => c / 255) as [number, number, number],
      near,
      far,
      strength: depthFadeStrength,
    };
    this.setShaderModuleProps({ depthFade: uniforms });
  }
}
