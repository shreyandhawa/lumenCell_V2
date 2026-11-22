/**
 * How each cell reads on a dense point cloud, on the GPU.
 *
 * - Real depth: cells are drawn opaque with the depth buffer on, so the nearest cell wins instead of
 *   the last one drawn. Colour alpha is kept as meaning, not transparency: a cell's colour is
 *   blended toward the stage by its alpha (missing values stay faint) before it is made opaque.
 * - Fog: far cells fade toward the stage colour, so orbiting never touches a per-cell JS buffer.
 * - Selection: while a selection exists, unselected cells fade toward the stage, shrink and drop
 *   behind every selected cell; selected cells skip the fog.
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
  float dim;
} depthFade;
`;

export type DepthFadeUniforms = {
  stage: [number, number, number];
  near: number;
  far: number;
  strength: number;
  dim: number;
};

// One float per cell: 1 = selected, 0 = not.
const VS = /* glsl */ `${UNIFORM_BLOCK}
in float depthFadeSelected;
bool depthFade_dimmed() { return depthFade.dim > 0.5 && depthFadeSelected < 0.5; }
bool depthFade_chosen() { return depthFade.dim > 0.5 && depthFadeSelected >= 0.5; }
`;

export const depthFadeModule = {
  name: "depthFade",
  dependencies: [project],
  vs: VS,
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
    // Selected cells keep their full colour: the selection is what is being looked at.
    if (!depthFade_chosen()) color.rgb = mix(color.rgb, depthFade.stage, depthFade.strength * depthFade_t);
    if (depthFade_dimmed()) color.rgb = mix(color.rgb, depthFade.stage, 0.8);
    `,
    // Dimmed cells go to the back of the depth range, so a selected cell is never covered by one
    // (in 2-D every cell has the same depth and the last drawn would win).
    "vs:DECKGL_FILTER_GL_POSITION": /* glsl */ `
    if (depthFade_dimmed()) position.z = position.w * 0.9999;
    `,
    "vs:DECKGL_FILTER_SIZE": /* glsl */ `
    if (depthFade_dimmed()) size *= 0.5;
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
    dim: "f32",
  },
} as const;

export interface DepthFadeProps {
  /** Stage colour, 0-255 RGB: what far and faint cells fade toward. */
  depthFadeStage?: readonly [number, number, number];
  /** 0 = no fog (2-D), up to 1 = the farthest cells reach the stage colour. */
  depthFadeStrength?: number;
  /** The box that holds every cell (min corner, max corner): the fog spans it. */
  depthFadeBounds?: readonly [Vec3, Vec3];
  /** A selection exists: cells whose `getSelected` is 0 are dimmed. */
  depthFadeDim?: boolean;
}

const UNIT_BOX: [Vec3, Vec3] = [
  [-0.5, -0.5, -0.5],
  [0.5, 0.5, 0.5],
];

export class DepthFadeExtension extends LayerExtension {
  static defaultProps = {
    depthFadeStage: [0, 0, 0],
    depthFadeStrength: 0,
    depthFadeBounds: UNIT_BOX,
    depthFadeDim: false,
    getSelected: { type: "accessor", value: 1 },
  };
  static extensionName = "DepthFadeExtension";

  getShaders() {
    return { modules: [depthFadeModule] };
  }

  initializeState(this: Layer<DepthFadeProps>) {
    this.getAttributeManager()?.add({
      // dynamic = one value per instance (per cell) on instanced layers; the default is per vertex.
      depthFadeSelected: { size: 1, stepMode: "dynamic", accessor: "getSelected", defaultValue: [1] },
    });
  }

  draw(this: Layer<DepthFadeProps>, params: { context: { viewport: { cameraPosition: number[] } } }) {
    const { depthFadeStage = [0, 0, 0], depthFadeStrength = 0, depthFadeBounds = UNIT_BOX, depthFadeDim = false } = this.props;
    const [near, far] = depthRange(params.context.viewport.cameraPosition, depthFadeBounds);
    const uniforms: DepthFadeUniforms = {
      stage: depthFadeStage.map((c) => c / 255) as [number, number, number],
      near,
      far,
      strength: depthFadeStrength,
      dim: depthFadeDim ? 1 : 0,
    };
    this.setShaderModuleProps({ depthFade: uniforms });
  }
}
