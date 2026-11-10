/**
 * Depth cue for the 3-D diffusion view: cells further from the camera are drawn dimmer
 * and more transparent, so a still frame of an orbited stage reads as depth.
 *
 * Runs in the shader (no per-frame CPU work at 100k cells). The camera-facing axis is
 * measured from the viewport itself: in an orthographic view, depth is linear in world
 * position, so projecting the unit axes gives the direction depth grows along.
 */
import { type Layer, LayerExtension, project, type Viewport } from "@deck.gl/core";
import type { ShaderModule } from "@luma.gl/shadertools";

export interface DepthCueProps {
  /** 0 turns the cue off (2-D embeddings); 1 is the full falloff. */
  depthCue?: number;
  /** Stage radius: the distance over which the cue runs from front to back. */
  depthCueRadius?: number;
}

type DepthCueUniforms = {
  axis: [number, number, number];
  radius: number;
  strength: number;
};

const uniformBlock = /* glsl */ `\
layout(std140) uniform depthCueUniforms {
  vec3 axis;
  float radius;
  float strength;
} depthCue;
`;

const depthCueModule: ShaderModule<DepthCueUniforms, DepthCueUniforms> = {
  name: "depthCue",
  dependencies: [project],
  vs: `${uniformBlock}\nout float depthCue_fade;\n`,
  fs: `${uniformBlock}\nin float depthCue_fade;\n`,
  inject: {
    "vs:DECKGL_FILTER_GL_POSITION": /* glsl */ `
      float depthCue_t = clamp(dot(geometry.worldPosition, depthCue.axis) / max(depthCue.radius, 1e-6), -1.0, 1.0);
      // t = 1 nearest the camera, -1 furthest; the furthest cells keep 30% of their light.
      depthCue_fade = 1.0 - depthCue.strength * 0.7 * (1.0 - (depthCue_t + 1.0) * 0.5);
    `,
    "fs:DECKGL_FILTER_COLOR": /* glsl */ `
      color.rgb *= depthCue_fade;
      color.a *= mix(1.0, depthCue_fade, 0.85);
    `,
  },
  uniformTypes: { axis: "vec3<f32>", radius: "f32", strength: "f32" },
};

/** The unit vector, in world space, that points from the stage towards the camera. */
export function cameraAxis(viewport: Viewport): [number, number, number] {
  const o = viewport.project([0, 0, 0]);
  const d = [
    viewport.project([1, 0, 0])[2]! - o[2]!,
    viewport.project([0, 1, 0])[2]! - o[2]!,
    viewport.project([0, 0, 1])[2]! - o[2]!,
  ];
  const len = Math.hypot(d[0]!, d[1]!, d[2]!) || 1;
  // Projected depth grows away from the camera, so the camera lies along -d.
  return [-d[0]! / len, -d[1]! / len, -d[2]! / len];
}

export class DepthCueExtension extends LayerExtension {
  static extensionName = "DepthCueExtension";
  static defaultProps = { depthCue: 0, depthCueRadius: 1 };

  getShaders() {
    return { modules: [depthCueModule] };
  }

  draw(this: Layer<DepthCueProps>, params: { context: { viewport: Viewport } }) {
    const strength = this.props.depthCue ?? 0;
    this.setShaderModuleProps({
      depthCue: {
        axis: strength > 0 ? cameraAxis(params.context.viewport) : [0, 0, 1],
        radius: this.props.depthCueRadius ?? 1,
        strength,
      },
    });
  }
}
