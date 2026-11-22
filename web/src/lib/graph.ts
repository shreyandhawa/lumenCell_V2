/**
 * Each edge's two ends at its cells' positions in one view (xyz is n x 3), as the flat arrays
 * deck.gl's LineLayer reads.
 */
export function edgeEnds(src: Uint32Array, dst: Uint32Array, xyz: Float32Array): { source: Float32Array; target: Float32Array } {
  const source = new Float32Array(src.length * 3);
  const target = new Float32Array(src.length * 3);
  for (let e = 0; e < src.length; e++) {
    for (let k = 0; k < 3; k++) {
      source[e * 3 + k] = xyz[src[e]! * 3 + k]!;
      target[e * 3 + k] = xyz[dst[e]! * 3 + k]!;
    }
  }
  return { source, target };
}
