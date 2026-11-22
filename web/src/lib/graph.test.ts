import { expect, it } from "vitest";
import { edgeEnds } from "@/lib/graph";

it("puts each edge's ends at its two cells", () => {
  const xyz = Float32Array.from([0, 0, 0, 1, 2, 3, 4, 5, 6]);
  const { source, target } = edgeEnds(Uint32Array.from([0, 1]), Uint32Array.from([2, 2]), xyz);
  expect([...source]).toEqual([0, 0, 0, 1, 2, 3]);
  expect([...target]).toEqual([4, 5, 6, 4, 5, 6]);
});
