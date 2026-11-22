import { describe, expect, it } from "vitest";
import { axesLabel, axesName, DEFAULT_AXES, setAxis } from "@/lib/axes";

describe("diffusion axes", () => {
  it("names the first three as the server's dc123 and any others as dc:i,j,k", () => {
    expect(axesName(DEFAULT_AXES)).toBe("dc123");
    expect(axesName([1, 2, 4])).toBe("dc:1,2,4");
    expect(axesName([2, 1, 3])).toBe("dc:2,1,3");
    expect(axesLabel(DEFAULT_AXES)).toBe("DC1–DC3");
    expect(axesLabel([1, 2, 4])).toBe("DC1, DC2, DC4");
  });

  it("keeps the three components different: choosing one already shown swaps the two axes", () => {
    expect(setAxis([1, 2, 3], 2, 4)).toEqual([1, 2, 4]);
    expect(setAxis([1, 2, 3], 0, 3)).toEqual([3, 2, 1]);
    expect(setAxis([1, 2, 3], 1, 2)).toEqual([1, 2, 3]);
  });
});
