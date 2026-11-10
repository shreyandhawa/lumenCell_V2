// @vitest-environment node
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { describe, expect, it } from "vitest";
import { decodeLcb, LcbError, row } from "@/api/lcb1";

// The same files the Python suite reads (tests/contracts/test_lcb1.py).
const FIXTURE_DIR = fileURLToPath(new URL("../../../contracts/lcb1/", import.meta.url));

function loadFixture(): ArrayBuffer {
  const bytes = readFileSync(`${FIXTURE_DIR}fixture.lcb`);
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer;
}

interface Expected {
  meta: Record<string, unknown>;
  arrays: Record<string, { dtype: string; shape: number[]; values: (number | string)[] }>;
}

const expected = JSON.parse(readFileSync(`${FIXTURE_DIR}fixture.json`, "utf-8")) as Expected;

function fromJson(v: number | string): number {
  if (v === "NaN") return Number.NaN;
  if (v === "Infinity") return Number.POSITIVE_INFINITY;
  if (v === "-Infinity") return Number.NEGATIVE_INFINITY;
  return v as number;
}

describe("decodeLcb with the shared fixture", () => {
  const decoded = decodeLcb(loadFixture());

  it("reads the meta block, including unicode", () => {
    expect(decoded.meta).toEqual(expected.meta);
  });

  it("keeps the manifest's array order", () => {
    expect([...decoded.arrays.keys()]).toEqual(Object.keys(expected.arrays));
  });

  it.each(Object.entries(expected.arrays))("decodes %s exactly", (name, spec) => {
    const arr = decoded.arrays.get(name);
    expect(arr).toBeDefined();
    expect(arr!.dtype).toBe(spec.dtype);
    expect(arr!.shape).toEqual(spec.shape);
    expect(Array.from(arr!.data)).toEqual(spec.values.map(fromJson));
  });

  it("returns zero-copy views into the response buffer", () => {
    const buffer = loadFixture();
    const pt = decodeLcb(buffer).get("pseudotime", "f32");
    expect(pt.data.buffer).toBe(buffer);
  });

  it("slices rows of a 2-D array without copying", () => {
    const m = decoded.get("matrix", "f32");
    expect(Array.from(row(m, 1))).toEqual([3.5, 4.5, 5.5]);
    expect(row(m, 1).buffer).toBe(m.data.buffer);
  });
});

describe("decodeLcb validation", () => {
  it("rejects a bad magic", () => {
    const buffer = loadFixture();
    new Uint8Array(buffer)[0] = 0x58;
    expect(() => decodeLcb(buffer)).toThrow(/magic/);
  });

  it("rejects an unknown version", () => {
    const buffer = loadFixture();
    new DataView(buffer).setUint32(4, 99, true);
    expect(() => decodeLcb(buffer)).toThrow(/version/);
  });

  it("rejects a truncated buffer", () => {
    const buffer = loadFixture().slice(0, -8);
    expect(() => decodeLcb(buffer)).toThrow(LcbError);
  });

  it("rejects a short buffer", () => {
    expect(() => decodeLcb(new ArrayBuffer(4))).toThrow(/shorter/);
  });

  it("type-checks names and dtypes on get()", () => {
    const fixture = decodeLcb(loadFixture());
    expect(() => fixture.get("missing", "f32")).toThrow(/not in this payload/);
    expect(() => fixture.get("pseudotime", "i32")).toThrow(/expected i32/);
  });
});
