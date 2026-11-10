/**
 * LCB1 decoder — the browser half of `src/lumencell/api/binary.py`.
 *
 * Every per-cell payload (embeddings, obs columns, expression, trajectory cells,
 * heatmaps) arrives in this one container. Arrays are returned as zero-copy
 * TypedArray views into the response buffer, so a 100k-cell embedding goes
 * from `fetch` to deck.gl without being copied.
 *
 * Layout: 16-byte header (`"LCB1"`, u32 version, u32 manifest length, u32 0),
 * a UTF-8 JSON manifest, zero padding to 8 bytes, then 8-byte-aligned blocks.
 * Array offsets in the manifest are relative to the start of the data section.
 */

export type LcbDtype = "u8" | "i8" | "u16" | "i16" | "u32" | "i32" | "f32";

export interface TypedArrayFor {
  u8: Uint8Array;
  i8: Int8Array;
  u16: Uint16Array;
  i16: Int16Array;
  u32: Uint32Array;
  i32: Int32Array;
  f32: Float32Array;
}

export interface LcbArray<D extends LcbDtype = LcbDtype> {
  name: string;
  dtype: D;
  shape: number[];
  data: TypedArrayFor[D];
}

export interface LcbContainer {
  version: number;
  meta: Record<string, unknown>;
  arrays: Map<string, LcbArray>;
  /** Fetch an array by name, asserting its dtype so the result is typed. */
  get<D extends LcbDtype>(name: string, dtype: D): LcbArray<D>;
}

export class LcbError extends Error {
  override name = "LcbError";
}

const MAGIC = "LCB1";
const SUPPORTED_VERSION = 1;
const HEADER_BYTES = 16;
const ALIGN = 8;

const CONSTRUCTORS = {
  u8: Uint8Array,
  i8: Int8Array,
  u16: Uint16Array,
  i16: Int16Array,
  u32: Uint32Array,
  i32: Int32Array,
  f32: Float32Array,
} as const;

// TypedArrays use host byte order and the wire format is little-endian. Every
// browser engine lumenCell targets is little-endian; fail loudly if that ever
// stops being true rather than rendering garbage.
const HOST_IS_LITTLE_ENDIAN = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;

interface ManifestEntry {
  name: string;
  dtype: string;
  shape: number[];
  offset: number;
  byteLength: number;
}

function isDtype(value: string): value is LcbDtype {
  return Object.hasOwn(CONSTRUCTORS, value);
}

function align(n: number): number {
  return Math.ceil(n / ALIGN) * ALIGN;
}

export function decodeLcb(buffer: ArrayBuffer): LcbContainer {
  if (!HOST_IS_LITTLE_ENDIAN) throw new LcbError("big-endian hosts are not supported");
  if (buffer.byteLength < HEADER_BYTES) throw new LcbError("buffer is shorter than the LCB1 header");

  const view = new DataView(buffer);
  const magic = String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3));
  if (magic !== MAGIC) throw new LcbError(`bad magic ${JSON.stringify(magic)}`);
  const version = view.getUint32(4, true);
  if (version !== SUPPORTED_VERSION) throw new LcbError(`unsupported LCB1 version ${version}`);
  const manifestLength = view.getUint32(8, true);
  const manifestEnd = HEADER_BYTES + manifestLength;
  if (buffer.byteLength < manifestEnd) throw new LcbError("buffer is truncated inside the manifest");

  let manifest: { arrays?: ManifestEntry[]; meta?: Record<string, unknown> };
  try {
    const text = new TextDecoder("utf-8", { fatal: true }).decode(
      new Uint8Array(buffer, HEADER_BYTES, manifestLength),
    );
    manifest = JSON.parse(text) as typeof manifest;
  } catch (err) {
    throw new LcbError(`manifest is not valid UTF-8 JSON: ${String(err)}`);
  }

  const dataStart = align(manifestEnd);
  const arrays = new Map<string, LcbArray>();
  for (const entry of manifest.arrays ?? []) {
    if (!isDtype(entry.dtype)) throw new LcbError(`array "${entry.name}": unknown dtype ${entry.dtype}`);
    const Ctor = CONSTRUCTORS[entry.dtype];
    const count = entry.shape.reduce((a, b) => a * b, 1);
    const start = dataStart + entry.offset;
    if (start % Ctor.BYTES_PER_ELEMENT !== 0) {
      throw new LcbError(`array "${entry.name}" is not aligned to its element size`);
    }
    if (entry.byteLength !== count * Ctor.BYTES_PER_ELEMENT) {
      throw new LcbError(`array "${entry.name}": byteLength does not match shape and dtype`);
    }
    if (start + entry.byteLength > buffer.byteLength) {
      throw new LcbError(`array "${entry.name}" runs past the end of the buffer`);
    }
    const data = new Ctor(buffer, start, count);
    arrays.set(entry.name, { name: entry.name, dtype: entry.dtype, shape: entry.shape, data });
  }

  return {
    version,
    meta: manifest.meta ?? {},
    arrays,
    get<D extends LcbDtype>(name: string, dtype: D): LcbArray<D> {
      const arr = arrays.get(name);
      if (!arr) throw new LcbError(`array "${name}" is not in this payload`);
      if (arr.dtype !== dtype) throw new LcbError(`array "${name}" is ${arr.dtype}, expected ${dtype}`);
      return arr as LcbArray<D>;
    },
  };
}

/** Zero-copy view of row `i` of a 2-D (row-major) array, e.g. one lineage or one gene. */
export function row<D extends LcbDtype>(arr: LcbArray<D>, i: number): TypedArrayFor[D] {
  if (arr.shape.length !== 2) throw new LcbError(`array "${arr.name}" is not 2-D`);
  const [rows = 0, cols = 0] = arr.shape;
  if (i < 0 || i >= rows) throw new LcbError(`row ${i} is out of range for "${arr.name}"`);
  return arr.data.subarray(i * cols, (i + 1) * cols) as TypedArrayFor[D];
}
