import { readFile } from 'node:fs/promises';

let h5wasm = null;

/** Lazy-init h5wasm */
async function ensureH5wasm() {
  if (!h5wasm) {
    h5wasm = await import('h5wasm/node');
    await h5wasm.ready;
  }
  return h5wasm;
}

// ── LRU file handle cache ──────────────────────────────────────────────
const MAX_CACHE = 5;
const cache = new Map(); // path → { file, lastUsed }

async function openFile(filePath) {
  if (cache.has(filePath)) {
    const entry = cache.get(filePath);
    entry.lastUsed = Date.now();
    return entry.file;
  }
  const hw = await ensureH5wasm();
  if (cache.size >= MAX_CACHE) {
    let oldest = null;
    for (const [p, e] of cache) {
      if (!oldest || e.lastUsed < oldest.lastUsed) oldest = { path: p, ...e };
    }
    if (oldest) {
      try { cache.get(oldest.path).file.close(); } catch (_) {}
      cache.delete(oldest.path);
    }
  }
  const file = new hw.File(filePath, 'r');
  cache.set(filePath, { file, lastUsed: Date.now() });
  return file;
}

// ── Metadata ───────────────────────────────────────────────────────────
export async function readMetadata(filePath) {
  const file = await openFile(filePath);
  const params = file.get('/Model/Params').value; // Float64Array
  return {
    time: params[0],
    width: params[1],
    height: params[2],
    nx: Math.round(params[3]),
    nz: Math.round(params[4]),
    dx: params[5],
    dz: params[6],
    dt: params[7],
  };
}

// ── Grid dimension table ───────────────────────────────────────────────
function gridInfo(nx, nz, gridType) {
  switch (gridType) {
    case 'center':
      return { rows: nx - 1, cols: nz - 1, xPath: '/Model/xc_coord', zPath: '/Model/zc_coord' };
    case 'vertex':
      return { rows: nx, cols: nz, xPath: '/Model/xg_coord', zPath: '/Model/zg_coord' };
    case 'vx':
      return { rows: nx, cols: nz + 1, xPath: '/Model/xg_coord', zPath: '/Model/zvx_coord' };
    case 'vz':
      return { rows: nx + 1, cols: nz, xPath: '/Model/xvz_coord', zPath: '/Model/zg_coord' };
    case 'vizgrid':
      return { rows: nx - 1, cols: nz - 1, xPath: '/VizGrid/xviz', zPath: '/VizGrid/zviz' };
    case 'vizgrid_hr':
      return { rows: 2 * (nx - 1), cols: 2 * (nz - 1), xPath: '/VizGrid/xviz_hr', zPath: '/VizGrid/zviz_hr' };
    default:
      throw new Error(`Unknown grid type: ${gridType}`);
  }
}

// ── Vertex → centre averaging ──────────────────────────────────────────
function vertexToCentre(vertexFlat, nx, nz) {
  // vertex flat is column-major: flat[ix + iz*nvx], nvx = nx
  // output flat is column-major: out[ix + iz*ncx]
  const ncx = nx - 1;
  const ncz = nz - 1;
  const out = new Float64Array(ncx * ncz);
  for (let ix = 0; ix < ncx; ix++) {
    for (let iz = 0; iz < ncz; iz++) {
      const v00 = vertexFlat[ix     + iz     * nx];
      const v10 = vertexFlat[(ix+1) + iz     * nx];
      const v01 = vertexFlat[ix     + (iz+1) * nx];
      const v11 = vertexFlat[(ix+1) + (iz+1) * nx];
      out[ix + iz * ncx] = 0.25 * (v00 + v10 + v01 + v11);
    }
  }
  return out;
}

// ── Derived field: stress second invariant ─────────────────────────────
function computeStressInvariant(file, nx, nz) {
  const ncx = nx - 1;
  const ncz = nz - 1;
  const sxxd = Float64Array.from(file.get('/Centers/sxxd').value);
  const szzd = Float64Array.from(file.get('/Centers/szzd').value);
  const sxz_v = Float64Array.from(file.get('/Vertices/sxz').value);
  const sxz_c = vertexToCentre(sxz_v, nx, nz);
  const n = ncx * ncz;
  const out = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const xx = sxxd[k];
    const zz = szzd[k];
    const yy = -(xx + zz);
    const xzc = sxz_c[k];
    out[k] = Math.sqrt(0.5 * (xx * xx + zz * zz + yy * yy) + xzc * xzc);
  }
  return out;
}

// ── Derived field: strain rate second invariant ────────────────────────
function computeStrainRateInvariant(file, nx, nz) {
  const ncx = nx - 1;
  const ncz = nz - 1;
  const exxd = Float64Array.from(file.get('/Centers/exxd').value);
  const ezzd = Float64Array.from(file.get('/Centers/ezzd').value);
  const exz_v = Float64Array.from(file.get('/Vertices/exz').value);
  const exz_c = vertexToCentre(exz_v, nx, nz);
  const n = ncx * ncz;
  const out = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const xx = exxd[k];
    const zz = ezzd[k];
    const yy = -(xx + zz);
    const xzc = exz_c[k];
    out[k] = Math.sqrt(0.5 * (xx * xx + zz * zz + yy * yy) + xzc * xzc);
  }
  return out;
}

// ── Air-cell masking ───────────────────────────────────────────────────
function maskAirCells(values, compoFlat) {
  for (let k = 0; k < values.length; k++) {
    if (compoFlat[k] === -1) values[k] = NaN;
  }
  return values;
}

// ── Min / max excluding NaN ────────────────────────────────────────────
export function computeMinMax(values) {
  let min = Infinity;
  let max = -Infinity;
  for (let k = 0; k < values.length; k++) {
    const v = values[k];
    if (Number.isFinite(v)) {
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  return { min: min === Infinity ? 0 : min, max: max === -Infinity ? 0 : max };
}

// ── Reshape flat → 2D (column-major: flat[ix + iz*ncx]) ───────────────
// Returns out[ix][iz] with ix=0..ncx-1, iz=0..ncz-1
function reshape2D(flat, ncx, ncz) {
  const out = new Array(ncx);
  for (let ix = 0; ix < ncx; ix++) {
    out[ix] = new Array(ncz);
    for (let iz = 0; iz < ncz; iz++) {
      out[ix][iz] = flat[ix + iz * ncx];
    }
  }
  return out;
}

// ── Main extraction ────────────────────────────────────────────────────
const deriveFunctions = {
  stress_invariant: computeStressInvariant,
  strainrate_invariant: computeStrainRateInvariant,
};

export async function extractField(filePath, fieldDef) {
  const file = await openFile(filePath);
  const meta = await readMetadata(filePath);
  const { nx, nz } = meta;
  const gt = fieldDef.gridType || 'center';
  const gi = gridInfo(nx, nz, gt);

  let flat;
  if (fieldDef.derive) {
    const fn = deriveFunctions[fieldDef.derive];
    if (!fn) throw new Error(`Unknown derive function: ${fieldDef.derive}`);
    flat = fn(file, nx, nz);
  } else {
    const ds = file.get(fieldDef.path);
    if (!ds) throw new Error(`Dataset not found: ${fieldDef.path}`);
    flat = Float64Array.from(ds.value);
  }

  // Apply offset (e.g., K → °C)
  if (fieldDef.offset) {
    for (let k = 0; k < flat.length; k++) flat[k] += fieldDef.offset;
  }

  // Mask air cells for centre-grid fields
  if (gt === 'center' && !fieldDef.discrete) {
    try {
      const compo = file.get('/VizGrid/compo').value;
      // compo might be Int8Array or similar
      const compoArr = Array.from(compo);
      maskAirCells(flat, compoArr);
    } catch (_) {
      // compo not available, skip masking
    }
  }

  // Read coordinate arrays
  const xCoords = Array.from(file.get(gi.xPath).value);
  const zCoords = Array.from(file.get(gi.zPath).value);

  const { min, max } = computeMinMax(flat);
  const values = reshape2D(flat, gi.rows, gi.cols);

  return {
    xCoords,
    zCoords,
    values,
    nx: gi.rows,
    nz: gi.cols,
    unit: fieldDef.unit || '',
    log: fieldDef.log || false,
    discrete: fieldDef.discrete || false,
    min,
    max,
  };
}

// ── Check if a dataset exists ──────────────────────────────────────────
export async function datasetExists(filePath, hdf5Path) {
  const file = await openFile(filePath);
  try {
    const ds = file.get(hdf5Path);
    return ds != null;
  } catch (_) {
    return false;
  }
}
