// ── overlays.mjs — client-side overlay math utilities ──────────────────

/**
 * Marching-squares ISO-contour for a single level.
 * Returns array of {x1,z1,x2,z2} line segments in physical (SI) coordinates.
 * @param {Float64Array|number[]} flat Column-major flat array (ncx × ncz)
 * @param {number[]} xCoords Centre x-coordinates (length ncx)
 * @param {number[]} zCoords Centre z-coordinates (length ncz)
 * @param {number} level Iso-value
 * @returns {{x1:number,z1:number,x2:number,z2:number}[]}
 */
export function marchingSquares(flat, xCoords, zCoords, level) {
  const ncx = xCoords.length;
  const ncz = zCoords.length;
  const segments = [];

  // Value at (ix, iz) — x-fastest layout: flat[ix + iz * ncx]
  const val = (ix, iz) => {
    const v = flat[ix + iz * ncx];
    return (v !== v) ? level - 1 : v;   // NaN → below threshold
  };

  // Lerp factor: where the isoline crosses between a and b
  const frac = (a, b) => {
    const d = b - a;
    return Math.abs(d) < 1e-30 ? 0.5 : (level - a) / d;
  };

  for (let ix = 0; ix < ncx - 1; ix++) {
    const x0 = xCoords[ix],  x1 = xCoords[ix + 1];
    for (let iz = 0; iz < ncz - 1; iz++) {
      const z0 = zCoords[iz], z1 = zCoords[iz + 1];
      // Corner values: BL=0, BR=1, TR=2, TL=3
      const v0 = val(ix, iz),     v1 = val(ix + 1, iz);
      const v2 = val(ix + 1, iz + 1), v3 = val(ix, iz + 1);

      let idx = 0;
      if (v0 >= level) idx |= 1;
      if (v1 >= level) idx |= 2;
      if (v2 >= level) idx |= 4;
      if (v3 >= level) idx |= 8;

      if (idx === 0 || idx === 15) continue;

      // Edge midpoints with linear interpolation
      // Bottom: v0→v1, Right: v1→v2, Top: v3→v2, Left: v0→v3
      const bx = x0 + frac(v0, v1) * (x1 - x0), bz = z0;
      const rx = x1,                              rz = z0 + frac(v1, v2) * (z1 - z0);
      const tx = x0 + frac(v3, v2) * (x1 - x0), tz = z1;
      const lx = x0,                              lz = z0 + frac(v0, v3) * (z1 - z0);

      // Lookup table for 16 cases
      switch (idx) {
        case 1: case 14: segments.push({ x1: bx, z1: bz, x2: lx, z2: lz }); break;
        case 2: case 13: segments.push({ x1: bx, z1: bz, x2: rx, z2: rz }); break;
        case 3: case 12: segments.push({ x1: lx, z1: lz, x2: rx, z2: rz }); break;
        case 4: case 11: segments.push({ x1: rx, z1: rz, x2: tx, z2: tz }); break;
        case 5:
          segments.push({ x1: bx, z1: bz, x2: rx, z2: rz });
          segments.push({ x1: lx, z1: lz, x2: tx, z2: tz });
          break;
        case 6: case 9: segments.push({ x1: bx, z1: bz, x2: tx, z2: tz }); break;
        case 7: case 8: segments.push({ x1: lx, z1: lz, x2: tx, z2: tz }); break;
        case 10:
          segments.push({ x1: bx, z1: bz, x2: lx, z2: lz });
          segments.push({ x1: rx, z1: rz, x2: tx, z2: tz });
          break;
      }
    }
  }
  return segments;
}

/**
 * Extract contour lines for multiple levels at once.
 * @returns {Map<number, {x1:number,z1:number,x2:number,z2:number}[]>}
 */
export function extractContours(flat, xCoords, zCoords, levels) {
  const map = new Map();
  for (const level of levels) {
    map.set(level, marchingSquares(flat, xCoords, zCoords, level));
  }
  return map;
}

/**
 * Closed-form 2×2 symmetric eigensolver.
 * Input: matrix [[a, b],[b, c]].
 * Returns { lambda1, lambda2, cos, sin } where (cos, sin) is the eigenvector
 * of the largest eigenvalue (lambda1 ≥ lambda2).
 */
export function eigen2x2(a, b, c) {
  const trace = a + c;
  const det   = a * c - b * b;
  const disc  = Math.sqrt(Math.max(trace * trace * 0.25 - det, 0));
  const lambda1 = trace * 0.5 + disc;
  const lambda2 = trace * 0.5 - disc;

  let ex, ez;
  if (Math.abs(b) > 1e-30) {
    ex = lambda1 - c;
    ez = b;
  } else if (a >= c) {
    ex = 1; ez = 0;
  } else {
    ex = 0; ez = 1;
  }
  const len = Math.sqrt(ex * ex + ez * ez) || 1;
  return { lambda1, lambda2, cos: ex / len, sin: ez / len };
}

/**
 * Interpolate staggered Vx (vxNx × vxNz) to cell centres (ncx × ncz).
 * x-fastest flat layout: flat[ix + iz * width].
 */
export function interpVxToCenter(Vx, vxNx, vxNz) {
  const ncx = vxNx - 1, ncz = vxNz - 2;
  const out = new Float64Array(ncx * ncz);
  for (let ix = 0; ix < ncx; ix++) {
    for (let iz = 0; iz < ncz; iz++) {
      const i00 = ix     + iz       * vxNx;
      const i01 = ix     + (iz + 1) * vxNx;
      out[ix + iz * ncx] = 0.25 * (Vx[i00] + Vx[i00 + 1] + Vx[i01] + Vx[i01 + 1]);
    }
  }
  return out;
}

/**
 * Interpolate staggered Vz (vzNx × vzNz) to cell centres (ncx × ncz).
 * x-fastest flat layout: flat[ix + iz * width].
 */
export function interpVzToCenter(Vz, vzNx, vzNz) {
  const ncx = vzNx - 2, ncz = vzNz - 1;
  const out = new Float64Array(ncx * ncz);
  for (let ix = 0; ix < ncx; ix++) {
    for (let iz = 0; iz < ncz; iz++) {
      const i00 = ix     + iz       * vzNx;
      const i01 = ix     + (iz + 1) * vzNx;
      out[ix + iz * ncx] = 0.25 * (Vz[i00] + Vz[i00 + 1] + Vz[i01] + Vz[i01 + 1]);
    }
  }
  return out;
}

/**
 * Place inline labels along contour segments at periodic pixel intervals.
 * @param {{x1:number,z1:number,x2:number,z2:number}[]} segments  Contour segments
 * @param {number} level  Contour level value (for label text)
 * @param {function} pixMapper  (x,z) → {px, pz} in canvas pixels
 * @param {number} intervalPx  Desired pixel spacing between labels
 * @returns {{x:number, z:number, angle:number, text:string}[]}
 */
export function inlineContourLabels(segments, level, pixMapper, intervalPx = 150) {
  if (segments.length === 0) return [];

  // Join segments into polylines by matching endpoints
  const chains = _joinSegments(segments);
  const text = String(Math.round(level));
  const labels = [];

  for (const chain of chains) {
    let accum = 0;
    // Place first label at half-interval offset
    let nextAt = intervalPx * 0.5;
    for (let k = 1; k < chain.length; k++) {
      const p0 = pixMapper(chain[k - 1].x, chain[k - 1].z);
      const p1 = pixMapper(chain[k].x, chain[k].z);
      const dx = p1.px - p0.px, dz = p1.pz - p0.pz;
      const segLen = Math.sqrt(dx * dx + dz * dz);
      if (segLen < 0.5) continue;

      while (accum + segLen >= nextAt) {
        const t = (nextAt - accum) / segLen;
        const mx = chain[k - 1].x + t * (chain[k].x - chain[k - 1].x);
        const mz = chain[k - 1].z + t * (chain[k].z - chain[k - 1].z);
        labels.push({ x: mx, z: mz, angle: Math.atan2(dz, dx), text });
        nextAt += intervalPx;
      }
      accum += segLen;
    }
  }
  return labels;
}

// ── Internal: join segments into polyline chains ───────────────────────
function _joinSegments(segments) {
  const EPS = 1e-12;
  const key = (x, z) => `${x.toFixed(10)},${z.toFixed(10)}`;

  // Build adjacency
  const adj = new Map();
  const addEdge = (k, seg, end) => {
    if (!adj.has(k)) adj.set(k, []);
    adj.get(k).push({ seg, end });
  };
  for (const s of segments) {
    const k1 = key(s.x1, s.z1), k2 = key(s.x2, s.z2);
    addEdge(k1, s, 1);
    addEdge(k2, s, 2);
  }

  const used = new Set();
  const chains = [];

  for (const s of segments) {
    if (used.has(s)) continue;
    used.add(s);
    const chain = [{ x: s.x1, z: s.z1 }, { x: s.x2, z: s.z2 }];

    // Extend forward from chain end
    let extended = true;
    while (extended) {
      extended = false;
      const last = chain[chain.length - 1];
      const k = key(last.x, last.z);
      const nbrs = adj.get(k);
      if (!nbrs) break;
      for (const { seg, end } of nbrs) {
        if (used.has(seg)) continue;
        used.add(seg);
        if (end === 1) chain.push({ x: seg.x2, z: seg.z2 });
        else chain.push({ x: seg.x1, z: seg.z1 });
        extended = true;
        break;
      }
    }

    // Extend backward from chain start
    extended = true;
    while (extended) {
      extended = false;
      const first = chain[0];
      const k = key(first.x, first.z);
      const nbrs = adj.get(k);
      if (!nbrs) break;
      for (const { seg, end } of nbrs) {
        if (used.has(seg)) continue;
        used.add(seg);
        if (end === 1) chain.unshift({ x: seg.x2, z: seg.z2 });
        else chain.unshift({ x: seg.x1, z: seg.z1 });
        extended = true;
        break;
      }
    }

    chains.push(chain);
  }
  return chains;
}
