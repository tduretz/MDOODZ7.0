// ── Title Render ──────────────────────────────────────────────────────
// Pure function: renderTitle(template, context) → string
// Replaces ${...} tokens with formatted values from context.

import { formatTime } from './time-display.mjs';

/**
 * Replace ${token} placeholders in a template string.
 * Unknown tokens are left as-is.
 */
export function renderTitle(template, context) {
  return template.replace(/\$\{(\w+)\}/g, (match, key) => {
    if (key in context && context[key] !== undefined && context[key] !== null) {
      return String(context[key]);
    }
    return match; // unknown token → leave as-is
  });
}

/**
 * Build a full title context from params + panel state + field defs.
 */
export function buildTitleContext(params, panelState, fieldDefs, timeUnit) {
  const ctx = {};

  // Scalar params
  if (params) {
    ctx.step = params.step ?? params.Step ?? '';
    const time = params.time ?? params.Time ?? 0;
    const { formatted } = formatTime(time, timeUnit);
    ctx.time = formatted;
    const nx = params.Nx ?? params.nx ?? '';
    const nz = params.Nz ?? params.nz ?? '';
    ctx.Nx = nx;
    ctx.Nz = nz;
    ctx.resolution = nx && nz ? `${nx} × ${nz}` : '';
  }

  // Field stats
  if (panelState && panelState.fieldData) {
    const data = panelState.fieldData;
    const def = panelState.fieldName && fieldDefs ? fieldDefs.get(panelState.fieldName) : null;
    const unit = def ? def.formattedUnit : (data.unit || '');

    ctx.field = def ? def.label : (panelState.fieldName || '');
    ctx.unit = unit;
    ctx.min = formatSciWithUnit(data.min, unit);
    ctx.max = formatSciWithUnit(data.max, unit);
    ctx.p2 = formatSciWithUnit(data.pMin != null ? data.pMin : data.min, unit);
    ctx.p98 = formatSciWithUnit(data.pMax != null ? data.pMax : data.max, unit);

    // Lazy mean/median — compute on demand, cache in panelState
    const stats = getOrComputeStats(panelState);
    ctx.mean = formatSciWithUnit(stats.mean, unit);
    ctx.median = formatSciWithUnit(stats.median, unit);
  }

  return ctx;
}

/**
 * Compute mean and median from flat values array. Cache in panelState._statsCache.
 */
export function getOrComputeStats(panelState) {
  if (panelState._statsCache) return panelState._statsCache;

  const data = panelState.fieldData;
  if (!data || !data.values) {
    return { mean: 0, median: 0 };
  }

  // Flatten 2D values[ix][iz] to a 1D array, skipping NaN/null
  const flat = [];
  const vals = data.values;
  for (let ix = 0; ix < vals.length; ix++) {
    const col = vals[ix];
    for (let iz = 0; iz < col.length; iz++) {
      const v = col[iz];
      if (v !== null && v !== undefined && !Number.isNaN(v)) {
        flat.push(v);
      }
    }
  }

  if (flat.length === 0) {
    panelState._statsCache = { mean: 0, median: 0 };
    return panelState._statsCache;
  }

  // Mean
  let sum = 0;
  for (let i = 0; i < flat.length; i++) sum += flat[i];
  const mean = sum / flat.length;

  // Median (sort a copy)
  flat.sort((a, b) => a - b);
  const mid = flat.length >> 1;
  const median = (flat.length & 1) ? flat[mid] : (flat[mid - 1] + flat[mid]) / 2;

  panelState._statsCache = { mean, median };
  return panelState._statsCache;
}

/**
 * Format a numeric value with scientific notation + unit suffix.
 */
function formatSciWithUnit(v, unit) {
  if (v === null || v === undefined) return '';
  const abs = Math.abs(v);
  let formatted;
  if (abs === 0) {
    formatted = '0.00e+00';
  } else if (abs >= 1e-1 && abs < 1e4) {
    formatted = v.toFixed(2);
  } else {
    formatted = v.toExponential(2);
  }
  return unit ? `${formatted} ${unit}` : formatted;
}
