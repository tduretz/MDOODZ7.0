// ── FieldCanvas View ──────────────────────────────────────────────────
// Renders 2D field data to <canvas> at display resolution (HiDPI-aware).
// Uses an offscreen canvas for field data (nearest-neighbour upscale),
// then draws title, axis labels, ticks, and integrated colour bar at
// full resolution.  Margins are computed dynamically from label widths so
// nothing is ever clipped.

import { renderTitle, buildTitleContext } from '../title-render.mjs';

// Fixed geometry (CSS px)
const CBAR_W    = 14;   // colour-bar strip width
const CBAR_PAD  = 6;    // gap between data rect and colour-bar strip
const LEGEND_BOX = 12;  // phase legend colour box size
const TICK_LEN  = 4;
const FONT_SIZE = 11;
const FONT      = `${FONT_SIZE}px monospace`;
const AXIS_FONT = `${FONT_SIZE}px sans-serif`;
const TITLE_FONT = 'bold 12px sans-serif';

// Minimum margins (always reserved, even when labels are tiny)
const MIN_MARGIN = { top: 22, right: 30, bottom: 30, left: 30 };

// Generous estimates exported for app.mjs aspect-ratio layout
export const MARGIN_ESTIMATE = { top: 24, right: 120, bottom: 38, left: 66 };

/** Read canvas theme colours from CSS custom properties (reactive to theme switch). */
function _getTheme() {
  const s = getComputedStyle(document.documentElement);
  return {
    bg:       s.getPropertyValue('--canvas-bg').trim()        || '#1e1e2e',
    nullCol:  s.getPropertyValue('--canvas-null').trim()      || '#1e1e2e',
    tick:     s.getPropertyValue('--canvas-tick').trim()       || '#aaa',
    tickLine: s.getPropertyValue('--canvas-tick-line').trim()  || '#666',
    label:    s.getPropertyValue('--canvas-label').trim()      || '#999',
    cbarText: s.getPropertyValue('--canvas-cbar-text').trim()  || '#ccc',
    title:    s.getPropertyValue('--canvas-title').trim()      || '#89b4fa',
    frame:    s.getPropertyValue('--canvas-frame').trim()      || '#000',
  };
}

/** Parse a hex colour string to [r, g, b]. */
function _hexToRgbArr(hex) {
  const m = hex.match(/^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i);
  return m ? [parseInt(m[1],16), parseInt(m[2],16), parseInt(m[3],16)] : [30,30,46];
}

export class FieldCanvas {
  constructor(model, canvasEl, colourMaps, panelState) {
    this.model = model;
    this.canvas = canvasEl;
    this.ctx = canvasEl.getContext('2d');
    this.colourMaps = colourMaps;
    this.panelState = panelState || null;
    this._offscreen = document.createElement('canvas');

    if (panelState) {
      const pid = panelState.id;
      this._onField   = e => { if (e.detail.panelId === pid) this.render(); };
      this._onCmap    = e => { if (e.detail.panelId === pid) this.render(); };
      this._onRange   = e => { if (e.detail.panelId === pid) this.render(); };
      this._onTitle   = e => { if (e.detail.panelId === pid) this.render(); };
      this._onParams  = () => this.render();
      this._onSpatial = () => this.render();
      this._onPhaseConfig = () => this.render();
      this._onTheme = () => this.render();
      this._onViewBounds = e => { if (e.detail.panelId === pid) this.render(); };
      model.addEventListener('panel:field-changed',      this._onField);
      model.addEventListener('panel:colourmap-changed',   this._onCmap);
      model.addEventListener('panel:range-changed',       this._onRange);
      model.addEventListener('panel:title-changed',       this._onTitle);
      model.addEventListener('params-loaded',             this._onParams);
      model.addEventListener('spatial-unit-changed',      this._onSpatial);
      model.addEventListener('phase-config-changed',      this._onPhaseConfig);
      model.addEventListener('theme-changed',             this._onTheme);
      model.addEventListener('panel:view-bounds-changed', this._onViewBounds);
    } else {
      this._onField   = () => this.render();
      this._onCmap    = () => this.render();
      this._onRange   = () => this.render();
      this._onTitle   = null;
      this._onParams  = null;
      this._onSpatial = null;
      model.addEventListener('field-loaded',       this._onField);
      model.addEventListener('colourmap-changed',   this._onCmap);
      model.addEventListener('range-changed',       this._onRange);
    }
  }

  // ── Main render ─────────────────────────────────────────────────────
  render() {
    const ps   = this.panelState;
    const data = ps ? ps.fieldData : this.model.fieldData;
    if (!data) return;

    const { values, nx, nz } = data;

    // HiDPI canvas sizing
    const displayW = this.canvas.clientWidth;
    const displayH = this.canvas.clientHeight;
    if (displayW <= 0 || displayH <= 0) return;

    const dpr = window.devicePixelRatio || 1;
    this.canvas.width  = Math.round(displayW * dpr);
    this.canvas.height = Math.round(displayH * dpr);
    const ctx = this.ctx;
    ctx.save();
    ctx.scale(dpr, dpr);

    // ── Theme colours ─────────────────────────────────────────────────
    const theme = _getTheme();
    const nullRgb = _hexToRgbArr(theme.nullCol);

    // ── Colour range & LUT ────────────────────────────────────────────
    const colourRange = ps ? ps.colourRange : this.model.colourRange;
    const { min, max } = colourRange;
    const isLog      = data.log;
    const isDiscrete = data.discrete;

    // ── View bounds (zoom/crop) ───────────────────────────────────────
    const vb = ps ? ps.viewBounds : null;
    let col0 = 0, col1 = nx - 1, row0 = 0, row1 = nz - 1;
    const hasCoords = data.xCoords && data.zCoords &&
                      data.xCoords.length > 0 && data.zCoords.length > 0;
    if (vb && hasCoords) {
      const xC = data.xCoords, zC = data.zCoords;
      if (vb.xMin != null) col0 = Math.max(0, xC.findIndex(x => x >= vb.xMin));
      if (vb.xMax != null) { for (let i = nx - 1; i >= 0; i--) { if (xC[i] <= vb.xMax) { col1 = i; break; } } }
      if (vb.zMin != null) row0 = Math.max(0, zC.findIndex(z => z >= vb.zMin));
      if (vb.zMax != null) { for (let i = nz - 1; i >= 0; i--) { if (zC[i] <= vb.zMax) { row1 = i; break; } } }
      if (col0 > col1) { col0 = 0; col1 = nx - 1; }
      if (row0 > row1) { row0 = 0; row1 = nz - 1; }
    }
    const cropNx = col1 - col0 + 1;
    const cropNz = row1 - row0 + 1;

    let lut, phasePalette;
    if (isDiscrete) {
      phasePalette = this.colourMaps.phases || [];
    } else {
      const cmapName = ps ? ps.colourMap : this.model.colourMap;
      lut = this.colourMaps[cmapName];
      if (!lut) { ctx.restore(); return; }
    }

    // ── Pre-compute all label strings for margin measurement ──────────
    const hasCbar  = (!isDiscrete && lut) || (isDiscrete && phasePalette);

    const spatialUnit = this.model.spatialUnit || 'km';
    const divisor   = spatialUnit === 'km' ? 1e3 : 1;
    const unitLabel = spatialUnit === 'km' ? 'km' : 'm';

    // Colour-bar labels
    const cbarStrings = [];
    if (!isDiscrete && lut) {
      const nTicks = 5;
      for (let i = 0; i < nTicks; i++) {
        const frac = i / (nTicks - 1);
        if (isLog) {
          const lMin = Math.log10(min > 0 ? min : 1);
          const lMax = Math.log10(max > 0 ? max : 1);
          cbarStrings.push(_fmtLog(lMin + frac * (lMax - lMin)));
        } else {
          cbarStrings.push(_fmtVal(min + frac * (max - min)));
        }
      }
    } else if (isDiscrete && phasePalette) {
      // Scan for visible phases in cropped area
      const seen = new Set();
      for (let col = col0; col <= col1; col++)
        for (let iz = row0; iz <= row1; iz++) {
          const v = values[col][iz];
          if (v !== null && v !== undefined && !Number.isNaN(v)) {
            const p = Math.round(v);
            if (p >= 0) seen.add(p);
          }
        }
      for (const p of Array.from(seen).sort((a, b) => a - b)) {
        const cfg = this.model.getPhaseConfig(p, phasePalette);
        cbarStrings.push(cfg.name);
      }
    }

    // Z-axis (left) labels — use cropped coordinate range
    const zStrings = [];
    let nTicksZ = 0;
    if (hasCoords) {
      const zCoords = data.zCoords;
      const czMin = zCoords[row0], czMax = zCoords[row1];
      const czExt = czMax - czMin;
      nTicksZ = Math.min(6, Math.max(3, Math.floor(200 / 60)));
      for (let i = 0; i < nTicksZ; i++) {
        const frac = i / (nTicksZ - 1);
        zStrings.push(_fmtAxis((czMin + frac * czExt) / divisor));
      }
    }

    // ── Measure text & compute dynamic margins ────────────────────────
    ctx.font = FONT;
    const maxCbarW = cbarStrings.length > 0
      ? Math.max(...cbarStrings.map(s => ctx.measureText(s).width)) : 0;
    const maxZTickW = zStrings.length > 0
      ? Math.max(...zStrings.map(s => ctx.measureText(s).width)) : 0;

    // Rotated colour-bar label (field [unit]) — adds ~FONT_SIZE to right margin
    const rotLabelExtra = (hasCbar && !isDiscrete) ? FONT_SIZE + 4 : 0;

    const m = {
      top:    MIN_MARGIN.top,                                          // title row
      right:  hasCbar
                ? (isDiscrete
                    ? CBAR_PAD + LEGEND_BOX + 4 + Math.ceil(maxCbarW) + 8
                    : CBAR_PAD + CBAR_W + 4 + Math.ceil(maxCbarW) + 6 + rotLabelExtra)
                : MIN_MARGIN.right,
      bottom: hasCoords
                ? TICK_LEN + 2 + FONT_SIZE + 4 + FONT_SIZE + 2       // ticks + "x [km]"
                : MIN_MARGIN.bottom,
      left:   hasCoords
                ? FONT_SIZE + 6 + Math.ceil(maxZTickW) + TICK_LEN + 3 // "z [km]" + ticks
                : MIN_MARGIN.left,
    };
    // Enforce minimums
    m.top    = Math.max(m.top,    MIN_MARGIN.top);
    m.right  = Math.max(m.right,  MIN_MARGIN.right);
    m.bottom = Math.max(m.bottom, MIN_MARGIN.bottom);
    m.left   = Math.max(m.left,   MIN_MARGIN.left);

    // ── Data area (aspect-ratio fit within remaining space) ───────────
    const availW = displayW - m.left - m.right;
    const availH = displayH - m.top  - m.bottom;
    if (availW <= 0 || availH <= 0) { ctx.restore(); return; }

    const dataAspect = cropNx / cropNz;
    let dW, dH;
    if (availW / availH > dataAspect) {
      dH = availH;
      dW = Math.round(dH * dataAspect);
    } else {
      dW = availW;
      dH = Math.round(dW / dataAspect);
    }
    // Centre data area in available space
    const dataX = m.left + (availW - dW) / 2;
    const dataY = m.top  + (availH - dH) / 2;

    // ── Clear background ──────────────────────────────────────────────
    ctx.fillStyle = theme.bg;
    ctx.fillRect(0, 0, displayW, displayH);

    // ── Field data via offscreen ImageData (cropped region) ───────────
    const off = this._offscreen;
    off.width  = cropNx;
    off.height = cropNz;
    const offCtx  = off.getContext('2d');
    const imgData = offCtx.createImageData(cropNx, cropNz);
    const pixels  = imgData.data;

    const logMin = isLog && min > 0 ? Math.log10(min) : 0;
    const logMax = isLog && max > 0 ? Math.log10(max) : 1;
    const range  = isLog ? (logMax - logMin) : (max - min);

    for (let r = 0; r < cropNz; r++) {
      const iz = row1 - r;  // top row of crop = row1 (highest z)
      for (let c = 0; c < cropNx; c++) {
        const col = col0 + c;
        const idx = (r * cropNx + c) * 4;
        const val = values[col][iz];

        if (val === null || val === undefined || Number.isNaN(val)) {
          pixels[idx] = nullRgb[0]; pixels[idx+1] = nullRgb[1]; pixels[idx+2] = nullRgb[2]; pixels[idx+3] = 255;
          continue;
        }
        if (isDiscrete) {
          const phase = Math.round(val);
          if (phase < 0) {
            pixels[idx] = nullRgb[0]; pixels[idx+1] = nullRgb[1]; pixels[idx+2] = nullRgb[2]; pixels[idx+3] = 255;
          } else {
            const cfg = this.model.getPhaseConfig(phase, phasePalette);
            const cc = cfg.color;
            pixels[idx] = cc[0]; pixels[idx+1] = cc[1]; pixels[idx+2] = cc[2]; pixels[idx+3] = 255;
          }
          continue;
        }
        let mapped = val;
        if (isLog) {
          if (val <= 0) { pixels[idx] = nullRgb[0]; pixels[idx+1] = nullRgb[1]; pixels[idx+2] = nullRgb[2]; pixels[idx+3] = 255; continue; }
          mapped = Math.log10(val);
        }
        const t  = range > 0 ? (mapped - (isLog ? logMin : min)) / range : 0.5;
        const ci = Math.max(0, Math.min(255, Math.round(t * 255)));
        const rgb = lut[ci];
        pixels[idx] = rgb[0]; pixels[idx+1] = rgb[1]; pixels[idx+2] = rgb[2]; pixels[idx+3] = 255;
      }
    }
    offCtx.putImageData(imgData, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, dataX, dataY, dW, dH);
    ctx.imageSmoothingEnabled = true;

    // ── Black frame around data area ──────────────────────────────────
    ctx.strokeStyle = theme.frame;
    ctx.lineWidth = 1;
    ctx.strokeRect(dataX - 0.5, dataY - 0.5, dW + 1, dH + 1);

    // ── Integrated colour bar ─────────────────────────────────────────
    if (!isDiscrete && lut) {
      this._drawCbar(ctx, lut, colourRange, data, dataX + dW, dataY, dH, theme);
    }
    if (isDiscrete && phasePalette) {
      this._drawPhaseLegend(ctx, phasePalette, values, col0, col1, row0, row1, dataX + dW, dataY, dH, theme);
    }

    // ── Axis ticks + labels ───────────────────────────────────────────
    if (hasCoords) {
      this._drawAxes(ctx, data, dataX, dataY, dW, dH, unitLabel, divisor, theme, col0, col1, row0, row1);
    }

    // ── Title ─────────────────────────────────────────────────────────
    this._drawTitle(ctx, displayW, theme);

    ctx.restore();
  }

  // ── Colour bar (continuous) ─────────────────────────────────────────
  _drawCbar(ctx, lut, colourRange, data, rightEdge, top, height, theme) {
    const cbarX = rightEdge + CBAR_PAD;

    // Gradient strip
    for (let y = 0; y < height; y++) {
      const t  = 1 - y / height;
      const ci = Math.max(0, Math.min(255, Math.round(t * 255)));
      const c  = lut[ci];
      ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
      ctx.fillRect(cbarX, top + y, CBAR_W, 1);
    }

    // Tick labels
    const { min: rMin, max: rMax } = colourRange;
    const isLog = data.log;
    const nTicks = 5;
    ctx.fillStyle = theme.cbarText;
    ctx.font = FONT;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';
    const labelX = cbarX + CBAR_W + 4;

    for (let i = 0; i < nTicks; i++) {
      const frac = i / (nTicks - 1);
      const y = top + height - frac * height;
      let label;
      if (isLog) {
        const lMin = Math.log10(rMin > 0 ? rMin : 1);
        const lMax = Math.log10(rMax > 0 ? rMax : 1);
        label = _fmtLog(lMin + frac * (lMax - lMin));
      } else {
        label = _fmtVal(rMin + frac * (rMax - rMin));
      }
      ctx.fillText(label, labelX, y);
    }

    // Rotated "field [unit]" label along right of tick labels
    const ps = this.panelState;
    const def = ps ? this.model.fieldDefs.get(ps.fieldName) : null;
    const fieldLabel = def ? def.label : '';
    const unitStr    = def ? def.formattedUnit : (data.unit || '');
    const cbarLabel  = unitStr ? `${fieldLabel} [${unitStr}]` : fieldLabel;
    if (cbarLabel) {
      ctx.save();
      ctx.fillStyle = theme.label;
      ctx.font = `${FONT_SIZE}px sans-serif`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      // Measure widest tick label to place rotated text beyond them
      ctx.font = FONT;
      let maxTickW = 0;
      for (let i = 0; i < nTicks; i++) {
        const frac = i / (nTicks - 1);
        let lbl;
        if (isLog) {
          const lMin = Math.log10(rMin > 0 ? rMin : 1);
          const lMax = Math.log10(rMax > 0 ? rMax : 1);
          lbl = _fmtLog(lMin + frac * (lMax - lMin));
        } else {
          lbl = _fmtVal(rMin + frac * (rMax - rMin));
        }
        maxTickW = Math.max(maxTickW, ctx.measureText(lbl).width);
      }
      ctx.font = `${FONT_SIZE}px sans-serif`;
      const rotX = labelX + maxTickW + FONT_SIZE + 6;
      ctx.translate(rotX, top + height / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(cbarLabel, 0, 0);
      ctx.restore();
    }
  }

  // ── Phase box legend (discrete phases — only visible ones) ────────────
  _drawPhaseLegend(ctx, palette, values, c0, c1, r0, r1, rightEdge, top, height, theme) {
    // Scan data for unique phase indices (skip air = -1 and NaN)
    const seen = new Set();
    for (let col = c0; col <= c1; col++) {
      for (let iz = r0; iz <= r1; iz++) {
        const v = values[col][iz];
        if (v !== null && v !== undefined && !Number.isNaN(v)) {
          const p = Math.round(v);
          if (p >= 0) seen.add(p);
        }
      }
    }
    const phases = Array.from(seen).sort((a, b) => a - b);
    if (phases.length === 0) return;

    const lineH = LEGEND_BOX + 4;  // box height + vertical gap
    const totalH = phases.length * lineH;
    const startY = top + (height - totalH) / 2;  // vertically centre
    const boxX   = rightEdge + CBAR_PAD;

    ctx.font = `${FONT_SIZE - 1}px sans-serif`;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'left';

    for (let i = 0; i < phases.length; i++) {
      const p = phases[i];
      const cfg = this.model.getPhaseConfig(p, palette);
      const y = startY + i * lineH;

      // Coloured box
      const c = cfg.color;
      ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
      ctx.fillRect(boxX, y, LEGEND_BOX, LEGEND_BOX);

      // Phase name
      ctx.fillStyle = theme.cbarText;
      ctx.fillText(cfg.name, boxX + LEGEND_BOX + 4, y + LEGEND_BOX / 2);
    }
  }

  // ── Axes: X on bottom, Z on left ──────────────────────────────────
  _drawAxes(ctx, data, dataX, dataY, dW, dH, unitLabel, divisor, theme, c0, c1, r0, r1) {
    const xCoords = data.xCoords;
    const zCoords = data.zCoords;
    const xMin = xCoords[c0], xMax = xCoords[c1], xExt = xMax - xMin;
    const zMin = zCoords[r0], zMax = zCoords[r1], zExt = zMax - zMin;

    ctx.strokeStyle = theme.tickLine;
    ctx.lineWidth = 1;

    // ── Bottom axis (X) ──────────────────────────────────────────────
    const nTicksX = Math.min(6, Math.max(3, Math.floor(dW / 80)));
    ctx.font = FONT;
    ctx.fillStyle = theme.tick;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const yBase = dataY + dH;

    for (let i = 0; i < nTicksX; i++) {
      const frac = i / (nTicksX - 1);
      const px   = dataX + frac * (dW - 1);
      const val  = (xMin + frac * xExt) / divisor;
      ctx.beginPath(); ctx.moveTo(px, yBase); ctx.lineTo(px, yBase + TICK_LEN); ctx.stroke();
      ctx.fillText(_fmtAxis(val), px, yBase + TICK_LEN + 1);
    }

    // "x [km]" axis label
    ctx.fillStyle = theme.label;
    ctx.font = AXIS_FONT;
    ctx.fillText(`x [${unitLabel}]`, dataX + dW / 2, yBase + TICK_LEN + FONT_SIZE + 4);

    // ── Left axis (Z) ────────────────────────────────────────────────
    const nTicksZ = Math.min(6, Math.max(3, Math.floor(dH / 60)));
    ctx.font = FONT;
    ctx.fillStyle = theme.tick;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';

    for (let i = 0; i < nTicksZ; i++) {
      const frac = i / (nTicksZ - 1);
      const py   = dataY + (1 - frac) * (dH - 1);
      const val  = (zMin + frac * zExt) / divisor;
      ctx.beginPath(); ctx.moveTo(dataX, py); ctx.lineTo(dataX - TICK_LEN, py); ctx.stroke();
      ctx.fillText(_fmtAxis(val), dataX - TICK_LEN - 2, py);
    }

    // Rotated "z [km]" axis label
    ctx.save();
    ctx.fillStyle = theme.label;
    ctx.font = AXIS_FONT;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    ctx.translate(dataX - TICK_LEN - 2, dataY + dH / 2);
    ctx.rotate(-Math.PI / 2);
    // Measure widest Z tick to know how far left we are
    ctx.font = FONT;
    const maxZW = Math.max(...Array.from({ length: nTicksZ }, (_, i) => {
      const frac = i / (nTicksZ - 1);
      return ctx.measureText(_fmtAxis((zMin + frac * zExt) / divisor)).width;
    }));
    ctx.font = AXIS_FONT;
    ctx.fillText(`z [${unitLabel}]`, 0, -maxZW - 2);
    ctx.restore();
  }

  // ── Title centred at top ────────────────────────────────────────────
  _drawTitle(ctx, displayW, theme) {
    const ps = this.panelState;
    if (!ps) return;
    const titleCtx = buildTitleContext(
      this.model.params, ps, this.model.fieldDefs, this.model.timeUnit
    );
    const text = renderTitle(ps.title, titleCtx);
    if (!text) return;

    ctx.fillStyle = theme.title;
    ctx.font = TITLE_FONT;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    ctx.fillText(text, displayW / 2, 4);
  }

  destroy() {
    if (this.panelState) {
      this.model.removeEventListener('panel:field-changed',      this._onField);
      this.model.removeEventListener('panel:colourmap-changed',   this._onCmap);
      this.model.removeEventListener('panel:range-changed',       this._onRange);
      if (this._onTitle)   this.model.removeEventListener('panel:title-changed', this._onTitle);
      if (this._onParams)  this.model.removeEventListener('params-loaded',       this._onParams);
      if (this._onSpatial) this.model.removeEventListener('spatial-unit-changed', this._onSpatial);
      if (this._onPhaseConfig) this.model.removeEventListener('phase-config-changed', this._onPhaseConfig);
      if (this._onTheme) this.model.removeEventListener('theme-changed', this._onTheme);
      if (this._onViewBounds) this.model.removeEventListener('panel:view-bounds-changed', this._onViewBounds);
    } else {
      this.model.removeEventListener('field-loaded',       this._onField);
      this.model.removeEventListener('colourmap-changed',   this._onCmap);
      this.model.removeEventListener('range-changed',       this._onRange);
    }
  }
}

// ── Formatting helpers ────────────────────────────────────────────────

const _SUP = { '0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹','-':'⁻','.':'·','+':'⁺' };

/** Format a log-scale value as "10" + Unicode superscript, e.g. "10⁻¹⁵·⁰". */
function _fmtLog(exponent) {
  const s = exponent.toFixed(1);
  const sup = s.replace(/./g, ch => _SUP[ch] || ch);
  return `10${sup}`;
}

function _fmtVal(v) {
  const abs = Math.abs(v);
  if (abs === 0) return '0';
  if (abs >= 1e6 || abs < 0.01) return v.toExponential(2);
  return v.toPrecision(4);
}

function _fmtAxis(v) {
  const abs = Math.abs(v);
  if (abs === 0) return '0';
  if (abs >= 1000) return v.toFixed(0);
  if (abs >= 1)    return v.toFixed(1);
  return v.toPrecision(3);
}
