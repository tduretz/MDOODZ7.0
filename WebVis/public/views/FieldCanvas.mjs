// ── FieldCanvas View ──────────────────────────────────────────────────
// Renders 2D field data to <canvas> via ImageData pixel-fill.
// Accepts a canvas element + panelState (or falls back to model for legacy).

import { renderTitle, buildTitleContext } from '../title-render.mjs';

export class FieldCanvas {
  constructor(model, canvasEl, colourMaps, panelState) {
    this.model = model;
    this.canvas = canvasEl;
    this.ctx = canvasEl.getContext('2d');
    this.colourMaps = colourMaps;
    this.panelState = panelState || null;

    if (panelState) {
      // Panel-scoped: listen to panel events, filter by panelId
      const pid = panelState.id;
      this._onField = e => { if (e.detail.panelId === pid) this.render(); };
      this._onCmap  = e => { if (e.detail.panelId === pid) this.render(); };
      this._onRange = e => { if (e.detail.panelId === pid) this.render(); };
      this._onTitle = e => { if (e.detail.panelId === pid) this.render(); };
      this._onParams = () => this.render();
      model.addEventListener('panel:field-changed',     this._onField);
      model.addEventListener('panel:colourmap-changed',  this._onCmap);
      model.addEventListener('panel:range-changed',      this._onRange);
      model.addEventListener('panel:title-changed',      this._onTitle);
      model.addEventListener('params-loaded',            this._onParams);
    } else {
      // Legacy global mode
      this._onField = () => this.render();
      this._onCmap  = () => this.render();
      this._onRange = () => this.render();
      this._onTitle = null;
      this._onParams = null;
      model.addEventListener('field-loaded',      this._onField);
      model.addEventListener('colourmap-changed',  this._onCmap);
      model.addEventListener('range-changed',      this._onRange);
    }
  }

  render() {
    const ps = this.panelState;
    const data = ps ? ps.fieldData : this.model.fieldData;
    if (!data) return;

    const { values, nx, nz } = data;
    const width  = nx;
    const height = nz;

    this.canvas.width  = width;
    this.canvas.height = height;

    const imgData = this.ctx.createImageData(width, height);
    const pixels  = imgData.data;

    const colourRange = ps ? ps.colourRange : this.model.colourRange;
    const { min, max } = colourRange;
    const isLog      = data.log;
    const isDiscrete = data.discrete;

    let lut, phasePalette;
    if (isDiscrete) {
      phasePalette = this.colourMaps.phases || [];
    } else {
      const cmapName = ps ? ps.colourMap : this.model.colourMap;
      lut = this.colourMaps[cmapName];
      if (!lut) return;
    }

    const logMin = isLog && min > 0 ? Math.log10(min) : 0;
    const logMax = isLog && max > 0 ? Math.log10(max) : 1;
    const range  = isLog ? (logMax - logMin) : (max - min);

    // Canvas row 0 = top of screen = highest z (iz = nz-1)
    for (let row = 0; row < height; row++) {
      const iz = height - 1 - row;   // flip z so surface is at top
      for (let col = 0; col < width; col++) {
        const ix = col;
        const idx = (row * width + col) * 4;
        const val = values[ix][iz];

        if (val === null || val === undefined || Number.isNaN(val)) {
          pixels[idx] = 0; pixels[idx+1] = 0; pixels[idx+2] = 0; pixels[idx+3] = 0;
          continue;
        }

        if (isDiscrete) {
          const phase = Math.round(val);
          if (phase < 0) {
            pixels[idx] = 0; pixels[idx+1] = 0; pixels[idx+2] = 0; pixels[idx+3] = 0;
          } else {
            const c = phasePalette[phase % phasePalette.length] || [128, 128, 128];
            pixels[idx] = c[0]; pixels[idx+1] = c[1]; pixels[idx+2] = c[2]; pixels[idx+3] = 255;
          }
          continue;
        }

        let mapped = val;
        if (isLog) {
          if (val <= 0) { pixels[idx] = 0; pixels[idx+1] = 0; pixels[idx+2] = 0; pixels[idx+3] = 0; continue; }
          mapped = Math.log10(val);
        }

        const t = range > 0 ? (mapped - (isLog ? logMin : min)) / range : 0.5;
        const ci = Math.max(0, Math.min(255, Math.round(t * 255)));
        const c = lut[ci];
        pixels[idx] = c[0]; pixels[idx+1] = c[1]; pixels[idx+2] = c[2]; pixels[idx+3] = 255;
      }
    }

    this.ctx.putImageData(imgData, 0, 0);

    // ── Title overlay ─────────────────────────────────────────────────
    if (this.panelState && this.panelState.title) {
      const titleCtx = buildTitleContext(
        this.model.params, this.panelState,
        this.model.fieldDefs, this.model.timeUnit
      );
      const text = renderTitle(this.panelState.title, titleCtx);
      if (text) {
        const ctx = this.ctx;
        const fontSize = Math.max(12, Math.round(height * 0.04));
        ctx.font = `${fontSize}px sans-serif`;
        const metrics = ctx.measureText(text);
        const pad = 4;
        const tx = pad;
        const ty = fontSize + pad;
        // Semi-transparent background
        ctx.fillStyle = 'rgba(0,0,0,0.55)';
        ctx.fillRect(tx - pad, 0, metrics.width + pad * 2, ty + pad);
        // White text
        ctx.fillStyle = '#fff';
        ctx.textBaseline = 'top';
        ctx.fillText(text, tx, pad);
      }
    }
  }

  destroy() {
    if (this.panelState) {
      this.model.removeEventListener('panel:field-changed',     this._onField);
      this.model.removeEventListener('panel:colourmap-changed',  this._onCmap);
      this.model.removeEventListener('panel:range-changed',      this._onRange);
      if (this._onTitle) this.model.removeEventListener('panel:title-changed', this._onTitle);
      if (this._onParams) this.model.removeEventListener('params-loaded', this._onParams);
    } else {
      this.model.removeEventListener('field-loaded',      this._onField);
      this.model.removeEventListener('colourmap-changed',  this._onCmap);
      this.model.removeEventListener('range-changed',      this._onRange);
    }
  }
}
