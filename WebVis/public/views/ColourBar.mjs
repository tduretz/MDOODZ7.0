// ── ColourBar View ────────────────────────────────────────────────────
// Draws vertical gradient bar + tick labels alongside the field canvas.
// Accepts a canvas element + panelState (or falls back to model for legacy).

export class ColourBar {
  constructor(model, canvasEl, colourMaps, panelState) {
    this.model = model;
    this.canvas = canvasEl;
    this.ctx = canvasEl.getContext('2d');
    this.colourMaps = colourMaps;
    this.panelState = panelState || null;

    if (panelState) {
      const pid = panelState.id;
      model.addEventListener('panel:field-changed',        e => { if (e.detail.panelId === pid) this.render(); });
      model.addEventListener('panel:colourmap-changed',     e => { if (e.detail.panelId === pid) this.render(); });
      model.addEventListener('panel:range-changed',         e => { if (e.detail.panelId === pid) this.render(); });
      model.addEventListener('panel:range-locked-changed',  e => { if (e.detail.panelId === pid) this.render(); });
    } else {
      model.addEventListener('field-loaded',         () => this.render());
      model.addEventListener('colourmap-changed',     () => this.render());
      model.addEventListener('range-changed',         () => this.render());
      model.addEventListener('range-locked-changed',  () => this.render());
    }
  }

  render() {
    const ps = this.panelState;
    const data = ps ? ps.fieldData : this.model.fieldData;
    if (!data) return;

    const { min, max, unit, log: isLog, discrete: isDiscrete } = data;
    const fieldName = ps ? ps.fieldName : this.model.currentField;
    const def = this.model.fieldDefs.get(fieldName);
    const label = def ? def.label : (fieldName || '');
    const fUnit = def ? def.formattedUnit : (unit || '');
    const rangeLocked = ps ? ps.rangeLocked : this.model.rangeLocked;
    const colourRange = ps ? ps.colourRange : this.model.colourRange;
    const cmapName = ps ? ps.colourMap : this.model.colourMap;

    const w = this.canvas.width;
    const h = this.canvas.height;
    this.ctx.clearRect(0, 0, w, h);

    const barX = 0;
    const barW = 20;
    const barTop = 10;
    const barBot = h - 30;
    const barH = barBot - barTop;

    if (isDiscrete) {
      const palette = this.colourMaps.phases || [];
      const n = palette.length;
      const cellH = barH / n;
      for (let i = 0; i < n; i++) {
        const y = barTop + (n - 1 - i) * cellH;
        const c = palette[i];
        this.ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
        this.ctx.fillRect(barX, y, barW, cellH);
      }
      this.ctx.fillStyle = '#ccc';
      this.ctx.font = '10px monospace';
      this.ctx.textBaseline = 'middle';
      for (let i = 0; i < n; i++) {
        const y = barTop + (n - 1 - i) * cellH + cellH / 2;
        this.ctx.fillText(String(i), barX + barW + 4, y);
      }
      return;
    }

    const lut = this.colourMaps[cmapName];
    if (!lut) return;

    // Draw gradient bar
    for (let y = 0; y < barH; y++) {
      const t = 1 - y / barH; // bottom=0, top=1
      const ci = Math.max(0, Math.min(255, Math.round(t * 255)));
      const c = lut[ci];
      this.ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
      this.ctx.fillRect(barX, barTop + y, barW, 1);
    }

    // Tick labels
    const { min: rMin, max: rMax } = colourRange;
    const nTicks = 5;
    this.ctx.fillStyle = '#ccc';
    this.ctx.font = '10px monospace';
    this.ctx.textBaseline = 'middle';

    for (let i = 0; i < nTicks; i++) {
      const frac = i / (nTicks - 1);
      const y = barBot - frac * barH;

      let tickLabel;
      if (isLog) {
        const logVal = Math.log10(rMin > 0 ? rMin : 1) + frac * (Math.log10(rMax > 0 ? rMax : 1) - Math.log10(rMin > 0 ? rMin : 1));
        tickLabel = `10^${logVal.toFixed(1)}`;
      } else {
        const val = rMin + frac * (rMax - rMin);
        tickLabel = formatVal(val);
      }
      this.ctx.fillText(tickLabel, barX + barW + 4, y);
    }

    // Unit label (scientific formatted)
    const unitStr = fUnit || unit || '';
    if (unitStr) {
      this.ctx.fillStyle = '#aaa';
      this.ctx.font = '10px sans-serif';
      this.ctx.textAlign = 'center';
      this.ctx.fillText(unitStr, barX + barW / 2, h - 14);
      this.ctx.textAlign = 'start';
    }

    // Field label
    if (label) {
      this.ctx.fillStyle = '#aaa';
      this.ctx.font = '10px sans-serif';
      this.ctx.textAlign = 'center';
      this.ctx.fillText(label, barX + barW / 2, h - 3);
      this.ctx.textAlign = 'start';
    }

    // Lock indicator
    if (rangeLocked) {
      this.ctx.fillStyle = '#ff9';
      this.ctx.font = '12px sans-serif';
      this.ctx.fillText('🔒', barX + barW + 4, barBot + 14);
    }
  }
}

function formatVal(v) {
  const abs = Math.abs(v);
  if (abs === 0) return '0';
  if (abs >= 1e6 || abs < 0.01) return v.toExponential(2);
  return v.toPrecision(4);
}
