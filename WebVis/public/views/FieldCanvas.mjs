// ── FieldCanvas View ──────────────────────────────────────────────────
// Renders 2D field data to <canvas> via ImageData pixel-fill.

export class FieldCanvas {
  constructor(model, canvasEl, colourMaps) {
    this.model = model;
    this.canvas = canvasEl;
    this.ctx = canvasEl.getContext('2d');
    this.colourMaps = colourMaps;

    model.addEventListener('field-loaded',      () => this.render());
    model.addEventListener('colourmap-changed',  () => this.render());
    model.addEventListener('range-changed',      () => this.render());
  }

  render() {
    const data = this.model.fieldData;
    if (!data) return;

    // values[ix][iz]: ix=x (horizontal, 0..nx-1), iz=z (vertical, 0..nz-1)
    const { values, nx, nz } = data;
    const width  = nx;   // x → canvas horizontal
    const height = nz;   // z → canvas vertical

    this.canvas.width  = width;
    this.canvas.height = height;

    const imgData = this.ctx.createImageData(width, height);
    const pixels  = imgData.data;

    const { min, max } = this.model.colourRange;
    const isLog      = data.log;
    const isDiscrete = data.discrete;

    let lut, phasePalette;
    if (isDiscrete) {
      phasePalette = this.colourMaps.phases || [];
    } else {
      lut = this.colourMaps[this.model.colourMap];
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
  }
}
