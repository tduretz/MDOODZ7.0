// ── ControlPanel View ─────────────────────────────────────────────────
// File dropdown, time-step slider, field dropdown, colour-map dropdown,
// min/max inputs with auto-reset.

export class ControlPanel {
  constructor(model, containerEl) {
    this.model = model;
    this.el = containerEl;
    this._build();

    model.addEventListener('file-list-loaded', () => this._updateFileList());
    model.addEventListener('file-selected',    () => this._syncFileUI());
    model.addEventListener('fields-loaded',    () => this._updateFieldList());
    model.addEventListener('field-loaded',     () => this._updateRange());
  }

  _build() {
    this.el.innerHTML = `
      <div class="ctrl-group">
        <label>File</label>
        <select id="file-select"></select>
        <input id="time-slider" type="range" min="0" max="0" value="0" step="1">
      </div>
      <div class="ctrl-group">
        <label>Field</label>
        <select id="field-select"></select>
      </div>
      <div class="ctrl-group">
        <label>Colour map</label>
        <select id="cmap-select">
          <option value="viridis" selected>viridis</option>
          <option value="turbo">turbo</option>
          <option value="inferno">inferno</option>
          <option value="plasma">plasma</option>
          <option value="coolwarm">coolwarm</option>
        </select>
      </div>
      <div class="ctrl-group">
        <label>Range</label>
        <input id="range-min" type="text" placeholder="min" size="10">
        <input id="range-max" type="text" placeholder="max" size="10">
        <button id="range-auto">Auto</button>
      </div>
    `;

    this.fileSelect  = this.el.querySelector('#file-select');
    this.timeSlider  = this.el.querySelector('#time-slider');
    this.fieldSelect = this.el.querySelector('#field-select');
    this.cmapSelect  = this.el.querySelector('#cmap-select');
    this.rangeMin    = this.el.querySelector('#range-min');
    this.rangeMax    = this.el.querySelector('#range-max');
    this.rangeAuto   = this.el.querySelector('#range-auto');

    // DOM events → dispatched as custom events on this.el
    this.fileSelect.addEventListener('change', () => {
      this.timeSlider.value = this.fileSelect.selectedIndex;
      this.el.dispatchEvent(new CustomEvent('ctrl:file-change', { detail: this.fileSelect.value }));
    });

    this.timeSlider.addEventListener('input', () => {
      this.fileSelect.selectedIndex = parseInt(this.timeSlider.value, 10);
      this.el.dispatchEvent(new CustomEvent('ctrl:file-change', { detail: this.fileSelect.value }));
    });

    this.fieldSelect.addEventListener('change', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:field-change', { detail: this.fieldSelect.value }));
    });

    this.cmapSelect.addEventListener('change', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:cmap-change', { detail: this.cmapSelect.value }));
    });

    const emitRange = () => {
      const min = parseFloat(this.rangeMin.value);
      const max = parseFloat(this.rangeMax.value);
      if (Number.isFinite(min) && Number.isFinite(max)) {
        this.el.dispatchEvent(new CustomEvent('ctrl:range-change', { detail: { min, max } }));
      }
    };
    this.rangeMin.addEventListener('change', emitRange);
    this.rangeMax.addEventListener('change', emitRange);

    this.rangeAuto.addEventListener('click', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:range-auto'));
    });
  }

  _updateFileList() {
    const files = this.model.files;
    this.fileSelect.innerHTML = files.map(f => `<option value="${f}">${f}</option>`).join('');
    this.timeSlider.max = String(Math.max(0, files.length - 1));
  }

  _syncFileUI() {
    const idx = this.model.files.indexOf(this.model.currentFile);
    if (idx >= 0) {
      this.fileSelect.selectedIndex = idx;
      this.timeSlider.value = String(idx);
    }
  }

  _updateFieldList() {
    const fields = this.model.fields;
    const prev = this.fieldSelect.value;
    this.fieldSelect.innerHTML = fields.map(f => `<option value="${f}">${f}</option>`).join('');
    // Preserve selection if available
    if (fields.includes(prev)) {
      this.fieldSelect.value = prev;
    }
  }

  _updateRange() {
    const { min, max } = this.model.colourRange;
    this.rangeMin.value = formatVal(min);
    this.rangeMax.value = formatVal(max);
  }
}

function formatVal(v) {
  const abs = Math.abs(v);
  if (abs === 0) return '0';
  if (abs >= 1e6 || abs < 0.01) return v.toExponential(3);
  return v.toPrecision(5);
}
