// ── Controller ────────────────────────────────────────────────────────
// Wires control panel DOM events to API fetch calls and model updates.

export class Controller {
  constructor(model, controlsEl) {
    this.model = model;

    // Listen to custom events from ControlPanel
    controlsEl.addEventListener('ctrl:file-change',  e => this.selectFile(e.detail));
    controlsEl.addEventListener('ctrl:field-change', e => this.selectField(e.detail));
    controlsEl.addEventListener('ctrl:cmap-change',  e => {
      model.fieldColourMaps.set(model.currentField, e.detail);
      model.colourMap = e.detail;
    });
    controlsEl.addEventListener('ctrl:range-change', e => { model.colourRange = e.detail; });
    controlsEl.addEventListener('ctrl:range-auto',   () => this._resetRange());
    controlsEl.addEventListener('ctrl:range-lock',   () => {
      model.rangeLocked = !model.rangeLocked;
      if (!model.rangeLocked) this._resetRange();
    });
    controlsEl.addEventListener('ctrl:time-unit',    e => { model.timeUnit = e.detail; });
  }

  async init() {
    this.model.loading = true;
    try {
      const res = await fetch('/api/files');
      const { files } = await res.json();
      // files is now [{ name, step, time }]
      this.model.fileList = files;
      this.model.files = files.map(f => f.name);
      if (files.length > 0) {
        await this.selectFile(files[0].name);
      }
    } finally {
      this.model.loading = false;
    }
  }

  async selectFile(filename) {
    this.model.loading = true;
    try {
      this.model.currentFile = filename;
      const res = await fetch(`/api/fields/${encodeURIComponent(filename)}`);
      const { fields, params } = await res.json();
      this.model.params = params;

      // fields is now [{ name, label, formattedUnit }]
      const fieldDefs = new Map();
      const fieldNames = fields.map(f => {
        fieldDefs.set(f.name, { label: f.label, formattedUnit: f.formattedUnit });
        return f.name;
      });
      this.model._fieldDefs = fieldDefs;
      this.model.fields = fieldNames;

      // Preserve field choice if available, else pick first
      const prev = this.model.currentField;
      const field = (prev && fieldNames.includes(prev)) ? prev : fieldNames[0];
      if (field) {
        await this.selectField(field);
      }
    } finally {
      this.model.loading = false;
    }
  }

  async selectField(fieldName) {
    this.model.loading = true;
    try {
      this.model.currentField = fieldName;
      // Restore per-field colour map if stored
      const savedMap = this.model.fieldColourMaps.get(fieldName);
      if (savedMap && savedMap !== this.model.colourMap) {
        this.model.colourMap = savedMap;
      }
      const file = this.model.currentFile;
      const res = await fetch(`/api/field-data/${encodeURIComponent(file)}/${encodeURIComponent(fieldName)}`);
      const data = await res.json();
      this.model.setFieldData(data);
    } finally {
      this.model.loading = false;
    }
  }

  _resetRange() {
    const data = this.model.fieldData;
    if (data) {
      const min = data.pMin != null ? data.pMin : data.min;
      const max = data.pMax != null ? data.pMax : data.max;
      this.model.colourRange = { min, max };
    }
  }
}
