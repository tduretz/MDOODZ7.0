// ── Controller ────────────────────────────────────────────────────────
// Wires control panel DOM events to API fetch calls and model updates.

export class Controller {
  constructor(model, controlsEl) {
    this.model = model;

    // Listen to custom events from ControlPanel
    controlsEl.addEventListener('ctrl:file-change',  e => this.selectFile(e.detail));
    controlsEl.addEventListener('ctrl:field-change', e => this.selectField(e.detail));
    controlsEl.addEventListener('ctrl:cmap-change',  e => { model.colourMap = e.detail; });
    controlsEl.addEventListener('ctrl:range-change', e => { model.colourRange = e.detail; });
    controlsEl.addEventListener('ctrl:range-auto',   () => this._resetRange());
  }

  async init() {
    this.model.loading = true;
    try {
      const res = await fetch('/api/files');
      const { files } = await res.json();
      this.model.files = files;
      if (files.length > 0) {
        await this.selectFile(files[0]);
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
      this.model.fields = fields;

      // Preserve field choice if available, else pick first
      const prev = this.model.currentField;
      const field = (prev && fields.includes(prev)) ? prev : fields[0];
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
      this.model.colourRange = { min: data.min, max: data.max };
    }
  }
}
