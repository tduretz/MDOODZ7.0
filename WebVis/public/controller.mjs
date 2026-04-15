// ── Controller ────────────────────────────────────────────────────────
// Wires control panel DOM events to API fetch calls and model updates.
// Supports both legacy single-panel events and panel-scoped events.

export class Controller {
  constructor(model, controlsEl) {
    this.model = model;
    this._abortCtrl = null;  // AbortController for in-flight fetches
    this._overlayCache = new Map();  // filename → overlayDataObject

    // Legacy global control events (from ControlPanel / global controls)
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
    controlsEl.addEventListener('ctrl:layout-change', e => this.changeLayout(e.detail));
    controlsEl.addEventListener('ctrl:dataset-change', e => this.selectDataset(e.detail));

    // Panel-scoped events bubble from #panels-grid, listen at document level
    document.addEventListener('ctrl:panel:field', e => {
      this.selectPanelField(e.detail.panelId, e.detail.fieldName);
    });
    document.addEventListener('ctrl:panel:cmap', e => {
      model.setPanelColourMap(e.detail.panelId, e.detail.colourMap);
    });
    document.addEventListener('ctrl:panel:range', e => {
      model.setPanelRange(e.detail.panelId, e.detail.range);
    });
    document.addEventListener('ctrl:panel:lock', e => {
      const panel = model.getPanel(e.detail.panelId);
      if (panel) {
        model.setPanelRangeLocked(e.detail.panelId, !panel.rangeLocked);
      }
    });
    document.addEventListener('ctrl:panel:range-auto', e => {
      const panel = model.getPanel(e.detail.panelId);
      if (panel && panel.fieldData) {
        const min = panel.fieldData.pMin != null ? panel.fieldData.pMin : panel.fieldData.min;
        const max = panel.fieldData.pMax != null ? panel.fieldData.pMax : panel.fieldData.max;
        model.setPanelRange(e.detail.panelId, { min, max });
      }
    });
    document.addEventListener('ctrl:panel:title', e => {
      model.setPanelTitle(e.detail.panelId, e.detail.template);
    });

    // Overlay events
    document.addEventListener('ctrl:panel:overlay-toggle', e => {
      const { panelId, layerType, enabled } = e.detail;
      model.setPanelOverlay(panelId, layerType, { enabled });
      this._ensureOverlayData(panelId);
    });
    document.addEventListener('ctrl:panel:overlay-config', e => {
      const { panelId, layerType, ...config } = e.detail;
      model.setPanelOverlay(panelId, layerType, config);
    });
  }

  async init() {
    this.model.loading = true;
    try {
      // Fetch available datasets (directories)
      const dsRes = await fetch('/api/datasets');
      const { datasets, active } = await dsRes.json();
      this.model.datasets = datasets;
      this.model.activeDataset = active;

      // Fetch files for the active dataset
      const res = await fetch('/api/files');
      const { files } = await res.json();
      this.model.fileList = files;
      this.model.files = files.map(f => f.name);
      if (files.length > 0) {
        await this.selectFile(files[0].name);
      }
    } finally {
      this.model.loading = false;
    }
  }

  _newAbort() {
    if (this._abortCtrl) this._abortCtrl.abort();
    this._abortCtrl = new AbortController();
    return this._abortCtrl.signal;
  }

  async fetchParams(filename, signal) {
    const res = await fetch(`/api/params/${encodeURIComponent(filename)}`, { signal });
    const params = await res.json();
    this.model.params = params;
    return params;
  }

  async selectPanelField(panelId, fieldName, signal) {
    const file = this.model.currentFile;
    if (!file) return;
    // Restore per-field colour map if stored
    const panel = this.model.getPanel(panelId);
    if (panel) {
      const savedMap = panel.fieldColourMaps.get(fieldName);
      if (savedMap && savedMap !== panel.colourMap) {
        this.model.setPanelColourMap(panelId, savedMap);
      }
    }
    const res = await fetch(`/api/field-data/${encodeURIComponent(file)}/${encodeURIComponent(fieldName)}`, { signal });
    const data = await res.json();
    this.model.setPanelField(panelId, fieldName, data);
  }

  async selectFile(filename) {
    const signal = this._newAbort();
    this.model.loading = true;
    this._overlayCache.clear();  // invalidate overlay cache on file change
    try {
      this.model.currentFile = filename;

      // Fetch params + fields in parallel
      const [, fieldsRes] = await Promise.all([
        this.fetchParams(filename, signal),
        fetch(`/api/fields/${encodeURIComponent(filename)}`, { signal }),
      ]);
      const { fields } = await fieldsRes.json();

      const fieldDefs = new Map();
      const fieldNames = fields.map(f => {
        fieldDefs.set(f.name, { label: f.label, formattedUnit: f.formattedUnit });
        return f.name;
      });
      this.model._fieldDefs = fieldDefs;
      this.model.fields = fieldNames;

      // Re-fetch each panel's field in parallel
      const fetches = this.model.panels.map(panel => {
        const prev = panel.fieldName;
        const field = (prev && fieldNames.includes(prev)) ? prev : fieldNames[0];
        if (field) return this.selectPanelField(panel.id, field, signal);
        return Promise.resolve();
      });
      await Promise.all(fetches);

      // Re-fetch overlay data for panels with enabled layers
      for (const panel of this.model.panels) {
        const hasEnabled = Object.values(panel.overlays).some(o => o.enabled);
        if (hasEnabled) this._ensureOverlayData(panel.id);
      }
    } catch (err) {
      if (err.name === 'AbortError') return; // superseded by newer request
      throw err;
    } finally {
      this.model.loading = false;
    }
  }

  async selectField(fieldName) {
    // Legacy: delegate to first panel
    await this.selectPanelField(this.model.panels[0].id, fieldName);
  }

  async selectDataset(name) {
    this.model.loading = true;
    try {
      // Tell server to switch active dataset
      await fetch('/api/dataset', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name }),
      });
      this.model.activeDataset = name;

      // Re-fetch file list for the new dataset
      const res = await fetch('/api/files');
      const { files } = await res.json();
      this.model.fileList = files;
      this.model.files = files.map(f => f.name);
      if (files.length > 0) {
        await this.selectFile(files[0].name);
      }
    } finally {
      this.model.loading = false;
    }
  }

  changeLayout(preset) {
    const prevCount = this.model.panels.length;
    this.model.setLayout(preset);

    // Trigger field fetches for newly added panels
    const file = this.model.currentFile;
    const fieldNames = this.model.fields;
    if (file && fieldNames.length > 0) {
      for (const panel of this.model.panels) {
        if (!panel.fieldName && fieldNames[0]) {
          this.selectPanelField(panel.id, fieldNames[0]);
        }
      }
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

  /** Fetch overlay data for all enabled layers on a panel that aren't cached yet. */
  async _ensureOverlayData(panelId) {
    const file = this.model.currentFile;
    if (!file) return;
    const panel = this.model.getPanel(panelId);
    if (!panel) return;

    // On first request for this file, fetch ALL layers to discover availability
    if (!this._overlayCache.has(file)) {
      const allLayers = ['phases','temperature','velocity','topo','director','sigma1','edot1','melt'];
      try {
        const res = await fetch(`/api/overlay-data/${encodeURIComponent(file)}?layers=${allLayers.join(',')}`);
        const data = await res.json();
        this._overlayCache.set(file, data);
      } catch (err) {
        if (err.name !== 'AbortError') console.error('Overlay fetch error:', err);
        return;
      }
    }

    const cached = this._overlayCache.get(file);
    const available = new Set(Object.keys(cached));
    this.model.setPanelOverlayAvailability(panelId, available);
    this._pushOverlayToPanel(panelId, file);
  }

  /** Push cached overlay data to model for a specific panel. */
  _pushOverlayToPanel(panelId, file) {
    const cached = this._overlayCache.get(file);
    if (cached) {
      this.model.setPanelOverlayData(panelId, cached);
    }
  }
}
