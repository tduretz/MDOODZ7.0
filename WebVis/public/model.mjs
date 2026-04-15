// ── Application Model ──────────────────────────────────────────────────
// Single source of truth. Extends EventTarget so views subscribe via
// addEventListener and receive CustomEvent on every state change.

let _nextPanelId = 1;

export function createDefaultPanel(id) {
  return {
    id: id ?? _nextPanelId++,
    fieldName: null,
    colourMap: 'viridis',
    colourRange: { min: 0, max: 1 },
    rangeLocked: false,
    fieldData: null,
    title: '${field} at t = ${time} — min = ${min}  max = ${max}',
    fieldColourMaps: new Map(),
    pMin: 0,
    pMax: 1,
    _statsCache: null,     // { mean, median } — invalidated on field change
  };
}

const LAYOUT_PANEL_COUNT = { '1x1': 1, '1x2': 2, '2x1': 2, '2x2': 4 };

export class Model extends EventTarget {
  constructor() {
    super();
    this._files       = [];
    this._currentFile = null;
    this._fields      = [];
    this._fieldDefs   = new Map();  // name → { label, formattedUnit }
    this._params      = null;
    this._loading     = false;
    this._fileList    = [];      // [{ name, step, time, nx, nz, dt }]
    this._timeUnit    = null;    // null=auto, 'yr', 'ka', 'Ma'

    // Panel state
    this._layout      = '1x1';
    this._panels      = [createDefaultPanel()];
    this._activePanelId = this._panels[0].id;
  }

  // ── Getters ────────────────────────────────────────────────────────
  get files()        { return this._files; }
  get currentFile()  { return this._currentFile; }
  get fields()       { return this._fields; }
  get fieldDefs()    { return this._fieldDefs; }
  get params()       { return this._params; }
  get loading()      { return this._loading; }
  get fileList()     { return this._fileList; }
  get timeUnit()     { return this._timeUnit; }
  get layout()       { return this._layout; }
  get panels()       { return this._panels; }
  get activePanelId(){ return this._activePanelId; }

  // ── Legacy single-panel compat getters (delegate to first panel) ───
  get currentField() { return this._panels[0]?.fieldName ?? null; }
  get fieldData()    { return this._panels[0]?.fieldData ?? null; }
  get colourMap()    { return this._panels[0]?.colourMap ?? 'viridis'; }
  get colourRange()  { return this._panels[0]?.colourRange ?? { min: 0, max: 1 }; }
  get rangeLocked()  { return this._panels[0]?.rangeLocked ?? false; }
  get fieldColourMaps() { return this._panels[0]?.fieldColourMaps ?? new Map(); }
  get pMin()         { return this._panels[0]?.pMin ?? 0; }
  get pMax()         { return this._panels[0]?.pMax ?? 1; }

  // ── Panel CRUD ─────────────────────────────────────────────────────
  getPanel(id) {
    return this._panels.find(p => p.id === id) || null;
  }

  addPanel() {
    const panel = createDefaultPanel();
    this._panels.push(panel);
    this._emit('panel-added', { panelId: panel.id });
    return panel;
  }

  removePanel(id) {
    const idx = this._panels.findIndex(p => p.id === id);
    if (idx === -1 || this._panels.length <= 1) return;
    this._panels.splice(idx, 1);
    if (this._activePanelId === id) {
      this._activePanelId = this._panels[0].id;
    }
    this._emit('panel-removed', { panelId: id });
  }

  setLayout(preset) {
    if (!LAYOUT_PANEL_COUNT[preset]) return;
    this._layout = preset;
    const target = LAYOUT_PANEL_COUNT[preset];
    // Add / remove panels to match target count
    while (this._panels.length < target) {
      this.addPanel();
    }
    while (this._panels.length > target) {
      this.removePanel(this._panels[this._panels.length - 1].id);
    }
    this._emit('layout-changed');
  }

  set activePanelId(id) {
    if (this.getPanel(id)) {
      this._activePanelId = id;
      this._emit('active-panel-changed', { panelId: id });
    }
  }

  // ── Panel-scoped setters ───────────────────────────────────────────
  setPanelField(panelId, fieldName, data) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.fieldName = fieldName;
    p.fieldData = data;
    p.pMin = data.pMin != null ? data.pMin : data.min;
    p.pMax = data.pMax != null ? data.pMax : data.max;
    p._statsCache = null;  // invalidate stats
    if (!p.rangeLocked) {
      p.colourRange = { min: p.pMin, max: p.pMax };
    }
    this._emit('panel:field-changed', { panelId });
  }

  setPanelColourMap(panelId, cmap) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.colourMap = cmap;
    p.fieldColourMaps.set(p.fieldName, cmap);
    this._emit('panel:colourmap-changed', { panelId });
  }

  setPanelRange(panelId, range) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.colourRange = range;
    this._emit('panel:range-changed', { panelId });
  }

  setPanelRangeLocked(panelId, locked) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.rangeLocked = locked;
    if (!locked && p.fieldData) {
      p.colourRange = { min: p.pMin, max: p.pMax };
    }
    this._emit('panel:range-locked-changed', { panelId });
  }

  setPanelTitle(panelId, template) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.title = template;
    this._emit('panel:title-changed', { panelId });
  }

  // ── Global setters (dispatch events) ────────────────────────────────
  set files(v) {
    this._files = v;
  }

  set currentFile(v) {
    this._currentFile = v;
    this._emit('file-selected');
  }

  set fields(v) {
    this._fields = v;
    this._emit('fields-loaded');
  }

  // Legacy single-panel setters — delegate to first panel
  set currentField(v) {
    if (this._panels[0]) this._panels[0].fieldName = v;
  }

  setFieldData(data) {
    this.setPanelField(this._panels[0].id, this._panels[0].fieldName, data);
    // Also emit legacy event for any non-panel-aware views
    this._emit('field-loaded');
  }

  set colourMap(v) {
    if (this._panels[0]) {
      this._panels[0].colourMap = v;
      this._panels[0].fieldColourMaps.set(this._panels[0].fieldName, v);
    }
    this._emit('colourmap-changed');
  }

  set colourRange(v) {
    if (this._panels[0]) this._panels[0].colourRange = v;
    this._emit('range-changed');
  }

  set params(v) {
    this._params = v;
    this._emit('params-loaded');
  }

  set loading(v) {
    this._loading = v;
    this._emit(v ? 'loading-start' : 'loading-end');
  }

  set rangeLocked(v) {
    if (this._panels[0]) this._panels[0].rangeLocked = v;
    this._emit('range-locked-changed');
  }

  set fileList(v) {
    this._fileList = v;
    this._emit('file-list-loaded');
  }

  set timeUnit(v) {
    this._timeUnit = v;
    this._emit('time-unit-changed');
  }

  // ── Internal ───────────────────────────────────────────────────────
  _emit(name, detail) {
    this.dispatchEvent(new CustomEvent(name, { detail: detail || this }));
  }
}
