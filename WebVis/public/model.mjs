// ── Application Model ──────────────────────────────────────────────────
// Single source of truth. Extends EventTarget so views subscribe via
// addEventListener and receive CustomEvent on every state change.

let _nextPanelId = 1;

export function createDefaultOverlays() {
  return {
    phases:      { enabled: false, color: '#ffffff', lineWidth: 1.5 },
    temperature: { enabled: false, color: '#000000', lineWidth: 1.5, dT: 200 },
    topo:        { enabled: false, color: '#000000', lineWidth: 2 },
    velocity:    { enabled: false, color: '#ffff00', density: 8 },
    director:    { enabled: false, color: '#ffffff', density: 6, lengthScale: 1.5 },
    sigma1:      { enabled: false, color: '#ffffff', density: 6, lengthScale: 1.5 },
    edot1:       { enabled: false, color: '#888888', density: 6, lengthScale: 1.5 },
    melt:        { enabled: false, color: '#ff4444', lineWidth: 2, levels: [0.01, 0.05, 0.1] },
  };
}

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
    viewBounds: null,      // null = full extent, or { xMin, xMax, zMin, zMax } in SI metres
    overlays: createDefaultOverlays(),
    overlayData: null,     // cached overlay data from server for this panel
    overlayAvailable: new Set(),  // layer types available for current file
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
    this._spatialUnit = 'km';   // 'm' or 'km'

    // Panel state
    this._layout      = '1x1';
    this._panels      = [createDefaultPanel()];
    this._activePanelId = this._panels[0].id;
    this._stashedPanels = [];  // panels preserved across layout round-trips

    // Phase config: Map<phaseId, { name: string, color: [r,g,b] }>
    this._phaseConfig = new Map();

    // Theme: 'dark' or 'light'
    this._theme = 'dark';

    // Dataset (directory) selection
    this._datasets      = [];   // [{ name, fileCount }]
    this._activeDataset = null; // dataset name string
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
  get spatialUnit()  { return this._spatialUnit; }
  get layout()       { return this._layout; }
  get panels()       { return this._panels; }
  get activePanelId(){ return this._activePanelId; }
  get theme()        { return this._theme; }

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

    // Shrink: move excess panels to stash (preserve state)
    while (this._panels.length > target) {
      const panel = this._panels.pop();
      if (this._activePanelId === panel.id) {
        this._activePanelId = this._panels[0].id;
      }
      this._stashedPanels.push(panel);
      this._emit('panel-removed', { panelId: panel.id });
    }

    // Grow: restore from stash first, then create new defaults
    while (this._panels.length < target) {
      if (this._stashedPanels.length > 0) {
        const restored = this._stashedPanels.pop();
        this._panels.push(restored);
        this._emit('panel-added', { panelId: restored.id, restored: true });
      } else {
        this.addPanel();
      }
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
    p.fieldData = null;  // release previous allocation for GC
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

  setPanelViewBounds(panelId, bounds) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.viewBounds = bounds;  // null or { xMin, xMax, zMin, zMax }
    this._emit('panel:view-bounds-changed', { panelId });
  }

  applyViewBoundsToAll(sourcePanelId) {
    const src = this.getPanel(sourcePanelId);
    if (!src) return;
    const bounds = src.viewBounds ? { ...src.viewBounds } : null;
    for (const p of this._panels) {
      if (p.id === sourcePanelId) continue;
      p.viewBounds = bounds ? { ...bounds } : null;
      this._emit('panel:view-bounds-changed', { panelId: p.id });
    }
  }

  /** Merge overlay config for a specific layer in a panel. */
  setPanelOverlay(panelId, layerType, config) {
    const p = this.getPanel(panelId);
    if (!p || !p.overlays[layerType]) return;
    Object.assign(p.overlays[layerType], config);
    this._emit('panel:overlay-changed', { panelId, layerType });
  }

  /** Store which overlay layer types are available for the current file. */
  setPanelOverlayAvailability(panelId, availableSet) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.overlayAvailable = availableSet instanceof Set ? availableSet : new Set(availableSet);
    this._emit('panel:overlay-availability-changed', { panelId });
  }

  /** Store fetched overlay data for a panel. */
  setPanelOverlayData(panelId, data) {
    const p = this.getPanel(panelId);
    if (!p) return;
    p.overlayData = data;
    this._emit('panel:overlay-data-changed', { panelId });
  }

  // ── Global setters (dispatch events) ────────────────────────────────

  // Dataset (directory) management
  get datasets()      { return this._datasets; }
  get activeDataset() { return this._activeDataset; }

  set datasets(v) {
    this._datasets = v;
    this._emit('datasets-loaded');
  }

  set activeDataset(v) {
    this._activeDataset = v;
    this._emit('dataset-changed');
  }

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

  set spatialUnit(v) {
    this._spatialUnit = v;
    this._emit('spatial-unit-changed');
  }

  set theme(v) {
    if (v !== 'dark' && v !== 'light') return;
    this._theme = v;
    if (v === 'light') {
      document.documentElement.setAttribute('data-theme', 'light');
    } else {
      document.documentElement.removeAttribute('data-theme');
    }
    this._emit('theme-changed');
  }

  // ── Phase config ───────────────────────────────────────────────────
  get phaseConfig() { return this._phaseConfig; }

  /** Get config for a single phase, with defaults if not customised. */
  getPhaseConfig(phaseId, defaultPalette) {
    if (this._phaseConfig.has(phaseId)) return this._phaseConfig.get(phaseId);
    const DEFAULTS = { 0: 'Crust', 1: 'Lithosphere', 2: 'Asthenosphere' };
    const name  = DEFAULTS[phaseId] || `Phase ${phaseId}`;
    const color = defaultPalette
      ? (defaultPalette[phaseId % defaultPalette.length] || [128,128,128])
      : [128,128,128];
    return { name, color };
  }

  setPhaseConfig(phaseId, cfg) {
    const prev = this._phaseConfig.get(phaseId) || {};
    this._phaseConfig.set(phaseId, { ...prev, ...cfg });
    this._emit('phase-config-changed', { phaseId });
  }

  // ── Internal ───────────────────────────────────────────────────────
  _emit(name, detail) {
    this.dispatchEvent(new CustomEvent(name, { detail: detail || this }));
  }
}
