// ── Application Model ──────────────────────────────────────────────────
// Single source of truth. Extends EventTarget so views subscribe via
// addEventListener and receive CustomEvent on every state change.

export class Model extends EventTarget {
  constructor() {
    super();
    this._files       = [];
    this._currentFile = null;
    this._fields      = [];
    this._fieldDefs   = new Map();  // name → { label, formattedUnit }
    this._currentField = null;
    this._fieldData   = null;
    this._colourMap   = 'viridis';
    this._colourRange = { min: 0, max: 1 };
    this._params      = null;
    this._loading     = false;
    this._rangeLocked = false;
    this._fieldColourMaps = new Map();
    this._pMin        = 0;
    this._pMax        = 1;
    this._fileList    = [];      // [{ name, step, time }]
    this._timeUnit    = null;    // null=auto, 'yr', 'ka', 'Ma'
  }

  // ── Getters ────────────────────────────────────────────────────────
  get files()        { return this._files; }
  get currentFile()  { return this._currentFile; }
  get fields()       { return this._fields; }
  get fieldDefs()    { return this._fieldDefs; }
  get currentField() { return this._currentField; }
  get fieldData()    { return this._fieldData; }
  get colourMap()    { return this._colourMap; }
  get colourRange()  { return this._colourRange; }
  get params()       { return this._params; }
  get loading()      { return this._loading; }
  get rangeLocked()  { return this._rangeLocked; }
  get fieldColourMaps() { return this._fieldColourMaps; }
  get pMin()         { return this._pMin; }
  get pMax()         { return this._pMax; }
  get fileList()     { return this._fileList; }
  get timeUnit()     { return this._timeUnit; }

  // ── Setters (dispatch events) ───────────────────────────────────────
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

  set currentField(v) {
    this._currentField = v;
  }

  setFieldData(data) {
    this._fieldData = data;
    this._pMin = data.pMin != null ? data.pMin : data.min;
    this._pMax = data.pMax != null ? data.pMax : data.max;
    if (!this._rangeLocked) {
      this._colourRange = { min: this._pMin, max: this._pMax };
    }
    this._emit('field-loaded');
  }

  set colourMap(v) {
    this._colourMap = v;
    this._emit('colourmap-changed');
  }

  set colourRange(v) {
    this._colourRange = v;
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
    this._rangeLocked = v;
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
  _emit(name) {
    this.dispatchEvent(new CustomEvent(name, { detail: this }));
  }
}
