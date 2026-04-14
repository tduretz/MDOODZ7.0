// ── Application Model ──────────────────────────────────────────────────
// Single source of truth. Extends EventTarget so views subscribe via
// addEventListener and receive CustomEvent on every state change.

export class Model extends EventTarget {
  constructor() {
    super();
    this._files       = [];
    this._currentFile = null;
    this._fields      = [];
    this._currentField = null;
    this._fieldData   = null;
    this._colourMap   = 'viridis';
    this._colourRange = { min: 0, max: 1 };
    this._params      = null;
    this._loading     = false;
  }

  // ── Getters ────────────────────────────────────────────────────────
  get files()        { return this._files; }
  get currentFile()  { return this._currentFile; }
  get fields()       { return this._fields; }
  get currentField() { return this._currentField; }
  get fieldData()    { return this._fieldData; }
  get colourMap()    { return this._colourMap; }
  get colourRange()  { return this._colourRange; }
  get params()       { return this._params; }
  get loading()      { return this._loading; }

  // ── Setters (dispatch events) ───────────────────────────────────────
  set files(v) {
    this._files = v;
    this._emit('file-list-loaded');
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
    this._colourRange = { min: data.min, max: data.max };
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

  // ── Internal ───────────────────────────────────────────────────────
  _emit(name) {
    this.dispatchEvent(new CustomEvent(name, { detail: this }));
  }
}
