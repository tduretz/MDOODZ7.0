// ── PanelView ─────────────────────────────────────────────────────────
// Per-panel DOM subtree: title area, field canvas + colour bar, mini controls.
// Each PanelView owns its FieldCanvas, ColourBar, and TitleInput instances.

import { FieldCanvas } from './FieldCanvas.mjs';
import { ColourBar }   from './ColourBar.mjs';
import { TitleInput }  from './TitleInput.mjs';

const CMAP_OPTIONS = ['viridis', 'turbo', 'inferno', 'plasma', 'coolwarm'];

export class PanelView {
  constructor(model, panelState, colourMaps, containerEl) {
    this.model = model;
    this.panelState = panelState;
    this.colourMaps = colourMaps;

    // Root element
    this.el = document.createElement('div');
    this.el.className = 'panel';
    this.el.dataset.panelId = panelState.id;
    if (model.activePanelId === panelState.id) {
      this.el.classList.add('active');
    }

    // ── Title area ────────────────────────────────────────────────────
    this.titleContainer = document.createElement('div');
    this.el.appendChild(this.titleContainer);
    this.titleInput = new TitleInput(model, this.titleContainer, panelState, this.el);

    // ── Canvas row (field canvas + colour bar) ────────────────────────
    this.canvasRow = document.createElement('div');
    this.canvasRow.className = 'panel-canvas-row';

    this.fieldCanvasEl = document.createElement('canvas');
    this.fieldCanvasEl.className = 'field-canvas';
    this.canvasRow.appendChild(this.fieldCanvasEl);

    this.colourBarEl = document.createElement('canvas');
    this.colourBarEl.className = 'colourbar-canvas';
    this.colourBarEl.width = 80;
    this.colourBarEl.height = 300;
    this.canvasRow.appendChild(this.colourBarEl);

    this.el.appendChild(this.canvasRow);

    // Instantiate FieldCanvas & ColourBar bound to this panel's state
    this.fieldCanvas = new FieldCanvas(model, this.fieldCanvasEl, colourMaps, panelState);
    this.colourBar   = new ColourBar(model, this.colourBarEl, colourMaps, panelState);

    // ── Mini control row ──────────────────────────────────────────────
    this.controls = document.createElement('div');
    this.controls.className = 'panel-controls';

    // Field select
    this.fieldSelect = document.createElement('select');
    this.fieldSelect.title = 'Field';
    this.controls.appendChild(this.fieldSelect);

    // Colour map select
    this.cmapSelect = document.createElement('select');
    this.cmapSelect.title = 'Colour map';
    for (const name of CMAP_OPTIONS) {
      const opt = document.createElement('option');
      opt.value = name;
      opt.textContent = name;
      this.cmapSelect.appendChild(opt);
    }
    this.cmapSelect.value = panelState.colourMap;
    this.controls.appendChild(this.cmapSelect);

    // Range inputs
    this.rangeMin = document.createElement('input');
    this.rangeMin.type = 'text';
    this.rangeMin.placeholder = 'min';
    this.rangeMin.size = 8;
    this.controls.appendChild(this.rangeMin);

    this.rangeMax = document.createElement('input');
    this.rangeMax.type = 'text';
    this.rangeMax.placeholder = 'max';
    this.rangeMax.size = 8;
    this.controls.appendChild(this.rangeMax);

    // Auto button
    this.autoBtn = document.createElement('button');
    this.autoBtn.textContent = 'Auto';
    this.controls.appendChild(this.autoBtn);

    // Lock button
    this.lockBtn = document.createElement('button');
    this.lockBtn.textContent = panelState.rangeLocked ? '🔒' : '🔓';
    this.lockBtn.title = 'Lock/unlock range';
    this.controls.appendChild(this.lockBtn);

    this.el.appendChild(this.controls);

    // Append to container
    containerEl.appendChild(this.el);

    // ── Wire DOM events → ctrl:panel:* ────────────────────────────────
    const pid = panelState.id;

    this.fieldSelect.addEventListener('change', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:panel:field', {
        bubbles: true,
        detail: { panelId: pid, fieldName: this.fieldSelect.value }
      }));
    });

    this.cmapSelect.addEventListener('change', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:panel:cmap', {
        bubbles: true,
        detail: { panelId: pid, colourMap: this.cmapSelect.value }
      }));
    });

    const emitRange = () => {
      const min = parseFloat(this.rangeMin.value);
      const max = parseFloat(this.rangeMax.value);
      if (Number.isFinite(min) && Number.isFinite(max)) {
        this.el.dispatchEvent(new CustomEvent('ctrl:panel:range', {
          bubbles: true,
          detail: { panelId: pid, range: { min, max } }
        }));
      }
    };
    this.rangeMin.addEventListener('change', emitRange);
    this.rangeMax.addEventListener('change', emitRange);

    this.autoBtn.addEventListener('click', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:panel:range-auto', {
        bubbles: true,
        detail: { panelId: pid }
      }));
    });

    this.lockBtn.addEventListener('click', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:panel:lock', {
        bubbles: true,
        detail: { panelId: pid }
      }));
    });

    // Click panel → set active
    this.el.addEventListener('click', () => {
      model.activePanelId = pid;
    });

    // ── Model events ──────────────────────────────────────────────────
    this._onFieldsLoaded = () => this._updateFieldList();
    this._onPanelField = (e) => {
      if (e.detail.panelId === pid) this._syncAfterFieldChange();
    };
    this._onPanelCmap = (e) => {
      if (e.detail.panelId === pid) this.cmapSelect.value = panelState.colourMap;
    };
    this._onPanelRange = (e) => {
      if (e.detail.panelId === pid) this._syncRange();
    };
    this._onPanelLock = (e) => {
      if (e.detail.panelId === pid) {
        this.lockBtn.textContent = panelState.rangeLocked ? '🔒' : '🔓';
        this._syncRange();
      }
    };
    this._onActiveChange = (e) => {
      this.el.classList.toggle('active', model.activePanelId === pid);
    };

    model.addEventListener('fields-loaded',              this._onFieldsLoaded);
    model.addEventListener('panel:field-changed',         this._onPanelField);
    model.addEventListener('panel:colourmap-changed',     this._onPanelCmap);
    model.addEventListener('panel:range-changed',         this._onPanelRange);
    model.addEventListener('panel:range-locked-changed',  this._onPanelLock);
    model.addEventListener('active-panel-changed',        this._onActiveChange);

    // Initial population of field list
    this._updateFieldList();
  }

  _updateFieldList() {
    const fields = this.model.fields;
    const defs = this.model.fieldDefs;
    const prev = this.fieldSelect.value;
    this.fieldSelect.innerHTML = '';
    for (const f of fields) {
      const opt = document.createElement('option');
      opt.value = f;
      const def = defs.get(f);
      opt.textContent = def ? def.label : f;
      this.fieldSelect.appendChild(opt);
    }
    // Restore selection
    if (fields.includes(prev)) {
      this.fieldSelect.value = prev;
    } else if (this.panelState.fieldName && fields.includes(this.panelState.fieldName)) {
      this.fieldSelect.value = this.panelState.fieldName;
    }
  }

  _syncAfterFieldChange() {
    // Sync field select
    if (this.panelState.fieldName) {
      this.fieldSelect.value = this.panelState.fieldName;
    }
    this._syncRange();
    this.cmapSelect.value = this.panelState.colourMap;
  }

  _syncRange() {
    const { min, max } = this.panelState.colourRange;
    this.rangeMin.value = _fmtVal(min);
    this.rangeMax.value = _fmtVal(max);
  }

  destroy() {
    this.titleInput.destroy();
    this.model.removeEventListener('fields-loaded',              this._onFieldsLoaded);
    this.model.removeEventListener('panel:field-changed',         this._onPanelField);
    this.model.removeEventListener('panel:colourmap-changed',     this._onPanelCmap);
    this.model.removeEventListener('panel:range-changed',         this._onPanelRange);
    this.model.removeEventListener('panel:range-locked-changed',  this._onPanelLock);
    this.model.removeEventListener('active-panel-changed',        this._onActiveChange);
    this.el.remove();
  }
}

function _fmtVal(v) {
  const abs = Math.abs(v);
  if (abs === 0) return '0';
  if (abs >= 1e6 || abs < 0.01) return v.toExponential(3);
  return v.toPrecision(5);
}
