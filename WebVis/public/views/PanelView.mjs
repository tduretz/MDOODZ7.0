// ── PanelView ─────────────────────────────────────────────────────────
// Per-panel DOM subtree: title (DOM, above grid), single field canvas
// (with integrated colour bar + axis ticks), and mini controls.

import { FieldCanvas } from './FieldCanvas.mjs';
import { TitleInput }  from './TitleInput.mjs';
// Title rendering moved into FieldCanvas (rendered on canvas for screenshots)

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

    // ── Title area (input + dropdown) ──────────────────────────────────
    this.titleContainer = document.createElement('div');
    this.el.appendChild(this.titleContainer);
    this.titleInput = new TitleInput(model, this.titleContainer, panelState, this.el);

    // Title is rendered inside the canvas by FieldCanvas._drawTitle()

    // ── Canvas area (single canvas with integrated cbar + ticks) ──────
    this.canvasRow = document.createElement('div');
    this.canvasRow.className = 'panel-canvas-row';

    this.fieldCanvasEl = document.createElement('canvas');
    this.fieldCanvasEl.className = 'field-canvas';
    this.canvasRow.appendChild(this.fieldCanvasEl);

    this.el.appendChild(this.canvasRow);

    // Instantiate FieldCanvas bound to this panel's state
    this.fieldCanvas = new FieldCanvas(model, this.fieldCanvasEl, colourMaps, panelState);

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

    // Phase config button + popover (shown only for discrete fields)
    this.phaseBtn = document.createElement('button');
    this.phaseBtn.textContent = 'Phases';
    this.phaseBtn.title = 'Edit phase names & colours';
    this.phaseBtn.className = 'phase-config-btn';
    this.phaseBtn.style.display = 'none';
    this.controls.appendChild(this.phaseBtn);

    this.phasePopover = document.createElement('div');
    this.phasePopover.className = 'phase-popover';
    this.phasePopover.style.display = 'none';
    this.el.appendChild(this.phasePopover);

    this.phaseBtn.addEventListener('click', (e) => {
      e.stopPropagation();
      const visible = this.phasePopover.style.display !== 'none';
      this.phasePopover.style.display = visible ? 'none' : 'block';
      if (!visible) this._buildPhaseEditor();
    });

    this.el.appendChild(this.controls);

    // ── Zoom controls row ─────────────────────────────────────────────
    this.zoomRow = document.createElement('div');
    this.zoomRow.className = 'panel-controls zoom-row';

    const mkLabel = (text) => {
      const s = document.createElement('span');
      s.className = 'zoom-label';
      s.textContent = text;
      return s;
    };

    this.zoomRow.appendChild(mkLabel('x:'));
    this.zoomXMin = document.createElement('input');
    this.zoomXMin.type = 'text'; this.zoomXMin.placeholder = 'min'; this.zoomXMin.size = 7;
    this.zoomRow.appendChild(this.zoomXMin);
    this.zoomXMax = document.createElement('input');
    this.zoomXMax.type = 'text'; this.zoomXMax.placeholder = 'max'; this.zoomXMax.size = 7;
    this.zoomRow.appendChild(this.zoomXMax);

    this.zoomRow.appendChild(mkLabel('z:'));
    this.zoomZMin = document.createElement('input');
    this.zoomZMin.type = 'text'; this.zoomZMin.placeholder = 'min'; this.zoomZMin.size = 7;
    this.zoomRow.appendChild(this.zoomZMin);
    this.zoomZMax = document.createElement('input');
    this.zoomZMax.type = 'text'; this.zoomZMax.placeholder = 'max'; this.zoomZMax.size = 7;
    this.zoomRow.appendChild(this.zoomZMax);

    this.zoomResetBtn = document.createElement('button');
    this.zoomResetBtn.textContent = 'Reset';
    this.zoomResetBtn.title = 'Reset to full extent';
    this.zoomResetBtn.style.display = 'none';
    this.zoomRow.appendChild(this.zoomResetBtn);

    this.zoomApplyAllBtn = document.createElement('button');
    this.zoomApplyAllBtn.textContent = '⇉ All';
    this.zoomApplyAllBtn.title = 'Apply zoom to all panels';
    this.zoomApplyAllBtn.className = 'zoom-apply-all';
    this.zoomApplyAllBtn.style.display = 'none';
    this.zoomRow.appendChild(this.zoomApplyAllBtn);

    this.el.appendChild(this.zoomRow);

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

    // ── Zoom input events ─────────────────────────────────────────────
    const emitZoom = () => {
      const d = this.panelState.fieldData;
      if (!d || !d.xCoords || !d.zCoords) return;
      const divisor = (this.model.spatialUnit || 'km') === 'km' ? 1e3 : 1;
      const xMin = this.zoomXMin.value !== '' ? parseFloat(this.zoomXMin.value) * divisor : null;
      const xMax = this.zoomXMax.value !== '' ? parseFloat(this.zoomXMax.value) * divisor : null;
      const zMin = this.zoomZMin.value !== '' ? parseFloat(this.zoomZMin.value) * divisor : null;
      const zMax = this.zoomZMax.value !== '' ? parseFloat(this.zoomZMax.value) * divisor : null;
      const anySet = xMin != null || xMax != null || zMin != null || zMax != null;
      const bounds = anySet ? { xMin, xMax, zMin, zMax } : null;
      model.setPanelViewBounds(pid, bounds);
      this._syncZoomBtns();
    };
    this.zoomXMin.addEventListener('change', emitZoom);
    this.zoomXMax.addEventListener('change', emitZoom);
    this.zoomZMin.addEventListener('change', emitZoom);
    this.zoomZMax.addEventListener('change', emitZoom);

    this.zoomResetBtn.addEventListener('click', () => {
      this.zoomXMin.value = '';
      this.zoomXMax.value = '';
      this.zoomZMin.value = '';
      this.zoomZMax.value = '';
      model.setPanelViewBounds(pid, null);
      this.zoomResetBtn.style.display = 'none';
      // Show apply-all so user can reset all panels at once
      this.zoomApplyAllBtn.style.display = '';
    });

    this.zoomApplyAllBtn.addEventListener('click', () => {
      model.applyViewBoundsToAll(pid);
      this.zoomApplyAllBtn.style.display = 'none';
    });

    // Click panel → set active
    this.el.addEventListener('click', () => {
      model.activePanelId = pid;
    });

    // ── Model events ──────────────────────────────────────────────────
    this._onFieldsLoaded = () => this._updateFieldList();
    this._onPanelField = (e) => {
      if (e.detail.panelId === pid) {
        this._syncAfterFieldChange();
      }
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
    this._onActiveChange = () => {
      this.el.classList.toggle('active', model.activePanelId === pid);
    };
    this._onParamsLoaded = () => {};
    this._onTitleChanged = () => {};
    this._onTimeUnit = () => {};
    this._onViewBounds = (e) => {
      if (e.detail.panelId === pid) {
        this._syncZoomInputs();
        // External sync: hide both buttons — this panel was just synced
        this.zoomApplyAllBtn.style.display = 'none';
        this.zoomResetBtn.style.display = this.panelState.viewBounds ? '' : 'none';
      }
    };
    this._onSpatialUnit = () => this._syncZoomInputs();

    model.addEventListener('fields-loaded',              this._onFieldsLoaded);
    model.addEventListener('panel:field-changed',         this._onPanelField);
    model.addEventListener('panel:colourmap-changed',     this._onPanelCmap);
    model.addEventListener('panel:range-changed',         this._onPanelRange);
    model.addEventListener('panel:range-locked-changed',  this._onPanelLock);
    model.addEventListener('active-panel-changed',        this._onActiveChange);
    model.addEventListener('params-loaded',               this._onParamsLoaded);
    model.addEventListener('panel:title-changed',         this._onTitleChanged);
    model.addEventListener('time-unit-changed',           this._onTimeUnit);
    model.addEventListener('panel:view-bounds-changed',   this._onViewBounds);
    model.addEventListener('spatial-unit-changed',        this._onSpatialUnit);

    // Initial population
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

    // Show/hide phase config button and continuous-only controls
    const isDiscrete = this.panelState.fieldData && this.panelState.fieldData.discrete;
    this.phaseBtn.style.display = isDiscrete ? '' : 'none';
    if (!isDiscrete) this.phasePopover.style.display = 'none';
    const cont = isDiscrete ? 'none' : '';
    this.cmapSelect.style.display = cont;
    this.rangeMin.style.display   = cont;
    this.rangeMax.style.display   = cont;
    this.autoBtn.style.display    = cont;
    this.lockBtn.style.display    = cont;

    // Reset zoom on field change
    this.zoomXMin.value = '';
    this.zoomXMax.value = '';
    this.zoomZMin.value = '';
    this.zoomZMax.value = '';
    this.zoomApplyAllBtn.style.display = 'none';
    this.zoomResetBtn.style.display = 'none';
  }

  _syncRange() {
    const { min, max } = this.panelState.colourRange;
    this.rangeMin.value = _fmtVal(min);
    this.rangeMax.value = _fmtVal(max);
  }

  _syncApplyAllBtn() {
    const vb = this.panelState.viewBounds;
    this.zoomApplyAllBtn.style.display = vb ? '' : 'none';
  }

  _syncZoomBtns() {
    const vb = this.panelState.viewBounds;
    this.zoomApplyAllBtn.style.display = vb ? '' : 'none';
    this.zoomResetBtn.style.display    = vb ? '' : 'none';
  }

  _syncZoomInputs() {
    const vb = this.panelState.viewBounds;
    const divisor = (this.model.spatialUnit || 'km') === 'km' ? 1e3 : 1;
    if (vb) {
      this.zoomXMin.value = vb.xMin != null ? _fmtVal(vb.xMin / divisor) : '';
      this.zoomXMax.value = vb.xMax != null ? _fmtVal(vb.xMax / divisor) : '';
      this.zoomZMin.value = vb.zMin != null ? _fmtVal(vb.zMin / divisor) : '';
      this.zoomZMax.value = vb.zMax != null ? _fmtVal(vb.zMax / divisor) : '';
    } else {
      this.zoomXMin.value = '';
      this.zoomXMax.value = '';
      this.zoomZMin.value = '';
      this.zoomZMax.value = '';
    }
    // Do NOT show/hide apply-all here — external sync should not toggle the button
    // The button is only shown/hidden by direct user actions (emitZoom, reset, apply-all)
  }

  _buildPhaseEditor() {
    const pop = this.phasePopover;
    pop.innerHTML = '';

    const data = this.panelState.fieldData;
    if (!data) return;
    const { values, nx, nz } = data;

    // Scan for visible phases
    const seen = new Set();
    for (let col = 0; col < nx; col++)
      for (let iz = 0; iz < nz; iz++) {
        const v = values[col][iz];
        if (v !== null && v !== undefined && !Number.isNaN(v)) {
          const p = Math.round(v);
          if (p >= 0) seen.add(p);
        }
      }
    const phases = Array.from(seen).sort((a, b) => a - b);
    const palette = this.colourMaps.phases || [];

    const heading = document.createElement('div');
    heading.className = 'phase-popover-heading';
    heading.textContent = 'Phase Legend';
    pop.appendChild(heading);

    for (const p of phases) {
      const cfg = this.model.getPhaseConfig(p, palette);
      const row = document.createElement('div');
      row.className = 'phase-row';

      const colorInput = document.createElement('input');
      colorInput.type = 'color';
      colorInput.value = _rgbToHex(cfg.color);
      colorInput.title = `Phase ${p} colour`;
      colorInput.addEventListener('input', () => {
        this.model.setPhaseConfig(p, { color: _hexToRgb(colorInput.value) });
      });
      row.appendChild(colorInput);

      const nameInput = document.createElement('input');
      nameInput.type = 'text';
      nameInput.value = cfg.name;
      nameInput.title = `Phase ${p} name`;
      nameInput.className = 'phase-name-input';
      nameInput.addEventListener('change', () => {
        this.model.setPhaseConfig(p, { name: nameInput.value });
      });
      row.appendChild(nameInput);

      pop.appendChild(row);
    }
  }

  destroy() {
    this.fieldCanvas.destroy();
    this.titleInput.destroy();
    this.model.removeEventListener('fields-loaded',              this._onFieldsLoaded);
    this.model.removeEventListener('panel:field-changed',         this._onPanelField);
    this.model.removeEventListener('panel:colourmap-changed',     this._onPanelCmap);
    this.model.removeEventListener('panel:range-changed',         this._onPanelRange);
    this.model.removeEventListener('panel:range-locked-changed',  this._onPanelLock);
    this.model.removeEventListener('active-panel-changed',        this._onActiveChange);
    this.model.removeEventListener('params-loaded',               this._onParamsLoaded);
    this.model.removeEventListener('panel:title-changed',         this._onTitleChanged);
    this.model.removeEventListener('time-unit-changed',           this._onTimeUnit);
    this.model.removeEventListener('panel:view-bounds-changed',   this._onViewBounds);
    this.model.removeEventListener('spatial-unit-changed',        this._onSpatialUnit);
    this.el.remove();
  }
}

function _fmtVal(v) {
  const abs = Math.abs(v);
  if (abs === 0) return '0';
  if (abs >= 1e6 || abs < 0.01) return v.toExponential(3);
  return v.toPrecision(5);
}

function _rgbToHex([r, g, b]) {
  return '#' + [r, g, b].map(c => c.toString(16).padStart(2, '0')).join('');
}

function _hexToRgb(hex) {
  const m = hex.match(/^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i);
  return m ? [parseInt(m[1], 16), parseInt(m[2], 16), parseInt(m[3], 16)] : [128, 128, 128];
}
