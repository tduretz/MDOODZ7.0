// ── StatusBar View ────────────────────────────────────────────────────
// Shows model time (adaptive), grid dimensions, panel count, active field.

import { formatTime } from '../time-display.mjs';

export class StatusBar {
  constructor(model, el) {
    this.model = model;
    this.el = el;

    model.addEventListener('panel:field-changed',  () => this.render());
    model.addEventListener('file-selected',        () => this.render());
    model.addEventListener('params-loaded',        () => this.render());
    model.addEventListener('time-unit-changed',     () => this.render());
    model.addEventListener('layout-changed',        () => this.render());
    model.addEventListener('active-panel-changed',  () => this.render());
  }

  render() {
    const parts = [];
    const p = this.model.params;
    if (p) {
      const { formatted } = formatTime(p.time, this.model.timeUnit);
      parts.push(`t = ${formatted}`);
      parts.push(`${p.Nx ?? p.nx ?? '?'} × ${p.Nz ?? p.nz ?? '?'}`);
    }
    const n = this.model.panels.length;
    parts.push(`${n} panel${n > 1 ? 's' : ''} (${this.model.layout})`);
    // Show active panel field
    const active = this.model.getPanel(this.model.activePanelId);
    if (active && active.fieldName) {
      const def = this.model.fieldDefs.get(active.fieldName);
      parts.push(def ? def.label : active.fieldName);
    }
    this.el.textContent = parts.join(' | ');
  }
}
