// ── StatusBar View ────────────────────────────────────────────────────
// Shows model time (adaptive), grid dimensions, current field label.

import { formatTime } from '../time-display.mjs';

export class StatusBar {
  constructor(model, el) {
    this.model = model;
    this.el = el;

    model.addEventListener('field-loaded',      () => this.render());
    model.addEventListener('file-selected',     () => this.render());
    model.addEventListener('params-loaded',     () => this.render());
    model.addEventListener('time-unit-changed',  () => this.render());
  }

  render() {
    const parts = [];
    const p = this.model.params;
    if (p) {
      const { formatted } = formatTime(p.time, this.model.timeUnit);
      parts.push(`t = ${formatted}`);
      parts.push(`${p.nx} × ${p.nz}`);
    }
    const field = this.model.currentField;
    if (field) {
      const def = this.model.fieldDefs.get(field);
      parts.push(def ? def.label : field);
    }
    this.el.textContent = parts.join(' | ');
  }
}
