// ── HeaderBanner View ─────────────────────────────────────────────────
// Shows: left=title, centre=field label+unit, right=time+unit selector.

import { formatTime } from '../time-display.mjs';

export class HeaderBanner {
  constructor(model, headerEl, controlsEl) {
    this.model = model;
    this.fieldEl = headerEl.querySelector('#header-field');
    this.timeEl  = headerEl.querySelector('#time-display');
    this.unitSel = headerEl.querySelector('#time-unit-select');

    model.addEventListener('field-loaded',      () => this.render());
    model.addEventListener('params-loaded',     () => this.render());
    model.addEventListener('time-unit-changed',  () => this.render());

    this.unitSel.addEventListener('change', () => {
      const val = this.unitSel.value || null;
      controlsEl.dispatchEvent(new CustomEvent('ctrl:time-unit', { detail: val }));
    });
  }

  render() {
    const field = this.model.currentField;
    const def = field && this.model.fieldDefs.get(field);
    if (field) {
      const label = def ? def.label : field;
      const unit  = def ? def.formattedUnit : '';
      this.fieldEl.textContent = unit ? `${label} [${unit}]` : label;
    }

    const p = this.model.params;
    if (p) {
      const { formatted } = formatTime(p.time, this.model.timeUnit);
      this.timeEl.textContent = formatted;
    }

    // Sync dropdown
    this.unitSel.value = this.model.timeUnit || '';
  }
}
