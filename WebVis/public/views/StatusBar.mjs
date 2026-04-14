// ── StatusBar View ────────────────────────────────────────────────────
// Shows model time (Myr), grid dimensions, current field name.

const SEC_TO_MYR = 1 / 3.1558e13;

export class StatusBar {
  constructor(model, el) {
    this.model = model;
    this.el = el;

    model.addEventListener('field-loaded',  () => this.render());
    model.addEventListener('file-selected', () => this.render());
    model.addEventListener('params-loaded', () => this.render());
  }

  render() {
    const parts = [];
    const p = this.model.params;
    if (p) {
      const myr = (p.time * SEC_TO_MYR).toFixed(3);
      parts.push(`t = ${myr} Myr`);
      parts.push(`${p.nx} × ${p.nz}`);
    }
    if (this.model.currentField) {
      parts.push(this.model.currentField);
    }
    this.el.textContent = parts.join(' | ');
  }
}
