// ── LoadingOverlay View ───────────────────────────────────────────────

export class LoadingOverlay {
  constructor(model, el) {
    this.el = el;
    model.addEventListener('loading-start', () => this.el.classList.remove('hidden'));
    model.addEventListener('loading-end',   () => this.el.classList.add('hidden'));
  }
}
