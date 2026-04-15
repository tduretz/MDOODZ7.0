// ── HeaderBanner View ─────────────────────────────────────────────────
// Shows: left=title, centre=layout selector, right=time+unit selector.

import { formatTime } from '../time-display.mjs';

const LAYOUT_PRESETS = ['1x1', '1x2', '2x1', '2x2'];

export class HeaderBanner {
  constructor(model, headerEl, controlsEl) {
    this.model = model;
    this.controlsEl = controlsEl;
    this.timeEl  = headerEl.querySelector('#time-display');
    this.unitSel = headerEl.querySelector('#time-unit-select');

    // Replace #header-field with layout selector
    const oldField = headerEl.querySelector('#header-field');
    this.layoutBar = document.createElement('span');
    this.layoutBar.className = 'header-layout';
    for (const preset of LAYOUT_PRESETS) {
      const btn = document.createElement('button');
      btn.className = 'layout-btn';
      btn.dataset.layout = preset;
      btn.textContent = preset;
      if (preset === model.layout) btn.classList.add('active');
      btn.addEventListener('click', () => {
        controlsEl.dispatchEvent(new CustomEvent('ctrl:layout-change', { detail: preset }));
      });
      this.layoutBar.appendChild(btn);
    }
    oldField.replaceWith(this.layoutBar);

    model.addEventListener('params-loaded',     () => this._renderTime());
    model.addEventListener('time-unit-changed',  () => this._renderTime());
    model.addEventListener('layout-changed',     () => this._syncLayoutBtns());

    this.unitSel.addEventListener('change', () => {
      const val = this.unitSel.value || null;
      controlsEl.dispatchEvent(new CustomEvent('ctrl:time-unit', { detail: val }));
    });

    this.spatialSel = headerEl.querySelector('#spatial-unit-select');
    this.spatialSel.value = model.spatialUnit || 'km';
    this.spatialSel.addEventListener('change', () => {
      model.spatialUnit = this.spatialSel.value;
    });

    // Theme toggle
    this.themeBtn = headerEl.querySelector('#theme-toggle');
    this.themeBtn.addEventListener('click', () => {
      model.theme = model.theme === 'dark' ? 'light' : 'dark';
    });
    model.addEventListener('theme-changed', () => {
      this.themeBtn.textContent = model.theme === 'dark' ? '🌙' : '☀️';
    });
  }

  _renderTime() {
    const p = this.model.params;
    if (p) {
      const { formatted } = formatTime(p.time, this.model.timeUnit);
      this.timeEl.textContent = formatted;
    }
    this.unitSel.value = this.model.timeUnit || '';
  }

  _syncLayoutBtns() {
    for (const btn of this.layoutBar.querySelectorAll('.layout-btn')) {
      btn.classList.toggle('active', btn.dataset.layout === this.model.layout);
    }
  }
}
