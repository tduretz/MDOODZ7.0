// ── ControlPanel View ─────────────────────────────────────────────────
// Global controls only: file selector, time-step slider.
// Per-field controls (colour map, range, lock) now live in each PanelView.

import { formatTime } from '../time-display.mjs';

export class ControlPanel {
  constructor(model, containerEl) {
    this.model = model;
    this.el = containerEl;
    this._build();

    model.addEventListener('file-list-loaded',  () => this._updateFileList());
    model.addEventListener('file-selected',     () => this._syncFileUI());
    model.addEventListener('time-unit-changed',  () => this._updateFileList());
  }

  _build() {
    this.el.innerHTML = `
      <div class="ctrl-group">
        <label>File</label>
        <select id="file-select"></select>
        <input id="time-slider" type="range" min="0" max="0" value="0" step="1">
        <span id="slider-label" class="slider-label"></span>
      </div>
    `;

    this.fileSelect  = this.el.querySelector('#file-select');
    this.timeSlider  = this.el.querySelector('#time-slider');
    this.sliderLabel = this.el.querySelector('#slider-label');

    // DOM events → dispatched as custom events on this.el
    this.fileSelect.addEventListener('change', () => {
      this.timeSlider.value = this.fileSelect.selectedIndex;
      this._updateSliderLabel();
      this.el.dispatchEvent(new CustomEvent('ctrl:file-change', { detail: this.fileSelect.value }));
    });

    this._sliderTimer = null;
    this.timeSlider.addEventListener('input', () => {
      this.fileSelect.selectedIndex = parseInt(this.timeSlider.value, 10);
      this._updateSliderLabel();  // immediate visual feedback
      clearTimeout(this._sliderTimer);
      this._sliderTimer = setTimeout(() => {
        this.el.dispatchEvent(new CustomEvent('ctrl:file-change', { detail: this.fileSelect.value }));
      }, 150);
    });
  }

  _updateFileList() {
    const fileList = this.model.fileList;
    const unit = this.model.timeUnit;
    this.fileSelect.innerHTML = fileList.map(f => {
      const { formatted } = formatTime(f.time, unit);
      return `<option value="${f.name}">Step ${f.step} — ${formatted}</option>`;
    }).join('');
    this.timeSlider.max = String(Math.max(0, fileList.length - 1));
    this._updateSliderLabel();
  }

  _syncFileUI() {
    const idx = this.model.files.indexOf(this.model.currentFile);
    if (idx >= 0) {
      this.fileSelect.selectedIndex = idx;
      this.timeSlider.value = String(idx);
      this._updateSliderLabel();
    }
  }

  _updateSliderLabel() {
    const idx = parseInt(this.timeSlider.value, 10);
    const fileList = this.model.fileList;
    if (fileList[idx]) {
      const { formatted } = formatTime(fileList[idx].time, this.model.timeUnit);
      this.sliderLabel.textContent = `Step ${fileList[idx].step} — ${formatted}`;
    }
  }
}
