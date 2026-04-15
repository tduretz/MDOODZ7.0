// ── ControlPanel View ─────────────────────────────────────────────────
// Global controls only: file selector, time-step slider.
// Per-field controls (colour map, range, lock) now live in each PanelView.

import { formatTime } from '../time-display.mjs';

export class ControlPanel {
  constructor(model, containerEl) {
    this.model = model;
    this.el = containerEl;
    this._build();

    model.addEventListener('datasets-loaded',    () => this._updateDatasetList());
    model.addEventListener('dataset-changed',    () => this._syncDatasetUI());
    model.addEventListener('file-list-loaded',   () => this._updateFileList());
    model.addEventListener('file-selected',      () => this._syncFileUI());
    model.addEventListener('time-unit-changed',  () => this._updateFileList());
  }

  _build() {
    this.el.innerHTML = `
      <div class="ctrl-group">
        <label>Dataset</label>
        <select id="dataset-select"></select>
        <button id="refresh-datasets" title="Rescan datasets and files">⟳</button>
      </div>
      <div class="ctrl-group">
        <label>File</label>
        <select id="file-select"></select>
        <input id="time-slider" type="range" min="0" max="0" value="0" step="1">
        <span id="slider-label" class="slider-label"></span>
      </div>
    `;

    this.datasetSelect = this.el.querySelector('#dataset-select');
    this.refreshBtn    = this.el.querySelector('#refresh-datasets');
    this.fileSelect    = this.el.querySelector('#file-select');
    this.timeSlider    = this.el.querySelector('#time-slider');
    this.sliderLabel   = this.el.querySelector('#slider-label');

    // Dataset change
    this.datasetSelect.addEventListener('change', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:dataset-change', { detail: this.datasetSelect.value }));
    });

    // Refresh datasets
    this.refreshBtn.addEventListener('click', () => {
      this.el.dispatchEvent(new CustomEvent('ctrl:refresh-datasets'));
    });

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

  _updateDatasetList() {
    const datasets = this.model.datasets;
    this.datasetSelect.innerHTML = datasets.map(d =>
      `<option value="${d.name}">${d.name} (${d.fileCount} files)</option>`
    ).join('');
    this._syncDatasetUI();
  }

  _syncDatasetUI() {
    if (this.model.activeDataset) {
      this.datasetSelect.value = this.model.activeDataset;
    }
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
