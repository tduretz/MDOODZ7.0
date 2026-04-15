// ── App Entry Point ───────────────────────────────────────────────────
import { Model }          from './model.mjs';
import { Controller }     from './controller.mjs';
import { PanelView }      from './views/PanelView.mjs';
import { MARGIN_ESTIMATE } from './views/FieldCanvas.mjs';
import { ControlPanel }   from './views/ControlPanel.mjs';
import { StatusBar }      from './views/StatusBar.mjs';
import { LoadingOverlay } from './views/LoadingOverlay.mjs';
import { HeaderBanner }   from './views/HeaderBanner.mjs';

async function main() {
  // Load colour maps
  const cmapRes = await fetch('/colormaps/maps.json');
  const colourMaps = await cmapRes.json();

  // Model
  const model = new Model();

  const controlsEl = document.getElementById('controls');
  const gridEl     = document.getElementById('panels-grid');

  // Views
  new ControlPanel(model, controlsEl);
  new StatusBar(model, document.getElementById('status-bar'));
  new LoadingOverlay(model, document.getElementById('loading-overlay'));
  new HeaderBanner(model, document.getElementById('header'), controlsEl);

  // Panel management
  const panelViews = new Map();  // panelId → PanelView

  function syncPanels() {
    // Update grid class
    gridEl.className = `grid-${model.layout}`;

    // Remove views for panels that no longer exist
    for (const [id, view] of panelViews) {
      if (!model.getPanel(id)) {
        view.destroy();
        panelViews.delete(id);
      }
    }

    // Create views for new panels (including restored from stash)
    for (const panel of model.panels) {
      if (!panelViews.has(panel.id)) {
        const view = new PanelView(model, panel, colourMaps, gridEl);
        panelViews.set(panel.id, view);
      }
    }

    // Trigger resize+render after layout settles (handles restored panels with existing data)
    requestAnimationFrame(() => requestAnimationFrame(resizePanels));
  }

  // Create initial panels
  syncPanels();

  // Sync on layout changes (setLayout dispatches layout-changed after all add/remove)
  model.addEventListener('layout-changed', () => syncPanels());

  // ── Responsive sizing ──────────────────────────────────────────────
  function resizePanels() {
    for (const [id, view] of panelViews) {
      const panel = model.getPanel(id);
      if (!panel || !panel.fieldData) continue;
      const { nx, nz } = panel.fieldData;
      const rowRect = view.canvasRow.getBoundingClientRect();
      const availW = rowRect.width;
      const availH = rowRect.height;
      if (availW <= 0 || availH <= 0) continue;

      // Budget for margins — FieldCanvas will compute exact values dynamically
      const MARGIN_LR = MARGIN_ESTIMATE.left + MARGIN_ESTIMATE.right;
      const MARGIN_TB = MARGIN_ESTIMATE.top  + MARGIN_ESTIMATE.bottom;
      const dataAspect = nx / nz;

      // Find best fit: data portion fills available space
      let dataW, dataH;
      const testW = availW - MARGIN_LR;
      const testH = availH - MARGIN_TB;
      if (testW <= 0 || testH <= 0) continue;

      if (testW / testH > dataAspect) {
        dataH = testH;
        dataW = Math.round(dataH * dataAspect);
      } else {
        dataW = testW;
        dataH = Math.round(dataW / dataAspect);
      }

      const totalW = dataW + MARGIN_LR;
      const totalH = dataH + MARGIN_TB;
      view.fieldCanvasEl.style.width  = totalW + 'px';
      view.fieldCanvasEl.style.height = totalH + 'px';
      view.fieldCanvas.render();
    }
  }

  const ro = new ResizeObserver(() => resizePanels());
  ro.observe(gridEl);

  // Also resize after field loads — ensures first render even if resize hasn't fired
  model.addEventListener('panel:field-changed', () => {
    // Double-rAF: wait for layout to settle, then size + render
    requestAnimationFrame(() => requestAnimationFrame(resizePanels));
  });

  // Controller
  const ctrl = new Controller(model, controlsEl);
  await ctrl.init();
}

main().catch(err => {
  console.error('App init failed:', err);
  document.getElementById('status-bar').textContent = 'Error: ' + err.message;
});
