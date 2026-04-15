// ── App Entry Point ───────────────────────────────────────────────────
import { Model }          from './model.mjs';
import { Controller }     from './controller.mjs';
import { PanelView }      from './views/PanelView.mjs';
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

  // Sync on layout changes
  model.addEventListener('layout-changed', () => syncPanels());
  model.addEventListener('panel-added',    () => syncPanels());
  model.addEventListener('panel-removed',  () => syncPanels());

  // ── Responsive sizing ──────────────────────────────────────────────
  function resizePanels() {
    for (const [id, view] of panelViews) {
      const panel = model.getPanel(id);
      if (!panel || !panel.fieldData) continue;
      const { nx, nz } = panel.fieldData;
      const rowRect = view.canvasRow.getBoundingClientRect();
      const cbarRect = view.colourBarEl.getBoundingClientRect();
      const availW = rowRect.width - cbarRect.width - 8; // 8px gap
      const availH = rowRect.height;
      if (availW <= 0 || availH <= 0) continue;
      const aspect = nx / nz;
      let w, h;
      if (availW / availH > aspect) {
        h = availH;
        w = Math.round(h * aspect);
      } else {
        w = availW;
        h = Math.round(w / aspect);
      }
      view.fieldCanvasEl.style.width  = w + 'px';
      view.fieldCanvasEl.style.height = h + 'px';
      view.colourBarEl.height = h;
      view.fieldCanvas.render();
      view.colourBar.render();
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
