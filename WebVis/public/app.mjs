// ── App Entry Point ───────────────────────────────────────────────────
import { Model }          from './model.mjs';
import { Controller }     from './controller.mjs';
import { FieldCanvas }    from './views/FieldCanvas.mjs';
import { ColourBar }      from './views/ColourBar.mjs';
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

  // Views
  new ControlPanel(model, controlsEl);
  new FieldCanvas(model, document.getElementById('field-canvas'), colourMaps);
  new ColourBar(model, document.getElementById('colourbar-canvas'), colourMaps);
  new StatusBar(model, document.getElementById('status-bar'));
  new LoadingOverlay(model, document.getElementById('loading-overlay'));
  new HeaderBanner(model, document.getElementById('header'), controlsEl);

  // Controller
  const ctrl = new Controller(model, controlsEl);
  await ctrl.init();
}

main().catch(err => {
  console.error('App init failed:', err);
  document.getElementById('status-bar').textContent = 'Error: ' + err.message;
});
