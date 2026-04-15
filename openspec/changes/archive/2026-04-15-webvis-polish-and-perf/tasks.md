# Tasks — webvis-polish-and-perf

## Group 1 — Panel stash for layout round-trips (D3)

- [x] 1.1  Add `_stashedPanels` array to Model constructor
- [x] 1.2  Rewrite `setLayout()` to stash excess panels instead of destroying, restore from stash before creating new
- [x] 1.3  Null `fieldData` on setPanelField before assigning new value

## Group 2 — Fix responsive sizing (D2)

- [x] 2.1  Replace hardcoded `row.clientWidth - 88` with `getBoundingClientRect()` minus colourbar actual width
- [x] 2.2  Fix tall-panel (2×1) sizing — use available height from grid cell

## Group 3 — Fix initial render timing (D1)

- [x] 3.1  In app.mjs `panel:field-changed` handler, schedule render via rAF and recompute sizes
- [x] 3.2  In `resizePanels()`, guard against zero-size canvas and render after size is set

## Group 4 — Canvas title overlay (D4)

- [x] 4.1  Add title-drawing pass at end of `FieldCanvas.render()` — use 2D context fillText with semi-transparent background
- [x] 4.2  Pass panelState + model refs so FieldCanvas can access title template and build context
- [x] 4.3  Remove `.panel-title-rendered` DOM element from TitleInput (keep input + dropdown only)

## Group 5 — Debounce slider + AbortController (D5)

- [x] 5.1  Add `_abortCtrl` to Controller, abort and replace on each `selectFile()` / `selectPanelField()`
- [x] 5.2  Pass `{ signal }` to all fetch calls in Controller
- [x] 5.3  Add 150ms debounce to slider `input` event in ControlPanel
- [x] 5.4  Update slider label immediately on input (before debounce fires)

## Group 6 — Memory cleanup (D6)

- [x] 6.1  Add `destroy()` to FieldCanvas — remove all model event listeners
- [x] 6.2  Add `destroy()` to ColourBar — remove all model event listeners
- [x] 6.3  Call `fieldCanvas.destroy()` and `colourBar.destroy()` from PanelView.destroy()
- [x] 6.4  In syncPanels, null fieldData for non-stashed destroyed panels
