## Context

WebVis is a zero-dependency browser HDF5 viewer (ES6+ MVC, Canvas ImageData, h5wasm). The `webvis-multi-panel-titles` change added multi-panel layouts (1×1/1×2/2×1/2×2), per-panel titles with `${...}` interpolation, and panel-scoped controls. Testing revealed six bugs and performance issues that need fixing before the feature is production-ready.

Current architecture: `app.mjs` creates PanelView instances via a `syncPanels()` function. Each PanelView owns a FieldCanvas, ColourBar, and TitleInput. A ResizeObserver on `#panels-grid` handles responsive sizing. The controller wires global events (file change, slider) and panel-scoped events (field, cmap, range) to model updates and API fetches.

## Goals / Non-Goals

**Goals:**
- Fix canvas not rendering on initial load (race between ResizeObserver and data arrival)
- Fix 2×1/1×2 distortion on window resize (wrong available-width calculation)
- Preserve panel states across layout round-trips (2×2 → 1×1 → 2×2)
- Render interpolated title text inside the canvas (for screenshots and exports)
- Debounce slider scrubbing and abort superseded fetch requests
- Eliminate memory leaks from accumulated fieldData and orphaned event listeners

**Non-Goals:**
- New features beyond what's listed (no new panel controls, no save/load state)
- Server-side changes
- Canvas download/export API (future change)
- Changing the panel maximum from 4

## Decisions

### D1: Fix initial render — explicit render after data arrival

**Problem**: `ResizeObserver` fires on layout change but field data hasn't arrived yet. `panel:field-changed` triggers `FieldCanvas.render()` but the canvas hasn't been sized.

**Decision**: After `selectPanelField()` resolves, call `resizePanels()` via `requestAnimationFrame`. Also add a `render()` call inside the `panel:field-changed` handler in `app.mjs` as a safety net.

**Alternative**: Could make FieldCanvas self-sizing (read `parentElement.getBoundingClientRect()`), but that couples the view to DOM layout concerns.

### D2: Fix resize math — use `getBoundingClientRect()`

**Problem**: Current code subtracts a hardcoded `88px` for the colour bar. This fails for narrow panels in 1×2 layouts.

**Decision**: Replace `row.clientWidth - 88` with measuring the actual available space: `canvasRow.getBoundingClientRect()` minus the colour bar element's actual width. This adapts correctly to all grid cell sizes.

### D3: Panel stash for layout round-trips

**Problem**: `model.setLayout()` calls `removePanel()` which destroys state. When layout grows back, new panels start blank.

**Decision**: Instead of destroying panels, stash them. Add `_stashedPanels` array to Model. When layout shrinks, move excess panels to stash. When layout grows, pop from stash before creating fresh defaults. `PanelView.destroy()` is deferred — views detach from DOM but panel data is preserved.

**Alternative**: Keep all 4 panels alive in model at all times, just hide DOM. Rejected because it wastes memory for data arrays of hidden panels and complicates event flow.

### D4: Render title in canvas

**Problem**: The `<div class="panel-title-rendered">` is a separate DOM element. Screenshots and canvas-to-blob exports don't include it.

**Decision**: Move title rendering into `FieldCanvas.render()`. After drawing the ImageData, use `ctx.fillText()` to draw the interpolated title at the top of the canvas with a semi-transparent background bar. Remove the DOM `<span>` from TitleInput. Keep the editable `<input>` and dropdown as DOM elements.

**Trade-off**: Canvas text is less crisp than DOM text at non-integer DPI. Mitigated by using `devicePixelRatio` scaling.

### D5: Debounce slider + AbortController

**Problem**: Each slider position change calls `selectFile()` → N parallel fetches (one per panel, plus params). Scrubbing through 29 files creates ~4×29 = 116 overlapping requests.

**Decision**: 
1. Add 150ms debounce on slider `input` event (fires `selectFile()` only after user pauses).
2. Wrap all fetch calls in controller with a shared `AbortController`. Each new `selectFile()` call aborts the previous controller. Catch `AbortError` silently.
3. The slider label updates immediately (no debounce) for visual responsiveness.

### D6: Memory leak fixes

**Problem 1**: `panelState.fieldData` holds multi-MB typed arrays. When new data arrives, old data is just overwritten. GC should collect it, but closures in event handlers may retain references.

**Problem 2**: `PanelView.destroy()` removes the DOM but FieldCanvas and ColourBar still hold event listener references on the model.

**Decision**:
1. In `model.setPanelField()`, explicitly set `p.fieldData = null` before assigning new data.
2. Add `destroy()` methods to FieldCanvas and ColourBar that call `model.removeEventListener()` for all handlers.
3. In `PanelView.destroy()`, call `this.fieldCanvas.destroy()`, `this.colourBar.destroy()`, `this.titleInput.destroy()`.
4. In `app.mjs`, when `syncPanels()` removes a view, ensure `view.destroy()` is called.

## Risks / Trade-offs

- **[Canvas title readability]** → Use white text with dark semi-transparent background bar; test at different zoom levels
- **[Debounce too aggressive]** → 150ms is a good balance; if users complain, make it configurable. Slider label updates instantly for feedback.
- **[Stash memory]** → Stashed panels retain fieldData arrays. Worst case: 4 panels × ~5MB = 20MB. Acceptable for modern browsers. Could add a `clearStash()` on file change if needed.
- **[AbortController browser compat]** → Supported in all modern browsers. h5wasm-based server uses Node 18+. No issue.
