## 1. Server — params endpoint

- [x] 1.1 Add `GET /api/params/:file` route to `server.mjs` returning `{ step, time, Nx, Nz, dt }` from the in-memory metadata cache (fallback: read HDF5 `/Model/Params`)
- [x] 1.2 Return HTTP 404 with JSON error for unknown filenames
- [x] 1.3 Smoke-test: `curl /api/params/Output00000.gzip.h5` returns correct JSON

## 2. Model — PanelState and panels array

- [x] 2.1 Define `PanelState` shape `{ id, fieldName, colourMap, colourRange, rangeLocked, fieldData, title }` and a `createDefaultPanel(id)` factory in `model.mjs`
- [x] 2.2 Add `panels` array and `activePanelId` properties to Model; initialise with one default panel on construction
- [x] 2.3 Add panel CRUD helpers: `addPanel()`, `removePanel(id)`, `getPanel(id)`, `setLayout(preset)`
- [x] 2.4 Add panel-scoped setters: `setPanelField(panelId, name, data)`, `setPanelColourMap(panelId, cmap)`, `setPanelRange(panelId, range)`, `setPanelRangeLocked(panelId, locked)`, `setPanelTitle(panelId, template)`
- [x] 2.5 Dispatch panel-scoped events with `{ panelId }` in detail (e.g. `panel:field-changed`, `panel:range-changed`, `panel:title-changed`)
- [x] 2.6 Add `layout` property (default `'1x1'`) and dispatch `layout-changed` event from `setLayout()`
- [x] 2.7 Migrate existing single-field state (`_colourMap`, `_colourRange`, `_rangeLocked`, `_fieldColourMaps`, etc.) into the first PanelState; keep global state (files, currentFile, fields, fieldDefs, timeUnit) on Model

## 3. Controller — panel-aware routing

- [x] 3.1 Add `selectPanelField(panelId, fieldName)` that fetches `/api/field-data/:file/:field` and calls `model.setPanelField()`
- [x] 3.2 Add `fetchParams()` that calls `/api/params/:file` and stores scalar metadata on Model (step, time, Nx, Nz, dt)
- [x] 3.3 Refactor `selectFile()` to call `fetchParams()` then re-fetch each panel's field in parallel via `selectPanelField()`
- [x] 3.4 Add `changeLayout(preset)` handler that calls `model.setLayout()`, adds/removes panels, and triggers field fetches for new panels
- [x] 3.5 Route per-panel control events (`ctrl:panel:field`, `ctrl:panel:cmap`, `ctrl:panel:lock`) to the correct panel setter

## 4. CSS Grid layout

- [x] 4.1 Add `#panels-grid` container to `index.html` replacing the single-canvas area
- [x] 4.2 Define CSS classes `.grid-1x1`, `.grid-1x2`, `.grid-2x1`, `.grid-2x2` in `style.css` with appropriate `grid-template-columns` / `grid-template-rows`
- [x] 4.3 Default class is `.grid-1x1`; switching layout swaps the class on `#panels-grid`

## 5. PanelView — per-panel DOM subtree

- [x] 5.1 Create `public/views/PanelView.mjs` that builds `<div class="panel">` containing: title area, `<canvas>`, colour bar, and mini control row (field select, cmap select, range inputs, lock button)
- [x] 5.2 PanelView instantiates its own `FieldCanvas` and `ColourBar` bound to its PanelState
- [x] 5.3 PanelView listens only to events matching its `panelId` and re-renders accordingly
- [x] 5.4 Mini control row dispatches `ctrl:panel:<event>` with `{ panelId }` in detail
- [x] 5.5 Clicking a panel sets `model.activePanelId`
- [x] 5.6 Style `.panel` with border/padding; highlight active panel with a distinct border colour

## 6. FieldCanvas & ColourBar — accept target element

- [x] 6.1 Refactor `FieldCanvas.mjs` to accept a canvas element as constructor parameter instead of querying the DOM for a single global canvas
- [x] 6.2 Refactor `ColourBar.mjs` to accept a container element as constructor parameter and render into it
- [x] 6.3 Ensure both read colour map, colour range, and range-lock from the PanelState passed to them (not from global Model properties)

## 7. Title interpolation engine

- [x] 7.1 Create `public/title-render.mjs` with `renderTitle(template, context) → string` that replaces `${...}` tokens; unknown tokens left as-is
- [x] 7.2 Build context object from scalar params (`step`, `time`, `Nx`, `Nz`, `resolution`) + field stats (`min`, `max`, `mean`, `median`, `p2`, `p98`, `field`, `unit`)
- [x] 7.3 Format stat values with scientific notation + unit from field-labels registry; format `time` with adaptive `formatTime()`; format `resolution` as `Nx × Nz`
- [x] 7.4 Implement lazy `mean` and `median` computation from the flat values array; cache results in PanelState; invalidate on field data change

## 8. TitleInput view — editable input + dropdown

- [x] 8.1 Create `public/views/TitleInput.mjs` with `<input type="text">` + rendered `<span>` below + ⊕ insert button
- [x] 8.2 Pre-fill input with default template `"${field} at t = ${time} — min = ${min}  max = ${max}"`
- [x] 8.3 On `$` keypress or ⊕ click, show floating dropdown with two groups: **Metadata** (Step number, Model time, Grid columns, Grid rows, Grid resolution) and **Field stats** (Min, Max, Mean, Median, Robust min (p2), Robust max (p98), Field name, Unit)
- [x] 8.4 Selecting a dropdown item inserts `${tokenName}` at cursor position and closes the dropdown
- [x] 8.5 Escape or outside click dismisses the dropdown
- [x] 8.6 On input change, call `renderTitle()` and update the rendered `<span>`; dispatch `ctrl:panel:title` event

## 9. Layout selector in header

- [x] 9.1 Add layout selector buttons (1×1, 1×2, 2×1, 2×2) to `HeaderBanner.mjs` centre section, replacing the single-field label
- [x] 9.2 Highlight the active layout button
- [x] 9.3 On click, dispatch layout-change event to controller

## 10. Web-UI integration — global vs per-panel controls

- [x] 10.1 Move file selector, time slider, and time-unit dropdown into a global controls bar above `#panels-grid`
- [x] 10.2 Remove per-field colour map / range / lock controls from the global `ControlPanel.mjs` (they now live in each PanelView's mini control row)
- [x] 10.3 Update `StatusBar.mjs` to show panel count instead of single-field name
- [x] 10.4 On file change / slider drag, trigger all panels to re-fetch their fields in parallel

## 11. Responsive sizing

- [x] 11.1 Add `ResizeObserver` on `#panels-grid` (or window resize listener) that recalculates each panel's canvas size to fill its grid cell while preserving the Nx:Nz aspect ratio
- [x] 11.2 On resize, re-render each panel's FieldCanvas and ColourBar at the new dimensions

## 12. Backward compatibility & smoke test

- [x] 12.1 Verify default 1×1 layout renders identically to current single-panel app (no regression)
- [x] 12.2 Test 1×2 layout: two panels side-by-side with independent fields, colour maps, ranges, and titles
- [x] 12.3 Test 2×2 layout: four panels with parallel field fetches on file switch
- [x] 12.4 Test title interpolation: edit template, insert token via dropdown, verify rendered output updates
- [x] 12.5 Test layout switching: grow 1→4, shrink 4→1, confirm panel 1 state preserved
