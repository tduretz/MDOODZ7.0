## 1. Model — overlay state

- [x] 1.1 Add default `overlays` config object to `createDefaultPanel()` in model.mjs with all 8 layer types (phases, temperature, topo, velocity, director, sigma1, edot1, melt), each with `enabled: false` and default colours/params
- [x] 1.2 Add `setPanelOverlay(panelId, layerType, config)` method that merges config into the panel's overlay entry and emits `panel:overlay-changed` event with `{ panelId, layerType }`
- [x] 1.3 Add `setPanelOverlayAvailability(panelId, availableSet)` method that stores which layers have data for the current file, for UI enable/disable

## 2. Server — batched overlay data endpoint

- [x] 2.1 Add `extractOverlayData(filePath, layers)` function in hdf5-reader.mjs that reads requested HDF5 datasets in one file-open, returns an object keyed by layer type with flat arrays + grid dimensions + coordinate arrays; pre-averages vertex sxz/exz to centres; skips layers whose datasets don't exist
- [x] 2.2 Add `GET /api/overlay-data/:filename?layers=...` route in server.mjs with filename validation, layer key parsing, dataset existence checks, gzip-compressed JSON response
- [x] 2.3 Add topography extraction (`/Topo/z_grid` + `/Model/xg_coord`) to the overlay data handler, included only when the `/Topo` group exists

## 3. Client — overlay math utilities (overlays.mjs)

- [x] 3.1 Create `WebVis/public/overlays.mjs` with `marchingSquares(grid, xCoords, zCoords, level)` returning `[{x1,z1,x2,z2}, ...]` segments in SI metres; handle NaN as below-threshold
- [x] 3.2 Add `extractContours(grid, xCoords, zCoords, levels)` that calls marchingSquares per level and returns `Map<level, segments[]>`
- [x] 3.3 Add `eigen2x2(a, b, c)` closed-form solver for symmetric `[[a,b],[b,c]]` matrix returning `{ lambda1, lambda2, cos, sin }` (eigenvector of largest eigenvalue)
- [x] 3.4 Add `interpVxToCenter(Vx, vxNx, vxNz)` and `interpVzToCenter(Vz, vzNx, vzNz)` staggered-grid averaging helpers returning flat Float64Arrays on cell-centre grid
- [x] 3.5 Add `inlineContourLabels(segments, level, pixelMapper, intervalPx)` that returns label positions `{x, z, angle, text}` at periodic arc-length intervals along joined polyline segments

## 4. Controller — overlay data fetching and caching

- [x] 4.1 Add overlay data cache in controller.mjs — `Map<filename, overlayDataObject>` — invalidated on file/dataset change
- [x] 4.2 Wire `ctrl:panel:overlay-toggle` event listener in controller: on enable, check cache, fetch `/api/overlay-data/:filename?layers=...` for all currently-enabled layers whose data is not yet cached, store result, call `model.setPanelOverlayAvailability()`, trigger re-render
- [x] 4.3 Wire `ctrl:panel:overlay-config` event listener in controller: call `model.setPanelOverlay()` with new config (no fetch needed — just model update triggers re-render)
- [x] 4.4 On file change (`selectFile`), clear overlay cache and re-fetch overlay data for all panels that have enabled layers

## 5. FieldCanvas — overlay rendering pipeline

- [x] 5.1 Add `_drawOverlays(ctx, dataX, dataY, dW, dH, col0, col1, row0, row1)` method in FieldCanvas.mjs that clips to data area and dispatches to per-layer draw helpers based on panelState.overlays config
- [x] 5.2 Call `_drawOverlays()` in `render()` after `ctx.strokeRect` (frame) and before `_drawCbar` / `_drawPhaseLegend`
- [x] 5.3 Listen for `panel:overlay-changed` event and trigger re-render
- [x] 5.4 Add `_coordToPixel(x, z, xCoords, zCoords, col0, col1, row0, row1, dataX, dataY, dW, dH)` helper mapping physical SI coordinates to canvas CSS pixel positions

## 6. FieldCanvas — contour layer renderers

- [x] 6.1 Add `_drawContourLayer(ctx, layerData, config, pixMapper)` that calls `extractContours()` and draws segments as `ctx.beginPath(); moveTo/lineTo; stroke()` with configured colour/lineWidth
- [x] 6.2 Implement phase boundary rendering: reshape compo flat array to 2D, generate levels at integer midpoints, call `_drawContourLayer`
- [x] 6.3 Implement temperature isotherm rendering: reshape T array to 2D, apply −273.15 offset, generate levels from ΔT, call `_drawContourLayer`, then draw inline labels via `inlineContourLabels()`
- [x] 6.4 Implement melt fraction rendering: reshape X array to 2D, use configured levels array, call `_drawContourLayer`

## 7. FieldCanvas — topography line renderer

- [x] 7.1 Add `_drawTopoLine(ctx, topoData, config, pixMapper)` that draws a polyline from `z_grid` / `x_grid` arrays using moveTo/lineTo, clipped to crop region

## 8. FieldCanvas — vector and tensor layer renderers

- [x] 8.1 Add `_drawArrow(ctx, x, z, dx, dz, headSize)` canvas arrow primitive (line + triangle head)
- [x] 8.2 Add `_drawVelocityLayer(ctx, velData, config, pixMapper)` that interpolates Vx/Vz to centres, sub-samples by density, draws arrows scaled to magnitude
- [x] 8.3 Add `_drawDirectorLayer(ctx, dirData, config, pixMapper)` that transforms (nx,nz) to fabric direction (−nz/nx, 1), normalises, draws crossing ticks at sub-sampled positions, skips where ani_fac < 0.01
- [x] 8.4 Add `_drawTensorLayer(ctx, tensorData, config, pixMapper, mode)` shared by sigma1 and edot1 — assembles tensor per sub-sampled cell, calls `eigen2x2()`, draws oriented arrows along eigenvector

## 9. PanelView — overlay controls UI

- [x] 9.1 Add collapsible "Overlays ▸" disclosure row below zoom row in PanelView constructor, collapsed by default
- [x] 9.2 Add 8 toggle checkboxes (Ph, T, Topo, V, Dir, σ₁, ε̇₁, Melt) inside the collapsible row, each dispatching `ctrl:panel:overlay-toggle` on click
- [x] 9.3 Add per-toggle gear icon that opens a settings popover with colour picker and layer-specific controls (ΔT for temperature, density/lengthScale sliders for vector layers, levels input for melt)
- [x] 9.4 Wire `ctrl:panel:overlay-config` event dispatch from popover inputs
- [x] 9.5 Add availability tracking: listen for overlay availability updates from model, disable/grey-out toggles for unavailable layers

## 10. Styling

- [x] 10.1 Add CSS for `.overlay-row`, `.overlay-toggle`, `.overlay-toggle:disabled`, `.overlay-gear`, `.overlay-popover` in style.css — compact layout matching existing panel-controls pattern, theme-aware colours

## 11. Integration and verification

- [x] 11.1 Start server with rifting dataset, enable phase boundary overlay — verify contour lines appear at phase transitions
- [x] 11.2 Enable temperature isotherm overlay, verify default ΔT=200 isotherms with inline °C labels, change ΔT and verify re-contour without server request
- [x] 11.3 Enable velocity overlay, verify arrows drawn at sub-sampled positions with correct direction
- [x] 11.4 Test conditional availability: load a file missing director data → verify Dir toggle is greyed out; load a file with director data → verify toggle becomes active
- [x] 11.5 Test overlay + crop/zoom: enable overlays, apply zoom crop → verify overlays clip correctly within cropped data area
- [x] 11.6 Test multi-panel: enable different overlays on different panels → verify independent per-panel config
