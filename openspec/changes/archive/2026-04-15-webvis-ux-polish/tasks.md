## 1. Field Labels & Registry

- [x] 1.1 Add `label` and `formattedUnit` properties to every entry in `fields.mjs`
- [x] 1.2 Export a helper `getFieldLabel(name)` and `getFieldUnit(name)` from `fields.mjs`

## 2. Time Display Module

- [x] 2.1 Create `public/time-display.mjs` with `formatTime(seconds, overrideUnit?)` returning `{ value, unit, formatted }`
- [x] 2.2 Add `timeUnit` property to model with `time-unit-changed` event

## 3. Percentile Range (hdf5-reader)

- [x] 3.1 Implement `quickselect(arr, k)` partial-sort helper in `hdf5-reader.mjs`
- [x] 3.2 Implement `computePercentileRange(flat, pLow, pHigh)` using quickselect
- [x] 3.3 Update `extractField` to include `pMin`/`pMax` in the returned response

## 4. Batch Metadata (hdf5-reader)

- [x] 4.1 Implement `readAllMetadata(dataDir, filenames)` returning `[{ name, step, time }]`
- [x] 4.2 Parse step number from filename regex `Output(\d+)\.gzip\.h5`

## 5. Disk Cache (server)

- [x] 5.1 Create `field-cache.mjs` with `getCached(dataDir, filename, field)` and `writeCache(dataDir, filename, field, json)`
- [x] 5.2 Implement sha256-based cache key generation (16 hex chars + `.json`)
- [x] 5.3 Implement mtime comparison for cache invalidation
- [x] 5.4 Auto-create `.cache/` directory on first write

## 6. Server API Changes

- [x] 6.1 Update `GET /api/files` to return `{ files: [{ name, step, time }] }` using `readAllMetadata`
- [x] 6.2 Cache file-list metadata in memory after first call
- [x] 6.3 Integrate disk cache into `GET /api/field-data/:file/:field` endpoint
- [x] 6.4 Include `pMin`/`pMax` in field-data response

## 7. Frontend Model Updates

- [x] 7.1 Add `rangeLocked` boolean to model with `range-locked-changed` event
- [x] 7.2 Add `fieldColourMaps` Map to model with save/restore logic
- [x] 7.3 Add `pMin`/`pMax` to model state from field-data response
- [x] 7.4 Add `fileList` array of `{ name, step, time }` objects to model

## 8. Controller Updates

- [x] 8.1 Update `init()` to parse new `/api/files` response shape (objects instead of strings)
- [x] 8.2 Update `selectField()` to skip range update when `rangeLocked` is true
- [x] 8.3 Update `selectFile()` to skip range update when `rangeLocked` is true
- [x] 8.4 Add colour-map save to `fieldColourMaps` on `ctrl:cmap-change`
- [x] 8.5 Add colour-map restore from `fieldColourMaps` on field switch
- [x] 8.6 Use `pMin`/`pMax` for auto-range instead of `min`/`max`

## 9. Header Banner

- [x] 9.1 Replace `<h1>` in `index.html` with three-part header layout (title / field+unit / time+unit-selector)
- [x] 9.2 Create `public/views/HeaderBanner.mjs` view subscribing to model events
- [x] 9.3 Wire time-unit dropdown to model `timeUnit` property

## 10. Control Panel Updates

- [x] 10.1 Add range lock toggle button (🔒/🔓) next to Auto button
- [x] 10.2 Update file selector to show `"Step N — X.XX Ma"` labels from `fileList`
- [x] 10.3 Update file selector labels on `time-unit-changed` event
- [x] 10.4 Update field dropdown to show scientific `label` instead of plain name
- [x] 10.5 Update range inputs to initialise from `pMin`/`pMax`
- [x] 10.6 Update slider label to show step + time

## 11. Colour Bar & Renderer Updates

- [x] 11.1 Update colour bar to render `formattedUnit` and scientific `label`
- [x] 11.2 Add locked-range visual indicator (🔒 icon) on colour bar when range is locked
- [x] 11.3 Update normalisation to use pMin/pMax for auto-range bounds

## 12. Status Bar Updates

- [x] 12.1 Update status bar to use `formatTime` for adaptive time display
- [x] 12.2 Update status bar to show scientific field label instead of plain name

## 13. Smoke Test

- [x] 13.1 Start server, load files, verify step+time labels and percentile coloring
- [x] 13.2 Test range lock: lock range, switch files, verify range persists
- [x] 13.3 Test per-field cmap: set turbo for one field, switch away and back
- [x] 13.4 Test disk cache: reload same field, confirm cache file in `.cache/`
- [x] 13.5 Test time-unit override: switch to ka, verify all displays update
