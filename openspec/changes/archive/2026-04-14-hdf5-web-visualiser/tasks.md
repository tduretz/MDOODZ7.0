## 1. Project Scaffold

- [x] 1.1 Create `WebVis/` directory with `package.json` (`"type": "module"`, dependency: `h5wasm`)
- [x] 1.2 Run `npm install` to install h5wasm
- [x] 1.3 Create directory structure: `public/`, `public/views/`, `public/colormaps/`
- [x] 1.4 Create `public/index.html` — single-page shell with `<script type="module" src="app.mjs">`, CSS link, container divs for controls, canvas, colour bar, status bar

## 2. HDF5 Reader Module

- [x] 2.1 Create `hdf5-reader.mjs` — init h5wasm, implement `openFile(path)` with LRU cache (max 5 handles)
- [x] 2.2 Implement `readMetadata(filePath)` — read `/Model/Params`, return `{ time, nx, nz, dx, dz, dt }`
- [x] 2.3 Implement grid dimension lookup — map grid types (`center`, `vertex`, `vx`, `vz`, `vizgrid`) to `(rows, cols)` and coordinate dataset paths
- [x] 2.4 Implement `extractField(filePath, fieldDef)` — read dataset, reshape flat typed array to 2D based on grid type, read coordinate arrays, return `{ xCoords, zCoords, values, nx, nz }`
- [x] 2.5 Implement air-cell masking — read `/VizGrid/compo`, set NaN at cells where phase == -1 for centre-grid fields
- [x] 2.6 Implement vertex-to-centre averaging — 2×2 block average helper for sxz/exz interpolation
- [x] 2.7 Implement derived field: stress second invariant (τ_II) from sxxd, szzd, sxz
- [x] 2.8 Implement derived field: strain rate second invariant (ε̇_II) from exxd, ezzd, exz
- [x] 2.9 Implement `computeMinMax(values)` — min/max excluding NaN
- [x] 2.10 Test reader against `benchmark-results/rifting-600x400-v2-h5/Output00000.gzip.h5` — verify shapes, coordinate lengths, min/max sanity

## 3. Field Registry

- [x] 3.1 Create `fields.mjs` — declarative registry mapping field names to HDF5 paths, grid types, units, log flag, offset, discrete flag, and derive function names
- [x] 3.2 Register all centre fields: Pressure, Temperature, Viscosity, Density, Grain size, Melt fraction, Porosity, Overpressure, Stress xx/zz, Strain rate xx/zz, accumulated strains (total, elastic, plastic, pwl, exp, lin, gbs, pl_vol), eII partitions, divergence fields, cohesion, friction
- [x] 3.3 Register vertex fields: Shear stress (sxz), Shear strain rate (exz), Vertex viscosity, Vertex density
- [x] 3.4 Register velocity fields: Vx, Vz
- [x] 3.5 Register derived fields: Stress II, Strain rate II
- [x] 3.6 Register VizGrid fields: Phases (compo, discrete), Phases HR (compo_hr, discrete)
- [x] 3.7 Implement `listAvailableFields(filePath)` — check which registered fields actually exist in a given file and return their names

## 4. Web Server

- [x] 4.1 Create `server.mjs` — `node:http` server, parse CLI arg for data directory, `PORT` env variable, default port 3000
- [x] 4.2 Implement static file serving for `public/` — correct MIME types for `.html`, `.css`, `.mjs`, `.json`, `.js`; `/` serves `index.html`
- [x] 4.3 Implement path traversal protection — reject `..` in static paths, validate API filename against `Output\d+\.gzip\.h5` regex
- [x] 4.4 Implement `GET /api/files` — glob `Output*.gzip.h5` in data directory, return sorted filename list
- [x] 4.5 Implement `GET /api/fields/:file` — call reader for metadata + listAvailableFields, return `{ fields, params }`
- [x] 4.6 Implement `GET /api/field-data/:file/:field` — look up field in registry, call reader, return full JSON response with xCoords, zCoords, values, nx, nz, unit, log, min, max
- [x] 4.7 Implement gzip compression — pipe JSON through `zlib.createGzip()` when `Accept-Encoding: gzip`
- [x] 4.8 Add CORS header (`Access-Control-Allow-Origin: *`) to all API responses
- [x] 4.9 Add error handling — 404 for missing files/fields, 400 for invalid filenames, JSON error bodies
- [x] 4.10 Manual test: start server, curl all 3 API endpoints, verify JSON responses

## 5. Colour Maps

- [x] 5.1 Create `public/colormaps/maps.json` — 256-entry LUTs for viridis, turbo, inferno, plasma, coolwarm
- [x] 5.2 Create discrete phase palette (10–15 distinct colours for phase indices 0–14, transparent for -1)

## 6. Frontend Model

- [x] 6.1 Create `public/model.mjs` — class extending `EventTarget` with state properties: `files`, `currentFile`, `fields`, `currentField`, `fieldData`, `colourMap`, `colourRange`, `params`
- [x] 6.2 Implement setter methods that dispatch `CustomEvent` on change: `file-list-loaded`, `file-selected`, `fields-loaded`, `field-loaded`, `colourmap-changed`, `range-changed`, `loading-start`, `loading-end`

## 7. Frontend Views

- [x] 7.1 Create `public/views/FieldCanvas.mjs` — subscribe to `field-loaded`, `colourmap-changed`, `range-changed`; render 2D values to `<canvas>` via `ImageData` pixel-fill using LUT; handle NaN as transparent; support log scaling; preserve aspect ratio
- [x] 7.2 Create `public/views/ColourBar.mjs` — draw vertical gradient bar + tick labels (min, max, 3 intermediates) + unit label; support log-scale exponent labels
- [x] 7.3 Create `public/views/ControlPanel.mjs` — file dropdown, time-step slider (range input), field dropdown, colour-map dropdown, min/max inputs with auto-reset button
- [x] 7.4 Create `public/views/FileList.mjs` — populate file dropdown from model, sync slider position with selected file index (merged into ControlPanel)
- [x] 7.5 Create `public/views/StatusBar.mjs` — display model time (Myr), grid dimensions, current field name; subscribe to `field-loaded` and `file-selected`

## 8. Frontend Controller

- [x] 8.1 Create `public/controller.mjs` — wire control panel DOM events to API fetch calls; on file select → fetch fields → update model; on field select → fetch field-data → update model; on colourmap/range change → update model (no API call, re-render only)
- [x] 8.2 Implement loading state management — set model `loading-start` before fetch, `loading-end` after render complete
- [x] 8.3 Implement field preservation across file changes — when file changes, re-request same field name if available

## 9. Frontend Entry Point and Styling

- [x] 9.1 Create `public/app.mjs` — import Model, Controller, all Views; instantiate and wire together; fetch initial file list on load; auto-select first file and first field
- [x] 9.2 Create `public/style.css` — dark theme, CSS Grid layout (sidebar controls | main canvas + colour bar | bottom status bar), responsive from 1024–2560px
- [x] 9.3 Style controls — grouped panels, consistent input sizing, hover/focus states, loading spinner/overlay

## 10. Integration Testing

- [x] 10.1 End-to-end test: start server pointing at `benchmark-results/rifting-600x400-v2-h5/`, open in browser, verify file list loads, select a file, select "Pressure", verify canvas renders
- [x] 10.2 Test all grid types: select a centre field (P), vertex field (sxz), velocity field (Vx), derived field (Stress II), phase map (Phases) — verify each renders without error
- [x] 10.3 Test colour map switching — change maps, verify re-render without refetch
- [x] 10.4 Test colour range adjustment — narrow range, verify clamping; click auto to reset
- [x] 10.5 Test time-step slider — drag across range, verify canvas updates for each file
- [x] 10.6 Test log-scale fields — select Viscosity, verify log₁₀ scaling in colour bar labels
- [x] 10.7 Verify air-cell masking — select Temperature, confirm air cells are transparent

## 11. Documentation

- [x] 11.1 Create `WebVis/README.md` — installation (`npm install`), usage (`node server.mjs /path/to/h5/dir`), supported fields, architecture overview, screenshots placeholder
