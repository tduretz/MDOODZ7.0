## Context

WebVis is a lightweight browser-based visualiser for MDOODZ HDF5 simulation output. The current MVP has:
- A Node.js server (`server.mjs`) with three API routes (`/api/files`, `/api/fields/:file`, `/api/field-data/:file/:field`)
- An h5wasm-based reader (`hdf5-reader.mjs`) that extracts fields on every request
- A declarative field registry (`fields.mjs`) with plain-text display names and units
- An ES6+ MVC frontend with `EventTarget`-based model, 5 views, and a controller

Every field extraction re-reads HDF5 from disk. The file selector shows raw filenames. The colour range resets on each file change. Field names and units use ASCII text.

## Goals / Non-Goals

**Goals:**
- Eliminate re-extraction overhead via server-side disk cache
- Make time-stepping exploration seamless with range locking and per-field colour-map memory
- Provide publication-quality field labels and unit display (Unicode symbols)
- Show human-readable model time everywhere (adaptive yr/ka/Ma)
- Keep zero additional npm dependencies

**Non-Goals:**
- Client-side (browser) caching or service workers — out of scope
- User accounts or persistent preferences across browser sessions — out of scope
- Support for non-MDOODZ HDF5 files — out of scope
- Animated playback / GIF export — future feature
- Real-time file watching / auto-reload on new output — future feature

## Decisions

### 1. Percentile-based auto-range via partial sort

Compute p2/p98 percentiles on the flat `Float64Array` before reshaping. Use a partial-sort (quickselect) algorithm on the non-NaN values to find the 2nd and 98th percentile values in O(n) average time. Return both the full min/max and the percentile range. The frontend auto-range uses percentiles; manual range still allows the full extent.

*Alternative: histogram binning.* Rejected — requires choosing bin count and doesn't give exact percentile values. Quickselect is simpler and precise.

### 2. Disk cache keyed by `filename:field` with mtime invalidation

Create a `.cache/` directory under the data directory. Cache key = `sha256(filename + ":" + fieldName)` truncated to 16 hex chars + `.json`. On each field-data request, check if the cache file exists and its mtime is newer than the source HDF5 file. If valid, stream the cached JSON directly. If stale or missing, extract, write cache, then respond.

*Alternative: in-memory LRU cache.* Rejected — memory usage would grow unbounded for 47 fields × N files. Disk cache survives server restarts and costs only disk I/O (which is fast for JSON reads).

### 3. Bulk metadata endpoint: `/api/files` returns step + time

Change the existing `/api/files` route to return `{ files: [{ name, step, time }] }` instead of just filenames. On startup (or first call), read `/Model/Params[0]` from every HDF5 file and cache the results in memory. The step number is parsed from the filename regex `Output(\d+)\.gzip\.h5`. This is a **breaking** change to the API shape.

*Alternative: separate `/api/file-list` endpoint.* Rejected — adds complexity for no benefit; the old endpoint was only used internally.

### 4. Scientific labels and formatted units in the field registry

Add two new properties to each registry entry: `label` (HTML/Unicode string for display, e.g., `"τ_II"`, `"ε̇_II"`, `"σ_xx"`) and `formattedUnit` (Unicode unit string, e.g., `"Pa·s"`, `"kg·m⁻³"`, `"s⁻¹"`). The existing `unit` field stays for data purposes; `formattedUnit` is for display only. Labels use Unicode characters (no HTML tags) so they work in both `<option>` elements and canvas `fillText`.

*Alternative: HTML `<span>` with sub/superscript tags.* Rejected — `<option>` elements don't render HTML; Unicode characters work everywhere.

### 5. Range lock as model state with toggle button

Add a boolean `rangeLocked` property to the model. When locked, `selectFile()` and `selectField()` skip updating `colourRange`. The control panel gets a lock/unlock toggle button (🔒/🔓) next to the Auto button. Locking captures the current range; unlocking triggers an auto-range recalculation.

### 6. Per-field colour-map memory in model

Add a `Map<string, string>` property `fieldColourMaps` to the model. When the user changes the colour map, store `fieldColourMaps.set(currentField, newMap)`. When switching fields, restore from the map if an entry exists, otherwise keep the current map. This is session-only (not persisted to disk).

### 7. Adaptive time formatting as a pure function

Create a `formatTime(seconds)` utility that returns `{ value, unit, formatted }`:
- `< 1000 yr` → `"X.XX yr"`
- `< 1e6 yr` → `"X.XX ka"`
- `≥ 1e6 yr` → `"X.XX Ma"`

Allow manual override via a model property `timeUnit` (`null` = auto, `'yr'`/`'ka'`/`'Ma'` = forced). Used in the header banner, file-selector labels, and status bar.

### 8. Header banner for time and field info

Replace the current `<h1>MDOODZ WebVis</h1>` header with a three-part layout: left = app name, center = scientific field label + formatted unit, right = adaptive time display with unit selector dropdown. This gives immediate visual context without looking at controls.

### 9. extractField returns percentile range alongside min/max

The `extractField` response gains two new fields: `pMin` and `pMax` (percentile-clipped min/max). The frontend uses `pMin`/`pMax` for auto-range and shows `min`/`max` as the absolute bounds for manual range inputs. Cached JSON includes both.

## Risks / Trade-offs

- **[Cache staleness]** → Mitigated by mtime comparison on every request. If the HDF5 file is overwritten mid-simulation, cache is automatically invalidated.
- **[Disk space for cache]** → Each cached field is ~2–5 MB JSON for a 600×400 grid. 47 fields × 29 files ≈ 3–7 GB worst case. Acceptable for local development. Could add a cache-size limit later.
- **[Breaking `/api/files` shape]** → Only the WebVis frontend consumes this API, so no external consumers are affected.
- **[Quickselect worst-case O(n²)]** → Mitigated by using randomized pivot. For 240K elements (600×400), even worst-case is fast.
- **[Unicode labels in `<select>` options]** → Tested: major browsers render Unicode in `<option>` correctly. Fallback: ASCII names still stored as registry keys.
