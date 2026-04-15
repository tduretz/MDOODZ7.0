## Why

The WebVis MVP renders MDOODZ HDF5 fields correctly but the user experience has several rough edges that make interactive exploration awkward: the colour range resets to auto on every file change (losing carefully chosen bounds), the auto-range algorithm includes statistical outliers (e.g., values near zero in strain-rate fields produce a collapsed range), file names are raw HDF5 filenames with no time information, field names and units lack scientific formatting, every field extraction re-reads HDF5 from scratch, and colour-map choices are not remembered per field.

## What Changes

- **Lock/save range toggle**: Add a toggle button that freezes the current colour range so it persists when switching files. When unlocked, revert to auto-range on each file change.
- **Robust auto-range (percentile clipping)**: Replace naive min/max with a percentile-based algorithm (e.g., 2nd–98th percentile) that ignores outlier tails (zeros, spikes). This eliminates the "all pink" problem seen with strain-rate fields.
- **Adaptive time display**: Show model time on top of the visualisation in human-readable units (yr, ka, Ma) that auto-adapt to magnitude; allow the user to override the unit.
- **Scientific field labels and unit formatting**: Display field names using proper scientific notation (τ\_II, ε̇\_II, σ\_xx, η, ρ, T, etc.) and format units with Unicode symbols (Pa·s, kg·m⁻³, s⁻¹, °C, etc.).
- **Disk cache for extracted field data**: Cache JSON responses to a local `.cache/` directory keyed by file+field so that revisiting a file/field serves from disk instead of re-reading HDF5.
- **Step + time in file selector**: Replace raw filenames in the dropdown and slider with "Step N — X.XX Ma" labels. Pre-extract time metadata from all files on startup.
- **Per-field colour-map memory**: Save the colour-map choice for each field name so switching between fields restores the previously chosen map.

## Capabilities

### New Capabilities
- `field-cache`: Server-side disk caching of extracted field JSON, keyed by filename + field name, with cache invalidation when source file changes.
- `field-labels`: Scientific display names and formatted units for each field in the registry — Unicode symbols, superscripts, subscripts for field titles and unit strings.
- `time-display`: Adaptive model-time formatting (yr / ka / Ma) with auto-unit selection and manual override, used in the header banner and file-selector labels.

### Modified Capabilities
- `web-ui`: Add range-lock toggle, per-field colour-map memory, step+time file-selector labels, header time display, and scientific field labels in dropdowns.
- `hdf5-reader`: Add percentile-based min/max computation (robust auto-range) and batch metadata extraction for all files.
- `web-server`: Add cache-aware field-data endpoint and bulk-metadata endpoint for file list with time/step info.
- `field-renderer`: Update colour bar to render scientific unit strings and handle locked vs auto range.

## Impact

- **Server** (`server.mjs`): New `/api/file-list` endpoint returning step+time for all files; field-data endpoint reads/writes disk cache; new cache directory `.cache/`.
- **Reader** (`hdf5-reader.mjs`): New `computePercentileRange()` function; new `readAllMetadata()` for batch time extraction.
- **Field registry** (`fields.mjs`): New `label` and `formattedUnit` properties per field entry.
- **Frontend** (`model.mjs`, `controller.mjs`, all views): Range-lock state, per-field colourmap map in model; ControlPanel gets lock toggle + scientific labels; StatusBar/header gets time display; ColourBar renders formatted units.
- **No new npm dependencies.** All changes are pure JS/CSS.
