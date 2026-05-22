## Why

The current WebVis renders a single field at a time. Geodynamic workflows routinely compare two or more fields side-by-side (e.g. strain rate + stress, viscosity + temperature) and annotate each plot with a custom title that embeds live data values — exactly like the MATLAB/Julia plots MDOODZ users already produce. Adding multi-panel layout and an inline data-interpolation title bar brings WebVis up to parity with existing desktop visualisation scripts.

## What Changes

- **Multi-panel canvas layout**: the viewer can display 1–N panels arranged in a configurable grid (e.g. 1×2, 2×2). Each panel has its own field, colour map, colour range, and colour bar.
- **Per-panel title input with data interpolation**: each panel gets an editable title text field. Typing a trigger character (`:` or `$`) opens a dropdown listing available data variables. Selecting a variable inserts an interpolation token that is rendered with its current value, scientific unit, and proper formatting.
- **Data variable catalogue**: variables are split into two kinds:
  - *Scalar* — single-valued metadata (time step index, model time, Nx, Nz, dt) — inserted directly, no aggregation needed.
  - *Array* — the panel's loaded field data — the user picks an aggregation: min, max, mean, median, p2 (robust min), p98 (robust max).
- **Scientific value rendering**: interpolated tokens display values with SI-prefixed scientific notation and the field's `formattedUnit` from the field-labels registry.
- **Panel state independence**: each panel tracks its own field selection, colour map, colour range, and range-lock state. Changing one panel does not affect others.

## Capabilities

### New Capabilities
- `multi-panel-layout`: grid-based panel container, panel add/remove, responsive sizing, per-panel canvas + colour bar
- `title-interpolation`: editable title input per panel, trigger-character dropdown, variable catalogue (scalar + array with aggregations), scientific token rendering

### Modified Capabilities
- `web-ui`: top-level layout changes from single-canvas to multi-panel container; controls become panel-scoped
- `field-renderer`: renderer must accept a target canvas element rather than assuming a single global canvas; colour bar draws per-panel
- `web-server`: `/api/field-data` response already includes all needed stats (min, max, pMin, pMax); may need a new `/api/params/:file` endpoint exposing scalar metadata (step, time, dt, Nx, Nz) without loading a full field

## Impact

- **Frontend (`public/`)**: `index.html` layout overhaul (panels container replaces single canvas); new `PanelManager` view; new `TitleInput` view with dropdown; `ControlPanel` becomes panel-scoped; `FieldCanvas`, `ColourBar` become per-panel instances; `model.mjs` gains a `panels[]` array of per-panel state; `controller.mjs` routes events by panel id.
- **Server (`server.mjs`)**: possible new `/api/params/:file` lightweight endpoint.
- **No breaking API changes**: existing single-panel UX remains the default (1-panel grid).
- **Dependencies**: no new external libraries; dropdown uses native HTML `<datalist>` or a lightweight custom element.
