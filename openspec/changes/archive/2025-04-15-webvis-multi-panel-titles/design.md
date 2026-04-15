## Context

WebVis is a zero-dependency browser-based HDF5 viewer for MDOODZ output. The current architecture is a single-panel MVC app:

- **Model** (`model.mjs`): single `EventTarget` holding one active field, one colour range, one colour map, one loaded file.
- **Controller** (`controller.mjs`): fetches from `/api/files`, `/api/fields/:file`, `/api/field-data/:file/:field`; writes to model.
- **Views**: `FieldCanvas`, `ColourBar`, `ControlPanel`, `StatusBar`, `HeaderBanner`, `LoadingOverlay` — all subscribe to model events and render from global state.

The screenshot reference shows the MATLAB-style layout: two field plots stacked vertically, each with its own colour bar and a title like `ε̇_II at t = 4.6928 Ma  - strain = 0 %  min = 0.00e+00 s⁻¹  max = 8.35e-14 s⁻¹`. The old MATLAB plots used ASCII abbreviations (eII, sII); WebVis already has a Unicode field-labels registry (ε̇_II, τ_II, η, etc.) so titles will use those scientific symbols instead.

## Goals / Non-Goals

**Goals:**
- Display 1–4 panels simultaneously in a CSS grid, each rendering an independent field with its own canvas, colour bar, and controls.
- Provide an editable title input per panel where the user types free text mixed with data tokens (e.g. `:time`, `:max`).
- Offer a trigger-character autocomplete dropdown that lets the user pick from scalar metadata or array field statistics without memorising token syntax.
- Render interpolated title tokens with scientific notation and units from the field-labels registry.
- Keep the default 1-panel experience identical to today (no regression).

**Non-Goals:**
- Arbitrary N-panel layouts (max 4 panels is sufficient).
- Drag-and-drop panel reordering or resizing.
- Cross-panel linked interactions (e.g. linked zoom, synchronised sliders). Each panel is independent.
- Server-side rendering or image export (future work).
- LaTeX rendering engine — titles use Unicode only (same as existing field-labels).

## Decisions

### D1: Panel state lives in Model as an array of PanelState objects

Each panel gets a `PanelState` plain object:
```js
{ id, fieldName, colourMap, colourRange, rangeLocked, fieldData, title }
```
`Model` gains a `panels` array and an `activePanelId`. Panel-scoped events carry `{ panelId }` in detail so views know which panel changed.

**Why not separate Model instances per panel?** A single Model keeps the file list, file selection, and params shared. Only the field-specific state fans out per panel.

**Alternative considered:** Make every panel a fully independent mini-app with its own model/controller. Rejected — too much duplication for shared state (file list, time unit, global metadata).

### D2: CSS Grid layout with 1×1, 1×2, 2×1, 2×2 presets

A `#panels-grid` container uses `display: grid` with `grid-template-columns` / `grid-template-rows` set via a CSS class (`.grid-1x1`, `.grid-1x2`, `.grid-2x1`, `.grid-2x2`). A layout selector in the header lets the user switch.

**Why CSS Grid over `flexbox`?** Grid gives explicit 2D placement; panels stay equal-sized without percentage hacks.

### D3: Each panel is a self-contained DOM subtree managed by a PanelView

`PanelView` creates: `<div class="panel">` → title input + canvas + colour bar + mini control row (field select, cmap select, lock button). Each `PanelView` instantiates its own `FieldCanvas` and `ColourBar` bound to one `PanelState`.

Panel-scoped events bubble as `ctrl:panel:<event>` with `detail: { panelId, ... }`.

### D4: Title interpolation uses `${...}` token syntax with trigger dropdown

The title input is a regular `<input type="text">`. When the user types `$` (or clicks an insert button), a floating dropdown appears listing available variables grouped by category:

**Scalars** (file-level metadata, no aggregation):
| Token | Dropdown label | Display example |
|-------|---------------|-----------------|
| `${step}` | Step number | `42` |
| `${time}` | Model time | `4.69 Ma` |
| `${Nx}` | Grid columns | `600` |
| `${Nz}` | Grid rows | `400` |
| `${resolution}` | Grid resolution | `600 × 400` |

**Field statistics** (computed from current panel's field data):
| Token | Dropdown label | Display example |
|-------|---------------|-----------------|
| `${min}` | Min | `0.00e+00 Pa` |
| `${max}` | Max | `9.35e+08 Pa` |
| `${mean}` | Mean | `4.21e+08 Pa` |
| `${median}` | Median | `3.98e+08 Pa` |
| `${p2}` | Robust min (p2) | `5.51e+07 Pa` |
| `${p98}` | Robust max (p98) | `3.09e+09 Pa` |
| `${field}` | Field name | `τ_II` |
| `${unit}` | Unit | `Pa` |

The dropdown shows human-readable names (e.g. "Model time", "Robust max (p98)") grouped under **Metadata** and **Field stats** headers. Selecting an item inserts the token text `${...}` at the cursor position.

**Why `${...}` over `:name`?** Unambiguous parsing — colons appear in titles naturally (e.g. "Stress: ..."). Dollar-brace is a standard interpolation syntax that won't conflict with prose.

**Why not a rich contentEditable div?** A plain text input with post-render interpolation is simpler, more accessible, and avoids contentEditable's cross-browser complexity. The rendered title appears as a separate read-only `<span>` below the input.

### D5: Title rendering is a pure function `renderTitle(template, context) → string`

Given a template string like `"${field} at t = ${time} — max = ${max}"` and a context object `{ step, time, Nx, Nz, min, max, mean, median, p2, p98, field, unit }`, the function returns the rendered string with formatted values.

Mean and median are computed lazily client-side from the flat values array only when the template references `${mean}` or `${median}`.

### D6: `mean` and `median` computed client-side on demand

They are cheap O(n) / O(n) operations on the already-loaded flat array. No server change needed. Computation happens once per field load in the panel, cached in `PanelState`.

### D7: File selection is still global (shared across all panels)

All panels view data from the same HDF5 file. The file selector and time slider remain global controls. When the file changes, each panel re-fetches its own field independently in parallel.

**Why not per-panel file selection?** Comparing fields at different time steps is a rare workflow; sharing the file keeps the slider/playback experience simple and avoids multiplying API calls.

### D8: New lightweight `/api/params/:file` endpoint

Returns `{ step, time, Nx, Nz, dt }` without extracting any field data. The server already reads `Params` in `readAllMetadata`; this endpoint reads it for a single file. Used for scalar title tokens.

## Risks / Trade-offs

**[Performance] Multi-panel = multiple field extractions per file switch** → Mitigation: disk cache already exists; parallel fetch calls; limit to 4 panels max.

**[Complexity] Model becomes more complex with panel array** → Mitigation: `PanelState` is a plain object, not a class. The Model merely dispatches panel-scoped events. Views only listen to events matching their panel ID.

**[UX] Title input syntax may not be discoverable** → Mitigation: default titles are pre-filled with a useful template (e.g. `"${field} at t = ${time} — min = ${min}  max = ${max}"`). The insert button (⊕ icon) next to the input makes the dropdown always discoverable without needing to know the trigger character.

**[Backwards compat] Existing 1-panel users see no change** → Default layout is 1×1 with a single panel. All current controls work as before. The layout selector is an opt-in addition.
