## Context

WebVis is a zero-dependency browser HDF5 visualiser for MDOODZ. Each panel renders a single scalar field on a HiDPI-aware `<canvas>` via `FieldCanvas`, with crop/zoom, axis ticks, and colour bar/phase legend. The rendering pipeline is: offscreen ImageData → main canvas blit → frame → colour bar → axes → title, all positioned relative to a computed `dataX, dataY, dW, dH` data rectangle in CSS pixels.

Fields are served by the Node.js server from HDF5 files using h5wasm, with a field registry (`fields.mjs`), caching layer (`field-cache.mjs`), and an LRU file-handle pool (`hdf5-reader.mjs`). The controller fetches field data per panel and stores it in the model's panel state.

The staggered grid means different fields live on different grids: cell centres (`Nx-1 × Nz-1`), vertices (`Nx × Nz`), Vx faces (`Nx × Nz+1`), Vz faces (`Nx+1 × Nz`), and vizgrid. Overlays that combine fields from different grids need interpolation to a common grid.

## Goals / Non-Goals

**Goals:**
- Enable up to 8 overlay layer types drawn on top of scalar fields: phase boundaries, temperature isotherms (with configurable ΔT and inline °C labels), topography line, velocity vectors, anisotropy directors, principal stress σ₁, principal strain-rate ε̇₁, melt fraction contours
- Per-panel overlay configuration: each layer independently toggled, coloured, and styled
- Overlay data fetched in a single batched API request per time step to avoid N round-trips
- Client-side contour extraction (marching squares) and eigenvalue decomposition — no server-side geometry computation
- Overlays respect crop/zoom bounds and use the same coordinate-to-pixel mapping as the field image

**Non-Goals:**
- Interactive overlay drawing or editing (overlays are read-only)
- Streamlines (line integration) — future follow-up
- 3D rendering or particle trajectory plots
- Overlay persistence across sessions (state is ephemeral)
- Server-side contour extraction or image rendering

## Decisions

### 1. Client-side contour extraction via marching squares

**Decision:** Implement marching-squares ISO-contour extraction in the browser (`overlays.mjs`), not on the server.

**Rationale:** Contour level parameters (ΔT, melt threshold) are user-adjustable per panel. Shipping raw 2D grids and contouring client-side avoids a server round-trip on every parameter change. The grids are small (600×400 ≈ 240K cells) so the computation is sub-millisecond in JS. Julia's `Mak.contour!` does the same — contour from grid data.

**Alternatives considered:**
- Server-side contouring: simpler client, but requires a request per parameter change and duplicates logic per contour type. Rejected.
- WebGL contour shader: overkill for line rendering at this grid size. Rejected.

### 2. Batched overlay data endpoint

**Decision:** New `GET /api/overlay-data/:filename?layers=...` endpoint. The `layers` query parameter is a comma-separated list of layer type keys (e.g., `phases,velocity,temperature,topo`). The server returns a JSON object keyed by layer type, each containing the raw field arrays and coordinate arrays needed for that layer.

**Rationale:** A single request avoids the latency of multiple sequential or parallel fetches for 3–5 separate fields. The server already has the HDF5 file open in its LRU cache, so reading multiple datasets from one open is cheap.

**Response shape:**
```json
{
  "phases":      { "compo": [...], "nx": N, "nz": M, "xCoords": [...], "zCoords": [...] },
  "temperature": { "T": [...], "nx": N, "nz": M, "xCoords": [...], "zCoords": [...] },
  "velocity":    { "Vx": [...], "Vz": [...], "vxNx": ..., "vxNz": ..., "vzNx": ..., "vzNz": ..., "xCoords": [...], "zCoords": [...] },
  "topo":        { "z_grid": [...], "x_grid": [...] },
  "director":    { "nx": [...], "nz": [...], "ani_fac": [...], "gridNx": N, "gridNz": M, "xCoords": [...], "zCoords": [...] },
  "sigma1":      { "sxxd": [...], "szzd": [...], "sxz_c": [...], "P": [...], "nx": N, "nz": M, "xCoords": [...], "zCoords": [...] },
  "edot1":       { "exxd": [...], "ezzd": [...], "exz_c": [...], "divu": [...], "nx": N, "nz": M, "xCoords": [...], "zCoords": [...] },
  "melt":        { "X": [...], "nx": N, "nz": M, "xCoords": [...], "zCoords": [...] }
}
```

Vertex fields (sxz, exz) are pre-averaged to cell centres on the server to keep the client simple. Vx/Vz stay on their native staggered grids — the client interpolates to cell centres for arrow placement (simple neighbour average).

**Alternatives considered:**
- Reuse existing `/api/field-data/` per field: works but 5–8 requests per overlay update. Rejected for latency.
- WebSocket streaming: complexity not justified for ~200 KB payloads. Rejected.

### 3. Overlay rendering order and insertion point in FieldCanvas

**Decision:** Overlays render in a fixed order after the frame `strokeRect` and before the colour bar, axes, and title. Order: topography line → filled/area overlays (none currently) → contour lines (phases, isotherms, melt) → vector/tensor overlays (velocity, σ₁, ε̇₁, directors).

**Rationale:** Contour lines on top of topography line, arrows on top of contour lines gives the best visual layering. Drawing before colour bar/axes ensures overlays don't extend into the margin area (clipped by the data rectangle).

**Implementation:** A new method `_drawOverlays(ctx, overlayData, overlayConfig, dataX, dataY, dW, dH, col0, col1, row0, row1)` in FieldCanvas, called between `ctx.strokeRect(...)` and `_drawCbar(...)`. Canvas clipping (`ctx.save(); ctx.beginPath(); ctx.rect(dataX, dataY, dW, dH); ctx.clip()`) ensures nothing leaks outside the data area.

### 4. Overlay configuration stored in panel state

**Decision:** Each panel gains an `overlays` map in its state:
```js
overlays: {
  phases:      { enabled: false, color: '#ffffff', lineWidth: 1.5 },
  temperature: { enabled: false, color: '#000000', lineWidth: 1.5, dT: 200 },
  topo:        { enabled: false, color: '#000000', lineWidth: 2 },
  velocity:    { enabled: false, color: '#ffff00', density: 8 },
  director:    { enabled: false, color: '#ffffff', density: 6, lengthScale: 1.5 },
  sigma1:      { enabled: false, color: '#ffffff', density: 6, lengthScale: 1.5 },
  edot1:       { enabled: false, color: '#888888', density: 6, lengthScale: 1.5 },
  melt:        { enabled: false, color: '#ff4444', lineWidth: 2, levels: [0.01, 0.05, 0.1] },
}
```

**Rationale:** Per-panel config lets users show different overlays on different panels (e.g., panel 1: phases + velocity, panel 2: isotherms + directors). The model emits `panel:overlay-changed` events so FieldCanvas re-renders reactively.

**Alternatives considered:**
- Global overlay config (same overlays on all panels): simpler but less useful for multi-panel comparison. Rejected.

### 5. Marching-squares contour algorithm

**Decision:** Implement a standard marching-squares algorithm in `overlays.mjs` that takes a 2D grid (column-major `values[ix][iz]`), coordinate arrays, and a threshold level, and returns an array of polyline segments `[{x1,z1,x2,z2}, ...]`. A higher-level `extractContours(grid, xCoords, zCoords, levels)` returns grouped polylines per level.

For phase boundaries, the algorithm runs on the phase-ID field with levels at `0.5, 1.5, 2.5, ...` (integer midpoints), producing contour lines at every phase transition.

For temperature isotherms, levels are generated from `dT`: `[dT, 2*dT, 3*dT, ...]` up to the max temperature in the grid. The raw Temperature field from the server is in K; the offset (−273.15) is applied before contouring so levels are in °C.

**Rationale:** Marching squares is O(N) in grid cells, trivially fast for 240K cells. Re-running on parameter change (ΔT slider) is instant.

### 6. Inline isotherm labels

**Decision:** Temperature contour lines have inline labels showing the temperature value in °C at periodic gaps along each polyline. Implementation: walk the polyline, at intervals of ~80 px insert a gap in the line and draw the label text (e.g., "600 °C") centred in the gap, rotated to match the local contour tangent.

**Rationale:** This matches standard geodynamic publication style (Duretz et al.) where isotherms are labelled at line breaks rather than in a legend. Users specifically requested "T in Celsius in a break in line".

### 7. 2×2 symmetric eigenvalue solver for principal stress/strain-rate

**Decision:** Implement a closed-form 2×2 symmetric eigenvalue decomposition in `overlays.mjs`:
$$\sigma = \begin{pmatrix} \sigma_{xx} & \sigma_{xz} \\ \sigma_{xz} & \sigma_{zz} \end{pmatrix}$$

Eigenvalues:
$$\lambda_{1,2} = \frac{\sigma_{xx} + \sigma_{zz}}{2} \pm \sqrt{\left(\frac{\sigma_{xx} - \sigma_{zz}}{2}\right)^2 + \sigma_{xz}^2}$$

Eigenvector of σ₁ (largest eigenvalue): `(cos θ, sin θ)` where $\theta = \frac{1}{2}\arctan\!\left(\frac{2\sigma_{xz}}{\sigma_{xx} - \sigma_{zz}}\right)$.

For stress: $\sigma_{xx} = -P + \sigma'_{xx}$, $\sigma_{zz} = -P + \sigma'_{zz}$, $\sigma_{xz}$ from vertices (averaged to centres on server).

For strain rate: same formula with $\epsilon'_{xx}, \epsilon'_{zz}, \epsilon_{xz}$, and $-\frac{1}{3}\nabla\cdot\mathbf{v}$ as the isotropic part.

**Rationale:** Closed-form is exact for 2×2 matrices and avoids iterative solvers. This matches Julia's `PrincipalStress()` function. Computed per grid cell at sub-sampled positions only (every `density`-th cell), so at density=6 on a 600×400 grid, only ~6700 eigensolves.

### 8. Overlay UI controls in PanelView

**Decision:** A collapsible "Overlays" row below the zoom row in each panel. Contains a compact grid of toggle checkboxes with layer icons/labels, plus a small settings popover per layer for colour, line width, density, and dT. The row is collapsed by default to avoid cluttering the interface.

**Rationale:** Overlays are secondary to the main field selection. A collapsible section keeps the UI clean for users who don't need overlays, while providing full control for those who do.

**Alternatives considered:**
- Sidebar overlay panel: breaks the per-panel paradigm and complicates the layout. Rejected.
- Context menu on right-click: poor discoverability. Rejected.

### 9. Overlay data caching strategy

**Decision:** The controller caches overlay data per `(filename, layerSet)` key. When the user toggles a new layer, the controller checks if its data is already cached; if not, it fetches a new batch including the missing layer. When the file/time step changes, all overlay caches for the previous file are discarded.

**Rationale:** Most users toggle overlays while viewing a single time step, so caching avoids re-fetching. The layerSet key (sorted list of enabled layer types) ensures cache hits when the same combination is reused.

### 10. Topography as a 1D polyline overlay

**Decision:** Topography is a special case: a 1D array `z_grid[x]` rather than a 2D field. The server reads `/Topo/z_grid` (length `Nx`) and `/Model/xg_coord` (vertex x-coordinates). The client draws it as a `ctx.beginPath(); moveTo/lineTo; stroke()` polyline.

**Rationale:** Topography is fundamentally different from 2D contour/vector overlays — it's a single curve. Handling it separately avoids forcing it through the marching-squares or arrow-drawing pipelines.

### 11. Conditional layer availability

**Decision:** The `/api/overlay-data/` endpoint checks dataset existence for each requested layer and only returns layers whose data is present. The client UI disables (grays out) toggles for layers whose data was not returned.

**Rationale:** Not all HDF5 files contain all datasets (e.g., director fields only exist when anisotropy is enabled, Topo group only exists with free surface). The UI must gracefully handle missing layers.

## Risks / Trade-offs

**[Marching-squares performance on large grids]** → Mitigation: For a 600×400 grid, marching squares takes <5 ms. For future 2000×1400 grids, it may reach 50 ms. If needed, offload to a Web Worker. Current grid sizes are safe.

**[Overlay data payload size]** → Mitigation: The batched endpoint returns raw Float64 arrays for up to 8 layers (~2 MB uncompressed for all layers on a 600×400 grid). With gzip compression (already implemented in the server), this drops to ~200–400 KB. Acceptable for local/LAN use.

**[Visual clutter with multiple overlays]** → Mitigation: Overlays are off by default. Per-layer colour and opacity controls let users manage visual density. The collapsible UI discourages enabling everything at once.

**[Staggered grid interpolation errors]** → Mitigation: Vertex-to-centre averaging is already implemented in `hdf5-reader.mjs` (`vertexToCentre`). Vx/Vz interpolation uses simple 2-point averaging along the staggered dimension. These are standard second-order accurate interpolations, matching the simulation's own internal stencils.

**[Canvas clipping edge artefacts]** → Mitigation: The `ctx.clip()` rectangle is set to `(dataX, dataY, dW, dH)` matching the field image. Overlay lines that cross the clip boundary are cleanly cut. Anti-aliased lines at the boundary may show faint sub-pixel artefacts but this is standard canvas behaviour.
