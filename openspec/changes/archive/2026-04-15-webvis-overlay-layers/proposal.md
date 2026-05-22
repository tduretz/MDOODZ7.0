## Why

MDOODZ geodynamic publications overlay additional information on top of scalar field plots — phase boundary contours, velocity arrows, temperature isotherms, topography lines, and anisotropy director ticks (see Duretz et al. rifting figures). WebVis currently renders only a single scalar field per panel with no mechanism to superimpose derived quantities. Users cannot visually correlate strain rate with phase geometry, thermal structure, velocity flow patterns, or fabric orientation without switching between fields.

## What Changes

- **Overlay layer system**: Each panel gains a stack of independently toggleable overlay layers drawn on top of the field data but below axes/colour bar. Layers are defined by type (contour, vector, director, line) and fetch their own field data from the server.
- **Phase boundary contour layer**: Draws contour lines at phase transitions (integer boundaries in the Phase field). Line colour and width are user-adjustable per panel.
- **Temperature isotherm layer**: Draws contour lines at constant temperature intervals. User-configurable ΔT spacing (default 200 °C). Labels show temperature in °C at line breaks (inline labels). Line colour and width adjustable.
- **Topography line layer**: Draws the free-surface elevation as a single polyline from `/Topo/z_grid`. Line colour and width are user-adjustable. Available when the Topo group exists in the HDF5 file.
- **Velocity vector layer**: Fetches Vx and Vz fields, sub-samples to a configurable grid density, and draws scaled arrows. Arrow colour is user-adjustable.
- **Anisotropy director layer**: Fetches director nx and nz fields, draws oriented line segments (ticks) at sub-sampled grid positions to show fabric direction. Only available when director datasets exist in the HDF5 file. Tick colour and length are user-adjustable.
- **Principal stress σ₁ layer**: Computes the largest principal stress eigenvector from stress tensor components (sxxd, szzd, sxz, P). Draws oriented arrows at sub-sampled grid positions. Available when all four stress datasets exist.
- **Principal strain-rate ε̇₁ layer**: Computes the largest principal strain-rate eigenvector from strain-rate tensor components (exxd, ezzd, exz, divu). Draws oriented arrows at sub-sampled grid positions.
- **Melt fraction contour layer**: Draws contour lines at melt-fraction iso-values from `/Centers/X`. User-configurable contour levels.
- **Overlay controls UI**: A compact control row per panel for toggling layers on/off, adjusting colours/widths, and setting vector/director density. Temperature isotherms expose a ΔT input field.
- **Server-side overlay data endpoint**: A new API route that returns pre-extracted overlay field data in a single request to avoid multiple round-trips.

## Capabilities

### New Capabilities
- `overlay-layers`: Overlay layer model, per-panel layer stack, layer types (contour, vector, director, line), layer configuration (colour, width, density, visibility, contour spacing), overlay rendering pipeline, coordinate-to-pixel mapping for overlay drawing, overlay data fetching, marching-squares contour extraction (for phase boundaries, isotherms, melt fraction), inline contour label rendering, 2×2 eigenvalue solver (for principal stress/strain-rate), and staggered-grid interpolation (Vx/Vz to cell centres, sxz vertices to centres).

### Modified Capabilities
- `field-renderer`: The render pipeline gains an overlay drawing phase between field image and colour bar/axes. FieldCanvas must accept and draw overlay data using the computed data-area coordinates.
- `web-ui`: New overlay toggle controls in PanelView's control area. Server gains `/api/overlay-data/` endpoint.

## Impact

- **WebVis/public/views/FieldCanvas.mjs**: New overlay drawing methods inserted after field image draw and frame, before colour bar and axes. Must map overlay data coordinates to canvas pixel positions using the crop-aware data area. Contour renderer (marching squares), arrow renderer (velocity/σ₁/ε̇₁), tick renderer (directors), line renderer (topography), and inline label renderer (isotherms with °C values at line breaks).
- **WebVis/public/views/PanelView.mjs**: New overlay control row (layer toggles, colour pickers, width/density sliders, ΔT input for isotherms).
- **WebVis/public/model.mjs**: Panel state extended with `overlays` array (per-layer config: type, enabled, colour, lineWidth, density, contourSpacing).
- **WebVis/public/controller.mjs**: Fetches overlay data when file/dataset changes, caches per time step.
- **WebVis/server.mjs**: New `/api/overlay-data/:filename` route returning requested overlay fields in one response.
- **WebVis/public/overlays.mjs** (new): Marching-squares contour extraction, staggered-grid interpolation helpers, 2×2 symmetric eigenvalue solver for principal stress/strain-rate, and inline contour label positioning logic.
- **WebVis/fields.mjs**: Overlay fields already defined (Vx, Vz, Phase, Temperature, director nx/nz, stress/strain-rate components) but need to be grouped for overlay extraction. Topography (`/Topo/z_grid`) added as a 1D overlay field.
- **WebVis/hdf5-reader.mjs**: May need a bulk-extract helper for multiple fields in one HDF5 open.
