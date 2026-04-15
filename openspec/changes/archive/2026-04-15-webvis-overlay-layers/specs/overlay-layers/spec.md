## ADDED Requirements

### Requirement: Overlay layer model

Each panel SHALL maintain an `overlays` configuration object in its PanelState. The object SHALL contain one entry per layer type: `phases`, `temperature`, `topo`, `velocity`, `director`, `sigma1`, `edot1`, `melt`. Each entry SHALL store: `enabled` (boolean), `color` (hex string), `lineWidth` (number), and type-specific parameters (density, dT, levels, lengthScale). All layers SHALL default to `enabled: false`. The model SHALL emit a `panel:overlay-changed` event when any overlay property changes, carrying `panelId` and `layerType` in the event detail.

#### Scenario: Default overlay state

- **WHEN** a new panel is created
- **THEN** all overlay layers SHALL have `enabled: false`
- **AND** default colours SHALL be: phases=#ffffff, temperature=#000000, topo=#000000, velocity=#ffff00, director=#ffffff, sigma1=#ffffff, edot1=#888888, melt=#ff4444

#### Scenario: Toggle overlay emits event

- **WHEN** `model.setPanelOverlay(panelId, 'temperature', { enabled: true })` is called
- **THEN** a `panel:overlay-changed` event SHALL fire with `detail.panelId` and `detail.layerType === 'temperature'`

#### Scenario: Overlay config survives file change

- **WHEN** the user changes to a different time step
- **THEN** all overlay `enabled` states and style settings SHALL be preserved

### Requirement: Marching-squares contour extraction

The client SHALL implement a marching-squares algorithm that takes a 2D grid (column-major `values[ix][iz]`), coordinate arrays (`xCoords`, `zCoords`), and one or more threshold levels, and returns an array of polyline segments per level. Each segment SHALL be an object `{ x1, z1, x2, z2 }` in physical coordinates (SI metres). The algorithm SHALL handle NaN values by treating them as below any threshold (no contour through air cells).

#### Scenario: Single contour level on uniform gradient

- **WHEN** a 3×3 grid has values linearly increasing from 0 to 8 and the threshold is 4
- **THEN** the algorithm SHALL return segments forming a contour line crossing cells where the value transitions through 4

#### Scenario: Multiple contour levels

- **WHEN** called with levels `[200, 400, 600]` on a temperature grid
- **THEN** the result SHALL contain three separate arrays of segments, one per level

#### Scenario: NaN cells produce no contour

- **WHEN** a cell value is NaN (air)
- **THEN** no contour segment SHALL pass through that cell

### Requirement: Phase boundary contour layer

When the `phases` overlay is enabled, the renderer SHALL draw contour lines at phase transition boundaries. Contour levels SHALL be at integer midpoints (`0.5, 1.5, 2.5, ...`) up to the maximum phase ID in the visible region. The layer SHALL use the Phase field from `/VizGrid/compo` (vizgrid, `Nx-1 × Nz-1`). Line colour and width SHALL be configurable per panel.

#### Scenario: Two-phase model

- **WHEN** the domain contains phases 0 and 1 and the phases overlay is enabled
- **THEN** a contour line SHALL be drawn at the boundary between phases 0 and 1 (level 0.5)

#### Scenario: Phase overlay colour change

- **WHEN** the user changes the phase overlay colour to red
- **THEN** the phase contour lines SHALL render in red on the next frame

### Requirement: Temperature isotherm layer

When the `temperature` overlay is enabled, the renderer SHALL draw contour lines at temperature intervals defined by the user-configurable ΔT parameter (default 200 °C). Contour levels SHALL be generated as `[ΔT, 2×ΔT, 3×ΔT, ...]` up to the maximum temperature in the grid. The temperature field SHALL be read from `/Centers/T` with the −273.15 K offset applied (converting to °C) before contouring. Line colour and width SHALL be configurable per panel.

#### Scenario: Default ΔT spacing

- **WHEN** the temperature overlay is enabled with default settings and the model has T from 0°C to 1400°C
- **THEN** isotherms SHALL be drawn at 200, 400, 600, 800, 1000, 1200, 1400 °C

#### Scenario: Custom ΔT

- **WHEN** the user sets ΔT to 100
- **THEN** isotherms SHALL be drawn at 100, 200, 300, ..., up to the maximum temperature

#### Scenario: ΔT change re-contours without server fetch

- **WHEN** the user changes ΔT from 200 to 50
- **THEN** the contours SHALL be re-extracted client-side from the cached temperature grid
- **AND** no additional server request SHALL be made

### Requirement: Inline isotherm labels

Temperature contour lines SHALL display inline labels showing the temperature value in °C. Labels SHALL appear at periodic intervals along each contour polyline (approximately every 80 CSS pixels of polyline arc length). At each label position, a gap SHALL be inserted in the contour line and the text (e.g., "600 °C") SHALL be drawn centred in the gap. The label text SHALL be rotated to match the local contour tangent direction.

#### Scenario: Label at line break

- **WHEN** an isotherm contour is long enough for a label
- **THEN** the line SHALL have a gap at the label position
- **AND** the temperature value with "°C" suffix SHALL be rendered in the gap

#### Scenario: Short contour has no label

- **WHEN** an isotherm contour is shorter than 80 CSS pixels
- **THEN** no inline label SHALL be drawn for that contour segment

#### Scenario: Label rotation matches contour

- **WHEN** an isotherm contour runs diagonally
- **THEN** the label text SHALL be rotated to align with the local tangent of the contour

### Requirement: Topography line layer

When the `topo` overlay is enabled, the renderer SHALL draw a polyline representing the free-surface elevation. The data SHALL be read from `/Topo/z_grid` (1D array of length `Nx`) with x-coordinates from `/Model/xg_coord`. The polyline SHALL be drawn using `ctx.moveTo`/`ctx.lineTo` in physical coordinates mapped to the data-area pixel space. The layer SHALL only be available when the `/Topo` group exists in the HDF5 file. Line colour and width SHALL be configurable per panel.

#### Scenario: Topography line drawn

- **WHEN** the topo overlay is enabled and `/Topo/z_grid` exists
- **THEN** a continuous polyline SHALL be drawn at the free-surface elevation across the domain width

#### Scenario: Topography respects crop/zoom

- **WHEN** the view is cropped to x: [−50 km, 50 km]
- **THEN** the topography line SHALL be drawn only within the visible x range

#### Scenario: Missing Topo group disables layer

- **WHEN** the HDF5 file does not contain a `/Topo` group
- **THEN** the topo overlay toggle SHALL be disabled (greyed out) in the UI

### Requirement: Velocity vector layer

When the `velocity` overlay is enabled, the renderer SHALL draw arrows representing the velocity field. Vx and Vz SHALL be read from `/VxNodes/Vx` (staggered, `Nx × Nz+1`) and `/VzNodes/Vz` (`Nx+1 × Nz`). The client SHALL interpolate both to cell centres by averaging adjacent staggered-grid values. Arrows SHALL be drawn at sub-sampled grid positions controlled by a `density` parameter (default 8, meaning every 8th cell in each direction). Arrow length SHALL be proportional to velocity magnitude. Arrow colour SHALL be configurable per panel.

#### Scenario: Velocity arrows drawn

- **WHEN** the velocity overlay is enabled
- **THEN** arrows SHALL be drawn at sub-sampled cell-centre positions with length proportional to speed

#### Scenario: Density parameter controls spacing

- **WHEN** density is set to 4 on a 600×400 grid
- **THEN** approximately 150×100 = 15000 arrows SHALL be drawn (every 4th cell)

#### Scenario: Zero velocity produces no arrow

- **WHEN** both Vx and Vz are zero at a sample point
- **THEN** no arrow SHALL be drawn at that position

### Requirement: Anisotropy director layer

When the `director` overlay is enabled, the renderer SHALL draw oriented line segments (crossing ticks) showing fabric direction. Director components SHALL be read from `/Centers/nx` and `/Centers/nz` (cell centres, `Nx-1 × Nz-1`). The director vector `(nx, nz)` SHALL be transformed to the perpendicular fabric direction `(-nz/nx, 1)` and normalised. Ticks SHALL be drawn at sub-sampled positions controlled by a `density` parameter (default 6). Tick length SHALL be controlled by a `lengthScale` parameter. The layer SHALL only be available when both `/Centers/nx` and `/Centers/nz` exist. Anisotropy factor (`/Centers/ani_fac`) SHALL optionally be used to skip ticks where `ani_fac ≈ 0`.

#### Scenario: Director ticks drawn

- **WHEN** the director overlay is enabled and director datasets exist
- **THEN** oriented line segments SHALL be drawn at sub-sampled positions showing fabric direction

#### Scenario: Missing director data disables layer

- **WHEN** the HDF5 file does not contain `/Centers/nx`
- **THEN** the director overlay toggle SHALL be disabled in the UI

#### Scenario: Skip where anisotropy is negligible

- **WHEN** `ani_fac` is below 0.01 at a sample point
- **THEN** no tick SHALL be drawn at that position

### Requirement: Principal stress σ₁ layer

When the `sigma1` overlay is enabled, the renderer SHALL compute the largest principal stress eigenvector at each sub-sampled cell-centre position and draw oriented arrows. The stress tensor SHALL be assembled from: `σ_xx = -P + sxxd`, `σ_zz = -P + szzd`, `σ_xz = sxz` (from `/Vertices/sxz` averaged to centres). The 2×2 symmetric eigenvalue problem SHALL be solved using the closed-form formula. The layer SHALL only be available when all four datasets (`/Centers/sxxd`, `/Centers/szzd`, `/Vertices/sxz`, `/Centers/P`) exist. Arrow colour, density, and length scale SHALL be configurable per panel.

#### Scenario: σ₁ arrows drawn

- **WHEN** the sigma1 overlay is enabled and all stress datasets exist
- **THEN** arrows aligned with the maximum principal stress direction SHALL be drawn at sub-sampled positions

#### Scenario: Missing stress data disables layer

- **WHEN** `/Centers/sxxd` does not exist in the HDF5 file
- **THEN** the sigma1 overlay toggle SHALL be disabled in the UI

#### Scenario: Eigenvector orientation

- **WHEN** σ_xx = 1, σ_zz = 0, σ_xz = 0 at a point
- **THEN** the σ₁ eigenvector SHALL be aligned horizontally (along x)

### Requirement: Principal strain-rate ε̇₁ layer

When the `edot1` overlay is enabled, the renderer SHALL compute the largest principal strain-rate eigenvector at each sub-sampled cell-centre position and draw oriented arrows. The strain-rate tensor SHALL be assembled from: `ε̇_xx = exxd`, `ε̇_zz = ezzd`, `ε̇_xz = exz` (from `/Vertices/exz` averaged to centres), with the isotropic part `−1/3 × divu`. The same 2×2 eigenvalue solver as σ₁ SHALL be used. The layer SHALL only be available when all four datasets (`/Centers/exxd`, `/Centers/ezzd`, `/Vertices/exz`, `/Centers/divu`) exist. Arrow colour, density, and length scale SHALL be configurable per panel.

#### Scenario: ε̇₁ arrows drawn

- **WHEN** the edot1 overlay is enabled and all strain-rate datasets exist
- **THEN** arrows aligned with the maximum principal strain-rate direction SHALL be drawn at sub-sampled positions

#### Scenario: Missing strain-rate data disables layer

- **WHEN** `/Centers/exxd` does not exist in the HDF5 file
- **THEN** the edot1 overlay toggle SHALL be disabled in the UI

### Requirement: Melt fraction contour layer

When the `melt` overlay is enabled, the renderer SHALL draw contour lines at user-configurable melt-fraction iso-values from `/Centers/X`. Default contour levels SHALL be `[0.01, 0.05, 0.1]`. The levels array SHALL be configurable per panel. Line colour and width SHALL be configurable per panel. The layer SHALL only be available when `/Centers/X` exists.

#### Scenario: Default melt contours

- **WHEN** the melt overlay is enabled with default levels
- **THEN** contour lines SHALL be drawn at melt fractions 0.01, 0.05, and 0.1

#### Scenario: Custom melt levels

- **WHEN** the user sets levels to `[0.02, 0.1, 0.2]`
- **THEN** contour lines SHALL redraw at those values without a server fetch

#### Scenario: Missing melt data disables layer

- **WHEN** `/Centers/X` does not exist in the HDF5 file
- **THEN** the melt overlay toggle SHALL be disabled in the UI

### Requirement: Batched overlay data endpoint

The server SHALL provide a `GET /api/overlay-data/:filename?layers=<comma-separated>` endpoint. The `layers` parameter SHALL accept a comma-separated list of layer type keys: `phases`, `temperature`, `velocity`, `topo`, `director`, `sigma1`, `edot1`, `melt`. The response SHALL be a JSON object keyed by layer type, containing the raw flat arrays and grid dimensions needed for client-side processing. Only layers whose required HDF5 datasets exist SHALL be included in the response. The response SHALL be gzip-compressed when the client accepts it. Vertex fields (sxz, exz) SHALL be pre-averaged to cell centres on the server.

#### Scenario: Request multiple layers

- **WHEN** the client requests `?layers=phases,temperature,velocity`
- **THEN** the response SHALL contain keys `phases`, `temperature`, and `velocity` with their respective data arrays

#### Scenario: Missing layer omitted

- **WHEN** `?layers=phases,director` is requested but `/Centers/nx` does not exist
- **THEN** the response SHALL contain `phases` but not `director`

#### Scenario: Invalid filename rejected

- **WHEN** the filename does not match `Output\d+\.gzip\.h5`
- **THEN** the server SHALL respond with 400 status

#### Scenario: Gzip compression

- **WHEN** the client sends `Accept-Encoding: gzip`
- **THEN** the response body SHALL be gzip-compressed

### Requirement: Overlay data caching

The controller SHALL cache overlay data per `(filename, layerSet)` key. When the user enables a new layer type whose data is not cached, the controller SHALL fetch a batch including the missing layer. When the file or time step changes, all overlay caches for the previous file SHALL be discarded. When the same overlay data is already cached for the current file, no additional fetch SHALL be made.

#### Scenario: Cache hit on re-enable

- **WHEN** the user disables and re-enables the temperature overlay on the same file
- **THEN** no new fetch SHALL be issued

#### Scenario: File change invalidates cache

- **WHEN** the user switches from Output00010 to Output00020
- **THEN** all cached overlay data for Output00010 SHALL be discarded

#### Scenario: New layer triggers fetch

- **WHEN** the user enables the director layer and it has not been fetched for the current file
- **THEN** the controller SHALL fetch overlay data including `director`

### Requirement: Coordinate-to-pixel mapping for overlays

All overlay renderers SHALL use the same coordinate-to-pixel mapping as the field image. Physical coordinates (SI metres) SHALL be mapped to pixel positions using the data-area rectangle (`dataX`, `dataY`, `dW`, `dH`) and the cropped coordinate range (`col0`, `col1`, `row0`, `row1`). Overlay drawing SHALL be clipped to the data-area rectangle using `ctx.clip()` to prevent lines or arrows from extending into margin areas.

#### Scenario: Contour line at known coordinate

- **WHEN** a contour passes through the physical midpoint of the domain
- **THEN** the line segment SHALL be drawn at the pixel centre of the data area

#### Scenario: Overlay clipped at crop boundary

- **WHEN** a contour extends beyond the cropped view bounds
- **THEN** the line SHALL be clipped at the data-area rectangle edge

### Requirement: 2×2 symmetric eigenvalue solver

The client SHALL implement a closed-form 2×2 symmetric matrix eigenvalue decomposition. Given a symmetric matrix `[[a, b], [b, c]]`, the solver SHALL return eigenvalues `λ₁ ≥ λ₂` and the unit eigenvector `(cos θ, sin θ)` of the largest eigenvalue, where `θ = 0.5 × atan2(2b, a − c)`. The solver SHALL handle the degenerate case `a == c` (θ = π/4 when b > 0, −π/4 when b < 0, 0 when b == 0).

#### Scenario: Diagonal matrix

- **WHEN** the matrix is `[[3, 0], [0, 1]]`
- **THEN** λ₁ SHALL be 3 and the eigenvector SHALL be `(1, 0)`

#### Scenario: Off-diagonal matrix

- **WHEN** the matrix is `[[0, 1], [1, 0]]`
- **THEN** λ₁ SHALL be 1 and the eigenvector SHALL be `(cos(π/4), sin(π/4))` ≈ `(0.707, 0.707)`

#### Scenario: Equal diagonal

- **WHEN** the matrix is `[[2, 1], [1, 2]]`
- **THEN** λ₁ SHALL be 3 and λ₂ SHALL be 1

### Requirement: Staggered-grid interpolation

The client SHALL interpolate staggered-grid velocity fields to cell centres for arrow placement. Vx (on vertical faces, `Nx × Nz+1`) SHALL be averaged as `Vxc[ix][iz] = 0.5 * (Vx[ix][iz] + Vx[ix][iz+1])` to produce `Nx × Nz` values. Vz (on horizontal faces, `Nx+1 × Nz`) SHALL be averaged as `Vzc[ix][iz] = 0.5 * (Vz[ix][iz] + Vz[ix+1][iz])` to produce `Nx × Nz` values.

#### Scenario: Vx interpolation

- **WHEN** Vx at face [3][5] = 2.0 and Vx at face [3][6] = 4.0
- **THEN** centre value Vxc[3][5] SHALL be 3.0

#### Scenario: Vz interpolation

- **WHEN** Vz at face [3][5] = 1.0 and Vz at face [4][5] = 3.0
- **THEN** centre value Vzc[3][5] SHALL be 2.0

### Requirement: Conditional layer availability

The client SHALL track which overlay layers have data available for the current file. After fetching overlay data, layers whose data was present in the response SHALL be marked as available. Layers whose data was absent SHALL be marked as unavailable. The UI SHALL disable (grey out) toggles for unavailable layers. Enabling an unavailable layer SHALL have no effect.

#### Scenario: Available layer is toggleable

- **WHEN** the overlay data response includes `phases`
- **THEN** the phases toggle SHALL be enabled and clickable

#### Scenario: Unavailable layer is greyed out

- **WHEN** the overlay data response does not include `director`
- **THEN** the director toggle SHALL be disabled and visually greyed out

#### Scenario: Availability refreshes on file change

- **WHEN** the user switches to a file that contains director data (previously missing)
- **THEN** the director toggle SHALL become enabled
