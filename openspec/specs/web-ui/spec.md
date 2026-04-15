## ADDED Requirements

### Requirement: Refresh datasets button

The UI SHALL display a refresh button (⟳ icon or "Refresh" label) near the dataset selector. Clicking the button SHALL send `POST /api/rescan` to the server, then re-fetch the dataset list via `GET /api/datasets`, and update the dataset selector dropdown. After updating datasets, the file list for the active dataset SHALL also be refreshed. A browser page refresh (F5 / Cmd+Shift+R) SHALL trigger the same rescan-then-fetch sequence on page load.

#### Scenario: Click refresh button

- **WHEN** the user clicks the refresh button
- **THEN** the client SHALL POST `/api/rescan`, then GET `/api/datasets` to update the dropdown
- **AND** the file list SHALL be refreshed for the active dataset

#### Scenario: Page refresh triggers rescan

- **WHEN** the user refreshes the browser page
- **THEN** on page load the client SHALL POST `/api/rescan` before fetching datasets and files
- **AND** newly appeared or removed datasets SHALL be reflected in the dropdown

#### Scenario: Active dataset preserved after refresh

- **WHEN** the user has selected "rifting-600x400-v2-h5" and clicks refresh
- **AND** that dataset directory still exists
- **THEN** the dropdown SHALL remain on "rifting-600x400-v2-h5" after the refresh completes

### Requirement: Overlay controls per panel

Each panel SHALL display a collapsible "Overlays" control row below the zoom controls row. The row SHALL contain a toggle checkbox for each overlay layer type, labelled with a short name: Ph (phases), T (temperature), Topo (topography), V (velocity), Dir (director), σ₁ (sigma1), ε̇₁ (edot1), Melt (melt). Clicking a toggle SHALL dispatch a `ctrl:panel:overlay-toggle` event with `panelId` and `layerType`. The row SHALL be collapsed by default and expand on click of an "Overlays ▸" disclosure button. Toggles for unavailable layers SHALL be disabled and visually greyed out.

#### Scenario: Overlay row collapsed by default

- **WHEN** a panel is created
- **THEN** the overlay control row SHALL be collapsed (only the "Overlays ▸" button visible)

#### Scenario: Expand overlay row

- **WHEN** the user clicks "Overlays ▸"
- **THEN** all layer toggles SHALL become visible
- **AND** the button label SHALL change to "Overlays ▾"

#### Scenario: Toggle dispatches event

- **WHEN** the user clicks the "T" (temperature) toggle
- **THEN** a `ctrl:panel:overlay-toggle` event SHALL bubble with `detail.panelId` and `detail.layerType === 'temperature'`

#### Scenario: Unavailable layer greyed out

- **WHEN** the director data is not available for the current file
- **THEN** the "Dir" toggle SHALL be disabled and styled with reduced opacity

### Requirement: Overlay layer settings popover

Each enabled overlay layer toggle SHALL have a small gear icon that opens a popover with layer-specific settings. Common settings: colour picker, line width slider. Temperature layer: ΔT number input (default 200, min 10, max 1000). Vector/director/tensor layers: density slider (range 2–20, default 6–8) and length scale slider. Melt layer: levels text input (comma-separated numbers). Changes SHALL dispatch `ctrl:panel:overlay-config` events with the updated configuration.

#### Scenario: Temperature ΔT input

- **WHEN** the user opens the temperature settings popover and changes ΔT to 100
- **THEN** a `ctrl:panel:overlay-config` event SHALL fire with `detail.config.dT === 100`

#### Scenario: Velocity density slider

- **WHEN** the user moves the velocity density slider to 4
- **THEN** a `ctrl:panel:overlay-config` event SHALL fire with `detail.config.density === 4`

#### Scenario: Colour picker change

- **WHEN** the user changes the phases overlay colour to #ff0000
- **THEN** a `ctrl:panel:overlay-config` event SHALL fire with `detail.config.color === '#ff0000'`

### Requirement: Overlay data API endpoint

The server SHALL expose `GET /api/overlay-data/:filename?layers=<comma-separated>` that returns a JSON object containing raw field data for the requested overlay layer types. The endpoint SHALL validate the filename against the `Output\d+\.gzip\.h5` pattern. The endpoint SHALL check dataset existence for each requested layer and only include layers whose data is present. The response SHALL be gzip-compressed when the client accepts it.

#### Scenario: Valid request returns data

- **WHEN** `GET /api/overlay-data/Output00010.gzip.h5?layers=phases,temperature` is requested
- **THEN** the server SHALL return 200 with JSON containing `phases` and `temperature` keys

#### Scenario: Invalid filename returns 400

- **WHEN** `GET /api/overlay-data/badfile.h5?layers=phases` is requested
- **THEN** the server SHALL return 400 with an error message

#### Scenario: Unknown layer key ignored

- **WHEN** `?layers=phases,unknown_layer` is requested
- **THEN** the response SHALL contain `phases` but not `unknown_layer`

#### Scenario: File not found returns 404

- **WHEN** the requested file does not exist on disk
- **THEN** the server SHALL return 404 with an error message
