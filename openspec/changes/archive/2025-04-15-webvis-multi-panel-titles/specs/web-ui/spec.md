## MODIFIED Requirements

### Requirement: File selector control

The UI SHALL display a dropdown or list of available HDF5 files fetched from `GET /api/files`. Each option SHALL display `"Step N — X.XX Ma"` (or ka/yr, per the current time unit setting) instead of the raw filename. The underlying value SHALL be the filename. Selecting a file SHALL trigger loading that file's available fields and metadata. The file selector SHALL remain a global control — selecting a file SHALL trigger each panel to re-fetch its own field data independently.

#### Scenario: File list populated on startup

- **WHEN** the application loads
- **THEN** the file selector SHALL be populated with labels like `"Step 0 — 0.00 Ma"`, `"Step 10 — 1.00 Ma"`, etc.

#### Scenario: Select a file

- **WHEN** the user selects `"Step 50 — 5.00 Ma"` from the file selector
- **THEN** each panel SHALL re-fetch its own field data for the new file
- **AND** all panel titles SHALL update with the new file's metadata

### Requirement: Time-step slider

The UI SHALL provide a range slider that maps to the available output files (by index). Dragging the slider SHALL update the file selector and trigger all panels to load the corresponding file's data for their respective fields. The slider label SHALL display the step and time for the current position.

#### Scenario: Slide to a different time step

- **WHEN** the user drags the slider to position corresponding to `Output00100.gzip.h5`
- **THEN** all panels SHALL re-fetch their respective field data for that file
- **AND** all panel titles SHALL update

#### Scenario: Slider range matches file count

- **WHEN** there are 29 output files
- **THEN** the slider SHALL have min=0, max=28, step=1

### Requirement: Header banner with field info and time

The UI SHALL display a header banner with: left = "MDOODZ WebVis", centre = layout selector (1×1, 1×2, 2×1, 2×2), right = adaptive time display with a unit selector dropdown (`Auto`/`yr`/`ka`/`Ma`). The centre section SHALL replace the single-field label with the layout selector, since field information is now shown per-panel in each panel's title and controls.

#### Scenario: Header shows layout selector

- **WHEN** the application is running with multiple panels
- **THEN** the header centre SHALL display the layout selector (not a single field label)

#### Scenario: Time display in header

- **WHEN** a file is loaded with time=3.1558e14 s
- **THEN** the header right SHALL show `"10.00 Ma"` with the time unit dropdown

### Requirement: Status bar showing metadata

The UI SHALL display a status bar showing: current model time (formatted with `formatTime`), grid dimensions (Nx × Nz), and the number of active panels. The single-field name display SHALL be removed since each panel now shows its own field label.

#### Scenario: Display metadata after file load

- **WHEN** a file is loaded with time=3.1558e13 s, Nx=601, Nz=401 and 2 panels are active
- **THEN** the status bar SHALL show `"t = 1.00 Ma | 601 × 401 | 2 panels"` (or similar)

## ADDED Requirements

### Requirement: Layout selector in header

The header banner SHALL contain a layout selector control (buttons or dropdown) allowing the user to choose between 1×1, 1×2, 2×1, and 2×2 grid presets. Changing the layout SHALL add or remove panels as needed and update the grid container CSS class.

#### Scenario: Select 2×2 layout

- **WHEN** the user clicks the 2×2 layout button
- **THEN** the grid SHALL switch to 2×2
- **AND** panels SHALL be created or removed to match (4 total)

#### Scenario: Layout selector shows current layout

- **WHEN** the current layout is 1×2
- **THEN** the 1×2 option SHALL be visually highlighted as active

### Requirement: Global controls remain shared

The file selector, time slider, and time unit dropdown SHALL remain global controls that apply to all panels. Per-panel controls (field selector, colour-map selector, range lock) SHALL be located within each panel's DOM subtree.

#### Scenario: File change affects all panels

- **WHEN** the user selects a new file from the global file selector
- **THEN** all panels SHALL load their respective fields from the new file

#### Scenario: Panel controls are independent

- **WHEN** the user changes the colour map in panel 1 to "turbo"
- **THEN** panel 2's colour map SHALL remain unchanged
