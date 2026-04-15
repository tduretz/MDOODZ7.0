## ADDED Requirements

### Requirement: Range lock toggle

The UI SHALL provide a lock/unlock toggle button (🔒/🔓) next to the colour range controls. When locked, the current colour range SHALL persist when switching files or fields. When unlocked, the range SHALL auto-reset to the percentile-clipped range on each data load.

#### Scenario: Lock the range

- **WHEN** the user clicks the lock button while it shows 🔓
- **THEN** the button SHALL change to 🔒
- **AND** `model.rangeLocked` SHALL be `true`
- **AND** the current min/max range SHALL be preserved on subsequent file/field changes

#### Scenario: Unlock the range

- **WHEN** the user clicks the lock button while it shows 🔒
- **THEN** the button SHALL change to 🔓
- **AND** `model.rangeLocked` SHALL be `false`
- **AND** the range SHALL immediately recalculate to the current data's pMin/pMax

#### Scenario: Locked range persists across file switch

- **WHEN** the range is locked at [1e20, 1e22] and the user switches to a different file
- **THEN** the colour range SHALL remain [1e20, 1e22]

### Requirement: Per-field colour-map memory

The model SHALL maintain a `Map<string, string>` named `fieldColourMaps` that stores the colour-map name for each field. When the user changes the colour map, the model SHALL store the mapping for the current field. When switching fields, the controller SHALL restore the stored colour map if one exists.

#### Scenario: Remember colour map for field

- **WHEN** the user selects "turbo" for "Viscosity" and then switches to "Temperature"
- **THEN** switching back to "Viscosity" SHALL restore the "turbo" colour map

#### Scenario: No stored colour map

- **WHEN** the user switches to a field that has never had a colour map set
- **THEN** the current colour map SHALL remain unchanged

#### Scenario: Session only

- **WHEN** the browser is refreshed
- **THEN** all per-field colour-map associations SHALL be lost (not persisted)

### Requirement: Header banner with field info and time

The UI SHALL display a header banner with: left = "MDOODZ WebVis", centre = layout selector (1×1, 1×2, 2×1, 2×2), right = adaptive time display with a unit selector dropdown (`Auto`/`yr`/`ka`/`Ma`). The centre section SHALL replace the single-field label with the layout selector, since field information is now shown per-panel in each panel's title and controls.

#### Scenario: Header shows layout selector

- **WHEN** the application is running with multiple panels
- **THEN** the header centre SHALL display the layout selector (not a single field label)

#### Scenario: Time display in header

- **WHEN** a file is loaded with time=3.1558e14 s
- **THEN** the header right SHALL show `"10.00 Ma"` with the time unit dropdown

### Requirement: Scientific field labels in dropdowns

The field selector dropdown SHALL display the scientific `label` from the field registry instead of the plain-text field name. The underlying value SHALL remain the field key for API calls.

#### Scenario: Dropdown shows scientific labels

- **WHEN** the field dropdown is populated
- **THEN** the option for "Viscosity" SHALL display `"η"` and the option for "Stress II" SHALL display `"τ_II"`

#### Scenario: Field key preserved for API calls

- **WHEN** the user selects the option showing `"η"`
- **THEN** the controller SHALL use `"Viscosity"` as the field name in API requests

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

### Requirement: Colour range controls

The UI SHALL provide min/max input fields (or a dual-range slider) allowing the user to adjust the colour mapping range. The initial range SHALL be set to the data's percentile-clipped pMin/pMax (instead of full min/max). Adjusting the range SHALL re-render the canvas without re-fetching data. An "Auto" button SHALL reset to pMin/pMax of the current data. Colour range controls SHALL be located within each panel's mini control row.

#### Scenario: Initial range uses percentile clipping

- **WHEN** a field is loaded with min=0, max=1e20, pMin=1e15, pMax=5e19
- **THEN** the colour range inputs SHALL be initialised to [1e15, 5e19]

#### Scenario: Narrow the colour range

- **WHEN** the user sets min=1e20 and max=1e22 for viscosity
- **THEN** values below 1e20 SHALL render as the minimum colour and values above 1e22 as the maximum colour

#### Scenario: Reset to auto range

- **WHEN** the user clicks "Auto" (or similar reset)
- **THEN** the range SHALL revert to the data's pMin and pMax

### Requirement: Status bar showing metadata

The UI SHALL display a status bar showing: current model time (formatted with `formatTime`), grid dimensions (Nx × Nz), and the number of active panels. The single-field name display SHALL be removed since each panel now shows its own field label.

#### Scenario: Display metadata after file load

- **WHEN** a file is loaded with time=3.1558e13 s, Nx=601, Nz=401 and 2 panels are active
- **THEN** the status bar SHALL show `"t = 1.00 Ma | 601 × 401 | 2 panels"` (or similar)
