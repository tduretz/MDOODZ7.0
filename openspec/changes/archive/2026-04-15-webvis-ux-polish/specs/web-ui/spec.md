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

The UI SHALL replace the simple `<h1>` header with a three-part banner: left = "MDOODZ WebVis", centre = scientific field label + formatted unit (from field-labels), right = adaptive time display with a unit selector dropdown (`Auto`/`yr`/`ka`/`Ma`).

#### Scenario: Header shows current field and time

- **WHEN** "Viscosity" is loaded from a file with time=3.1558e14 s
- **THEN** the centre SHALL show `"η [Pa·s]"` and the right SHALL show `"10.00 Ma"`

#### Scenario: Time unit selector

- **WHEN** the user selects "ka" from the time unit dropdown
- **THEN** the time display SHALL update to show ka units
- **AND** all file-selector labels SHALL also update

### Requirement: Scientific field labels in dropdowns

The field selector dropdown SHALL display the scientific `label` from the field registry instead of the plain-text field name. The underlying value SHALL remain the field key for API calls.

#### Scenario: Dropdown shows scientific labels

- **WHEN** the field dropdown is populated
- **THEN** the option for "Viscosity" SHALL display `"η"` and the option for "Stress II" SHALL display `"τ_II"`

#### Scenario: Field key preserved for API calls

- **WHEN** the user selects the option showing `"η"`
- **THEN** the controller SHALL use `"Viscosity"` as the field name in API requests

## MODIFIED Requirements

### Requirement: File selector control

The UI SHALL display a dropdown or list of available HDF5 files fetched from `GET /api/files`. Each option SHALL display `"Step N — X.XX Ma"` (or ka/yr, per the current time unit setting) instead of the raw filename. The underlying value SHALL be the filename. Selecting a file SHALL trigger loading that file's available fields and metadata.

#### Scenario: File list populated on startup

- **WHEN** the application loads
- **THEN** the file selector SHALL be populated with labels like `"Step 0 — 0.00 Ma"`, `"Step 10 — 1.00 Ma"`, etc.

#### Scenario: Select a file

- **WHEN** the user selects `"Step 50 — 5.00 Ma"` from the file selector
- **THEN** the field dropdown SHALL update with the fields available in that file
- **AND** the header banner SHALL show the model time for that file

### Requirement: Time-step slider

The UI SHALL provide a range slider that maps to the available output files (by index). Dragging the slider SHALL update the file selector and load the corresponding file's data for the currently selected field. The slider label SHALL display the step and time for the current position.

#### Scenario: Slide to a different time step

- **WHEN** the user drags the slider to position corresponding to `Output00100.gzip.h5`
- **THEN** the file selector SHALL update to the matching step entry
- **AND** the canvas SHALL re-render with data from that file

#### Scenario: Slider range matches file count

- **WHEN** there are 29 output files
- **THEN** the slider SHALL have min=0, max=28, step=1

### Requirement: Colour range controls

The UI SHALL provide min/max input fields (or a dual-range slider) allowing the user to adjust the colour mapping range. The initial range SHALL be set to the data's percentile-clipped pMin/pMax (instead of full min/max). Adjusting the range SHALL re-render the canvas without re-fetching data. An "Auto" button SHALL reset to pMin/pMax of the current data.

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

The UI SHALL display a status bar showing: current model time (formatted with `formatTime`), grid dimensions (Nx × Nz), and the currently loaded field name (scientific label). The time SHALL use the adaptive formatting from the time-display module.

#### Scenario: Display metadata after file load

- **WHEN** a file is loaded with time=3.1558e13 s, Nx=601, Nz=401
- **THEN** the status bar SHALL show `"t = 1.00 Ma | 601 × 401"` (or similar, using adaptive units)
