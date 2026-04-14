## ADDED Requirements

### Requirement: MVC architecture with observer pattern

The frontend SHALL use a Model-View-Controller architecture. The Model SHALL extend `EventTarget` and dispatch `CustomEvent` instances when state changes. Views SHALL subscribe to model events and re-render in response. The Controller SHALL handle user interactions, call APIs, and update the Model. Views SHALL NOT call APIs directly.

#### Scenario: Model emits event on field data load

- **WHEN** the controller sets new field data on the model
- **THEN** the model SHALL dispatch a `field-loaded` event
- **AND** all subscribed views SHALL re-render with the new data

#### Scenario: Views do not call APIs

- **WHEN** a view needs new data (e.g., user clicks a field)
- **THEN** the view SHALL emit a DOM event that the controller handles
- **AND** the controller SHALL fetch data and update the model

### Requirement: File selector control

The UI SHALL display a dropdown or list of available HDF5 files fetched from `GET /api/files`. Selecting a file SHALL trigger loading that file's available fields and metadata.

#### Scenario: File list populated on startup

- **WHEN** the application loads
- **THEN** the file selector SHALL be populated with all available output files

#### Scenario: Select a file

- **WHEN** the user selects `Output00050.gzip.h5` from the file selector
- **THEN** the field dropdown SHALL update with the fields available in that file
- **AND** the status bar SHALL show the model time for that file

### Requirement: Time-step slider

The UI SHALL provide a range slider that maps to the available output files (by index). Dragging the slider SHALL update the file selector and load the corresponding file's data for the currently selected field.

#### Scenario: Slide to a different time step

- **WHEN** the user drags the slider to position corresponding to `Output00100.gzip.h5`
- **THEN** the file selector SHALL update to `Output00100.gzip.h5`
- **AND** the canvas SHALL re-render with data from that file

#### Scenario: Slider range matches file count

- **WHEN** there are 29 output files
- **THEN** the slider SHALL have min=0, max=28, step=1

### Requirement: Field selector dropdown

The UI SHALL display a dropdown listing all available field names for the current file (as returned by `GET /api/fields/:file`). Selecting a field SHALL trigger loading that field's data and re-rendering the canvas.

#### Scenario: Select a field

- **WHEN** the user selects "Temperature" from the field dropdown
- **THEN** the controller SHALL fetch `/api/field-data/:file/Temperature`
- **AND** the canvas SHALL re-render with the temperature data

#### Scenario: Preserve field selection across files

- **WHEN** the user has "Pressure" selected and switches to a different file
- **THEN** "Pressure" SHALL remain selected (if available in the new file)
- **AND** the canvas SHALL show pressure data from the new file

### Requirement: Colour map selector

The UI SHALL provide a dropdown to select the active colour map from the available maps (viridis, turbo, inferno, plasma, coolwarm). Changing the colour map SHALL re-render the canvas without re-fetching data.

#### Scenario: Change colour map

- **WHEN** the user switches from viridis to turbo
- **THEN** the canvas SHALL re-render using the turbo LUT with the existing data (no API call)

### Requirement: Colour range controls

The UI SHALL provide min/max input fields (or a dual-range slider) allowing the user to adjust the colour mapping range. The initial range SHALL be set to the data's actual min/max. Adjusting the range SHALL re-render the canvas without re-fetching data.

#### Scenario: Narrow the colour range

- **WHEN** the user sets min=1e20 and max=1e22 for viscosity
- **THEN** values below 1e20 SHALL render as the minimum colour and values above 1e22 as the maximum colour

#### Scenario: Reset to auto range

- **WHEN** the user clicks "Auto" (or similar reset)
- **THEN** the range SHALL revert to the data's actual min and max

### Requirement: Status bar showing metadata

The UI SHALL display a status bar showing: current model time (in Myr), grid dimensions (Nx × Nz), and the currently loaded field name. The time SHALL be computed as `Params[0] / 3.1558e13` (seconds to Myr).

#### Scenario: Display metadata after file load

- **WHEN** a file is loaded with time=3.1558e13 s, Nx=601, Nz=401
- **THEN** the status bar SHALL show "t = 1.000 Myr | 601 × 401" (or similar)

### Requirement: Loading indicator

The UI SHALL show a loading indicator (e.g., spinner or progress text) while an API request is in progress. The indicator SHALL be hidden when the response is received and rendering is complete.

#### Scenario: Show loading during fetch

- **WHEN** the controller initiates a field-data fetch
- **THEN** a loading indicator SHALL be visible
- **AND** it SHALL be hidden after the canvas finishes rendering

### Requirement: No build step, no framework

The frontend SHALL consist of plain ES6+ modules (`.mjs` files) served directly by the Node.js server. There SHALL be no bundler, no transpiler, and no framework dependency. The HTML file SHALL use `<script type="module">` to load the entry point.

#### Scenario: Browser loads frontend

- **WHEN** a browser navigates to `http://localhost:3000/`
- **THEN** the page SHALL load and render without any build or compilation step

### Requirement: Responsive layout

The UI SHALL use CSS Grid or Flexbox to arrange the control panel, canvas, and colour bar. The canvas SHALL resize to fill available space while preserving the data's physical aspect ratio. The layout SHALL work in viewports from 1024px to 2560px width.

#### Scenario: Resize browser window

- **WHEN** the browser window is resized
- **THEN** the canvas SHALL resize proportionally while maintaining the domain aspect ratio

### Requirement: Clean visual design

The UI SHALL use a dark background with light text (suitable for scientific visualisation). Controls SHALL be visually grouped. Font choices SHALL be legible at standard sizes. The overall appearance SHALL be professional and uncluttered.

#### Scenario: Visual appearance

- **WHEN** the application loads
- **THEN** the background SHALL be dark (e.g., #1a1a2e or similar)
- **AND** text and controls SHALL have sufficient contrast for readability
