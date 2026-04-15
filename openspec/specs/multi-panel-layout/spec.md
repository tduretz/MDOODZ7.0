## Purpose

Multi-panel layout system for WebVis — CSS Grid presets, panel CRUD, responsive sizing, per-panel DOM subtrees, and panel state management in the Model.

## Requirements

### Requirement: Grid container with layout presets

The application SHALL provide a `#panels-grid` container using CSS Grid that arranges panels in one of four presets: 1×1 (single), 1×2 (side-by-side), 2×1 (stacked), 2×2 (quad). The active preset SHALL be controlled by a layout selector in the header. The default layout SHALL be 1×1, preserving the existing single-panel experience.

#### Scenario: Default layout is 1×1

- **WHEN** the application loads for the first time
- **THEN** the grid container SHALL use the `.grid-1x1` class
- **AND** exactly one panel SHALL be visible

#### Scenario: Switch to 1×2 layout

- **WHEN** the user selects the 1×2 layout from the layout selector
- **THEN** the grid container SHALL switch to `.grid-1x2`
- **AND** two panels SHALL be displayed side-by-side (two columns, one row)

#### Scenario: Switch to 2×1 layout

- **WHEN** the user selects the 2×1 layout from the layout selector
- **THEN** the grid container SHALL switch to `.grid-2x1`
- **AND** two panels SHALL be displayed stacked vertically (one column, two rows)

#### Scenario: Switch to 2×2 layout

- **WHEN** the user selects the 2×2 layout from the layout selector
- **THEN** the grid container SHALL switch to `.grid-2x2`
- **AND** four panels SHALL be displayed in a 2×2 grid

### Requirement: Panel add and remove

When switching to a layout with more panels than currently exist, the application SHALL automatically create new panels initialised with default state (no field selected, default colour map, empty title template). When switching to a layout with fewer panels, the application SHALL remove excess panels starting from the highest panel index. The minimum number of panels SHALL be 1.

#### Scenario: Grow from 1 to 2 panels

- **WHEN** the user switches from 1×1 to 1×2 layout
- **THEN** a second panel SHALL be created with default state
- **AND** the first panel SHALL retain its current field, colour map, colour range, and title

#### Scenario: Shrink from 4 to 1 panel

- **WHEN** the user switches from 2×2 to 1×1 layout
- **THEN** panels 2, 3, and 4 SHALL be removed
- **AND** panel 1 SHALL retain all its state

#### Scenario: Cannot remove last panel

- **WHEN** exactly one panel exists
- **THEN** no panel removal action SHALL be available

### Requirement: Responsive panel sizing

Each panel SHALL fill its grid cell completely. When the browser window is resized, all panel canvases and colour bars SHALL resize proportionally to fit their grid cells. The canvas aspect ratio SHALL match the data grid dimensions (Nx × Nz).

#### Scenario: Window resize updates all panels

- **WHEN** the browser window is resized while 4 panels are displayed
- **THEN** all four canvases SHALL resize to fit their new grid cell dimensions

#### Scenario: Canvas aspect ratio matches data

- **WHEN** a panel displays a 600×400 field
- **THEN** the canvas aspect ratio SHALL be 3:2, letterboxing within the grid cell if necessary

### Requirement: Per-panel DOM subtree

Each panel SHALL be a self-contained DOM subtree containing: a title display area, a canvas element for the field visualisation, a colour bar, and a mini control row (field selector, colour-map selector, lock button). Each panel's views SHALL be bound to its own PanelState and SHALL NOT affect other panels.

#### Scenario: Panel contains all required elements

- **WHEN** a panel is created
- **THEN** it SHALL contain a title area, a `<canvas>` element, a colour bar, and control widgets for field selection, colour map, and range lock

#### Scenario: Panel independence

- **WHEN** the user changes the field in panel 1
- **THEN** panel 2 SHALL remain unchanged (same field, colour map, colour range, title)

### Requirement: Panel state in Model

The Model SHALL maintain a `panels` array of PanelState objects, each containing `{ id, fieldName, colourMap, colourRange, rangeLocked, fieldData, title }`. The Model SHALL also track an `activePanelId`. Panel-scoped events SHALL carry `{ panelId }` in their detail so views can filter events for their own panel.

#### Scenario: Initial panels array

- **WHEN** the application starts with 1×1 layout
- **THEN** `model.panels` SHALL contain exactly one PanelState object with a unique `id`

#### Scenario: Panel event carries panelId

- **WHEN** a field is loaded for panel 2
- **THEN** the dispatched event SHALL include `detail.panelId` matching the id of panel 2

#### Scenario: Active panel tracking

- **WHEN** the user clicks on panel 2
- **THEN** `model.activePanelId` SHALL update to panel 2's id
