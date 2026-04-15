## MODIFIED Requirements

### Requirement: Panel add and remove

When switching to a layout with more panels than currently exist, the application SHALL first restore panels from the stash (most-recently-stashed first) before creating new defaults. When switching to a layout with fewer panels, the application SHALL move excess panels to a stash array (preserving their full PanelState including fieldData, colourMap, colourRange, title, and rangeLocked) instead of destroying them. The stash SHALL enable round-trip layout changes (e.g. 2×2 → 1×1 → 2×2) without losing panel state. The minimum number of active panels SHALL be 1.

#### Scenario: Grow from 1 to 2 panels — restore from stash

- **WHEN** the user switches from 1×1 to 1×2 layout and a stashed panel exists
- **THEN** the stashed panel SHALL be restored with its previous field, colour map, colour range, and title
- **AND** the first panel SHALL retain its current state

#### Scenario: Grow from 1 to 2 panels — no stash

- **WHEN** the user switches from 1×1 to 1×2 layout and no stashed panels exist
- **THEN** a new panel SHALL be created with default state

#### Scenario: Shrink from 4 to 1 panel — stash preserved

- **WHEN** the user switches from 2×2 to 1×1 layout
- **THEN** panels 2, 3, and 4 SHALL be moved to the stash (not destroyed)
- **AND** panel 1 SHALL retain all its state
- **AND** expanding back to 2×2 SHALL restore panels 2, 3, and 4 from stash

#### Scenario: Cannot remove last panel

- **WHEN** exactly one panel exists
- **THEN** no panel removal action SHALL be available

### Requirement: Responsive panel sizing

Each panel SHALL fill its grid cell completely. When the browser window is resized, all panel canvases and colour bars SHALL resize proportionally to fit their grid cells. The canvas aspect ratio SHALL match the data grid dimensions (Nx × Nz). The available width SHALL be computed from the panel element's actual bounding rect minus the colour bar's actual rendered width, not from a hardcoded offset. The sizing logic SHALL handle all grid presets (1×1, 1×2, 2×1, 2×2) correctly without distortion.

#### Scenario: Window resize updates all panels

- **WHEN** the browser window is resized while 4 panels are displayed
- **THEN** all four canvases SHALL resize to fit their new grid cell dimensions

#### Scenario: Canvas aspect ratio matches data

- **WHEN** a panel displays a 600×400 field
- **THEN** the canvas aspect ratio SHALL be 3:2, letterboxing within the grid cell if necessary

#### Scenario: Narrow panels in 1×2 layout

- **WHEN** the browser window is narrow and two panels are displayed side-by-side (1×2)
- **THEN** each canvas SHALL shrink proportionally without distortion or overflow

#### Scenario: Tall panels in 2×1 layout

- **WHEN** the browser window is tall and two panels are stacked vertically (2×1)
- **THEN** each canvas SHALL fill its cell height correctly without distortion

### Requirement: Per-panel DOM subtree

Each panel SHALL be a self-contained DOM subtree containing: a title display area, a canvas element for the field visualisation, a colour bar, and a mini control row (field selector, colour-map selector, lock button). Each panel's views SHALL be bound to its own PanelState and SHALL NOT affect other panels. When a PanelView is destroyed, it SHALL call `destroy()` on all owned sub-views (FieldCanvas, ColourBar, TitleInput) to remove event listeners and prevent memory leaks.

#### Scenario: Panel contains all required elements

- **WHEN** a panel is created
- **THEN** it SHALL contain a title area, a `<canvas>` element, a colour bar, and control widgets for field selection, colour map, and range lock

#### Scenario: Panel independence

- **WHEN** the user changes the field in panel 1
- **THEN** panel 2 SHALL remain unchanged (same field, colour map, colour range, title)

#### Scenario: Panel destroy cleans up listeners

- **WHEN** a PanelView is destroyed (e.g. layout shrink)
- **THEN** all event listeners on the model registered by FieldCanvas, ColourBar, and TitleInput SHALL be removed
- **AND** the panel DOM subtree SHALL be removed from the document

### Requirement: Panel state in Model

The Model SHALL maintain a `panels` array of PanelState objects, each containing `{ id, fieldName, colourMap, colourRange, rangeLocked, fieldData, title }`. The Model SHALL also track an `activePanelId` and a `_stashedPanels` array for layout round-trips. Panel-scoped events SHALL carry `{ panelId }` in their detail so views can filter events for their own panel. When new field data is assigned to a panel, the previous `fieldData` reference SHALL be explicitly nulled before the new value is set, to allow garbage collection.

#### Scenario: Initial panels array

- **WHEN** the application starts with 1×1 layout
- **THEN** `model.panels` SHALL contain exactly one PanelState object with a unique `id`

#### Scenario: Panel event carries panelId

- **WHEN** a field is loaded for panel 2
- **THEN** the dispatched event SHALL include `detail.panelId` matching the id of panel 2

#### Scenario: Active panel tracking

- **WHEN** the user clicks on panel 2
- **THEN** `model.activePanelId` SHALL update to panel 2's id

#### Scenario: Previous fieldData nulled on update

- **WHEN** new field data is loaded for a panel that already has data
- **THEN** the previous `fieldData` reference SHALL be set to `null` before the new data is assigned
