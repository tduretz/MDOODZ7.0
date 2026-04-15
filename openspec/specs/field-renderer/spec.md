## MODIFIED Requirements

### Requirement: Map values to colours using a 256-entry LUT

The renderer SHALL normalise data values to the range [0, 1] using the provided min/max bounds, multiply by 255, clamp to [0, 255], and use the integer index to look up an RGB triplet from the active colour map LUT. When the range is auto (unlocked), the bounds SHALL be pMin/pMax. When the range is locked or manually set, the bounds SHALL be the user-specified values. The renderer SHALL accept a target canvas element as a parameter rather than assuming a single global canvas, enabling per-panel rendering. After drawing the field data via ImageData, the renderer SHALL draw the panel's interpolated title text as an overlay at the top of the canvas with a semi-transparent dark background bar. The title SHALL be rendered using `ctx.fillText()` so it is included in any canvas-to-blob or screenshot export. The renderer SHALL use `devicePixelRatio` scaling for crisp text at non-integer DPI.

#### Scenario: Value at minimum maps to first colour

- **WHEN** a cell value equals the minimum
- **THEN** it SHALL render with the colour at LUT index 0

#### Scenario: Value at maximum maps to last colour

- **WHEN** a cell value equals the maximum
- **THEN** it SHALL render with the colour at LUT index 255

#### Scenario: Value outside range is clamped

- **WHEN** a cell value exceeds the maximum
- **THEN** it SHALL render with the colour at LUT index 255 (clamped)

#### Scenario: Auto-range uses percentile bounds

- **WHEN** the range is not locked and not manually set
- **THEN** the normalisation bounds SHALL be pMin and pMax from the field data

#### Scenario: Render to specific canvas element

- **WHEN** the renderer is called with a target canvas element for panel 2
- **THEN** the field image SHALL be drawn on that specific canvas, not on any other panel's canvas

#### Scenario: Title rendered inside canvas

- **WHEN** a panel has a title template and field data is loaded
- **THEN** the interpolated title text SHALL be drawn at the top of the canvas over a semi-transparent dark background bar
- **AND** the title SHALL be included in any `canvas.toBlob()` or `canvas.toDataURL()` export

#### Scenario: Canvas render triggered by data arrival

- **WHEN** field data arrives for a panel (via `panel:field-changed` event)
- **THEN** the canvas SHALL re-render immediately, even if the ResizeObserver has not yet fired

### Requirement: Draw colour bar

The renderer SHALL draw a vertical colour bar alongside the field image showing the mapping from value to colour. The bar SHALL display the min and max values, the scientific `formattedUnit` string (from field-labels), and 3–5 intermediate tick labels. For log-scale fields, tick labels SHALL show the exponents (e.g., 10¹⁸). The colour bar title SHALL use the field's scientific `label`. Each panel SHALL have its own independent colour bar instance bound to its PanelState's colour range and colour map. The colour bar view SHALL provide a `destroy()` method that removes all event listeners it registered on the model.

#### Scenario: Colour bar for pressure field

- **WHEN** the field "Pressure" is rendered in a panel with min=0, max=1e9
- **THEN** that panel's colour bar SHALL appear with tick labels in Pa, the title `"P"`, and the unit `"Pa"` rendered using the `formattedUnit` string

#### Scenario: Colour bar for log-scale viscosity

- **WHEN** the field "Viscosity" is rendered in a panel with log=true
- **THEN** that panel's colour bar tick labels SHALL display as powers of 10, the title SHALL be `"η"`, and the unit SHALL be `"Pa·s"`

#### Scenario: Colour bar with locked range indicator

- **WHEN** a panel's range is locked (rangeLocked === true)
- **THEN** that panel's colour bar SHALL display a visual indicator (e.g., 🔒 icon or border) showing the range is locked

#### Scenario: Independent colour bars across panels

- **WHEN** panel 1 shows Viscosity [1e18, 1e24] and panel 2 shows Pressure [0, 1e9]
- **THEN** each panel's colour bar SHALL display its own independent range and colour map

#### Scenario: Destroy removes event listeners

- **WHEN** a ColourBar's `destroy()` method is called
- **THEN** all event listeners registered on the model SHALL be removed
- **AND** no further model events SHALL trigger re-renders on this instance

### Requirement: Phase box legend for discrete fields

When the active field is discrete (e.g. Phase), the renderer SHALL draw a box legend instead of a continuous colour bar. Each visible phase SHALL be represented by a coloured square and a label. Phase names and colours SHALL be customizable via a popover editor with colour pickers and text inputs. Default phase names SHALL be "Phase 0", "Phase 1", etc. Custom names and colours SHALL be stored in the model and persist across file/time-step changes. Negative phase values (air) SHALL be rendered with the null-pixel colour.

#### Scenario: Phase legend shows only visible phases

- **WHEN** the cropped view contains phases 0, 1, 3
- **THEN** the legend SHALL show exactly three entries (phases 0, 1, 3)

#### Scenario: Custom phase name

- **WHEN** the user renames phase 0 to "Crust" via the popover
- **THEN** the legend SHALL display "Crust" next to the phase-0 colour square

#### Scenario: Custom phase colour

- **WHEN** the user changes phase 1's colour to red via the colour picker
- **THEN** the field image and legend SHALL both use the new red colour for phase 1

#### Scenario: Discrete field hides continuous controls

- **WHEN** a discrete field is selected
- **THEN** the colour-map selector, range inputs, Auto, and Lock controls SHALL be hidden
- **AND** the Phases config button SHALL be shown

### Requirement: HiDPI rendering with dynamic margins

The canvas SHALL render at `devicePixelRatio` resolution for crisp output on Retina/HiDPI displays. Margins (top, right, bottom, left) SHALL be computed dynamically based on actual text measurements of axis tick labels, colour-bar labels, and axis titles. The data area SHALL be aspect-ratio-fitted into the remaining space after margins. A black frame SHALL be drawn around the data area.

#### Scenario: 2x DPI display

- **WHEN** `devicePixelRatio` is 2
- **THEN** the canvas backing store SHALL be 2× the CSS pixel dimensions
- **AND** all drawing SHALL be scaled by 2×

#### Scenario: Dynamic left margin

- **WHEN** Z-axis tick labels are wide (e.g. "-109.9")
- **THEN** the left margin SHALL expand to accommodate the widest label plus the rotated "z [km]" axis title

### Requirement: Axis ticks and labels

The renderer SHALL draw bottom (X) and left (Z) axes with tick marks and numeric labels. For log-scale fields, tick labels SHALL use Unicode superscripts (e.g. 10⁻¹⁸). Axis titles SHALL show the coordinate name and current spatial unit (e.g. "x [km]"). The Z-axis title SHALL be rotated 90° counter-clockwise. Tick count SHALL adapt to the available space (3–6 ticks depending on axis length). When zoom is active, ticks SHALL reflect the cropped coordinate range.

#### Scenario: Axis labels update with spatial unit

- **WHEN** the user switches from km to m
- **THEN** axis labels SHALL show "x [m]" and "z [m]"
- **AND** tick values SHALL be in metres

#### Scenario: Log-scale tick superscripts

- **WHEN** a log-scale field has tick value 1e-18
- **THEN** the tick label SHALL render as "10⁻¹⁸" using Unicode superscript characters
