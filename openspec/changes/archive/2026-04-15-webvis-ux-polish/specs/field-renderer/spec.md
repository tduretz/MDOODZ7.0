## MODIFIED Requirements

### Requirement: Draw colour bar

The renderer SHALL draw a vertical colour bar alongside the field image showing the mapping from value to colour. The bar SHALL display the min and max values, the scientific `formattedUnit` string (from field-labels), and 3–5 intermediate tick labels. For log-scale fields, tick labels SHALL show the exponents (e.g., 10¹⁸). The colour bar title SHALL use the field's scientific `label`.

#### Scenario: Colour bar for pressure field

- **WHEN** the field "Pressure" is rendered with min=0, max=1e9
- **THEN** a colour bar SHALL appear with tick labels in Pa, the title `"P"`, and the unit `"Pa"` rendered using the `formattedUnit` string

#### Scenario: Colour bar for log-scale viscosity

- **WHEN** the field "Viscosity" is rendered with log=true
- **THEN** tick labels SHALL display as powers of 10, the title SHALL be `"η"`, and the unit SHALL be `"Pa·s"`

#### Scenario: Colour bar with locked range indicator

- **WHEN** the range is locked (model.rangeLocked === true)
- **THEN** the colour bar SHALL display a visual indicator (e.g., 🔒 icon or border) showing the range is locked

### Requirement: Map values to colours using a 256-entry LUT

The renderer SHALL normalise data values to the range [0, 1] using the provided min/max bounds, multiply by 255, clamp to [0, 255], and use the integer index to look up an RGB triplet from the active colour map LUT. When the range is auto (unlocked), the bounds SHALL be pMin/pMax. When the range is locked or manually set, the bounds SHALL be the user-specified values.

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
