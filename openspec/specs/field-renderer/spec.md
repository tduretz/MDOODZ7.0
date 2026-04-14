## ADDED Requirements

### Requirement: Render 2D field data to canvas as colour-mapped image

The renderer SHALL take a 2D array of values, 1D x-coordinate and z-coordinate arrays, and a colour map LUT, and draw the corresponding colour-mapped image onto an HTML `<canvas>` element using `ImageData` pixel writing.

#### Scenario: Render a centre field

- **WHEN** the renderer receives a 600×400 float array with coordinate arrays
- **THEN** it SHALL draw a colour-mapped image onto the canvas where each data cell maps to one or more pixels

#### Scenario: Canvas matches data aspect ratio

- **WHEN** the domain width and height differ (e.g., 300 km wide × 200 km tall)
- **THEN** the canvas SHALL preserve the physical aspect ratio (width/height)

### Requirement: Map values to colours using a 256-entry LUT

The renderer SHALL normalise data values to the range [0, 1] using the provided min/max bounds, multiply by 255, clamp to [0, 255], and use the integer index to look up an RGB triplet from the active colour map LUT.

#### Scenario: Value at minimum maps to first colour

- **WHEN** a cell value equals the minimum
- **THEN** it SHALL render with the colour at LUT index 0

#### Scenario: Value at maximum maps to last colour

- **WHEN** a cell value equals the maximum
- **THEN** it SHALL render with the colour at LUT index 255

#### Scenario: Value outside range is clamped

- **WHEN** a cell value exceeds the maximum
- **THEN** it SHALL render with the colour at LUT index 255 (clamped)

### Requirement: Support logarithmic scaling

When the field's `log` property is true, the renderer SHALL apply `log10()` to each value before normalisation. Values ≤ 0 SHALL be treated as NaN (transparent).

#### Scenario: Log-scale viscosity

- **WHEN** the field is "Viscosity" with `log: true` and values range from 1e18 to 1e24
- **THEN** the colour normalisation SHALL operate on log₁₀ values (18 to 24)

#### Scenario: Zero or negative values in log mode

- **WHEN** a cell value is ≤ 0 and log mode is active
- **THEN** the cell SHALL render as transparent (NaN treatment)

### Requirement: Render NaN cells as transparent

Cells with NaN values (air cells, invalid data) SHALL be rendered with alpha = 0 (fully transparent), allowing a background colour or pattern to show through.

#### Scenario: Air cells appear transparent

- **WHEN** the data contains NaN values from air-cell masking
- **THEN** those pixels SHALL have RGBA alpha = 0

### Requirement: Support discrete colour maps for phase fields

For fields marked as `discrete: true`, the renderer SHALL use a discrete colour palette where each integer phase value maps to a specific colour. Phase -1 (air) SHALL be transparent.

#### Scenario: Phase map rendering

- **WHEN** the field is "Phases" with integer values 0, 1, 2, 3, -1
- **THEN** each phase SHALL render with its assigned palette colour and phase -1 SHALL be transparent

### Requirement: Provide multiple scientific colour maps

The renderer SHALL support at least 5 colour maps: viridis, turbo, inferno, plasma, and coolwarm. Each SHALL be a 256-entry array of [R, G, B] triplets (0–255). The colour maps SHALL be loaded from a bundled JSON file.

#### Scenario: Switch colour map

- **WHEN** the active colour map is changed from viridis to turbo
- **THEN** the canvas SHALL re-render with the turbo colour map applied to the same data

### Requirement: Draw colour bar

The renderer SHALL draw a vertical colour bar alongside the field image showing the mapping from value to colour. The bar SHALL display the min and max values, the field unit, and 3–5 intermediate tick labels. For log-scale fields, tick labels SHALL show the exponents (e.g., 10¹⁸).

#### Scenario: Colour bar for pressure field

- **WHEN** the field "Pressure" is rendered with min=0, max=1e9, unit="Pa"
- **THEN** a colour bar SHALL appear with tick labels in Pa and the colour gradient matching the field image

#### Scenario: Colour bar for log-scale field

- **WHEN** the field "Viscosity" is rendered with log=true
- **THEN** tick labels SHALL display as powers of 10

### Requirement: Draw axis labels with physical coordinates

The renderer SHALL draw x and z axis tick labels using the physical coordinate arrays (in metres or kilometres). The axes SHALL show at least 3 ticks on each axis.

#### Scenario: Axis labels in kilometres

- **WHEN** the domain width is 300,000 m
- **THEN** x-axis labels SHALL be displayed in km (e.g., -150, -100, -50, 0, 50, 100, 150)

### Requirement: Render performance

The renderer SHALL render a 600×400 field in under 50 ms on a modern browser, measured from receiving the data array to completing the canvas draw.

#### Scenario: Benchmark render time

- **WHEN** a 600×400 float array is rendered
- **THEN** the time from `ImageData` creation to `putImageData` completion SHALL be under 50 ms
