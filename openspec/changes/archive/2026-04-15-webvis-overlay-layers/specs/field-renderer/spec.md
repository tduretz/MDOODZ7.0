## MODIFIED Requirements

### Requirement: Map values to colours using a 256-entry LUT

The renderer SHALL normalise data values to the range [0, 1] using the provided min/max bounds, multiply by 255, clamp to [0, 255], and use the integer index to look up an RGB triplet from the active colour map LUT. When the range is auto (unlocked), the bounds SHALL be pMin/pMax. When the range is locked or manually set, the bounds SHALL be the user-specified values. The renderer SHALL accept a target canvas element as a parameter rather than assuming a single global canvas, enabling per-panel rendering. After drawing the field data via ImageData, the renderer SHALL draw the panel's interpolated title text as an overlay at the top of the canvas with a semi-transparent dark background bar. The title SHALL be rendered using `ctx.fillText()` so it is included in any canvas-to-blob or screenshot export. The renderer SHALL use `devicePixelRatio` scaling for crisp text at non-integer DPI. After drawing the field image and black frame, and before drawing the colour bar, axes, and title, the renderer SHALL call the overlay drawing phase. The overlay phase SHALL draw all enabled overlay layers using the computed data-area coordinates (`dataX`, `dataY`, `dW`, `dH`) and crop indices (`col0`, `col1`, `row0`, `row1`). The overlay drawing region SHALL be clipped to the data-area rectangle.

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

#### Scenario: Overlay phase executes between frame and colour bar

- **WHEN** the render method is called and at least one overlay layer is enabled
- **THEN** the overlay drawing SHALL occur after the black frame `strokeRect` and before `_drawCbar` / `_drawPhaseLegend`
- **AND** the overlay SHALL be clipped to the data-area rectangle via `ctx.clip()`

#### Scenario: Overlay re-renders on config change

- **WHEN** a `panel:overlay-changed` event fires for this panel
- **THEN** the canvas SHALL re-render including the updated overlay

#### Scenario: No overlays does not affect pipeline

- **WHEN** no overlay layers are enabled
- **THEN** the render pipeline SHALL behave identically to the pre-overlay version (no performance impact)
