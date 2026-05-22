## Purpose

Title interpolation system for WebVis panels — editable title templates with `${...}` token replacement, a token dropdown catalogue, and lazy statistics caching.

## Requirements

### Requirement: Editable title input per panel

Each panel SHALL display an editable `<input type="text">` field where the user can type a title template containing free text and interpolation tokens. Below the input, a read-only rendered title `<span>` SHALL display the interpolated result with live data values. The title input SHALL be pre-filled with a default template: `"${field} at t = ${time} — min = ${min}  max = ${max}"`.

#### Scenario: Default title on panel creation

- **WHEN** a new panel is created
- **THEN** the title input SHALL contain `"${field} at t = ${time} — min = ${min}  max = ${max}"`
- **AND** the rendered title SHALL show the interpolated values (or empty placeholders if no field is loaded)

#### Scenario: User edits the title template

- **WHEN** the user changes the title input to `"${field} — ${mean} (${unit})"`
- **THEN** the rendered title SHALL update to show the current field name, mean value, and unit

#### Scenario: Title persists across field changes

- **WHEN** the user has set a custom title template and switches to a different field
- **THEN** the title template text SHALL remain unchanged
- **AND** the rendered title SHALL update with the new field's data values

### Requirement: Trigger-character dropdown for token insertion

When the user types `$` in the title input (or clicks an insert button ⊕ next to the input), a floating dropdown SHALL appear listing all available interpolation variables grouped under two headers: **Metadata** (scalar tokens) and **Field stats** (array-derived tokens). Selecting an item SHALL insert the token text `${...}` at the current cursor position in the input and close the dropdown. Pressing Escape or clicking outside SHALL dismiss the dropdown without inserting.

#### Scenario: Trigger dropdown with $ character

- **WHEN** the user types `$` in the title input
- **THEN** a dropdown SHALL appear showing all available tokens grouped by category

#### Scenario: Trigger dropdown with insert button

- **WHEN** the user clicks the ⊕ insert button next to the title input
- **THEN** the same dropdown SHALL appear

#### Scenario: Select a token from dropdown

- **WHEN** the dropdown is open and the user clicks "Model time" under the Metadata heading
- **THEN** `${time}` SHALL be inserted at the cursor position in the title input
- **AND** the dropdown SHALL close
- **AND** the rendered title SHALL update immediately

#### Scenario: Dismiss dropdown without selection

- **WHEN** the dropdown is open and the user presses Escape
- **THEN** the dropdown SHALL close without modifying the title input

### Requirement: Variable catalogue — scalar tokens

The following scalar tokens SHALL be available in the dropdown under the **Metadata** heading. Each token SHALL resolve to a single value from the current file's parameters:

| Token | Dropdown label | Description |
|-------|---------------|-------------|
| `${step}` | Step number | Integer time-step index |
| `${time}` | Model time | Formatted with adaptive units (yr/ka/Ma) |
| `${Nx}` | Grid columns | Number of horizontal grid cells |
| `${Nz}` | Grid rows | Number of vertical grid cells |
| `${resolution}` | Grid resolution | Composite display: `Nx × Nz` (e.g. `600 × 400`) |

#### Scenario: Step token renders integer

- **WHEN** the loaded file is `Output00050.gzip.h5` with step 50
- **THEN** `${step}` SHALL render as `50`

#### Scenario: Time token renders with units

- **WHEN** the loaded file has time = 4.6928e13 s and the time unit is Ma
- **THEN** `${time}` SHALL render as `1.49 Ma`

#### Scenario: Resolution token renders composite

- **WHEN** the file has Nx=600, Nz=400
- **THEN** `${resolution}` SHALL render as `600 × 400`

### Requirement: Variable catalogue — field statistics tokens

The following field statistics tokens SHALL be available in the dropdown under the **Field stats** heading. Each token SHALL resolve using the current panel's loaded field data:

| Token | Dropdown label | Description |
|-------|---------------|-------------|
| `${min}` | Min | Absolute minimum of the field values |
| `${max}` | Max | Absolute maximum of the field values |
| `${mean}` | Mean | Arithmetic mean of the field values |
| `${median}` | Median | Median of the field values |
| `${p2}` | Robust min (p2) | 2nd percentile of the field values |
| `${p98}` | Robust max (p98) | 98th percentile of the field values |
| `${field}` | Field name | Scientific label from the field-labels registry (e.g. `ε̇_II`, `τ_II`, `η`) |
| `${unit}` | Unit | Formatted unit string from the field-labels registry (e.g. `Pa`, `Pa·s`, `s⁻¹`) |

#### Scenario: Min and max tokens render with scientific notation

- **WHEN** the panel shows "Viscosity" with min=1e18, max=1e24
- **THEN** `${min}` SHALL render as `1.00e+18 Pa·s` and `${max}` SHALL render as `1.00e+24 Pa·s`

#### Scenario: Mean token computed on demand

- **WHEN** the title template contains `${mean}` and the field data is loaded
- **THEN** the mean SHALL be computed from the flat values array and rendered with scientific notation and the field's unit

#### Scenario: Median token computed on demand

- **WHEN** the title template contains `${median}` and the field data is loaded
- **THEN** the median SHALL be computed from the flat values array and rendered with scientific notation and the field's unit

#### Scenario: Field name token uses scientific label

- **WHEN** the panel's field is "Strain rate II"
- **THEN** `${field}` SHALL render as `ε̇_II`

#### Scenario: Unit token uses formatted unit

- **WHEN** the panel's field is "Stress II"
- **THEN** `${unit}` SHALL render as `Pa`

### Requirement: Title rendering function

A pure function `renderTitle(template, context)` SHALL accept a template string and a context object containing all token values, and SHALL return the fully interpolated string. The function SHALL replace all `${...}` tokens with their corresponding formatted values from the context. Unrecognised tokens SHALL be rendered as-is (left in the template unchanged).

#### Scenario: Full interpolation

- **WHEN** `renderTitle("${field} at t = ${time}", { field: "η", time: "4.69 Ma", ... })` is called
- **THEN** the return value SHALL be `"η at t = 4.69 Ma"`

#### Scenario: Mixed text and tokens

- **WHEN** `renderTitle("Title: ${field} — max = ${max}", context)` is called
- **THEN** tokens SHALL be replaced and surrounding text SHALL be preserved verbatim

#### Scenario: Unknown token left as-is

- **WHEN** the template contains `${unknown}`
- **THEN** `${unknown}` SHALL remain literally in the rendered output

#### Scenario: No tokens in template

- **WHEN** the template is `"My custom title"` with no tokens
- **THEN** the rendered output SHALL be `"My custom title"` unchanged

### Requirement: Statistics caching in PanelState

When `${mean}` or `${median}` tokens are computed for a panel's current field data, the computed values SHALL be cached in the PanelState. The cache SHALL be invalidated when the panel's field data changes (field switch or file switch). Subsequent renders of the same title with the same field data SHALL use the cached values without recomputation.

#### Scenario: Cache populated on first render

- **WHEN** a title with `${mean}` is rendered for the first time after loading a field
- **THEN** the mean SHALL be computed and stored in PanelState

#### Scenario: Cache reused on re-render

- **WHEN** the title is re-rendered (e.g. template edited) without changing the field data
- **THEN** the cached mean value SHALL be used without recomputation

#### Scenario: Cache invalidated on field change

- **WHEN** the user switches the panel to a different field
- **THEN** the cached mean and median values SHALL be cleared
