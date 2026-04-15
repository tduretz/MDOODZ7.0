## ADDED Requirements

### Requirement: Scientific display label per field

Each entry in the field registry SHALL have a `label` property containing a Unicode scientific display name. Labels SHALL use Unicode characters for Greek letters, subscripts, and superscripts (no HTML tags). Example labels: `"τ_II"`, `"ε̇_II"`, `"σ_xx"`, `"η"`, `"ρ"`, `"T"`, `"P"`.

#### Scenario: Greek letter labels

- **WHEN** the field registry entry for "Viscosity" is queried
- **THEN** `label` SHALL be `"η"` (Unicode eta)

#### Scenario: Subscripted invariant labels

- **WHEN** the field registry entry for "Stress II" is queried
- **THEN** `label` SHALL be `"τ_II"` (Unicode tau with subscript notation)

#### Scenario: All fields have labels

- **WHEN** any field in the registry is queried
- **THEN** it SHALL have a non-empty `label` string

### Requirement: Formatted unit string per field

Each entry in the field registry SHALL have a `formattedUnit` property containing a Unicode-formatted unit string for display. The existing `unit` property SHALL be preserved unchanged for data purposes. Example formatted units: `"Pa·s"`, `"kg·m⁻³"`, `"s⁻¹"`, `"°C"`, `"W·m⁻³"`.

#### Scenario: Viscosity unit formatting

- **WHEN** the field registry entry for "Viscosity" is queried
- **THEN** `formattedUnit` SHALL be `"Pa·s"`

#### Scenario: Density unit formatting

- **WHEN** the field registry entry for "Density" is queried
- **THEN** `formattedUnit` SHALL be `"kg·m⁻³"` (using Unicode superscript minus-three)

#### Scenario: Dimensionless fields

- **WHEN** a field has no unit (e.g., phase map, strain)
- **THEN** `formattedUnit` SHALL be `""` (empty string)

### Requirement: Labels work in option elements and canvas

Labels and formatted units SHALL use only Unicode characters (no HTML entities or tags) so they render correctly in `<option>` elements, canvas `fillText()`, and plain text contexts.

#### Scenario: Label in select option

- **WHEN** the field dropdown is populated
- **THEN** each `<option>` SHALL display the `label` text correctly without escaped HTML

#### Scenario: Unit in canvas fillText

- **WHEN** the colour bar renders the unit string using `ctx.fillText(formattedUnit, x, y)`
- **THEN** the Unicode characters SHALL render correctly on the canvas
