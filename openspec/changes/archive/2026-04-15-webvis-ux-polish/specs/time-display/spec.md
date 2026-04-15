## ADDED Requirements

### Requirement: Adaptive time formatting function

The module SHALL export a `formatTime(seconds)` function that converts model time in seconds to a human-readable string with an auto-selected unit. The unit thresholds SHALL be:
- `< 1000 yr` → display as years (`"X.XX yr"`)
- `< 1e6 yr` → display as kiloannum (`"X.XX ka"`)
- `≥ 1e6 yr` → display as megaannum (`"X.XX Ma"`)

Conversion factor: 1 yr = 3.1558e7 s.

#### Scenario: Time in years

- **WHEN** `formatTime(6.3116e8)` is called (≈ 20 yr)
- **THEN** it SHALL return an object with `value ≈ 20.0`, `unit = "yr"`, `formatted = "20.00 yr"`

#### Scenario: Time in kiloannum

- **WHEN** `formatTime(3.1558e11)` is called (≈ 10 ka)
- **THEN** it SHALL return `value ≈ 10.0`, `unit = "ka"`, `formatted = "10.00 ka"`

#### Scenario: Time in megaannum

- **WHEN** `formatTime(3.1558e14)` is called (≈ 10 Ma)
- **THEN** it SHALL return `value ≈ 10.0`, `unit = "Ma"`, `formatted = "10.00 Ma"`

#### Scenario: Zero time

- **WHEN** `formatTime(0)` is called
- **THEN** it SHALL return `value = 0`, `unit = "yr"`, `formatted = "0.00 yr"`

### Requirement: Manual unit override

The `formatTime` function SHALL accept an optional second argument specifying the desired unit (`'yr'`, `'ka'`, `'Ma'`). When provided, the function SHALL format with that unit regardless of magnitude.

#### Scenario: Force megaannum display for small time

- **WHEN** `formatTime(3.1558e10, 'Ma')` is called (≈ 1 ka)
- **THEN** it SHALL return `formatted = "0.00 Ma"` (forced unit)

#### Scenario: Force years display for large time

- **WHEN** `formatTime(3.1558e14, 'yr')` is called (≈ 10 Ma)
- **THEN** it SHALL return `formatted = "10000000.00 yr"` (forced unit)

### Requirement: Model stores time unit preference

The model SHALL have a `timeUnit` property that is `null` (auto) by default and can be set to `'yr'`, `'ka'`, or `'Ma'` for manual override. Changing this property SHALL dispatch a `time-unit-changed` event.

#### Scenario: Default is auto

- **WHEN** the model is initialised
- **THEN** `timeUnit` SHALL be `null`

#### Scenario: User overrides unit

- **WHEN** the user selects "Ma" from the time unit control
- **THEN** `model.timeUnit` SHALL be `'Ma'`
- **AND** a `time-unit-changed` event SHALL be dispatched

### Requirement: File selector labels include step and formatted time

Each entry in the file selector SHALL display `"Step N — X.XX Ma"` (or ka/yr) instead of the raw filename. The step number SHALL be parsed from the filename. The time SHALL be formatted using `formatTime` with the model's current `timeUnit` setting.

#### Scenario: File selector with time labels

- **WHEN** the file list is loaded with files having times 0 s, 3.1558e13 s, 6.3116e13 s
- **THEN** the file selector options SHALL show `"Step 0 — 0.00 Ma"`, `"Step 10 — 1.00 Ma"`, `"Step 20 — 2.00 Ma"` (if auto selects Ma)

#### Scenario: Time unit change updates file labels

- **WHEN** the user changes timeUnit from auto to "ka"
- **THEN** all file selector labels SHALL update to use ka units
