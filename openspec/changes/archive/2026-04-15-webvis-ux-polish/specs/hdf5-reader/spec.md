## ADDED Requirements

### Requirement: Compute percentile-clipped range

The reader SHALL export a `computePercentileRange(flat, pLow, pHigh)` function that computes the pLow-th and pHigh-th percentile values from a flat `Float64Array`, ignoring NaN values. The default percentiles SHALL be 2 and 98. The function SHALL use a randomized quickselect (partial sort) algorithm with O(n) average time complexity.

#### Scenario: Percentile range excludes outliers

- **WHEN** `computePercentileRange` is called on an array where 1% of values are 0 and 1% are extremely large
- **THEN** `pMin` SHALL be above 0 and `pMax` SHALL be below the extreme values

#### Scenario: All values identical

- **WHEN** `computePercentileRange` is called on an array where all non-NaN values are equal
- **THEN** `pMin` SHALL equal `pMax` SHALL equal that value

#### Scenario: NaN values excluded from percentile computation

- **WHEN** the array contains NaN values (air cells)
- **THEN** NaN values SHALL be excluded before computing rank positions

### Requirement: Batch metadata extraction

The reader SHALL export a `readAllMetadata(dataDir, filenames)` function that reads `/Model/Params[0]` (time in seconds) from every file in the list and returns an array of `{ name, step, time }` objects. The step number SHALL be parsed from the filename regex `Output(\d+)\.gzip\.h5`.

#### Scenario: Read metadata from all files

- **WHEN** `readAllMetadata` is called with 29 filenames
- **THEN** it SHALL return 29 objects each with `name` (string), `step` (integer), and `time` (number in seconds)

#### Scenario: Step number extraction

- **WHEN** the filename is `Output00050.gzip.h5`
- **THEN** `step` SHALL be `50`

## MODIFIED Requirements

### Requirement: Return structured field response

The reader SHALL return field data in a structured format containing: `xCoords` (1D array), `zCoords` (1D array), `values` (2D array, row-major), `nx` (number of columns), `nz` (number of rows), `unit` (string), `log` (boolean — whether log₁₀ display is recommended), `min` (numeric, excluding NaN), `max` (numeric, excluding NaN), `pMin` (numeric, 2nd percentile excluding NaN), `pMax` (numeric, 98th percentile excluding NaN).

#### Scenario: Complete response structure

- **WHEN** the reader extracts any field
- **THEN** the response SHALL contain all required properties including `pMin` and `pMax`
- **AND** `min` and `max` SHALL exclude NaN values
- **AND** `pMin` SHALL be >= `min` and `pMax` SHALL be <= `max`

#### Scenario: Percentile range used for auto-range

- **WHEN** the frontend receives a field response
- **THEN** it SHALL use `pMin` and `pMax` for the default colour range instead of `min` and `max`
