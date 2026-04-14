## ADDED Requirements

### Requirement: PCG residual CSV export parameter
The system SHALL read an `export_pcg_residuals` integer parameter (default 0) from the `.txt` input file. When set to 1 and `thermal_solver = 1`, the PCG solver SHALL write per-iteration residual norms to CSV files.

#### Scenario: Parameter read from .txt
- **WHEN** the `.txt` file contains `export_pcg_residuals = 1`
- **THEN** the system SHALL set `model.export_pcg_residuals = 1` and enable residual CSV output for all PCG thermal solves

#### Scenario: Default is disabled
- **WHEN** the `.txt` file does not contain `export_pcg_residuals`
- **THEN** the system SHALL default to 0 and no residual files SHALL be written

### Requirement: Per-step residual CSV file
When residual export is enabled, the PCG solver SHALL write one CSV file per thermal solve invocation, named `pcg_residuals_step<NNNNN>.csv` in the simulation's output directory (determined by `writer_subfolder`).

#### Scenario: CSV file format
- **WHEN** a PCG thermal solve completes with export enabled
- **THEN** the CSV file SHALL contain a header line `iteration,abs_residual,rel_residual` followed by one row per iteration, where `abs_residual` is $\|r_k\|$ and `rel_residual` is $\|r_k\|/\|r_0\|$

#### Scenario: Step numbering matches simulation step
- **WHEN** the thermal solve occurs at simulation step 5
- **THEN** the output file SHALL be named `pcg_residuals_step00005.csv`

#### Scenario: Cold start vs warm start visible in data
- **WHEN** residual files from step 0 (cold start) and step N (warm start) are compared
- **THEN** step 0 SHALL have more iterations than step N (assuming temporal coherence), and this SHALL be visible in the iteration count and residual curves
