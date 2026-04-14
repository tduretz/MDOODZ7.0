## ADDED Requirements

### Requirement: PCG convergence plot
The system SHALL include a C++ program and gnuplot script that produce a PNG figure showing PCG residual convergence (log₁₀ of relative residual vs iteration number) with two curves: step 0 (cold start) and a later step (warm start).

#### Scenario: Convergence plot generated
- **WHEN** the visual test runs with `pcg_residuals_step00000.csv` and `pcg_residuals_step00004.csv` available
- **THEN** it SHALL produce `VISUAL_TESTS/img/pcg_convergence.png` with two labelled curves, log-scale y-axis, and iteration count on x-axis

#### Scenario: Plot shows cold vs warm start difference
- **WHEN** the plot is generated
- **THEN** the step 0 curve SHALL show more iterations to reach tolerance than the later step, visually demonstrating temporal warm-starting

### Requirement: L2 error bar chart
The system SHALL include a C++ program and gnuplot script that produce a PNG bar chart comparing L2 errors: analytic vs CHOLMOD and analytic vs PCG, for both the geotherm and Gaussian diffusion tests.

#### Scenario: Bar chart generated
- **WHEN** the visual test runs after both CHOLMOD and PCG test HDF5 outputs are available
- **THEN** it SHALL produce `VISUAL_TESTS/img/thermal_accuracy.png` with grouped bars showing L2 errors for each solver and test case

#### Scenario: Bars for both test cases
- **WHEN** the chart is generated
- **THEN** it SHALL contain bars for SteadyStateGeotherm (CHOLMOD, PCG) and GaussianDiffusion (CHOLMOD, PCG), with error values printed above each bar

### Requirement: Interp difference heatmap
The system SHALL include a C++ program and gnuplot script that produce a PNG heatmap of the absolute temperature difference $|T_{mode0} - T_{mode3}|$ across the grid.

#### Scenario: Heatmap generated
- **WHEN** the visual test runs after both interp_mode=0 and interp_mode=3 HDF5 outputs are available
- **THEN** it SHALL produce `VISUAL_TESTS/img/interp_diff.png` with a 2D colour map, a colour bar showing the difference magnitude, and axis labels in physical coordinates

#### Scenario: Heatmap shows near-zero differences
- **WHEN** the heatmap is generated for a correct implementation
- **THEN** the maximum value in the colour bar SHALL be less than 1e-6, visually confirming equivalence
