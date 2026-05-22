## ADDED Requirements

### Requirement: L2 error comparison between CHOLMOD and PCG against analytic geotherm
The system SHALL include a GTest case that runs the SteadyStateGeotherm problem with both `thermal_solver=0` (CHOLMOD) and `thermal_solver=1` (PCG) and computes the relative L2 error of each against the analytic linear geotherm $T(z) = T_{top} + (T_{bot} - T_{top})(z_{top} - z)/H$.

#### Scenario: Both solvers achieve the same L2 error
- **WHEN** both SteadyStateGeotherm and SteadyStateGeothermPCG tests complete on a 51×51 grid with 20 steps
- **THEN** the L2 error of each against the analytic solution SHALL be less than 1e-4, and the absolute difference between the two L2 errors SHALL be less than 1e-6

#### Scenario: L2 errors are printed for inspection
- **WHEN** the comparison test runs
- **THEN** the test SHALL print both L2 errors and their difference to stdout (e.g. `CHOLMOD L2: 2.4e-6, PCG L2: 2.4e-6, diff: 1.2e-12`)

### Requirement: Transient Gaussian diffusion L2 test with analytic solution
The system SHALL include GTest cases (`GaussianDiffusionL2` and `GaussianDiffusionL2PCG`) that compute the L2 error of the numerical temperature field against the transient analytic solution $T(r,t) = T_{bg} + \Delta T \cdot R^2/(R^2+4\kappa t) \cdot \exp(-r^2/(R^2+4\kappa t))$ where $\kappa = k/(\rho C_p)$.

#### Scenario: CHOLMOD Gaussian diffusion L2
- **WHEN** the GaussianDiffusionL2 test runs with `thermal_solver=0` on a 51×51 grid for 5 steps
- **THEN** the L2 error of the final temperature field against the analytic solution SHALL be less than 1e-2

#### Scenario: PCG Gaussian diffusion L2
- **WHEN** the GaussianDiffusionL2PCG test runs with `thermal_solver=1` on a 51×51 grid for 5 steps
- **THEN** the L2 error SHALL be less than 1e-2, and the absolute difference between the CHOLMOD and PCG L2 errors SHALL be less than 1e-6

### Requirement: Dedicated parameter files for L2 comparison tests
The system SHALL include `.txt` parameter files in `TESTS/Thermal/` for each new test case, using 51×51 grid, realistic scaling, and distinct `writer_subfolder` values.

#### Scenario: GaussianDiffusionL2 parameter file
- **WHEN** the test `GaussianDiffusionL2` is run
- **THEN** it SHALL read `Thermal/GaussianDiffusionL2.txt` with `thermal_solver = 0` and `writer_subfolder = GaussianDiffusionL2`

#### Scenario: GaussianDiffusionL2PCG parameter file
- **WHEN** the test `GaussianDiffusionL2PCG` is run
- **THEN** it SHALL read `Thermal/GaussianDiffusionL2PCG.txt` with `thermal_solver = 1` and `writer_subfolder = GaussianDiffusionL2PCG`
