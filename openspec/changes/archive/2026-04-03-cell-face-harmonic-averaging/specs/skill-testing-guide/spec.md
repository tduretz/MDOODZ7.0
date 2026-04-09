## ADDED Requirements

### Requirement: Document cell_avg experimental comparison
The skill-testing-guide SHALL include a "Cell-Face Averaging" subsection after the existing "Viscosity Averaging — eta_average" section. It SHALL contain a comparison table of L2 errors and convergence orders for `cell_avg = 0` vs `cell_avg = 1` on the SolVi benchmark, and a recommendation for when to use each mode.

#### Scenario: User setting up a benchmark test
- **WHEN** a user is configuring a SolVi-type benchmark and wants best accuracy
- **THEN** the skill SHALL recommend `cell_avg = 1` for pressure accuracy and document the measured improvement
