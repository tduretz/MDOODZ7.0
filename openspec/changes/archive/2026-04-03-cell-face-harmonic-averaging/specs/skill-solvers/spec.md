## ADDED Requirements

### Requirement: Document cell_avg parameter in solver accuracy section
The "Numerical Accuracy & Convergence Orders" section of skill-solvers SHALL document the `cell_avg` parameter, its effect on convergence orders, and the measured L2 improvements from the SolVi benchmark.

#### Scenario: User asks about improving pressure convergence
- **WHEN** a user asks about improving convergence orders or pressure accuracy
- **THEN** the skill SHALL mention `cell_avg = 1` as an option alongside the existing `eta_average` discussion
