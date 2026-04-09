## MODIFIED Section: SolVi Convergence

### Requirement: Document L1 convergence findings
The SolVi convergence subsection of skill-solvers SHALL document that L1 norm at resolutions >80 cells recovers expected first-order pressure convergence (order ~1.0), while L2 norm at coarse resolutions (≤81) shows a misleadingly low ~0.75 due to the inclusion interface.

#### Scenario: Skills reflects L1 vs L2 insight
- **WHEN** a developer reads about SolVi convergence in skill-solvers
- **THEN** the section SHALL explain that L1 norm is recommended for evaluating P convergence and list measured orders at 101→201
