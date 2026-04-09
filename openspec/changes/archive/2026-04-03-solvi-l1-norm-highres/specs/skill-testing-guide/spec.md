## MODIFIED Section: SolVi Benchmark Tests

### Requirement: Document L1 norm findings
The SolVi section of the testing guide skill SHALL document that L1 norm (`mean(|num - ana|)`) was added alongside L2 to better measure pressure convergence. It SHALL record the measured L1 convergence orders at high resolution and note how they compare to L2 orders.

#### Scenario: Skills reflects current test inventory
- **WHEN** a developer reads skill-testing-guide
- **THEN** the SolVi benchmark section SHALL mention both L2 and L1 norms and list the HighResL1Convergence test alongside existing CI tests
