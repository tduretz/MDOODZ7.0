## MODIFIED Section: Benchmark Accuracy Findings

### Requirement: Document per-test accuracy expectations and parameter recommendations
The skill-testing-guide SHALL include a section documenting the measured accuracy for each analytical benchmark test, which parameters matter most, and recommended settings.

#### Scenario: Developer reads about benchmark accuracy
- **WHEN** a developer reads skill-testing-guide for benchmark information
- **THEN** the skill SHALL list for each test: measured L2 and L1 errors, recommended .txt parameters, and the key finding (e.g., "penalty dominates P accuracy" or "Nt must reach steady state")

### Requirement: Document L1 vs L2 findings across non-SolVi tests
The skill-testing-guide SHALL document whether L1 norm provides meaningfully different accuracy measurements compared to L2 for the non-SolVi benchmarks.

#### Scenario: Developer asks about norm choice
- **WHEN** a developer asks which norm to use for a new benchmark
- **THEN** the skill SHALL provide guidance summarizing L1 vs L2 findings from both SolVi and non-SolVi experiments
