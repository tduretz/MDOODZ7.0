# gmg-performance-harness — delta spec

This change removes the EC2/local performance harness that implemented repeatable GMG vs CHOLMOD sweeps. Benchmark results already captured (for example in `PERFORMANCE_REPORT.md` and CSVs) remain valid historical records outside this spec.

## REMOVED Requirements

### Requirement: Harness lives under `benchmarks/ec2/` as reusable infrastructure

**Reason**: The entire `benchmarks/ec2/` tree (scripts, Python helpers, YAML grid configs) is deleted as part of GMG closure.

**Migration**: Recover harness from git history at pre-removal commits if sweeps must be reproduced.

### Requirement: EC2 instances are pinned to `c7i.8xlarge` in `eu-west-1`

**Reason**: EC2 provisioning scripts are deleted with the harness.

**Migration**: N/A.

### Requirement: Every configuration is measured 5 times with median + IQR reporting and throttle rejection

**Reason**: Automated measurement harness is removed.

**Migration**: Manual or ad-hoc benchmarking if needed; historical methodology remains in archived defence materials.

### Requirement: Results format is one-row-per-run CSV with full provenance

**Reason**: Harness-generated CSV pipeline is removed.

**Migration**: Retain existing CSV files on disk or in version control where already committed.

### Requirement: Harness enforces a cost safeguard and guarantees teardown

**Reason**: `provision.sh` / `teardown.sh` and related scripts are deleted.

**Migration**: N/A.

### Requirement: Grid sweep covers the regression envelope declared in the previous change

**Reason**: Grid sweep automation is removed.

**Migration**: Historical sweep definitions live in archived commits and `grids*.yaml` history.

## ADDED Requirements

### Requirement: Published `gmg-performance-harness` spec is a retirement tombstone

The file `openspec/specs/gmg-performance-harness/spec.md` SHALL be replaced with a short **RETIRED** notice stating that the `benchmarks/ec2/` harness is deleted, pointing to archived changes that defined it, and naming `PERFORMANCE_REPORT.md` (where available on disk) as the consolidated results record.

#### Scenario: Reader checks harness status

- **WHEN** a reader opens `openspec/specs/gmg-performance-harness/spec.md` after apply
- **THEN** the document states the harness no longer exists in-tree and where to find historical methodology

### Requirement: Harness tree is absent from the repository

The paths `benchmarks/ec2/` and `benchmarks/vcycle_vis/` SHALL NOT exist in the working tree after apply (except possibly empty directories or gitignored result subtrees the user keeps locally).

#### Scenario: Repository layout after apply

- **WHEN** a developer lists `benchmarks/` at the repository root
- **THEN** neither `ec2/` nor `vcycle_vis/` contains tracked harness sources from this capability
