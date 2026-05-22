# gmg-defence-materials — delta spec

This change retires ongoing maintenance of GMG defence materials as a live capability. The defence write-ups, plots, and build scripts produced during `add-gmg-stokes-defence` remain discoverable via git history and archived OpenSpec artifacts; the capability to extend or rebuild them as part of the main repo workflow is removed.

## REMOVED Requirements

### Requirement: Defence document follows the fixed nine-section structure

**Reason**: No new GMG defence documents are produced or maintained in-repo after GMG removal.

**Migration**: Consult archived change `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/` and `benchmarks/ec2/results/PERFORMANCE_REPORT.md` (on disk where present) for the historical narrative and numbers.

### Requirement: Measurements section carries defensible numbers

**Reason**: Defence document maintenance is retired with the GMG experiment closure.

**Migration**: Use `PERFORMANCE_REPORT.md` and post-fix CSVs under `benchmarks/ec2/results/post_fixF/` as the canonical measurement record.

### Requirement: V-cycle animation is produced and embedded

**Reason**: V-cycle tooling (`benchmarks/vcycle_vis/`) is deleted; no new animations are produced.

**Migration**: Recover scripts and outputs from git history at pre-removal commits if needed.

### Requirement: Process diagrams are clean, single-concept, print-friendly

**Reason**: Defence diagram maintenance is retired.

**Migration**: Historical diagrams remain in archived change attachments / git history.

### Requirement: Engine fixes surfaced during benchmarking are itemised in the defence

**Reason**: The defence workflow is closed; future engine fixes are tracked in their own changes.

**Migration**: Reference archived `design.md` / `proposal.md` under `2026-04-21-add-gmg-stokes-defence`.

## ADDED Requirements

### Requirement: Published `gmg-defence-materials` spec is a retirement tombstone

The file `openspec/specs/gmg-defence-materials/spec.md` SHALL be replaced with a short **RETIRED** notice that points to archived change `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/` for the historical defence scope, names `PERFORMANCE_REPORT.md` as the external measurement record where available, and states that no further defence deliverables are maintained in the main workflow.

#### Scenario: Reader checks defence capability status

- **WHEN** a reader opens `openspec/specs/gmg-defence-materials/spec.md` after this change is applied
- **THEN** the first sections make clear the capability is retired and reference archives rather than implying active maintenance

### Requirement: Non-canonical defence working files are removed from the archived change directory

Files named `STATUS.md`, `defence.md`, `defence.html`, and `build_pdf.py` under `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/` SHALL be deleted from the working tree. Canonical artifacts `proposal.md`, `design.md`, `tasks.md`, `specs/`, and `.openspec.yaml` SHALL remain.

#### Scenario: Archive directory inspection

- **WHEN** a maintainer lists `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/`
- **THEN** only spec-driven canonical files plus `.openspec.yaml` remain; process-only deliverables are absent
