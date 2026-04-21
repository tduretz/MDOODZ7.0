## Why

The GMG experiment is complete with a **negative outcome**: `lin_solver = 3` is not faster than CHOLMOD on any MDOODZ workload and cannot be made faster without architectural changes the user has explicitly ruled out (Galerkin coarsening → bandwidth widening, AMG → new dependency). Post-Fix F (`ad39d7f`) measurements in `benchmarks/ec2/results/PERFORMANCE_REPORT.md §§11–13` establish the envelope precisely:

- **2D SolVi, 81²–801²**: GMG is consistently 1.8–3.4× slower than CHOLMOD (at best — after tuning `gmg_levels = 4, tol = 1e-7` at 801²; auto-level is 29× slower with fallback).
- **2D Rifting, any developed state**: GMG stagnates at ρ ≈ 1.000 after one order-of-magnitude residual reduction, hits the 600-iter budget, falls back to CHOLMOD every solve — 67 s wasted per linear solve at 151×101.
- **Memory**: not a constraint for the user; CHOLMOD fits all production grids in available RAM.

No regime exists where current GMG beats CHOLMOD on speed, and speed is the user's sole remaining optimisation axis. The ~4,900 LOC of GMG code in `MDLIB/` plus its ~20 test fixtures, its `benchmarks/ec2/` harness (2,200 LOC of Python + shell), and its `benchmarks/vcycle_vis/` diagnostic tooling impose a continuous maintenance tax (CMake, HDF5-dump paths, parameter parsing, equivalence fixtures) for zero production value. This change **removes the implementation cleanly**, revert-style, while preserving the full experimental record in archived OpenSpec changes and git history so any future GMG investigation has a direct reference.

The closure is deliberate and final: **no further GMG work is planned**. If CHOLMOD ever becomes an actual speed bottleneck, the successor architecture is `add-recycled-cholmod-preconditioner` (cached factorisation reused across time steps) or a block-Schur preconditioner — *not* more multigrid. This change does not reserve scope for either; those are new proposals if and when the need arises.

## What Changes

### Code removals (revert MDLIB to pre-GMG baseline)

- **REMOVE** newly-added MDLIB files (6 files, ~4,880 LOC introduced post-`cab45077`):
  - `MDLIB/MultigridStokes.{c,h}`
  - `MDLIB/MultigridLevels.{c,h}`
  - `MDLIB/StokesAssemblyGMG.{c,h}`
- **REVERT** modified MDLIB files to their content at `cab45077c7f6ccaf447e3422b820a96d3c94169e` (commit "add new reseeding mode", 2026-04-17):
  - `MDLIB/CMakeLists.txt` — drop GMG sources from the build target.
  - `MDLIB/InputOutput.c` — drop `gmg_*` parameter parsing.
  - `MDLIB/Main_DOODZ.c` — drop GMG dispatch cleanup.
  - `MDLIB/StokesRoutines.c` — drop `lin_solver = 3` branch in `SolveStokes`; user attempts to set `lin_solver = 3` SHALL produce a clear `LOG_ERR + abort` at startup (not a silent fallback, so nobody discovers the removal mid-simulation).
  - `MDLIB/ChemicalRoutines.c` — revert incidental GMG-era edit.
  - `MDLIB/include/mdoodz.h` — drop `gmg_*` fields from the `params` struct.
  - `MDLIB/mdoodz-private.h` — drop GMG function declarations.

### Test removals

- **REMOVE** all GMG-specific test fixtures under `TESTS/` (equivalence, dispatch, multigrid-internal tests):
  - `TESTS/GmgStokesEquivalence.cpp`
  - `TESTS/MultigridStokesTests.cpp`
  - `TESTS/MultigridTests.cpp`
  - `TESTS/DisabledBenchmarks/*Gmg*.cpp`, `UplegAmplificationBound.cpp`, `GmgLevelsSweepProbe.cpp`
  - `TESTS/AnalyticalBench/SolKz_{chol,gmg}.txt`
  - `TESTS/ShearBand/{Conjugate,Single}Bands*{chol,gmg}.txt`
  - `TESTS/SolViBenchmark/SolVi*{gmg,chol}*.txt` and `SolViUpleg161.txt`
- **REVERT** `TESTS/CMakeLists.txt` to the `cab45077` baseline (drop GMG-test registrations).
- **PRESERVE** `TESTS/RotationAdvectionTests.cpp` with its `887e8bc9` fix intact (this commit touches `TESTS/` only, not `MDLIB`, and is independent of GMG).

### Benchmark and driver removals

- **REMOVE** `benchmarks/ec2/` entirely — harness (`run_local.sh`, `run_perf_sweep.sh`, `provision.sh`, `teardown.sh`, `_common.sh`, `_parse_grids.py`, `analyse.py`, `_aggregate_perf_detail.py`, `_analyse_vcycle_dump.py`), yaml configs, and the `README.md`. All were created for GMG performance measurement.
- **REMOVE** `benchmarks/vcycle_vis/` entirely — V-cycle residual visualisation tooling, single-purpose for GMG diagnostic work.
- **REMOVE** the `SolViPerf` performance driver: `SETS/SolViPerf.c`, `SETS/SolViPerf.txt`, `cmake-exec/SolViPerf/` (CMake target). Built solely for the GMG perf sweep.

### OpenSpec artifact cleanup

- **RETIRE** specs by replacing their content with a one-page "capability retired" notice pointing to the archived changes for historical requirements (keeps the spec file present so cross-references from archived deltas still resolve, but makes it unmistakable the capability is gone):
  - `openspec/specs/gmg-stokes-solver/spec.md`
  - `openspec/specs/gmg-defence-materials/spec.md`
  - `openspec/specs/gmg-performance-harness/spec.md`
- **CLEAN** non-canonical files from archived GMG changes (they describe process/status/retrospective, not specification content, and don't belong in the permanent record per the user's "clean up artifacts that are not spec, proposal, design etc" directive):
  - `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/{STATUS,RETROSPECTIVE,DECISIONS}.md`
  - `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/{STATUS.md, defence.{md,html}, build_pdf.py}`
  - `openspec/changes/archive/2026-04-18-add-gmg-upleg-fix/{STATUS.md,*.md}` (any non-canonical files)
- **PRESERVE** the canonical artifacts in each archive (`proposal.md`, `design.md`, `tasks.md`, `specs/`). These are the experimental record the user wants to keep.

### What this change does NOT do

- **NOT** implement Galerkin coarsening, AMG, or any other multigrid variant. Ruled out by user.
- **NOT** implement `add-recycled-cholmod-preconditioner`. Separate future proposal if needed.
- **NOT** touch `lin_solver ∈ {-1, 0, 1, 2}` paths, the Powell-Hestenes iteration, the direct solver, or the thermal PCG path.
- **NOT** remove git history. The current branch head and every archived change preserve the experimental code; `git show cab45077..HEAD -- MDLIB` reproduces the full implementation.
- **NOT** remove archived change `proposal.md` / `design.md` / `tasks.md` / `specs/` content. These are the permanent record of what was tried and why it did not work.

## Capabilities

### New Capabilities

None. This change is a net removal.

### Modified Capabilities

- `gmg-stokes-solver`: all requirements MARKED REMOVED. The capability is retired; the spec file is retained as a tombstone with a pointer to the archived `add-gmg-stokes-solver`, `add-gmg-stokes-defence`, and `add-gmg-upleg-fix` changes that historically defined and measured it. No active requirements remain.
- `gmg-defence-materials`: all requirements MARKED REMOVED. The defence materials (analysis write-ups, plots, figures) described by this capability were produced, archived with `2026-04-21-add-gmg-stokes-defence`, and are no longer maintained going forward. Spec tombstoned with a pointer to the archive.
- `gmg-performance-harness`: all requirements MARKED REMOVED. The `benchmarks/ec2/` harness that implemented this capability is deleted in this change. Spec tombstoned with a pointer to the archived `add-gmg-stokes-defence` change where the harness was defined and the `PERFORMANCE_REPORT.md` where it was used.

## Impact

**Removed files** (git-tracked, ~5,500 LOC in `MDLIB` + tests + benchmarks):

- MDLIB: 6 files deleted, 7 files reverted (diff ~95 lines against `cab45077`).
- TESTS: ~20 fixture files + registrations.
- benchmarks/: 2 subdirectories, ~20 files.
- SETS / cmake-exec: 3 files + 1 target.
- openspec archives: ~7 non-canonical process files.

**Modified files**:

- `openspec/specs/gmg-{stokes-solver, defence-materials, performance-harness}/spec.md` — content replaced with retirement notices.
- None in MDLIB besides the revert (the revert itself is the modification).

**APIs / user-facing**:

- `lin_solver = 3` is REMOVED as a valid input. Users who still have it in `.txt` files will receive a clear startup error listing the valid `lin_solver` range and pointing at the archived GMG changes.
- All `gmg_*` parameters (`gmg_levels`, `gmg_nu_pre`, `gmg_nu_post`, `gmg_fgmres_restart`, `gmg_fgmres_max_restarts`, `gmg_fgmres_tol`, `gmg_standalone`, `gmg_dump_vcycle`) are REMOVED. Their presence in `.txt` files is ignored with a single WARN line at startup naming the removed parameter family (so legacy config files do not crash but the unused parameters surface visibly).
- No HDF5 schema change.
- No CLI flag change.

**Dependencies**:

- UMFPACK: the GMG coarse solver used it. Check whether any other MDLIB path still links UMFPACK; if not, the dependency can also be dropped in `MDLIB/CMakeLists.txt`. Scoped under D2 in `design.md`.
- No other dependency changes.

**Build system**:

- `MDLIB/CMakeLists.txt` and `TESTS/CMakeLists.txt` shrink.
- `cmake-exec/SolViPerf/CMakeLists.txt` deleted entirely; root CMake loses one `add_subdirectory` or equivalent.

**Risks**:

- **Revert introduces a conflict with non-GMG edits made in the GMG era**: `MDLIB/ChemicalRoutines.c` was touched incidentally (4 lines per the diff stat). If any of those 4 lines fix a real bug and we revert them, we'd regress a non-GMG feature. Mitigation: `design.md` D3 requires a line-by-line audit of the `cab45077..HEAD` diff for `ChemicalRoutines.c`, `InputOutput.c`, `Main_DOODZ.c`, `StokesRoutines.c`, `include/mdoodz.h`, `mdoodz-private.h` before applying the revert. Any non-GMG change gets re-applied as a separate fixup commit on top of the revert.
- **`887e8bc9` rotation-test fix is lost**: `887e8bc9` only touches `TESTS/RotationAdvectionTests.cpp`, not MDLIB, so reverting MDLIB cannot lose it. Verified via `git show 887e8bc9 --stat`. Mitigation not required; documented for transparency.
- **Archived changes reference files that no longer exist**: after the removal, archived `proposal.md` / `design.md` / `tasks.md` will reference `MDLIB/MultigridStokes.c` etc. that are deleted in main. This is acceptable and intentional — archives describe history, not current state. Mitigation: a single top-of-spec pointer in each retirement notice redirects forward-looking readers; the archive itself is left untouched.
- **User has `lin_solver = 3` in a pinned `.txt`**: the defence recorded zero production users of `lin_solver = 3`; the user is the sole operator and is the one requesting removal. No third-party impact.
- **Test suite shrinks, coverage metrics drop**: expected and accepted. The removed tests covered the removed code; nothing else regresses.

**Rollback**: the change is a single coherent revert+cleanup on a feature branch. If the user changes their mind later, `git revert <this-change-commit>` restores everything, or (more likely) the user re-checks out the current branch (`add-performance-metrics-clean`) which carries the full GMG history as of `c59ce62`.
