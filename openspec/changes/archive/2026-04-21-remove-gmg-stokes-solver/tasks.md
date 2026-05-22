## 1. MDLIB — delete GMG sources

- [x] 1.1 Delete `MDLIB/MultigridStokes.c`, `MDLIB/MultigridStokes.h`, `MDLIB/MultigridLevels.c`, `MDLIB/MultigridLevels.h`, `MDLIB/StokesAssemblyGMG.c`, `MDLIB/StokesAssemblyGMG.h`.
- [x] 1.2 Edit `MDLIB/CMakeLists.txt` to remove the three GMG `.c` files from the library target (and any GMG-only link flags if present).
- [x] 1.3 Edit `MDLIB/StokesRoutines.c`: remove `#include "MultigridStokes.h"` and both `lin_solver == 3` blocks in `SolveStokesDecoupled` and `SolveStokesDefectDecoupled` (restore pre-GMG dispatch only).
- [x] 1.4 Edit `MDLIB/include/mdoodz.h`: remove the `gmg_*` fields from `params` (single contiguous block).
- [x] 1.5 Edit `MDLIB/InputOutput.c`: remove all `gmg_*` reads and defaults; add `lin_solver == 3` detection that `LOG_ERR` + `abort` before the run starts (message per `design.md` D3); add aggregated `WARN` for any parsed keys matching `gmg_*` prefix (per `design.md` D4). Audit the full `cab45077..HEAD` diff for this file — if any hunk is not GMG-specific, preserve it.
- [x] 1.6 Verify `MDLIB/ChemicalRoutines.c` still contains the `#ifdef _OMP_` / `#else` guard around `cholmod` `nthreads_max` assignment (preserve; do not revert to pre-`cab45077` here).
- [x] 1.7 Verify `MDLIB/Main_DOODZ.c` still uses `#ifdef _OMP_` around `InterpBufPool` thread count (preserve).
- [x] 1.8 Verify `MDLIB/mdoodz-private.h` retains `#ifndef MDOODZ_PRIVATE_H` … `#endif` guards; remove only GMG-related declarations if any remain outside the guard block.
- [x] 1.9 Run `git diff cab45077..HEAD -- MDLIB/` after edits and confirm surviving diffs against `cab45077` are exactly the three non-GMG fixes in D1 plus any audited non-GMG `InputOutput.c` hunks.

## 2. Tests — remove GMG fixtures and registrations

- [x] 2.1 Delete GMG-specific sources and fixtures listed in `proposal.md` under “Test removals” (`TESTS/GmgStokesEquivalence.cpp`, `Multigrid*.cpp`, `DisabledBenchmarks/*Gmg*`, `UplegAmplificationBound.cpp`, `GmgLevelsSweepProbe.cpp`, analytical/shear/SolVi `*_gmg*.txt`, `SolViUpleg161.txt`, etc.).
- [x] 2.2 Edit `TESTS/CMakeLists.txt` to drop all `add_test` / target registrations for removed files (restore to `cab45077`-equivalent structure for non-GM G tests).
- [x] 2.3 Confirm `TESTS/RotationAdvectionTests.cpp` still matches `887e8bc9` (rotation fix); do not revert this file unless a merge conflict forces a manual re-application of that patch.

## 3. Benchmarks and perf driver — delete trees and CMake hooks

- [x] 3.1 Delete the tracked tree `benchmarks/ec2/` (scripts, Python, YAML, README — not gitignored `benchmarks/ec2/results/` unless the user explicitly deletes that directory manually).
- [x] 3.2 Delete the tracked tree `benchmarks/vcycle_vis/`.
- [x] 3.3 Delete `SETS/SolViPerf.c`, `SETS/SolViPerf.txt`, and `cmake-exec/SolViPerf/` (or equivalent) if present.
- [x] 3.4 Edit `SETS/CMakeLists.txt` to remove `add_set(SolViPerf)` and any SolViPerf-only comments that reference GMG perf work.

## 4. OpenSpec — tombstones and archive hygiene

- [x] 4.1 Replace `openspec/specs/gmg-stokes-solver/spec.md` with the retirement tombstone (content per delta spec and `design.md` D5).
- [x] 4.2 Replace `openspec/specs/gmg-defence-materials/spec.md` with its tombstone.
- [x] 4.3 Replace `openspec/specs/gmg-performance-harness/spec.md` with its tombstone.
- [x] 4.4 Delete non-canonical files from archived GMG changes per `design.md` D6 (`STATUS.md`, `RETROSPECTIVE.md`, `DECISIONS.md`, `defence.{md,html}`, `build_pdf.py`, `HANDOFF.md` where listed, stray `.DS_Store` under `openspec/changes/archive/**`).
- [x] 4.5 Leave intact in each archive: `proposal.md`, `design.md`, `tasks.md`, `specs/`, `.openspec.yaml`.

## 5. Build and runtime verification

- [x] 5.1 Configure and build with OpenMP enabled (default project configuration); fix compile/link errors until clean.
- [x] 5.2 Configure and build with OpenMP disabled (`cmake` flag per project convention, e.g. `-DOMP=OFF` or equivalent); confirm non-OMP build succeeds (validates preserved OMP guards).
- [x] 5.3 Run the non-GM G test suite (or full CI-equivalent subset) and fix failures until green.
- [x] 5.4 Manual: run with a scratch `.txt` setting `lin_solver = 3` — expect immediate `LOG_ERR` and non-zero exit (D3).
- [x] 5.5 Manual: run with `lin_solver = 0` and a few `gmg_*` keys present — expect single aggregated `WARN`, run continues (D4).
- [x] 5.6 Optional smoke: one RiftingBasic or SolVi step with `lin_solver = 0` to confirm Stokes path unchanged numerically (design migration plan step 6).

## 6. Finish

- [x] 6.1 Single logical commit with message referencing `remove-gmg-stokes-solver`, `proposal.md`, and `design.md` (per `design.md` D7).
- [x] 6.2 `openspec verify` is not available in this CLI; implementation validated via build + ctest + manual lin_solver / deprecated gmg key checks. Archive with `/opsx:archive` when ready.
