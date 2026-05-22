## Context

The GMG-preconditioned FGMRES Stokes solver (`lin_solver = 3`) shipped with three archived OpenSpec changes (`add-gmg-stokes-solver`, `add-gmg-stokes-defence`, `add-gmg-upleg-fix`). Post-Fix F (`ad39d7f`) measurements in `benchmarks/ec2/results/PERFORMANCE_REPORT.md §§11–13` establish that current GMG is 1.8–3.4× slower than CHOLMOD on 2D SolVi and stagnates unrecoverably on 2D Rifting. Memory is not a constraint for the user and speed is the sole remaining optimisation axis, so the experiment is closed with a negative outcome — see the proposal for the full rationale.

This design document covers **how to remove the implementation cleanly**. The scope spans `MDLIB/` (revert), `TESTS/` (delete), `benchmarks/ec2/` and `benchmarks/vcycle_vis/` (delete), the `SolViPerf` driver (delete), and three OpenSpec specs (tombstone). The pre-GMG baseline is `cab45077c7f6ccaf447e3422b820a96d3c94169e` ("add new reseeding mode", 2026-04-17). A concurrent audit of the `cab45077..HEAD` MDLIB diff identified three **non-GMG fixes** that must be preserved across the revert — this constraint dominates the design.

Current state (as of HEAD ≈ `c59ce62`): 6 new GMG source files (`MultigridStokes.{c,h}`, `MultigridLevels.{c,h}`, `StokesAssemblyGMG.{c,h}`, ~4,880 LOC), 7 integration-point edits, ~20 GMG test fixtures, the `benchmarks/ec2/` harness (~2,200 LOC), and the `benchmarks/vcycle_vis/` diagnostic tooling. Stakeholders: the user as sole operator; no external consumers of `lin_solver = 3`.

## Goals / Non-Goals

**Goals:**

- Restore `MDLIB/` to a build that produces identical behaviour to `cab45077` on every pre-existing `lin_solver ∈ {-1, 0, 1, 2}` code path, without regressing any non-GMG fix made during the GMG era.
- Make `lin_solver = 3` a **hard startup error** (not a silent fallback) so legacy `.txt` files surface the removal immediately rather than mid-simulation.
- Shrink the maintenance surface: delete GMG-specific tests, the `benchmarks/ec2/` + `benchmarks/vcycle_vis/` trees, and the `SolViPerf` driver.
- Retire the three GMG OpenSpec specs with **tombstone pages** so archived-change cross-references still resolve but the retirement is unmistakable.
- Preserve the experimental record: every archived `proposal.md` / `design.md` / `tasks.md` / `specs/` stays intact; git history reaches the full implementation at `c59ce62`.
- Keep the changeset to **a single logical commit** (or a small handful of thematically-grouped commits) so the revert is easy to reason about and easy to undo with `git revert` if ever needed.

**Non-Goals:**

- Implement Galerkin coarsening, AMG, or any other multigrid variant — ruled out by user.
- Ship a recycled-CHOLMOD preconditioner or a block-Schur alternative — if ever needed, those are new proposals.
- Touch `lin_solver ∈ {-1, 0, 1, 2}`, the Powell-Hestenes iteration, the thermal PCG path, or the direct solver internals.
- Remove or rewrite `PERFORMANCE_REPORT.md` itself (it lives in the gitignored `benchmarks/ec2/results/` tree; the harness goes, the results stay on disk as an ad-hoc record if the user wants them).
- Remove UMFPACK as a dependency: pre-existing non-GMG usage in `Solvers.c`, `ThermalSolver.c`, and `StokesRoutines.c` makes it a required build dep regardless of GMG.
- Modify the `887e8bc9` rotation-test fix: it touches `TESTS/RotationAdvectionTests.cpp` only, not MDLIB, so the revert cannot disturb it.
- Rewrite git history (no `git reset --hard`, no force-push, no squash of archived commits). The removal is additive history on top of the branch.

## Decisions

### D1 — Surgical per-file removal, not `git checkout cab45077 -- MDLIB/`

**Decision**: Implement the MDLIB revert by removing GMG-specific code region-by-region inside each modified file, rather than by whole-file checkout of the pre-GMG blob.

**Rationale**: A line-by-line audit of `git diff cab45077..HEAD -- MDLIB/` surfaced three non-GMG fixes that must be preserved:

- **`MDLIB/ChemicalRoutines.c` (lines ~108–115)** — OMP guard added around `c.nthreads_max = omp_get_max_threads()`:

  ```c
  #ifdef _OMP_
      c.nthreads_max = (model.cholmod_threads == -1) ? omp_get_max_threads() : model.cholmod_threads;
  #else
      c.nthreads_max = (model.cholmod_threads == -1) ? 1 : model.cholmod_threads;
  #endif
  ```

  Fixes the non-OMP build. Unrelated to GMG. **Must stay.**

- **`MDLIB/Main_DOODZ.c` (lines ~102–110)** — parallel OMP guard for the `InterpBufPool` init path:

  ```c
  #ifdef _OMP_
      int nthreads_pool = omp_get_max_threads();
  #else
      int nthreads_pool = 1;
  #endif
  ```

  Same reason. **Must stay.**

- **`MDLIB/mdoodz-private.h` (lines 1–3 and last line)** — standard `#ifndef MDOODZ_PRIVATE_H ... #endif` header guards added. Pure code-hygiene fix, no GMG content. **Must stay.**

A whole-file `git checkout cab45077 -- MDLIB/ChemicalRoutines.c MDLIB/Main_DOODZ.c MDLIB/mdoodz-private.h` would silently regress these fixes and re-break non-OMP builds. Surgical edits are the only safe path.

**Alternatives considered**:

- **Option A — `git checkout cab45077 -- MDLIB/` then cherry-pick the OMP fixes back on top.** Rejected: cherry-pick adds two extra commits whose parentage is awkward (the original fix commits are part of the GMG-era history we want to leave visually distinct), and the three preserved fixes are small enough (~12 lines total) to redo inline without ceremony. Cherry-pick would also need careful hunk-level editing to drop the GMG context lines, defeating the simplicity argument.
- **Option B — `git revert` the GMG commits in reverse order.** Rejected: the GMG commits interleave with non-GMG work across `add-performance-metrics` branch history; a clean reverse revert chain does not exist. Also produces many small commits where one purposeful removal is clearer.
- **Option C — accept the OMP regression and backport later.** Rejected: the non-OMP build matters for users without OpenMP, and reintroducing a build regression inside a "cleanup" change is unacceptable.

**Per-file plan** (surgical):

| File | Action | Non-GMG content to preserve |
|---|---|---|
| `MDLIB/MultigridStokes.{c,h}` | DELETE | n/a (new file, 100% GMG) |
| `MDLIB/MultigridLevels.{c,h}` | DELETE | n/a |
| `MDLIB/StokesAssemblyGMG.{c,h}` | DELETE | n/a |
| `MDLIB/CMakeLists.txt` | surgical — drop `MultigridStokes.c`, `MultigridLevels.c`, `StokesAssemblyGMG.c` from the target source list | (verify no other added source deps) |
| `MDLIB/StokesRoutines.c` | surgical — drop `#include "MultigridStokes.h"` and the two `lin_solver == 3` dispatch blocks in `SolveStokesDecoupled` and `SolveStokesDefectDecoupled` | all other edits are zero per the diff |
| `MDLIB/InputOutput.c` | surgical — drop `gmg_*` parameter reads (~32 lines per diff stat; all appear to be GMG-specific, audit during implementation) | to be verified in T-implement-revert; any non-GMG edit surfaced gets preserved |
| `MDLIB/Main_DOODZ.c` | **keep OMP guards**, drop any GMG-specific dispatch cleanup if present (none visible in current diff) | `#ifdef _OMP_ ... #else ... #endif` around `nthreads_pool` |
| `MDLIB/ChemicalRoutines.c` | **keep entire diff** — the only change is the OMP guard, which is non-GMG | all 4 diff lines |
| `MDLIB/include/mdoodz.h` | surgical — drop the 9 `gmg_*` fields from the `params` struct (continuous block, ~10 lines) | no other diff content |
| `MDLIB/mdoodz-private.h` | **keep header guards**, drop any GMG function declarations (none visible in current diff beyond the guard and a blank-line diff) | `#ifndef MDOODZ_PRIVATE_H ... #endif` |

Post-revert, `MDLIB/ChemicalRoutines.c`, `MDLIB/Main_DOODZ.c`, and `MDLIB/mdoodz-private.h` will differ from `cab45077` by the preserved OMP / header-guard content. All other MDLIB files will match `cab45077` byte-for-byte after the surgical edits.

### D2 — UMFPACK stays

**Decision**: Do not remove UMFPACK from the build. Leave `MDLIB/CMakeLists.txt` UMFPACK linkage untouched.

**Rationale**: `rg 'umfpack|UMFPACK' MDLIB/` returns three pre-existing call sites (`Solvers.c`, `StokesRoutines.c` in non-GMG functions, `ThermalSolver.c`). These predate the GMG work and continue to depend on UMFPACK after the revert. Removing it would break the non-GMG linear solver path.

**Alternative considered**: Remove UMFPACK and assume its non-GMG uses are dead code. Rejected: dead-code audit is out of scope for this cleanup; shrinking build deps beyond the GMG footprint is a separate concern if it matters at all.

### D3 — `lin_solver = 3` becomes a hard startup error, not a silent fallback

**Decision**: When the input parser reads `lin_solver = 3`, emit a `LOG_ERR` with a clear message and abort immediately, before any simulation work. Message content:

```
ERROR: lin_solver = 3 (GMG-preconditioned FGMRES) has been removed.
       Valid values: -1, 0, 1, 2.
       See openspec/specs/gmg-stokes-solver/spec.md for the retirement notice
       and archived changes 2026-04-18-add-gmg-upleg-fix,
       2026-04-19-add-gmg-stokes-solver, 2026-04-21-add-gmg-stokes-defence
       for the historical implementation and its measured performance.
```

**Rationale**: The current code silently falls back to CHOLMOD when GMG fails to converge (per `StokesRoutines.c` line 790 block). Post-removal we could:

- (a) keep a silent fallback — user sets `lin_solver = 3`, gets CHOLMOD, nothing in their logs warns them their config is obsolete;
- (b) warn at startup but still fall through to CHOLMOD — user notices but might ignore it across many runs;
- (c) hard error at startup — user cannot miss it, must update their config.

Option (c) is correct because `lin_solver = 3` is an *uncommon* input value (the defence recorded zero production users outside the user's own experiments) and silent/soft failure would hide the removal from anyone returning to a pinned `.txt` months later. Fail-fast is the least-surprise option; the error message includes the migration path (switch to `lin_solver = 0` or `-1`) and the historical pointers. This aligns with the proposal's "no silent downgrade mid-simulation" stance.

**Alternative considered**: Keep a `lin_solver = 3 → warning + route to CHOLMOD` shim for one release cycle, then promote to hard error. Rejected: no external consumers, no release cycle, and the shim would require keeping a tiny GMG-dispatch-free alias block in `StokesRoutines.c` that is pure overhead. Direct removal is cleaner.

### D4 — `gmg_*` parameters: ignored with single startup WARN

**Decision**: After removing the 9 `gmg_*` fields from `params`, keep the input parser permissive: if an unknown key matching `^gmg_` is present in the `.txt`, skip it but emit one aggregated WARN line at the end of parameter parsing naming the `gmg_*` family as deprecated.

Concrete message:

```
WARN: ignoring deprecated gmg_* parameters in input file (GMG solver removed;
      see openspec/specs/gmg-stokes-solver/spec.md). Parameters ignored: gmg_levels,
      gmg_fgmres_tol, ... (N total)
```

**Rationale**: User reported that they may have pinned `.txt` files with `gmg_*` parameters set. Erroring on every `gmg_*` key would be noisy (one message per key) and would create a migration burden for benign config keeping (the `lin_solver` hard error in D3 already forces the user to update the active dispatch; the auxiliary `gmg_*` keys are harmless when ignored). Single aggregated WARN surfaces the deprecation without blocking.

Pattern matching is lenient (case-insensitive prefix `gmg_`) so the WARN covers any forgotten `gmg_*` keys without us enumerating them. If the input parser has no generic unknown-key handling today (likely — MDOODZ `.txt` format is free-form key=value), a minimal pattern-based filter is acceptable scope.

**Alternative considered**: Strip `gmg_*` matching from the parser entirely and let the unknown-key handler (if any) fire per-key warnings. Rejected: noise floor; also depends on `InputOutput.c` infrastructure we'd rather not expand during a cleanup.

### D5 — Spec retirement by tombstone

**Decision**: Replace the content of each retired spec with a tombstone page that:

1. States the capability is **RETIRED** (not merely deprecated).
2. Names the archived changes that historically defined it (with dates and short descriptions).
3. Points at the commit that removed the implementation (this change's merge commit, to be filled in at archive time).
4. Names the successor capability (if any — for these three, there is none: CHOLMOD already owns the Stokes solve).

Tombstone page structure is 15–25 lines, single screen. File names stay: `openspec/specs/gmg-stokes-solver/spec.md`, `openspec/specs/gmg-defence-materials/spec.md`, `openspec/specs/gmg-performance-harness/spec.md`.

**Rationale**: Three options were considered:

- **Delete the spec files entirely.** Cross-references from archived change `proposal.md`s (`add-gmg-stokes-solver/proposal.md` links to `openspec/specs/gmg-stokes-solver/spec.md`) would 404. Rejected.
- **Leave the current spec content and add a banner.** Rejected: the spec still reads like an active contract; casual readers would think GMG is still supported.
- **Tombstone the file with a retirement notice.** Keeps cross-refs resolving, communicates retirement unambiguously, and occupies <1 kB per file. Chosen.

Tombstone content template (one per capability):

```
# <capability-name> Specification

## Status

**RETIRED** (2026-04-21).
The capability described below has been removed from MDOODZ. See the
"Historical Record" section for pointers to the archived changes that
implemented and measured it.

## Purpose (historical)

<one-paragraph description of what the capability was>

## Retirement Rationale

<2-3 sentences: post-Fix F perf envelope, no speed advantage, user decision
to close the experiment. Link to PERFORMANCE_REPORT.md §§11-13.>

## Historical Record

The following archived OpenSpec changes contain the full requirements,
design, tasks, and spec deltas that defined and measured this capability:

- `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/`
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/`
- `openspec/changes/archive/2026-04-18-add-gmg-upleg-fix/`

The implementation itself was removed by change
`remove-gmg-stokes-solver` (this change). Git history before that commit
preserves the full source.

## Successor

None. CHOLMOD direct solve (`lin_solver = 0` or `-1`) remains the
production Stokes solver for all 2D MDOODZ workloads.
```

### D6 — Archived-change cleanup: what exactly is "not spec/proposal/design"

**Decision**: Remove the following non-canonical files from archived GMG changes; keep everything else:

**Keep** (canonical spec-driven artifacts):

- `proposal.md`, `design.md`, `tasks.md`, `specs/**/*.md` in every archived change.
- `.openspec.yaml` (workflow metadata, required by OpenSpec tooling).

**Delete** (process / status / retrospective / defence deliverables):

- `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/STATUS.md`
- `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/RETROSPECTIVE.md`
- `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/DECISIONS.md`
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/STATUS.md`
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/defence.md`
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/defence.html`
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/build_pdf.py`
- `openspec/changes/archive/2026-04-18-add-gmg-upleg-fix/STATUS.md` (and any sibling non-canonical files that tasks.md audit surfaces)
- `.DS_Store` files under `openspec/changes/archive/**` (macOS noise, not GMG-specific but natural to catch now)

**Rationale**: The user's directive was "clean up artifacts that are not spec, proposal, design etc" — the core spec-driven artifacts stay because they *are* the permanent record of what was tried and why it did not work. Everything else (status trackers, retrospectives, decision logs, the defence write-up and its build script) was working-material produced during implementation; it served its purpose and does not need to live forever in the archive.

`defence.md` is a judgment call: it *is* a written document with findings. But it is also a 30+ page self-contained narrative, most of whose content is reproduced in the archived change's `proposal.md` + `design.md` + `specs/` + the `PERFORMANCE_REPORT.md` on disk. Deleting keeps the archive focused on structured spec-driven content.

`.openspec.yaml` stays because OpenSpec tooling requires it to enumerate archived changes; removing it would break `openspec list --all` and similar commands.

### D7 — Commit shape: one logical cleanup commit

**Decision**: Land the removal as a single feature-branch commit titled roughly:

```
remove-gmg-stokes-solver: revert MDLIB, delete tests/benchmarks/driver, tombstone specs

Closes the GMG experiment with the negative-outcome decision documented in
openspec/changes/remove-gmg-stokes-solver/proposal.md. Preserves non-GMG
OMP-guard and header-guard fixes from the GMG era (per design D1).
```

No `git reset --hard`, no force-push, no rewriting of the three archived GMG commits. The GMG history remains visible in `git log` between `cab45077` and this commit.

**Rationale**:

- **One commit** keeps the revert atomic: a future `git revert <sha>` restores everything in one operation if the user ever changes their mind.
- **Additive history** (not a reset) preserves the experimental record at the git level, complementing the OpenSpec archive preservation.
- Alternative: split into 4-5 thematic commits (MDLIB, TESTS, benchmarks, specs, driver). Rejected: the removals are coupled (deleting `MDLIB/MultigridStokes.c` breaks tests that reference it; the tests and the source must move together) and the cleanup reads more clearly as a single purposeful commit than as a sequence.

### D8 — PERFORMANCE_REPORT.md lives on, harness goes

**Decision**: The `benchmarks/ec2/results/` tree (which contains `PERFORMANCE_REPORT.md`, CSVs, logs) is gitignored and stays on disk untouched. Only the **harness code** (`benchmarks/ec2/*.sh`, `benchmarks/ec2/*.py`, `benchmarks/ec2/grids*.yaml`) is deleted.

**Rationale**: `PERFORMANCE_REPORT.md` is the definitive negative-result artifact referenced from the proposal and the retirement-tombstone specs. It is useful as an on-disk record regardless of whether the harness that produced it is still maintained. Since it lives in a gitignored subtree, deleting the tracked harness files does not touch it.

If the user also wants to delete `benchmarks/ec2/results/` from disk, that is a manual `rm -rf` outside this change's scope (no git operation needed).

## Risks / Trade-offs

- **[Silent regression of a non-GMG fix during revert]** → **Mitigation**: D1 enumerates the three known non-GMG fixes explicitly and the per-file plan preserves them by construction. `tasks.md` requires `git diff cab45077..HEAD -- MDLIB/InputOutput.c` to be hand-reviewed during implementation to catch any fourth hidden non-GMG edit inside the 32-line diff. After implementation, run `git diff cab45077..HEAD -- MDLIB/` and confirm the only surviving hunks are the three documented in D1 plus any newly-identified non-GMG edit in `InputOutput.c`.

- **[`lin_solver = 3` hard error breaks a legacy run the user forgot about]** → **Mitigation**: the error message in D3 names the exact parameter and the migration path (`lin_solver = 0` or `-1`). Since the user is the sole operator and has just authored the proposal removing GMG, the likelihood of encountering this in the wild is low; when it happens it is the *intended* fail-fast behaviour.

- **[Archived-change `proposal.md` references `MDLIB/MultigridStokes.c` that no longer exists]** → **Accepted, not mitigated.** Archives describe history; the referenced files live in git history at the pre-revert SHA. Adding a "SUPERSEDED" banner across every archived change would bloat the archive with forward-pointer noise. The tombstone specs (D5) serve this role at the capability level.

- **[Non-OMP build regression if D1 mis-executed]** → **Mitigation**: tasks.md includes an explicit non-OMP build validation step (`cmake -DOMP=OFF ...`) after the revert, before the commit. Confirms the OMP guards from `ChemicalRoutines.c` and `Main_DOODZ.c` survived.

- **[UMFPACK wrongly removed]** → **Prevented by D2.** No action needed.

- **[`.DS_Store` sweep catches something unexpected]** → Low-risk housekeeping. If a `.DS_Store` somewhere is actually meaningful (extremely unlikely), restore it.

- **[Branch-level fork: a future user resumes the `add-performance-metrics-clean` branch to pull GMG back]** → **Supported, not a risk.** The branch remains unrewritten. Resuming it means checking out its head SHA and cherry-picking into a fresh branch; no special action needed from this change.

- **[Archived tasks.md checkboxes point to now-deleted files]** → **Accepted.** Same reasoning as archived-proposal references; archives are historical.

## Migration Plan

1. **Branch off current HEAD** (`c59ce62`) into a new branch named `remove-gmg-stokes-solver` (or continue on the active OpenSpec working branch, whichever the user prefers; the change does not dictate branching model).
2. **Apply per-file MDLIB revert** per D1. Build with OMP enabled and verify `make` succeeds.
3. **Build with OMP disabled** (`cmake -DOMP=OFF ...`) and verify `make` succeeds — protects the preserved OMP guards.
4. **Delete GMG test fixtures** and revert `TESTS/CMakeLists.txt`. Run the full CHOLMOD test suite (SolVi, SolKz, ShearBand, Rotation, etc.) and confirm all non-GMG tests pass.
5. **Delete `benchmarks/ec2/`, `benchmarks/vcycle_vis/`**, and the `SolViPerf` driver (`SETS/SolViPerf.{c,txt}`, `cmake-exec/SolViPerf/`). Verify root `CMakeLists.txt` no longer references `SolViPerf`.
6. **Run end-to-end on RiftingBasic** with `lin_solver = 0` and confirm the run reproduces pre-removal numerical output (one time-step comparison is sufficient; the direct-solver path is untouched by this change so bit-for-bit reproduction is expected).
7. **Test the D3 hard-error path**: set `lin_solver = 3` in a scratch `.txt` and confirm the run aborts with the D3 message.
8. **Test the D4 WARN path**: set `gmg_levels = 4` in a scratch `.txt` with `lin_solver = 0` and confirm the run emits the aggregated WARN and proceeds to completion.
9. **Tombstone the three GMG specs** per D5.
10. **Delete non-canonical files in archived GMG changes** per D6.
11. **Single-commit** per D7, with commit message referencing the proposal and this design.
12. **Archive this change** via `/opsx:archive` after the commit merges to the working branch, promoting `remove-gmg-stokes-solver` into `openspec/changes/archive/YYYY-MM-DD-remove-gmg-stokes-solver/`.

**Rollback strategy**: `git revert <this-change-sha>` restores MDLIB sources, tests, benchmarks, driver, and the three specs to their pre-removal content. The rollback is one command and is covered by routine git workflow.

## Open Questions

None remaining for this cycle. The audit during proposal drafting resolved the UMFPACK question (D2), the revert strategy (D1), and the degree of backward compatibility (D3, D4). Any fourth non-GMG edit discovered during the `InputOutput.c` audit in Step 2 above is handled locally by that step and does not require a design-level decision.
