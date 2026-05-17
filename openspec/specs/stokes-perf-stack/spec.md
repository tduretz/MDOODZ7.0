# stokes-perf-stack Specification

## Purpose

This capability describes a coherent set of four source-resident performance interventions on the MDOODZ7.0 Stokes non-linear hot path, plus two configuration recommendations, that together reduce wall-time on the `RiftingRoman` 400×210 `ani_fstrain=3` production workload by **57.8 %** vs the 4-thread reference baseline (10.72 → 4.525 s/step over a 60-step measurement window restarted from `Breakpoint00840.dat`).

The four interventions are independent in source but **stack non-linearly**:
1. Anisotropy bisection cap reduction (H06) — `MDLIB/AnisotropyRoutines.c`
2. SparseMat CSR sparsity-pattern cache (H02) — `MDLIB/mdoodz-private.h`, `MDLIB/MemoryAllocFree.c`, `MDLIB/StokesAssemblyDecoupled.c`
3. Predictor-corrector line search in `LineSearchDecoupled` (PC-LS) — `MDLIB/StokesRoutines.c`
4. Accelerate single-threaded BLAS via `dlsym` (BLASSetThreading) — `MDLIB/Main_DOODZ.c`

All four preserve physics within the solver's documented nonlinear tolerance band (mean relative diff `1e-4` on velocity, `1e-3` on `eta_n` at step 845 of the verification window) and have been validated on two distinct rifting states (Breakpoint 840 and 850).

Research provenance, measurement protocol, and full physics-diff reports are archived at `~/.claude/research/perf_engine_2026-05-16/`.

## Requirements

### Requirement: Anisotropy inner-solver bisection cap

The `ViscosityConciseAniso` inner non-linear visco-elastic solver SHALL bracket the effective viscosity `eta_ve` with **at most 10 geometric bisection halvings** before handing the root to Newton's method. The bracket SHALL initially span `[min_eta * 10^-10, max_eta * 10^2]` (the existing range; ~30 dex). After 10 halvings the bracket SHALL be ~3e-3 dex wide, comfortably inside Newton's quadratic-convergence basin for every cell exercised by the rifting workloads validated to date.

#### Scenario: Steady-state rifting cell with cached eta_ve nearby

- **WHEN** `ViscosityConciseAniso` is called on a cell whose previous-iter `eta_ve` is within ~3 dex of the current converged answer
- **THEN** the bisection loop SHALL terminate after at most 10 iterations and Newton's method SHALL converge to machine epsilon in ≤ 4 iterations

#### Scenario: Cold-lithosphere stiff cell

- **WHEN** `ViscosityConciseAniso` is called on a cold-lithosphere cell where the residual function has a steep curvature
- **THEN** the bisection loop SHALL still produce a bracket inside Newton's basin within 10 halvings; verified by 60-step long-run convergence on `RiftingRoman` ani_fstrain=3 with no cell-level Newton failure

#### Scenario: Bisection floor below 10

- **WHEN** `bisection_steps` is reduced below 10 (tested values 5 and 3)
- **THEN** the inner Newton SHALL fail to converge on a non-empty subset of cells and the simulation SHALL abort within the first 1–2 timesteps; therefore 10 is the documented safe floor

### Requirement: SparseMat sparsity-pattern cache

The `SparseMat` struct SHALL carry allocated-capacity counters `nnz_alloc` (entries in `J` and `A`) and `neq_alloc` (entries in `Ic`, `bbc`) that persist across non-linear iterations within a timestep. `AllocMat` SHALL reuse existing buffers via `memset(0)` when the requested size is ≤ allocated capacity, and SHALL only call `DoodzCalloc` to grow capacity. `FreeMat` SHALL be a no-op for cached buffers. `SFree` (called once at simulation teardown) SHALL release the underlying memory.

#### Scenario: Within-timestep non-linear iteration

- **WHEN** `AllocMat` is called for a `SparseMat` whose `nnz_alloc ≥ requested_nnz` and `neq_alloc ≥ requested_neq`
- **THEN** the existing `Ic`, `J`, `A`, `bbc` buffers SHALL be zeroed via `memset` and the buffers SHALL be reused; no `DoodzCalloc` / `DoodzFree` cycle SHALL occur

#### Scenario: Free-surface advection grows the pattern

- **WHEN** marker advection across cell boundaries causes a new `AllocMat` call with `requested_nnz > nnz_alloc`
- **THEN** the existing buffers SHALL be freed and re-allocated at the new (larger) size, and `nnz_alloc` / `neq_alloc` SHALL be updated

#### Scenario: `StokesAssembly_Decoupled` post-assembly shrink

- **WHEN** `BuildStokesOperatorDecoupled` completes assembly and the actual non-zero count is less than the conservative initial estimate
- **THEN** the over-allocated tail SHALL NOT be `Realloc`'d to shrink; `nnz_alloc` remains the truth and the over-allocation is harmless

#### Scenario: End-of-simulation teardown

- **WHEN** `SFree` is called on a `SparseMat` at simulation end
- **THEN** `Ic`, `J`, `A`, `bbc` SHALL be freed and `nnz_alloc` / `neq_alloc` SHALL be reset to 0

#### Scenario: Numerical equivalence with the pre-cache codebase

- **WHEN** any solver step is run with the cache active versus a baseline build with the cache disabled
- **THEN** the assembled `SparseMat` contents SHALL be byte-identical, and all downstream solver outputs (`/Centers/P`, `/Centers/eta_n`, `/VxNodes/Vx`) SHALL match to bit precision

### Requirement: Predictor-corrector line search in `LineSearchDecoupled`

The Newton-branch of `LineSearchDecoupled` SHALL test the full Newton step (`alpha = -1.0`) first via a single residual evaluation and accept it when the resulting combined residual `||(rx,rz,rp)||` is strictly less than `frac * ||(prev_rx, prev_rz, prev_rp)||`, where `frac` is the existing Armijo coefficient (`Nmodel->LineSearch_armijo`, default `1.0`). On acceptance the function SHALL skip the original 6-trial geometric bracket scan. On rejection the function SHALL fall through to the existing 6-trial scan with no degradation versus the pre-patch behaviour.

The patch SHALL be confined to the Newton branch; the Picard branch (`Newton == 0`) is unchanged.

#### Scenario: Newton step satisfies Armijo

- **WHEN** the full Newton step at `alpha = -1.0` produces a residual ratio `||f_new|| / ||f_prev|| < frac`
- **THEN** the line search SHALL accept `alpha = -1.0`, set `success = 1`, log "PC-LS ACCEPTED on first trial", and return without scanning the 6 trial alphas

#### Scenario: Newton step fails Armijo

- **WHEN** the full Newton step at `alpha = -1.0` does NOT satisfy `||f_new|| / ||f_prev|| < frac`
- **THEN** the line search SHALL set `pc_accepted = 0` and execute the original 6-trial bracket scan over `alpha ∈ {-1.0, -0.5, -0.25, -0.125, -0.0625, -0.03125}`, selecting the best alpha by minimum residual; semantically equivalent to the pre-patch behaviour for this branch

#### Scenario: Picard outer iteration (Newton == 0)

- **WHEN** `LineSearchDecoupled` is called with `Newton == 0` (Picard mode)
- **THEN** the predictor-corrector block SHALL NOT execute and the function SHALL run the original Picard line-search code path unchanged

#### Scenario: Aggressive Armijo coefficient

- **WHEN** `Nmodel->LineSearch_armijo` is reduced from the default `1.0` (e.g. to `0.9` or below)
- **THEN** the predictor SHALL reject more often and the function SHALL fall through to the 6-trial scan; net wall-time effect SHALL be neutral relative to the pre-patch code

#### Scenario: Convergence properties

- **WHEN** the predictor-corrector line search is active on the `RiftingRoman` ani_fstrain=3 workload from step 840
- **THEN** the mean outer Newton iteration count over a 60-step run SHALL decrease compared to the pre-patch code (measured: `nit` 3.68 → 2.35), because Armijo's longest-step rule yields a more effective per-iteration descent than the bracket-scan's argmin-of-trial selection

#### Scenario: Physics agreement with pre-patch run

- **WHEN** the patched solver is run against a baseline run with the predictor disabled at the same breakpoint
- **THEN** the difference in `/Centers/P` and `/VxNodes/Vx` SHALL have mean relative error ≤ `1e-3` and the difference in `/Centers/eta_n` SHALL have mean relative error ≤ `1e-3` (single shear-band-transition cells MAY show higher pointwise differences due to nonlinear branching, but the mean across the field SHALL remain within solver tolerance)

### Requirement: Accelerate BLAS single-threading via `dlsym`

On Apple platforms (`__APPLE__` defined at compile time) the simulation entry point `RunMDOODZ` SHALL attempt to resolve `BLASSetThreading` from `RTLD_DEFAULT` at runtime via `dlsym`. When the symbol resolves (macOS 15 "Sequoia" and later) the function SHALL be called with `BLAS_THREADING_SINGLE_THREADED` (value `1`) and the return code SHALL be logged via `LOG_INFO`. When the symbol does NOT resolve (older macOS, non-Apple platforms when `__APPLE__` is defined for compatibility reasons, or future symbol removal) the call site SHALL silently skip and SHALL NOT introduce a link-time dependency on `BLASSetThreading`.

#### Scenario: macOS 15 or later

- **WHEN** the binary is run on macOS 15+ (Accelerate.framework exports `BLASSetThreading`)
- **THEN** `dlsym` SHALL return a non-NULL function pointer, the call `BLASSetThreading(1)` SHALL succeed, and Accelerate vecLib SHALL execute subsequent BLAS calls (e.g. CHOLMOD's internal `dgemm_`) in single-threaded mode

#### Scenario: macOS 14 or earlier

- **WHEN** the binary is run on a macOS version that does not expose `BLASSetThreading`
- **THEN** `dlsym` SHALL return `NULL`, the call site SHALL skip without error, and BLAS threading SHALL fall back to its default behaviour (controlled by deprecated `VECLIB_MAXIMUM_THREADS` or implicit)

#### Scenario: Symbol absence does not break the build

- **WHEN** the codebase is built on any platform
- **THEN** linking SHALL succeed without requiring `BLASSetThreading` to be present in the linker's symbol search path

### Requirement: Recommended runtime and workload configuration for the full perf stack

To achieve the documented `4.525 s/step` headline on the validated workload, runtime and workload configuration SHALL be set per the following table. These knobs live outside the source tree but are an integral part of the validated performance result.

| Knob | Recommended value | Lives in |
|---|---|---|
| `OMP_NUM_THREADS` | `12` | environment (export before `RiftingRoman` invocation) |
| Powell-Hestenes `penalty` | `1e1` | workload `.txt` (e.g. `cmake-exec/RiftingRoman/run_aniso3_v3/RiftingRomanAniso3_v3.txt:28`) |
| Clang PGO | `-fprofile-instr-use=<profdata>` | CMake build flag (3-stage instrumented build + training run + use rebuild) |

#### Scenario: Production rifting run (Nx·Nz ≥ 10⁴)

- **WHEN** the user invokes `RiftingRoman` (or any large-grid Stokes scenario) for production
- **THEN** the documentation SHALL recommend `OMP_NUM_THREADS=12`, `penalty = 1e1`, and a PGO-built binary; deviating from any of these is permitted but SHALL be expected to forfeit a documented fraction of the win (OMP=4 forfeits ~38 %, `penalty = 1e2` forfeits ~0.9 %, non-PGO build forfeits ~1.1 %)

#### Scenario: Micro-scenario test (Nx·Nz < 10³)

- **WHEN** the user runs a sub-1000-cell test (e.g. lab-validation, MWE, `AniFstrainXxx` torsion suite)
- **THEN** the documentation SHALL recommend `OMP_NUM_THREADS=4` and the default `penalty`; the per-NL-iter overhead of H02's CSR-cache bookkeeping and the line-search predictor's residual probe become measurable on sub-millisecond per-step workloads and may regress wall-time by a few percent

#### Scenario: Workload `.txt` carries an explicit `penalty` value

- **WHEN** a workload `.txt` specifies `penalty` (any value)
- **THEN** that value SHALL be honoured by the solver as-is; the documentation SHALL note that `1e1` is the empirically tuned value for `RiftingRoman` ani_fstrain=3 but may not transfer to other rheologies (lower γ destabilises some configurations; higher γ regresses wall by ~3-5 %)

### Requirement: Reproducibility provenance

The performance research that produced this stack SHALL be retained at `~/.claude/research/perf_engine_2026-05-16/` (outside the repository, in the shared Claude library) and SHALL include:

- 60-step `perf.csv` measurement files for the four headline configurations (`vS3_now2`, `v26_pgo_g10`, `v28_pcls`, `v29_pcls_strict`)
- The measurement harness `run_variant.sh`
- The summary documents `SUMMARY.md`, `FINAL_REPORT.md`, and per-session `notes/0X_*.md`

#### Scenario: Re-verifying the headline number

- **WHEN** a future developer wants to confirm that the stack still delivers `4.525 s/step` on the validation workload
- **THEN** they SHALL invoke `~/.claude/research/perf_engine_2026-05-16/run_variant.sh <name> <build_dir> 900` against a build of the patched MDLIB and compare the resulting `cmake-exec/perf_eng/<name>/perf.csv` to the archived `v28_pcls/perf.csv`

#### Scenario: Auditing a single intervention

- **WHEN** a future developer wants to isolate the wall-time contribution of one of the four interventions
- **THEN** the archived `hypotheses/H*.md` files plus the per-session notes SHALL provide the experimental protocol and the per-phase measurement breakdown to enable a controlled re-test
