# Comprehensive Performance Test — Handoff Document

## Git State

- **Branch**: `performance-tests`
- **Last commit**: `44e5793 add benchmarks`
- **Uncommitted changes** in 5 files (all changes are unstaged/uncommitted):
  - `MDLIB/include/mdoodz.h` — added `t_end` and `L0` restored on struct field line, added `time_unit` field
  - `MDLIB/InputOutput.c` — added `ReadDou2(fin, "t_end", 0.0)/scaling.t` and `ReadInt2(fin, "time_unit", 0)`
  - `MDLIB/Main_DOODZ.c` — time-based loop, 8 subsystem timing accumulators, expanded perf.csv
  - `MDLIB/MdoodzLog.c` — added `model_time`/`time_unit` fields, `mdoodz_log_set_model_time()`, metadata prefix with `|T<value><unit>`
  - `MDLIB/include/mdoodz-log.h` — declared `mdoodz_log_set_model_time()`

**IMPORTANT**: Commit these changes before continuing. Run:
```bash
git add -A && git commit -m "perf.csv expansion: subsystem timing, time_ma, log metadata"
```

## Build Status

The last build had a **stale cache issue** — a previous edit left `dt_anisotropy` in the pre-loop init section, which was reverted in the source file but the build system didn't pick up the change. Fix:
```bash
make clean && make build
```
This should compile cleanly. If it still shows `'dt_anisotropy' undeclared` at line ~344, verify that `MDLIB/Main_DOODZ.c` line ~344 does NOT have `dt_anisotropy` (it should be plain `UpdateAnisoFactor` without timing wrapper in the init section).

## Completed Tasks (Group 1: 1.1–1.15)

All core C instrumentation is done:

| Task | Description | Status |
|------|-------------|--------|
| 1.1 | `t_end` field in model struct + ReadDou2 in InputOutput.c | DONE |
| 1.2 | `time_unit` field (0=Ma,1=Ka,2=yr) + ReadInt2 | DONE |
| 1.3 | Loop condition: `t_end > 0 ? time < t_end : step <= Nt` | DONE |
| 1.4 | 8 timing accumulators declared at top of loop body | DONE |
| 1.5 | `dt_thermal` from existing `t_omp` around EnergyDirectSolve | DONE |
| 1.6 | `dt_advection` captured after advection LOG_TIME | DONE |
| 1.7 | `dt_free_surface` wrapping both free_surface blocks in advection substep | DONE |
| 1.8 | `dt_reseeding` captured after CountPartCell LOG_TIME | DONE |
| 1.9 | `dt_melting` wrapping both MeltFractionGrid+UpdateAlphaCp/UpdateDensity blocks | DONE |
| 1.10 | `dt_anisotropy` wrapping UpdateAnisoFactor (line ~606) and NonNewtonianViscosityGridAniso in chemical section (~1054). **NOT** in pre-loop init section (lines ~344-349) | DONE |
| 1.11 | `dt_gse` wrapping both UpdateParticleGrainSize calls | DONE |
| 1.12 | `dt_output` wrapping HDF5 write + breakpoint block | DONE |
| 1.13 | perf.csv header: `step,wall_s,time_ma,rheology_s,assembly_s,solve_s,thermal_s,advection_s,free_surface_s,reseeding_s,melting_s,anisotropy_s,gse_s,output_s,nit,n_particles,neq_mom,neq_cont,peak_rss_mb,user_cpu_s,sys_cpu_s` | DONE |
| 1.14 | `time_ma = model.time * scaling.t / (1e6 * 365.25 * 24 * 3600)` computed in perf.csv fprintf | DONE |
| 1.15 | `mdoodz_log_set_model_time(time_s, time_unit)` — new function, struct fields, metadata prefix `\|T%.2f<unit>` | DONE |

## Next Task: 1.16 — Build and Verify

1. `make clean && make build` — must compile with zero errors
2. Run RiftingBasic for 3 steps to confirm perf.csv has the new 21-column header
3. Verify disabled subsystems (melting, anisotropy, GSE) show 0.0000 in their columns

## Remaining Tasks (Groups 2–8)

See `openspec/changes/comprehensive-perf-test/tasks.md` for the full checklist. Summary:

### Group 2: RiftingComprehensive Scenario (tasks 2.1–2.10)
- Merge `origin/add-meltin-vp-model` branch
- Fork RiftingCombinedYield → RiftingComprehensive (.c + .txt)
- Add anisotropy on mantle phases (IDs 1,2,3,5), GSE on phase 1
- Create 4 resolution .txt files: `_lowres` (100×80), default (200×160), `_medres` (501×401), `_highres` (1000×800)
- Register in SETS/CMakeLists.txt

### Group 3: Model Stability Validation (tasks 3.1–3.6)
- Run each resolution for 3 steps, verify no crash/NaN
- Diagnose and fix convergence issues if any
- **This is a gate**: do not proceed to groups 4+ until all resolutions pass

### Group 4: Benchmark Pipeline (tasks 4.1–4.7)
- Add `--scenario`, `--resolutions`, `--validate` flags to `misc/benchmark.sh`
- Update `misc/benchmark-report.sh` for new perf.csv columns

### Group 5: AWS Preflight (tasks 5.1–5.3)
- Create `misc/benchmark-preflight.sh`
- **This is a gate**: AWS tasks blocked until preflight passes

### Group 6: AWS Automation (tasks 6.1–6.8)
- `misc/benchmark-aws.sh` + `misc/benchmark-ec2-setup.sh`
- EC2 start/stop lifecycle, rsync, remote exec, S3 upload

### Group 7: Skill Documentation (tasks 7.1–7.3)
### Group 8: E2E Verification (tasks 8.1–8.4)

## Key Design Decisions (Summary)

- **perf.csv**: 21 columns total (was 12). Header in design.md Decision 3
- **t_end**: Read in seconds (like dt), non-dimensionalized on read. Default 0 = use Nt
- **time_unit**: 0=Ma (default), 1=Ka, 2=yr. Affects log metadata prefix and perf.csv column
- **Anisotropy timing**: Only instrument calls INSIDE the time loop. The pre-loop "Initialize viscosity" section (~lines 340-349) has the same functions but must NOT be instrumented (they run before timing accumulators are declared)
- **Resolution tiers**: Separate .txt files, not runtime patching
- **Validation before benchmarking**: 3 steps per resolution, lowest→highest, bail on first failure
- **AWS preflight as hard gate**: Don't implement AWS scripts until connectivity verified

## Key Files

| File | Role |
|------|------|
| `openspec/changes/comprehensive-perf-test/design.md` | Full design document with 8 decisions |
| `openspec/changes/comprehensive-perf-test/tasks.md` | 58-task checklist with completion status |
| `openspec/changes/comprehensive-perf-test/proposal.md` | Original proposal |
| `openspec/changes/comprehensive-perf-test/specs/` | 5 spec files |
| `MDLIB/Main_DOODZ.c` | Main simulation loop — all timing instrumentation lives here |
| `MDLIB/MdoodzLog.c` | Logger — model time metadata |
| `MDLIB/include/mdoodz.h` | Model struct — `t_end`, `time_unit` fields |
| `MDLIB/InputOutput.c` | Parameter reading — `t_end`, `time_unit` |
| `misc/benchmark.sh` | Benchmark pipeline (needs --scenario, --resolutions, --validate) |
| `misc/benchmark-report.sh` | Report generator (needs update for new columns) |

## Build Environment

- **macOS**: Use `make build` from repo root. Dependencies: SuiteSparse, HDF5, BLAS/LAPACK via Homebrew
- CMake build directory: `cmake-build/`
- Default scenario: controlled by `-DSET=<name>` in cmake (set in `makefile`)
- The `env.cmake` file may need adjusting for macOS paths (see `env.cmake.example`)
