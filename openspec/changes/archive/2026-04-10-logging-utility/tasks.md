## 1. Logger Core Implementation

- [x] 1.1 Create `MDLIB/include/mdoodz-log.h` with log level enum, `MdoodzTimestampMode` enum, `MdoodzLogConfig` struct, `LOG_INFO`/`LOG_WARN`/`LOG_ERR`/`LOG_DBG`/`LOG_TIME` macros, and function declarations (`mdoodz_log_init`, `mdoodz_log_reconfigure`, `mdoodz_log_shutdown`, `mdoodz_log_flush`, `mdoodz_log_emit`, `mdoodz_log_set_step`, `mdoodz_log_set_iteration`, `mdoodz_log_clear_iteration`)
- [x] 1.2 Create `MDLIB/MdoodzLog.c` with global `g_logger` static struct, `mdoodz_log_init()` (open file, store t0, set defaults), `mdoodz_log_shutdown()` (flush + fclose), `mdoodz_log_flush()`
- [x] 1.3 Implement `mdoodz_log_emit()` — format prefix (timestamp + metadata + level), write to console with ANSI colors, write to file with ANSI stripped, auto-append newline
- [x] 1.4 Implement ANSI escape stripping state machine for file output (skip `\033[...m` sequences)
- [x] 1.5 Implement timestamp formatting — relative mode (`omp_get_wtime() - t0`), absolute mode (`localtime_r` + ms), both mode (absolute + `+relative`)
- [x] 1.6 Implement metadata prefix — `|S%04d` for step, `|N%02d` for iteration, only when `show_metadata=1` and values are set
- [x] 1.7 Implement `mdoodz_log_reconfigure()` — update g_logger fields from config, reopen log file if path changed
- [x] 1.8 Implement metadata setters: `mdoodz_log_set_step()`, `mdoodz_log_set_iteration()`, `mdoodz_log_clear_iteration()`

## 2. Build Integration

- [x] 2.1 Add `MdoodzLog.c` to the source list in `MDLIB/CMakeLists.txt`
- [x] 2.2 Add `mdoodz-log.h` as a public header in `MDLIB/CMakeLists.txt`
- [x] 2.3 Verify the project builds with the new files (no existing code changed yet)

## 3. Configuration via .txt Parameter File

- [x] 3.1 Add `MdoodzLogConfig log` field to the `params` struct in `MDLIB/include/mdoodz.h`
- [x] 3.2 Add `ReadInt2()` calls in `ReadInputFile()` in `InputOutput.c` for `log_dest` (default 2), `log_level` (default 2), `log_timestamp` (default 1), `log_ts_mode` (default 0), `log_metadata` (default 0)
- [x] 3.3 Populate the `model.log` config struct from the parsed values

## 4. Logger Lifecycle in Main Loop

- [x] 4.1 In `RunMDOODZ()` (`Main_DOODZ.c`): call `mdoodz_log_init(NULL)` before `ReadInputFile()`, call `mdoodz_log_reconfigure(&input.model.log)` after
- [x] 4.2 Call `mdoodz_log_shutdown()` at the end of `RunMDOODZ()`
- [x] 4.3 Add `mdoodz_log_set_step(model.step)` at the top of the timestep loop
- [x] 4.4 Add `mdoodz_log_set_iteration(Nmodel.nit)` at the top of the nonlinear iteration loop
- [x] 4.5 Add `mdoodz_log_clear_iteration()` after the nonlinear iteration loop ends
- [x] 4.6 Add `mdoodz_log_flush()` at the end of each timestep

## 5. Migrate printf Calls — Core Files

- [x] 5.1 Migrate `Main_DOODZ.c` — replace all `printf()` with appropriate `LOG_*` macros, remove trailing `\n` from format strings
- [x] 5.2 Migrate `InputOutput.c` — replace `printf()` and `fprintf(stderr, ...)` calls
- [x] 5.3 Migrate `Solvers.c` — replace `printf()` calls, convert existing timing prints to `LOG_TIME`
- [x] 5.4 Migrate `StokesRoutines.c` — replace `printf()` calls, convert timing prints to `LOG_TIME`
- [x] 5.5 Migrate `StokesAssemblyDecoupled.c` — replace `printf()` calls
- [x] 5.6 Migrate `ThermalRoutines.c` and `ThermalSolver.c` — replace `printf()` calls, convert timing prints to `LOG_TIME`

## 6. Migrate printf Calls — Rheology & Particles

- [x] 6.1 Migrate `RheologyDensity.c` — replace `printf()` calls
- [x] 6.2 Migrate `RheologyParticles.c` — replace `printf()` calls
- [x] 6.3 Migrate `ParticleRoutines.c` — replace `printf()` calls
- [x] 6.4 Migrate `ParticleReseeding.c` — replace `printf()` calls
- [x] 6.5 Migrate `FlowLaws.c` — replace `printf()` calls

## 7. Migrate printf Calls — Remaining Files

- [x] 7.1 Migrate `MiscFunctions.c` — replace `printf()` calls (including inside `MinMaxArray`/`MinMaxArrayTag`)
- [x] 7.2 Migrate `AdvectionRoutines.c` — replace `printf()` calls, convert `if (quiet==0)` guards to `LOG_DBG`
- [x] 7.3 Migrate `FreeSurface.c` — replace `printf()` calls
- [x] 7.4 Migrate `GridRoutines.c` — replace `printf()` calls
- [x] 7.5 Migrate `AnisotropyRoutines.c` — replace `printf()` calls
- [x] 7.6 Migrate `ChemicalRoutines.c` — replace `printf()` calls
- [x] 7.7 Migrate `MeltingRoutines.c` — replace `printf()` calls
- [x] 7.8 Migrate `MemoryAllocFree.c` — replace `printf()` calls
- [x] 7.9 Migrate `SparseTools.c` — replace `printf()` calls
- [x] 7.10 Migrate `Setup.c` — replace `printf()` calls
- [x] 7.11 Migrate `HDF5Output.c` — replace `printf()` calls
- [x] 7.12 Migrate `FD_Jacobian.c` — N/A (no printf calls in file)

## 8. Phase Timing Instrumentation

- [x] 8.1 Add timing around particle interpolation phase (`P2Mastah` block) in `Main_DOODZ.c`
- [x] 8.2 Add timing around `UpdateNonLinearity` + `RheologicalOperators` in the nonlinear iteration loop
- [x] 8.3 Add timing around `BuildStokesOperatorDecoupled` + `BuildJacobianOperatorDecoupled`
- [x] 8.4 Add timing around direct solver (CHOLMOD/UMFPACK solve block)
- [x] 8.5 Add timing around thermal solver phase
- [x] 8.6 Add timing around advection phase (`RogerGunther` and related calls)
- [x] 8.7 Add timing around HDF5 output (`WriteOutputHDF5`, `WriteOutputHDF5Particles`)
- [x] 8.8 Add per-iteration timing summary at end of nonlinear iteration
- [x] 8.9 Add per-timestep timing summary at end of timestep (total + breakdown)

## 9. Verification

- [x] 9.1 Verify build compiles cleanly with all changes (WSL Ubuntu, cmake 3.22)
- [x] 9.2 Run `grep -r "printf(" MDLIB/*.c` and confirm zero active matches outside `MdoodzLog.c` and `omptest.c`
- [x] 9.3 Run RiftingBasic scenario and verify `.log` file is created with correct format
- [x] 9.4 Verify log output contains timing summaries for all instrumented phases
- [x] 9.5 Test with `log_timestamp = 0` — verify no timestamps in output
- [x] 9.6 Test with `log_metadata = 1` — verify step/iteration in prefix
- [x] 9.7 Test with `log_dest = 0` (console only) — verify no log file created
- [ ] 9.8 Verify existing unit tests still pass (deferred to CI — WSL build too slow)

## 10. Polish

- [x] 10.1 Add colored output per log level (ERROR=bold red, WARN=yellow, DEBUG=cyan, TIMING=green)
- [x] 10.2 Fix absolute timestamps to use `gettimeofday()` for proper millisecond precision
- [x] 10.3 Audit all LOG_INFO calls — reclassify to LOG_ERR (fatal/exit) or LOG_WARN (warnings)
- [x] 10.4 Remove empty log file on disk when `log_dest=0` (console only) via `remove()` in reconfigure
