## Why

MDOODZ currently uses raw `printf()` calls (~500+ across MDLIB) with no timestamps, no log levels, and no file output. This makes it impossible to systematically profile where wall-clock time is spent during a simulation. Before optimizing OpenMP parallelization (fork-join overhead, memory bandwidth), we need instrumented logging that can produce `.log` files with timing data for each major computation phase. The existing partial `omp_get_wtime()` timing covers only a few phases and prints to stdout interleaved with diagnostic output.

## What Changes

- Add a lightweight logging utility (`mdoodz_log`) that wraps `fprintf` with wall-clock timestamps and log-level prefixes
- Timestamps are optional but enabled by default; can be toggled off at initialization
- Logging output is configurable: file only, console only, or both (default: both file and console)
- Log file is created automatically per simulation run (`.log` extension)
- **BREAKING**: All existing `printf()` calls across all MDLIB source files are migrated to the new logging API — no raw `printf()` remains in the codebase. Console output format changes (gains timestamp + level prefix when timestamps are enabled)
- Add structured timing instrumentation around each major phase in the main time-stepping loop: particle interpolation, rheology evaluation, Stokes assembly, direct solver, thermal solver, advection, I/O
- Color codes (ANSI escapes) are written to console output only, stripped from log file output

## Capabilities

### New Capabilities
- `log-api`: Core logging API — initialization, shutdown, log-level filtering, fprintf wrapper with optional timestamps (default: on) supporting relative, absolute, or both modes, optional metadata context (timestep/iteration) in prefix (default: off), configurable output destination (file/console/both), ANSI color stripping for file output. All existing `printf()` calls across MDLIB are migrated to this API
- `phase-timing`: Wall-clock timing instrumentation for major computation phases in the main loop, per-timestep and per-iteration timing summaries

### Modified Capabilities

_(none — no existing spec-level requirements change)_

## Impact

- **Code**: All `.c` files in `MDLIB/` — every `printf()` call is replaced with `mdoodz_log()` or a level-specific macro
- **Build**: New source files added to `MDLIB/CMakeLists.txt`
- **Public API**: New header `mdoodz-log.h` exposed via `include/`; `mdoodz.h` unchanged
- **Dependencies**: None — uses only `stdio.h`, `time.h`, `omp_get_wtime()`
- **Output**: Simulations now produce a `.log` file alongside HDF5 output
- **Performance**: Negligible — buffered file I/O adds ~microseconds per write vs existing terminal printf
