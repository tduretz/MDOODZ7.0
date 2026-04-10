---
name: skill-logging
description: MDOODZ logging system — LOG_INFO/WARN/ERR/DBG/TIME macros, .txt configuration parameters (log_dest, log_level, log_timestamp, log_ts_mode, log_metadata), colored console output, ANSI-stripped log files, phase timing instrumentation, and metadata prefixes.
---

# MDOODZ Logging System

## Overview

All console/file output in MDLIB uses structured logging macros defined in `MDLIB/include/mdoodz-log.h`, implemented in `MDLIB/MdoodzLog.c`. The logger supports configurable destinations, log levels, timestamps, and per-step/iteration metadata.

## Log Macros

Use these instead of `printf()`:

| Macro | Level | Console Color | When to use |
|-------|-------|---------------|-------------|
| `LOG_ERR(fmt, ...)` | ERROR (0) | Bold red | Fatal errors before `exit()` |
| `LOG_WARN(fmt, ...)` | WARN (1) | Yellow | Non-fatal warnings, missing params |
| `LOG_TIME(fmt, ...)` | TIMING (4) | Green | Phase timing measurements |
| `LOG_INFO(fmt, ...)` | INFO (2) | Default | Normal progress messages |
| `LOG_DBG(fmt, ...)` | DEBUG (3) | Cyan | Verbose diagnostics |

**Important:** Do not append `\n` to format strings — the logger adds newlines automatically.

## .txt Configuration Parameters

Add these to the `OUTPUT FILES` section of any scenario `.txt` file:

```
log_dest      = 2   / 0=console, 1=file, 2=both (default: 2)
log_level     = 2   / 0=error, 1=warn, 2=info, 3=debug (default: 2)
log_timestamp = 1   / 0=off, 1=on (default: 1)
log_ts_mode   = 0   / 0=relative, 1=absolute, 2=both (default: 0)
log_metadata  = 0   / 0=off, 1=on (default: 0)
```

All parameters are optional — the defaults above are used when omitted.

## Output Format

The log line format is:

```
[<timestamp>|<metadata>] <LEVEL> | <message>
```

### Timestamp Modes

| `log_ts_mode` | Example prefix |
|---------------|----------------|
| 0 (relative) | `[    5.620]` |
| 1 (absolute) | `[04:55:10.857]` |
| 2 (both) | `[04:55:10.857 +    5.620]` |

Relative = seconds since simulation start (`omp_get_wtime()`). Absolute = wall clock via `gettimeofday()` with millisecond precision.

### Metadata

When `log_metadata = 1`, the step and nonlinear iteration are appended:

```
[    6.565|S0001|N00] INFO  | *** Picard it. 00 of 01 ...
```

- `S0001` = timestep number
- `N00` = nonlinear (Newton/Picard) iteration number
- Metadata only appears when values have been set (not during init)

## Log File

- Default path: `mdoodz.log` in the execution directory
- ANSI escape codes are automatically stripped from file output
- When `log_dest = 0` (console only), no log file is created

## Lifecycle (for developers)

The logger is managed in `Main_DOODZ.c`:

```c
mdoodz_log_init(NULL);                        // before ReadInputFile
mdoodz_log_reconfigure(&inputFile.model.log); // after ReadInputFile
// ... simulation loop ...
mdoodz_log_set_step(model.step);              // top of timestep loop
mdoodz_log_set_iteration(nit);                // top of nonlinear iteration
mdoodz_log_clear_iteration();                 // after nonlinear loop
mdoodz_log_flush();                           // end of each timestep
mdoodz_log_shutdown();                        // end of RunMDOODZ
```

## Choosing the Right Level

- **LOG_ERR**: Condition makes it impossible to continue → always followed by `exit()`
- **LOG_WARN**: Something unexpected but recoverable (missing parameter, convergence difficulty)
- **LOG_TIME**: Timing measurements (`omp_get_wtime()` deltas for phases)
- **LOG_INFO**: Normal progress (step start, residuals, solver status)
- **LOG_DBG**: Verbose output only needed during development/debugging

## Phase Timing

Key simulation phases are instrumented with `LOG_TIME` in the main loop:

- Particle interpolation
- Rheology update (`UpdateNonLinearity`)
- Stokes assembly
- Direct solve (CHOLMOD factorization + back-substitution)
- Thermal solver
- Advection
- HDF5 output

Per-timestep breakdown summaries are emitted at the end of each step.

## Performance CSV (`perf.csv`)

Every simulation automatically writes a `perf.csv` file in the execution directory with one row per timestep:

```csv
step,wall_s,rheology_s,assembly_s,solve_s,nit,n_particles,neq_mom,neq_cont,peak_rss_mb,user_cpu_s,sys_cpu_s
1,2.1240,0.0024,0.0229,0.0000,0,224692,27918,14006,220.8,40.9,2.4
```

| Column | Description |
|--------|-------------|
| `step` | Timestep number |
| `wall_s` | Wall-clock time for the timestep (seconds) |
| `rheology_s` | Cumulative rheology update time |
| `assembly_s` | Cumulative Stokes assembly time |
| `solve_s` | Cumulative direct solve time |
| `nit` | Number of nonlinear iterations |
| `n_particles` | Particle count |
| `neq_mom` | Momentum equation count |
| `neq_cont` | Continuity equation count |
| `peak_rss_mb` | Peak resident memory (MB, from `getrusage`) |
| `user_cpu_s` | Cumulative user CPU time (seconds) |
| `sys_cpu_s` | Cumulative system CPU time (seconds) |

This file is always written (no config parameter needed). Use it for scaling studies, regression checks, and capacity planning.

## Memory and CPU Monitoring

When `log_level = 3` (debug), per-timestep memory and CPU usage is also emitted to the log:

```
DEBUG | Memory: peak RSS = 268.3 MB, user CPU = 141.3 s, sys CPU = 5.2 s
```

This uses `getrusage(RUSAGE_SELF)` — zero overhead, no external tools needed.

## Source Files

- `MDLIB/include/mdoodz-log.h` — Public API, enums, macros
- `MDLIB/MdoodzLog.c` — Implementation (emit, init, reconfigure, shutdown)
- `MDLIB/InputOutput.c` — Config parsing (`ReadInt2` calls for `log_*` params)
- `MDLIB/Main_DOODZ.c` — Lifecycle calls and phase timing instrumentation
