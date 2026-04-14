## Context

MDOODZ uses ~200+ raw `printf()` calls across 24 `.c` files in MDLIB. Output goes exclusively to stdout (3 `fprintf(stderr, ...)` for fatal errors). Some calls use ANSI color macros (`GREEN`, `RED`, `RESET`), some are gated behind runtime flags (`noisy`, `quiet`), and a handful of timing measurements use `omp_get_wtime()` with ad-hoc `printf` of elapsed time. There is no log file, no timestamps, and no unified output control.

The codebase is C11, builds with CMake, and conditionally uses OpenMP via `#ifdef _OMP_`. Every source file that uses OpenMP already includes a fallback block defining `omp_get_wtime()` as `clock()/CLOCKS_PER_SEC`.

## Goals / Non-Goals

**Goals:**
- Single logging API that replaces all `printf()`/`fprintf(stderr, ...)` calls in MDLIB
- Configurable output destination: file, console, or both (default: both)
- Optional timestamps on every line (default: on, togglable at init) with configurable mode: absolute, relative, or both
- Optional metadata context in log prefix (timestep, iteration) — off by default, togglable at init
- Log levels (ERROR, WARN, INFO, DEBUG, TIMING) with runtime filtering
- ANSI color codes rendered on console, stripped in file output
- Thread-safe for OpenMP parallel regions
- Zero external dependencies beyond libc and OpenMP

**Non-Goals:**
- Structured/machine-parseable log format (JSON, CSV) — plain text is sufficient
- Log rotation or size limits — simulations produce finite output
- Remote logging, syslog, or network output
- Changing the semantic content of any existing message — only the delivery mechanism changes
- Async logging / ring buffers — the volume (~200 calls per timestep) doesn't warrant it

## Decisions

### 1. API shape: variadic macro wrapping a single emit function

```c
// Public macros (mdoodz-log.h)
#define LOG_INFO(fmt, ...)    mdoodz_log_emit(MDOODZ_LOG_INFO,  __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_WARN(fmt, ...)    mdoodz_log_emit(MDOODZ_LOG_WARN,  __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_ERR(fmt, ...)     mdoodz_log_emit(MDOODZ_LOG_ERROR, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_DBG(fmt, ...)     mdoodz_log_emit(MDOODZ_LOG_DEBUG, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#define LOG_TIME(fmt, ...)    mdoodz_log_emit(MDOODZ_LOG_TIMING, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
```

**Why macros, not functions?** Macros capture `__FILE__`/`__LINE__` at call site for free. The `##__VA_ARGS__` GNU extension (already used in the codebase via `_GNU_SOURCE`) handles zero-argument calls. The backing function `mdoodz_log_emit` does all the work.

**Alternative considered:** Direct function calls `mdoodz_log(level, fmt, ...)` — rejected because it loses call-site info without extra manual parameters.

### 2. Global logger state via static struct, not a passed-around context

```c
// Internal (MdoodzLog.c)
static struct {
    FILE *log_file;
    int   min_level;       // MDOODZ_LOG_ERROR=0 .. MDOODZ_LOG_TIMING=4
    int   dest;            // MDOODZ_LOG_CONSOLE | MDOODZ_LOG_FILE
    int   show_timestamp;  // 1 (default) or 0
    int   timestamp_mode;  // MDOODZ_TS_RELATIVE (default), MDOODZ_TS_ABSOLUTE, MDOODZ_TS_BOTH
    int   show_metadata;   // 0 (default) or 1 — show step/iteration in prefix
    double t0;             // omp_get_wtime() at init — for relative timestamps
    // Metadata context — updated by setter calls, costs ~12 bytes
    int   step;            // current timestep (-1 = not set)
    int   iteration;       // current nonlinear iteration (-1 = not set)
} g_logger;
```

**Why global?** Every printf call in MDLIB would need a logger pointer threaded through function signatures — that's 100+ function signature changes. A file-scoped static with init/shutdown functions is the lightest touch.

**Alternative considered:** Logger pointer in the `model` struct — rejected because it would require changing every function that calls printf to accept model, and many helper functions don't have access to it.

### 3. Thread safety via independent writes, no mutex

Each `mdoodz_log_emit` call does one `fprintf` per destination (console, file). `fprintf` to different `FILE*` handles is thread-safe by POSIX. Writing to the same `FILE*` from multiple threads may interleave lines but won't corrupt — this matches current `printf` behavior in OpenMP parallel regions.

**No mutex needed.** A mutex would serialize all logging and create contention in parallel regions. The current code already accepts interleaved output from `printf` in `#pragma omp parallel for` sections.

### 4. ANSI color handling: strip at write time for file output

When writing to the log file, `mdoodz_log_emit` scans the formatted string for `\033[...m` sequences and skips them. This is a simple state machine (~10 lines) applied only to file writes.

**Why not separate format strings?** The existing color macros (`GREEN "text" RESET`) concatenate at compile time — the format string already contains escape codes. Stripping is simpler than maintaining parallel color/no-color format strings.

### 5. Timestamp format: configurable mode with three options

The timestamp mode is set at init and controls the prefix format:

**Relative (default)** — seconds since `mdoodz_log_init()`, best for profiling:
```
[  12.345] INFO  | Selected dt = 1.23e+10
```

**Absolute** — wall-clock time-of-day, useful for correlating with external events:
```
[14:23:01.123] INFO  | Selected dt = 1.23e+10
```

**Both** — absolute first, then relative, for maximum information:
```
[14:23:01.123 +12.345] INFO  | Selected dt = 1.23e+10
```

When timestamps are disabled entirely, the prefix is just the level:
```
INFO  | Selected dt = 1.23e+10
```

Absolute timestamps use `localtime` + `gettimeofday`/`omp_get_wtime` fractional part. Cost is ~50ns extra for the `localtime_r` call — negligible.

### 6. Optional metadata context in log prefix

When `show_metadata` is enabled, the log prefix includes timestep and iteration context:

```
[14:23:01.123 +12.345|S0042|N03] INFO  | Selected dt = 1.23e+10
[14:23:01.124 +12.346|S0042|N03] WARN  | CFL violation detected
```

Fields appear only when set (iteration absent between loops):
```
[  12.345|S0042] INFO  | Starting thermal solve
```

Metadata is off by default — no prefix change for users who don't opt in. The context is updated via setter calls at natural boundaries in the main loop:

```c
void mdoodz_log_set_step(int step);          // top of timestep loop
void mdoodz_log_set_iteration(int nit);      // top of nonlinear iteration
void mdoodz_log_clear_iteration(void);       // after iteration loop ends
```

**Cost**: 2 extra ints in the global struct (~8 bytes, always in L1 cache). One extra `snprintf` field per enabled tag (~5-10 ns). Zero heap allocation.

**Risk**: Forgetting to update context → stale step/iteration in log. Mitigation: setter calls are placed at the top of well-defined loops in Main_DOODZ.c — only 2-3 call sites.

### 7. Log file naming and location

The log file is created in the current working directory (where HDF5 output goes) with the name `mdoodz.log`. Overwritten each run — no appending. Created at `mdoodz_log_init()`, flushed at end of each timestep, closed at `mdoodz_log_shutdown()`.

**Why not timestamped filenames?** Simulations are run one at a time in dedicated directories. A single `mdoodz.log` is simpler to find and doesn't accumulate stale files.

### 8. Migration strategy for existing printf calls

Mapping:
| Current pattern | Replacement |
|---|---|
| `printf("text\n", ...)` | `LOG_INFO("text", ...)` |
| `printf(GREEN "text" RESET, ...)` | `LOG_INFO(GREEN "text" RESET, ...)` — colors preserved, stripped in file |
| `printf(RED "WARNING..." RESET)` | `LOG_WARN("WARNING...")` — level conveys severity, color added by logger for WARN |
| `fprintf(stderr, "error...")` | `LOG_ERR("error...")` |
| `if (noisy) printf(...)` | `LOG_DBG(...)` — filtered by log level instead of ad-hoc flags |
| `printf("** Time for X = %lf sec\n", elapsed)` | `LOG_TIME("Time for X = %lf sec", elapsed)` |

The `\n` at end of format strings is dropped — the logger appends newlines automatically.

### 9. Init/shutdown API and configuration source

Log configuration is read from the `.txt` parameter file via `ReadInputFile()` in `InputOutput.c`, using the same `ReadInt2()` pattern as all other parameters. Default values are set in `InputOutput.c` so that logging works even if the `.txt` file omits the logging section entirely.

**`.txt` parameters** (all optional, defaults in parentheses):
```
/**** LOGGING ****/
log_dest       = 2       / 0=console, 1=file, 2=both (default: 2)
log_level      = 2       / 0=error, 1=warn, 2=info, 3=debug, 4=timing (default: 2)
log_timestamp  = 1       / 0=off, 1=on (default: 1)
log_ts_mode    = 0       / 0=relative, 1=absolute, 2=both (default: 0)
log_metadata   = 0       / 0=off, 1=show step/iteration in prefix (default: 0)
```

**Startup sequence**: The logger is initialized with defaults *before* `ReadInputFile()`, so that early startup messages (file reading, banner) are captured. After `ReadInputFile()` returns, the logger is reconfigured with the parsed values. The first ~10 lines use defaults — this is acceptable.

```c
// In RunMDOODZ():
mdoodz_log_init(NULL);                          // 1. defaults: both, info, timestamps on
Input input = ReadInputFile(inputFileName);      // 2. reads .txt including log_* params
mdoodz_log_reconfigure(&input.model.log);        // 3. apply user settings from .txt
```

**Config struct and API:**

```c
typedef enum {
    MDOODZ_TS_RELATIVE = 0,  // seconds since init (default)
    MDOODZ_TS_ABSOLUTE = 1,  // wall-clock HH:MM:SS.mmm
    MDOODZ_TS_BOTH     = 2,  // absolute first, then +relative
} MdoodzTimestampMode;

typedef struct {
    int dest;                    // MDOODZ_LOG_CONSOLE | MDOODZ_LOG_FILE (bitfield, default: both)
    int min_level;               // Minimum level to emit (default: MDOODZ_LOG_INFO)
    int show_timestamp;          // 1 = show timestamps (default), 0 = hide
    MdoodzTimestampMode ts_mode; // default: MDOODZ_TS_RELATIVE
    int show_metadata;           // 1 = show step/iteration in prefix, 0 = hide (default)
    const char *log_path;        // NULL = "mdoodz.log" in cwd
} MdoodzLogConfig;

void mdoodz_log_init(const MdoodzLogConfig *config);   // NULL = all defaults
void mdoodz_log_reconfigure(const MdoodzLogConfig *config); // apply new settings (e.g. after .txt parse)
void mdoodz_log_shutdown(void);
void mdoodz_log_flush(void);  // explicit flush at timestep boundaries

// Metadata context setters
void mdoodz_log_set_step(int step);
void mdoodz_log_set_iteration(int nit);
void mdoodz_log_clear_iteration(void);
```

The `model` struct in `mdoodz.h` gains a `MdoodzLogConfig log;` field, populated by `ReadInputFile()` with defaults from `ReadInt2()`'s default-value parameter.

## Risks / Trade-offs

**[Interleaved output in parallel regions]** → Accepted. Same behavior as current `printf`. Logging from inside `#pragma omp parallel for` may interleave lines. Mitigation: timing instrumentation is placed outside parallel regions; diagnostic output inside loops is DEBUG level.

**[Large diff touching every MDLIB file]** → Accepted. This is unavoidable for a complete migration. Mitigation: mechanical replacement (search/replace with verification), reviewed per-file. Existing unit/visual tests validate output correctness.

**[Trailing newline inconsistency]** → The logger auto-appends `\n`. Some existing printf calls have `\n` in the middle of format strings (multi-line output). These will need case-by-case handling — either split into multiple LOG_ calls or use a raw `mdoodz_log_emit_raw()` variant that doesn't append newlines.

**[Performance in tight loops]** → Some printf calls are inside particle loops (200k+ iterations). These should be removed or guarded with DEBUG level. The logger itself adds ~100ns overhead per call (timestamp + fprintf), but calls inside tight loops were already problematic with printf.
