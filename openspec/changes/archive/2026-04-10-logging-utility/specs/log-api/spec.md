## ADDED Requirements

### Requirement: Logger initialization with defaults
The system SHALL provide a `mdoodz_log_init()` function that initializes the logger with default settings when called with NULL. Default settings SHALL be: destination=both (file+console), level=INFO, timestamps=on, timestamp mode=relative, metadata=off, log path="mdoodz.log".

#### Scenario: Init with NULL config
- **WHEN** `mdoodz_log_init(NULL)` is called
- **THEN** the logger SHALL write to both console and file, at INFO level, with relative timestamps enabled, metadata disabled, and log file "mdoodz.log" in the current working directory

#### Scenario: Init with custom config
- **WHEN** `mdoodz_log_init()` is called with a `MdoodzLogConfig` struct specifying `dest=MDOODZ_LOG_FILE`, `min_level=MDOODZ_LOG_DEBUG`, `show_timestamp=0`
- **THEN** the logger SHALL write to file only, at DEBUG level, with no timestamp prefix

### Requirement: Logger reconfiguration after .txt parse
The system SHALL provide a `mdoodz_log_reconfigure()` function that updates logger settings at runtime without restarting. This SHALL be called after `ReadInputFile()` parses the `.txt` parameter file.

#### Scenario: Reconfigure from .txt parameters
- **WHEN** `ReadInputFile()` reads `log_level = 3` and `log_metadata = 1` from the `.txt` file, and `mdoodz_log_reconfigure()` is called with the parsed config
- **THEN** subsequent log calls SHALL use DEBUG level filtering and include metadata (step/iteration) in the prefix

#### Scenario: .txt omits logging parameters
- **WHEN** the `.txt` file contains no `log_*` parameters
- **THEN** `ReadInputFile()` SHALL use default values (dest=2, level=2, timestamp=1, ts_mode=0, metadata=0) and the logger SHALL continue with those defaults after reconfiguration

### Requirement: .txt parameter file logging section
The system SHALL read the following parameters from the `.txt` file via `ReadInt2()` in `ReadInputFile()`: `log_dest`, `log_level`, `log_timestamp`, `log_ts_mode`, `log_metadata`. Each parameter SHALL have a default value set in `InputOutput.c`.

#### Scenario: All logging parameters specified
- **WHEN** the `.txt` file contains `log_dest = 1`, `log_level = 0`, `log_timestamp = 0`, `log_ts_mode = 0`, `log_metadata = 0`
- **THEN** the logger SHALL write to file only, at ERROR level, with no timestamps and no metadata

### Requirement: Log level filtering
The system SHALL define five log levels: ERROR (0), WARN (1), INFO (2), DEBUG (3), TIMING (4). Messages below the configured minimum level SHALL be suppressed. The default minimum level SHALL be INFO.

#### Scenario: INFO level filters DEBUG
- **WHEN** the logger is configured with `min_level=MDOODZ_LOG_INFO` and `LOG_DBG("debug message")` is called
- **THEN** no output SHALL be produced

#### Scenario: DEBUG level passes INFO
- **WHEN** the logger is configured with `min_level=MDOODZ_LOG_DEBUG` and `LOG_INFO("info message")` is called
- **THEN** the message SHALL be emitted to all configured destinations

#### Scenario: ERROR always passes
- **WHEN** the logger is configured with any `min_level` and `LOG_ERR("error")` is called
- **THEN** the message SHALL be emitted (ERROR is the lowest level, always passes)

### Requirement: Log level macros with call-site info
The system SHALL provide variadic macros `LOG_INFO`, `LOG_WARN`, `LOG_ERR`, `LOG_DBG`, and `LOG_TIME` that capture `__FILE__` and `__LINE__` at the call site and delegate to `mdoodz_log_emit()`.

#### Scenario: Macro captures file and line
- **WHEN** `LOG_INFO("test")` is called from `Main_DOODZ.c` line 100
- **THEN** the internal call to `mdoodz_log_emit` SHALL receive `__FILE__="Main_DOODZ.c"` and `__LINE__=100`

### Requirement: Output destination configuration
The system SHALL support three output destinations: console only (0), file only (1), or both (2, default). When file output is enabled, a log file SHALL be created at the configured path (default: "mdoodz.log" in cwd). The file SHALL be overwritten on each run, not appended.

#### Scenario: Console only
- **WHEN** destination is configured as console only
- **THEN** all log output SHALL go to stdout and no log file SHALL be created

#### Scenario: File only
- **WHEN** destination is configured as file only
- **THEN** all log output SHALL go to the log file and nothing SHALL be printed to stdout

#### Scenario: Both destinations
- **WHEN** destination is configured as both
- **THEN** each log message SHALL be written to both stdout and the log file

#### Scenario: File overwritten per run
- **WHEN** a log file "mdoodz.log" already exists and a new simulation starts
- **THEN** the existing file SHALL be overwritten (not appended to)

### Requirement: Timestamp modes
The system SHALL support three timestamp modes: RELATIVE (0, default), ABSOLUTE (1), and BOTH (2). Timestamps SHALL be optional (configurable via `show_timestamp`).

#### Scenario: Relative timestamp
- **WHEN** timestamp mode is RELATIVE and timestamps are enabled
- **THEN** each log line prefix SHALL contain elapsed seconds since `mdoodz_log_init()` with millisecond precision, formatted as `[  12.345]`

#### Scenario: Absolute timestamp
- **WHEN** timestamp mode is ABSOLUTE and timestamps are enabled
- **THEN** each log line prefix SHALL contain wall-clock time formatted as `[HH:MM:SS.mmm]`

#### Scenario: Both timestamps
- **WHEN** timestamp mode is BOTH and timestamps are enabled
- **THEN** each log line prefix SHALL contain absolute time followed by relative time, formatted as `[HH:MM:SS.mmm +12.345]`

#### Scenario: Timestamps disabled
- **WHEN** `show_timestamp` is 0
- **THEN** log lines SHALL have no timestamp in the prefix, only the level label

### Requirement: Metadata context in prefix
The system SHALL support optional metadata fields (timestep, iteration) in the log prefix, controlled by `show_metadata`. Context SHALL be updated via `mdoodz_log_set_step()`, `mdoodz_log_set_iteration()`, and `mdoodz_log_clear_iteration()`.

#### Scenario: Metadata enabled with step and iteration
- **WHEN** `show_metadata=1`, step is set to 42, and iteration is set to 3
- **THEN** log prefix SHALL include `|S0042|N03` after the timestamp

#### Scenario: Metadata enabled, iteration not set
- **WHEN** `show_metadata=1`, step is set to 42, and iteration has been cleared
- **THEN** log prefix SHALL include `|S0042` without iteration field

#### Scenario: Metadata disabled
- **WHEN** `show_metadata=0`
- **THEN** no step or iteration fields SHALL appear in the log prefix regardless of setter calls

### Requirement: ANSI color stripping for file output
The system SHALL strip ANSI escape sequences (`\033[...m`) from log messages when writing to the log file. Console output SHALL preserve ANSI codes.

#### Scenario: Color in console, stripped in file
- **WHEN** `LOG_INFO(GREEN "Starting phase" RESET)` is called with both destinations enabled
- **THEN** console output SHALL contain the ANSI escape codes and file output SHALL contain only "Starting phase" without escape sequences

### Requirement: Automatic newline
The system SHALL append a newline (`\n`) to every log message. Callers SHALL NOT include trailing `\n` in format strings.

#### Scenario: Single line output
- **WHEN** `LOG_INFO("message")` is called
- **THEN** the output SHALL be the formatted prefix + "message" + newline

### Requirement: Thread safety in OpenMP parallel regions
The system SHALL be thread-safe for concurrent calls from OpenMP parallel regions without using mutexes. Line interleaving between threads is acceptable.

#### Scenario: Concurrent logging from parallel for
- **WHEN** multiple threads call `LOG_DBG(...)` inside a `#pragma omp parallel for` region
- **THEN** all messages SHALL be emitted without data corruption, though lines MAY interleave

### Requirement: Logger shutdown
The system SHALL provide `mdoodz_log_shutdown()` that flushes and closes the log file. After shutdown, log calls SHALL write to console only.

#### Scenario: Clean shutdown
- **WHEN** `mdoodz_log_shutdown()` is called
- **THEN** the log file SHALL be flushed and closed

### Requirement: Explicit flush
The system SHALL provide `mdoodz_log_flush()` that flushes the log file buffer to disk. This SHALL be called at timestep boundaries.

#### Scenario: Flush at timestep end
- **WHEN** `mdoodz_log_flush()` is called at the end of a timestep
- **THEN** all buffered log data SHALL be written to disk

### Requirement: All existing printf calls migrated
All existing `printf()` and `fprintf(stderr, ...)` calls across all MDLIB `.c` source files SHALL be replaced with the appropriate `LOG_*` macro. No raw `printf()` calls SHALL remain in MDLIB after migration.

#### Scenario: No remaining raw printf
- **WHEN** the migration is complete
- **THEN** `grep -r "printf(" MDLIB/*.c` SHALL return zero matches (excluding the logger implementation file itself)
