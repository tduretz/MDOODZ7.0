## ADDED Requirements

### Requirement: Cache extracted field JSON to disk

The field-cache module SHALL write the JSON response of each field extraction to a `.cache/` directory under the data directory. The cache key SHALL be `sha256(filename + ":" + fieldName)` truncated to 16 hex characters, with a `.json` extension.

#### Scenario: First extraction creates cache file

- **WHEN** field "Pressure" is extracted from `Output00010.gzip.h5` for the first time
- **THEN** the cache directory SHALL contain a file named `<sha256("Output00010.gzip.h5:Pressure")[0:16]>.json`
- **AND** the file SHALL contain the complete JSON response (values, coords, min, max, pMin, pMax, unit, log, nx, nz)

#### Scenario: Cache directory does not exist

- **WHEN** the `.cache/` directory does not exist on the first extraction
- **THEN** the module SHALL create it automatically

### Requirement: Serve from cache on subsequent requests

When a field-data request is received, the module SHALL check if a valid cache file exists before performing HDF5 extraction. If valid, the cached JSON SHALL be streamed directly as the response.

#### Scenario: Cache hit

- **WHEN** field "Pressure" for `Output00010.gzip.h5` is requested and a valid cache file exists
- **THEN** the response SHALL be served from the cache file without opening the HDF5 file
- **AND** the response content SHALL be identical to a fresh extraction

#### Scenario: Cache miss

- **WHEN** no cache file exists for the requested file+field combination
- **THEN** the module SHALL extract from HDF5, write the cache file, and then respond

### Requirement: Invalidate cache when source file changes

The module SHALL compare the mtime of the cache file against the mtime of the source HDF5 file on every request. If the source file's mtime is newer than the cache file's mtime, the cache entry SHALL be considered stale and re-extracted.

#### Scenario: Source file updated after caching

- **WHEN** `Output00010.gzip.h5` is overwritten (new mtime) after its "Pressure" field was cached
- **THEN** the next request for "Pressure" SHALL re-extract from HDF5 and overwrite the cache file

#### Scenario: Source file unchanged

- **WHEN** `Output00010.gzip.h5` has not been modified since caching
- **THEN** the cached JSON SHALL be served directly

### Requirement: Cache file format includes percentile range

Cached JSON files SHALL include both the full min/max and the percentile-clipped pMin/pMax values so that the frontend receives all range information from the cache.

#### Scenario: Cached response contains percentile bounds

- **WHEN** a cache file is read and served
- **THEN** the JSON SHALL contain `min`, `max`, `pMin`, and `pMax` numeric properties

### Requirement: No cache size limit in initial implementation

The cache SHALL NOT enforce a maximum size or entry count in this initial implementation. All extracted fields for all files SHALL be cached indefinitely until invalidated by source file change or manual deletion.

#### Scenario: Many fields cached

- **WHEN** 47 fields × 29 files are all extracted and cached
- **THEN** all 1363 cache files SHALL coexist without eviction
