## MODIFIED Requirements

### Requirement: List available HDF5 files

The server SHALL expose `GET /api/files` which returns a JSON object with a `files` array. Each entry SHALL be an object `{ name, step, time }` where `name` is the filename, `step` is the integer time-step number parsed from the filename, and `time` is the model time in seconds read from `/Model/Params[0]`. The files SHALL be sorted by step number. On the first call, the server SHALL read metadata from all files using `readAllMetadata` and cache the result in memory.

#### Scenario: Directory with output files

- **WHEN** the data directory contains `Output00000.gzip.h5` (t=0), `Output00010.gzip.h5` (t=3.1558e13), `Output00020.gzip.h5` (t=6.3116e13)
- **THEN** `GET /api/files` SHALL return `{ "files": [{ "name": "Output00000.gzip.h5", "step": 0, "time": 0 }, { "name": "Output00010.gzip.h5", "step": 10, "time": 3.1558e13 }, { "name": "Output00020.gzip.h5", "step": 20, "time": 6.3116e13 }] }`

#### Scenario: Empty directory

- **WHEN** the data directory contains no matching files
- **THEN** `GET /api/files` SHALL return `{ "files": [] }`

### Requirement: Params endpoint for scalar metadata

The server SHALL expose `GET /api/params/:file` which returns a JSON object `{ step, time, Nx, Nz, dt }` with scalar metadata from the specified HDF5 file. The endpoint SHALL read from the `/Model/Params` dataset without extracting any field data. The response SHALL be served from the in-memory metadata cache when available, falling back to reading the HDF5 file if necessary.

#### Scenario: Get params for a valid file

- **WHEN** `GET /api/params/Output00010.gzip.h5` is requested
- **THEN** the server SHALL return HTTP 200 with `{ "step": 10, "time": 3.1558e13, "Nx": 601, "Nz": 401, "dt": ... }`

#### Scenario: File not found

- **WHEN** `GET /api/params/NonExistent.h5` is requested
- **THEN** the server SHALL return HTTP 404 with a JSON error message

#### Scenario: Params served from metadata cache

- **WHEN** `GET /api/params/Output00010.gzip.h5` is requested after `/api/files` has already been called
- **THEN** the server SHALL return the params from the in-memory metadata cache without reopening the HDF5 file

### Requirement: Return field data

The server SHALL expose `GET /api/field-data/:file/:field` which uses the hdf5-reader to extract the named field and return JSON with `xCoords`, `zCoords`, `values`, `nx`, `nz`, `unit`, `log`, `min`, `max`, `pMin`, `pMax`. The endpoint SHALL check the disk cache (via field-cache module) before performing extraction. If the cache contains a valid entry, the cached JSON SHALL be served directly. If the cache misses or is stale, the server SHALL extract the field, write the result to cache, and then respond.

#### Scenario: Extract a simple centre field (cache miss)

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/Pressure` is requested and no cache exists
- **THEN** the server SHALL extract from HDF5, write a cache file, and return HTTP 200 with a JSON body including `pMin` and `pMax`

#### Scenario: Serve from cache (cache hit)

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/Pressure` is requested and a valid cache file exists
- **THEN** the server SHALL stream the cached JSON directly without opening the HDF5 file

#### Scenario: Cache invalidated by source update

- **WHEN** `Output00010.gzip.h5` has been modified since the cache was written
- **THEN** the server SHALL re-extract from HDF5, overwrite the cache file, and respond with fresh data

#### Scenario: Extract a derived field

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/Stress%20II` is requested
- **THEN** the server SHALL compute the stress invariant and return the result with percentile range

#### Scenario: Unknown field name

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/UnknownField` is requested
- **THEN** the server SHALL respond with HTTP 400 and a JSON error message
