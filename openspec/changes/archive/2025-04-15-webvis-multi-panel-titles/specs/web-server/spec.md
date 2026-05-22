## ADDED Requirements

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
