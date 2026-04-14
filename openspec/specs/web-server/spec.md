## ADDED Requirements

### Requirement: Serve static frontend files

The server SHALL serve all files under the `public/` directory for any request path not matching `/api/*`. It SHALL set appropriate `Content-Type` headers based on file extension (`.html`, `.css`, `.mjs`, `.json`, `.js`). Requests to `/` SHALL serve `public/index.html`.

#### Scenario: Serve index page

- **WHEN** a browser requests `GET /`
- **THEN** the server SHALL respond with `public/index.html` and `Content-Type: text/html`

#### Scenario: Serve ES6 module

- **WHEN** a browser requests `GET /app.mjs`
- **THEN** the server SHALL respond with the file contents and `Content-Type: application/javascript`

#### Scenario: File not found

- **WHEN** a browser requests a path that does not match any static file or API route
- **THEN** the server SHALL respond with HTTP 404

### Requirement: List available HDF5 files

The server SHALL expose `GET /api/files` which returns a JSON array of `Output*.gzip.h5` filenames found in the configured data directory, sorted by filename (which corresponds to time-step order).

#### Scenario: Directory with output files

- **WHEN** the data directory contains `Output00000.gzip.h5`, `Output00010.gzip.h5`, `Output00020.gzip.h5`
- **THEN** `GET /api/files` SHALL return `{ "files": ["Output00000.gzip.h5", "Output00010.gzip.h5", "Output00020.gzip.h5"] }`

#### Scenario: Empty directory

- **WHEN** the data directory contains no matching files
- **THEN** `GET /api/files` SHALL return `{ "files": [] }`

### Requirement: List available fields for a file

The server SHALL expose `GET /api/fields/:file` which returns the list of field names available for that file (from the field registry) along with model metadata (time, Nx, Nz, dx, dz, dt). The server SHALL only include fields whose HDF5 datasets actually exist in the file.

#### Scenario: List fields for a valid file

- **WHEN** `GET /api/fields/Output00010.gzip.h5` is requested
- **THEN** the response SHALL contain `fields` (array of field name strings) and `params` (object with `time`, `nx`, `nz`, `dx`, `dz`, `dt`)

#### Scenario: File does not exist

- **WHEN** `GET /api/fields/nonexistent.h5` is requested
- **THEN** the server SHALL respond with HTTP 404 and a JSON error message

### Requirement: Return field data

The server SHALL expose `GET /api/field-data/:file/:field` which uses the hdf5-reader to extract the named field and return JSON with `xCoords`, `zCoords`, `values`, `nx`, `nz`, `unit`, `log`, `min`, `max`.

#### Scenario: Extract a simple centre field

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/Pressure` is requested
- **THEN** the server SHALL return HTTP 200 with a JSON body containing all required field properties

#### Scenario: Extract a derived field

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/Stress%20II` is requested
- **THEN** the server SHALL compute the stress invariant and return the result

#### Scenario: Unknown field name

- **WHEN** `GET /api/field-data/Output00010.gzip.h5/UnknownField` is requested
- **THEN** the server SHALL respond with HTTP 400 and a JSON error message

### Requirement: Configure data directory via CLI argument

The server SHALL accept the data directory path as the first command-line argument: `node server.mjs /path/to/hdf5/files`. If no argument is provided, it SHALL default to the current working directory.

#### Scenario: Custom data directory

- **WHEN** the server is started with `node server.mjs ./benchmark-results/rifting-600x400-v2-h5`
- **THEN** `GET /api/files` SHALL list files from that directory

#### Scenario: Default to cwd

- **WHEN** the server is started with `node server.mjs` (no argument)
- **THEN** the data directory SHALL be the current working directory

### Requirement: Configure port

The server SHALL listen on port 3000 by default. The port SHALL be configurable via the `PORT` environment variable.

#### Scenario: Default port

- **WHEN** the server starts without `PORT` set
- **THEN** it SHALL listen on port 3000

#### Scenario: Custom port

- **WHEN** the server starts with `PORT=8080`
- **THEN** it SHALL listen on port 8080

### Requirement: Enable gzip compression for API responses

The server SHALL gzip-compress JSON API responses when the request includes `Accept-Encoding: gzip`.

#### Scenario: Client accepts gzip

- **WHEN** a request to `/api/field-data/...` includes `Accept-Encoding: gzip`
- **THEN** the response SHALL have `Content-Encoding: gzip` and the body SHALL be gzip-compressed

#### Scenario: Client does not accept gzip

- **WHEN** a request omits `Accept-Encoding` or does not include `gzip`
- **THEN** the response SHALL be uncompressed JSON

### Requirement: CORS headers for local development

The server SHALL include `Access-Control-Allow-Origin: *` on all API responses to allow cross-origin requests during development.

#### Scenario: API response includes CORS header

- **WHEN** any `/api/*` request is made
- **THEN** the response SHALL include the `Access-Control-Allow-Origin: *` header

### Requirement: Path traversal protection

The server SHALL reject any request path containing `..` or attempting to access files outside the `public/` directory (for static files) or outside the data directory (for API file parameters). The filename parameter in API routes SHALL be validated to match the pattern `Output\d+\.gzip\.h5`.

#### Scenario: Path traversal attempt in static files

- **WHEN** a request for `GET /../../../etc/passwd` is made
- **THEN** the server SHALL respond with HTTP 400

#### Scenario: Path traversal attempt in API file parameter

- **WHEN** a request for `GET /api/fields/../../secret.txt` is made
- **THEN** the server SHALL respond with HTTP 400

#### Scenario: Valid filename accepted

- **WHEN** a request for `GET /api/fields/Output00010.gzip.h5` is made
- **THEN** the server SHALL process the request normally
