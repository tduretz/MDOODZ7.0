## ADDED Requirements

### Requirement: Open MDOODZ HDF5 output files via h5wasm

The hdf5-reader module SHALL open `Output*.gzip.h5` files using h5wasm in read-only mode. It SHALL support gzip-compressed datasets. It SHALL read the `/Model/Params` array and extract grid dimensions: Nx = `Params[3]`, Nz = `Params[4]`.

#### Scenario: Open a valid output file

- **WHEN** the reader opens `Output00010.gzip.h5`
- **THEN** it SHALL return a file handle without error and Nx, Nz SHALL be positive integers

#### Scenario: Open a non-existent file

- **WHEN** the reader attempts to open a file that does not exist
- **THEN** it SHALL throw an error with a descriptive message including the file path

#### Scenario: Open a corrupt or non-HDF5 file

- **WHEN** the reader attempts to open a file that is not valid HDF5
- **THEN** it SHALL throw an error rather than returning partial data

### Requirement: Extract cell-centre fields

The reader SHALL extract any dataset from the `/Centers/` group and reshape the flat 1D array into a 2D array of dimensions `(Nx-1) × (Nz-1)`. The corresponding coordinate arrays SHALL be `xc_coord` (length Nx-1) and `zc_coord` (length Nz-1) read from `/Model/`.

#### Scenario: Read pressure field

- **WHEN** the reader extracts field `/Centers/P` from a file with Nx=601, Nz=401
- **THEN** the returned values array SHALL have shape 600 × 400
- **AND** xCoords SHALL have length 600 and zCoords SHALL have length 400

#### Scenario: Read all centre fields listed in the field registry

- **WHEN** the reader extracts any centre field (`P`, `T`, `eta_n`, `rho_n`, `sxxd`, `szzd`, `exxd`, `ezzd`, `d`, `X`, `phi`, `strain`, `strain_el`, `strain_pl`, `strain_pwl`, `cohesion`, `friction`, `divu`, `eII_el`, `eII_pl`, `eII_pwl`, `eII_exp`, `eII_lin`, `eII_gbs`)
- **THEN** each SHALL reshape to `(Nx-1) × (Nz-1)` without error

### Requirement: Extract vertex fields

The reader SHALL extract datasets from `/Vertices/` and reshape to `Nx × Nz`. Coordinate arrays SHALL be `xg_coord` (length Nx) and `zg_coord` (length Nz).

#### Scenario: Read vertex shear stress

- **WHEN** the reader extracts `/Vertices/sxz` from a file with Nx=601, Nz=401
- **THEN** the returned values array SHALL have shape 601 × 401

### Requirement: Extract velocity fields on staggered nodes

The reader SHALL extract `/VxNodes/Vx` and reshape to `Nx × (Nz+1)` with coordinate arrays `xg_coord` and `zvx_coord`. It SHALL extract `/VzNodes/Vz` and reshape to `(Nx+1) × Nz` with coordinate arrays `xvz_coord` and `zg_coord`.

#### Scenario: Read Vx field

- **WHEN** the reader extracts `/VxNodes/Vx` from a file with Nx=601, Nz=401
- **THEN** the returned values array SHALL have shape 601 × 402

#### Scenario: Read Vz field

- **WHEN** the reader extracts `/VzNodes/Vz` from a file with Nx=601, Nz=401
- **THEN** the returned values array SHALL have shape 602 × 401

### Requirement: Extract VizGrid phase maps

The reader SHALL extract `/VizGrid/compo` and reshape to `(Nx-1) × (Nz-1)`. It SHALL also support `compo_hr` and `compo_dual_hr` using the hi-res coordinate arrays `xviz_hr` and `zviz_hr` from `/VizGrid/`.

#### Scenario: Read normal-resolution phase map

- **WHEN** the reader extracts `/VizGrid/compo` from a file with Nx=601, Nz=401
- **THEN** the returned values array SHALL have shape 600 × 400 with integer phase values

#### Scenario: Read hi-res phase map

- **WHEN** the reader extracts `/VizGrid/compo_hr`
- **THEN** the returned dimensions SHALL match `(2*(Nx-1)) × (2*(Nz-1))`

### Requirement: Compute stress second invariant

The reader SHALL compute τ_II at cell centres from `sxxd`, `szzd` (centres) and `sxz` (vertices). The vertex shear stress SHALL be averaged to centres using a 2×2 block average: `sxz_c[i,j] = 0.25 * (sxz[i,j] + sxz[i+1,j] + sxz[i,j+1] + sxz[i+1,j+1])`. The invariant SHALL be: `τ_II = sqrt(0.5 * (sxxd² + szzd² + (sxxd + szzd)²) + sxz_c²)`.

#### Scenario: Compute stress invariant

- **WHEN** the reader derives "Stress II" for a valid output file
- **THEN** all τ_II values SHALL be non-negative
- **AND** the result shape SHALL be `(Nx-1) × (Nz-1)`

### Requirement: Compute strain rate second invariant

The reader SHALL compute ε̇_II at cell centres from `exxd`, `ezzd` (centres) and `exz` (vertices). The vertex-to-centre averaging SHALL use the same 2×2 block average as for stress. The invariant SHALL be: `ε̇_II = sqrt(0.5 * (exxd² + ezzd² + (exxd + ezzd)²) + exz_c²)`.

#### Scenario: Compute strain rate invariant

- **WHEN** the reader derives "Strain rate II" for a valid output file
- **THEN** all ε̇_II values SHALL be non-negative
- **AND** the result shape SHALL be `(Nx-1) × (Nz-1)`

### Requirement: Mask air cells

The reader SHALL identify air cells by reading `/VizGrid/compo` and finding cells where phase == -1. For any centre-grid field, air cells SHALL have their values replaced with `NaN` before returning data to the caller.

#### Scenario: Air cells masked in temperature field

- **WHEN** the reader extracts temperature and the phase map contains air cells (phase == -1)
- **THEN** the temperature values at air cell positions SHALL be `NaN`

#### Scenario: No air cells present

- **WHEN** the phase map contains no air cells
- **THEN** no values SHALL be replaced with `NaN`

### Requirement: Read model metadata

The reader SHALL extract model metadata from `/Model/Params` and return: time (seconds), domain width (metres), domain height (metres), Nx, Nz, dx, dz, dt. All values are already in SI units in HDF5.

#### Scenario: Extract metadata

- **WHEN** the reader reads metadata from a valid output file
- **THEN** it SHALL return time >= 0, Nx > 0, Nz > 0, dx > 0, dz > 0

### Requirement: File handle caching

The reader SHALL cache open h5wasm file handles using an LRU strategy with a configurable maximum (default: 5). When the cache is full and a new file is opened, the least recently used handle SHALL be closed.

#### Scenario: Cache hit

- **WHEN** the same file is requested twice consecutively
- **THEN** the second request SHALL reuse the cached handle without reopening

#### Scenario: Cache eviction

- **WHEN** 6 different files are opened with a cache max of 5
- **THEN** the first file's handle SHALL be closed
- **AND** re-requesting the first file SHALL reopen it

### Requirement: Return structured field response

The reader SHALL return field data in a structured format containing: `xCoords` (1D array), `zCoords` (1D array), `values` (2D array, row-major), `nx` (number of columns), `nz` (number of rows), `unit` (string), `log` (boolean — whether log₁₀ display is recommended), `min` (numeric, excluding NaN), `max` (numeric, excluding NaN).

#### Scenario: Complete response structure

- **WHEN** the reader extracts any field
- **THEN** the response SHALL contain all required properties
- **AND** `min` and `max` SHALL exclude NaN values
