---
name: skill-visualisation
description: Visualising MDOODZ results — HDF5 output structure, Julia/Makie plotting scripts, field names, overlay options (contours, director, velocity vectors, principal stress), multi-panel layouts, marker trajectory analysis, batch PNG/video export.
---

# Visualisation Toolkit

## HDF5 Output Structure

Each time step produces a file `Output<NNNNN>.gzip.h5` containing groups:

### Groups and Key Datasets

| Group | Datasets | Grid location |
|-------|----------|---------------|
| **Model** | `Params` (array of model parameters), `xc_coord`, `zc_coord`, `xg_coord`, `zg_coord`, `xvz_coord`, `zvx_coord` | Metadata |
| **Centers** | `P`, `sxxd`, `szzd`, `T`, `eta_n`, `rho_n`, `d`, `phi`, `X`, `strain`, `strain_el`, `strain_pl`, `strain_pwl`, `strain_exp`, `strain_lin`, `strain_gbs`, `exxd`, `ezzd`, `eII_el`, `eII_pl`, `eII_pwl`, `eII_exp`, `eII_lin`, `eII_gbs`, `divu`, `divu_el`, `divu_pl`, `divu_th`, `divu_r`, `cohesion`, `friction` | Cell centres |
| **Vertices** | `sxz`, `exz`, `eta_s`, `rho_s` | Cell vertices |
| **VxNodes** | `Vx` | Vertical cell edges |
| **VzNodes** | `Vz` | Horizontal cell edges |
| **VizGrid** | `compo`, `compo_hr`, `compo_dual` | Phase visualisation (normal + 2× resolution) |
| **Topo** | `z_grid`, `z_mark`, `x_mark` | Topography and surface markers |
| **Flags** | `tag_n`, `tag_s`, `tag_u`, `tag_v` | BC flags |
| **Iterations** | Convergence residuals per nonlinear iteration | Solver diagnostics |
| **TimeSeries** | 14 time-averaged quantities (Tmax, Vmax, τ_II max, etc.) | Evolution tracking |

### Optional Datasets (enabled by model switches)

| Dataset | Condition | Group |
|---------|-----------|-------|
| `nx`, `nz`, `ani_fac` | `anisotropy = 1` | Centers |
| `Fxx`, `Fxz`, `Fzx`, `Fzz` | Finite strain tracking | Centers |
| `strain_pl_vol` | Volumetric plastic strain | Centers |

### Particle Files

Separate files `Particles<NNNNN>.gzip.h5` contain raw marker data under `/Particles/`:
`x`, `z`, `phase`, `generation`, `T`, `P`, `d`, `phi`, `X`, `Fxx`, `Fxz`, `Fzx`, `Fzz`

## Params Array (Model Metadata)

The `/Model/Params` dataset is a 1D array of 8 doubles:

| Index | Content | Units |
|-------|---------|-------|
| 0 | Current model time | s (SI, already scaled) |
| 1 | Domain width (xmax − xmin) | m |
| 2 | Domain height (zmax − zmin) | m |
| 3 | Nx (number of **vertex** nodes in x) | – |
| 4 | Nz (number of **vertex** nodes in z) | – |
| 5 | dx (cell width) | m |
| 6 | dz (cell height) | m |
| 7 | dt (current time step) | s |

Derive grid sizes from Params:
```
nvx, nvz = Int(Params[4]), Int(Params[5])   # vertex counts
ncx, ncz = nvx - 1, nvz - 1                 # cell-centre counts
```

## Array Dimensions and Staggered Grid Indexing

All HDF5 datasets are stored as **flat 1D arrays**. You must reshape them to 2D using the correct dimensions for each grid location.

### Dimension Table

| Grid location | Reshape size | Coordinate arrays | Which datasets |
|---------------|-------------|-------------------|----------------|
| Cell centres | `(ncx, ncz)` = `(Nx-1, Nz-1)` | `xc_coord` (Nx-1), `zc_coord` (Nz-1) | P, T, sxxd, szzd, eta_n, rho_n, strain*, exxd, ezzd, eII_*, d, X, phi, cohesion, friction, divu*, nx, nz, ani_fac, Fxx/Fxz/Fzx/Fzz |
| Cell vertices | `(nvx, nvz)` = `(Nx, Nz)` | `xg_coord` (Nx), `zg_coord` (Nz) | sxz, exz, eta_s, rho_s |
| Vx nodes | `(nvx, nvz+1)` = `(Nx, Nz+1)` | `xg_coord` (Nx), `zvx_coord` (Nz+1) | Vx |
| Vz nodes | `(nvx+1, nvz)` = `(Nx+1, Nz)` | `xvz_coord` (Nx+1), `zg_coord` (Nz) | Vz |
| VizGrid normal | `(ncx, ncz)` | same as centres | compo |
| VizGrid hi-res | `(ncx_hr, ncz_hr)` | `xviz_hr`, `zviz_hr` | compo_hr, compo_dual_hr |
| Topography | `(Nx,)` 1D | `xg_coord` | z_grid, Vx_grid |
| BC flags | same as their grid | — | tag_n, tag_s, tag_u, tag_v |

> **Warning — coordinate mapping**: For Vx nodes, the z-coordinates include ghost rows at iz=0 and iz=Nz with half-cell positions *outside* the domain. Interior nodes (iz=1..Nz-1) sit at cell-center z positions. If you reconstruct coordinates from `xmin`/`zmin` instead of reading the HDF5 coordinate arrays, the correct formulas are: `z = zmin + (iz - 0.5)*dz` for Vx, `x = xmin + (ix - 0.5)*dx` for Vz. Using `(iz + 0.5)*dz` is a common error — it shifts by one full cell width. Always prefer reading the dedicated coordinate arrays (`zvx_coord`, `xvz_coord`) from HDF5 when available.

### Extraction Pattern (Julia)

```julia
using HDF5

filename = "Output00010.gzip.h5"

# 1. Read metadata
Params = h5read(filename, "/Model/Params")
nvx, nvz = Int(Params[4]), Int(Params[5])
ncx, ncz = nvx - 1, nvz - 1

# 2. Read coordinates
xc = h5read(filename, "/Model/xc_coord")   # length ncx
zc = h5read(filename, "/Model/zc_coord")   # length ncz
xv = h5read(filename, "/Model/xg_coord")   # length nvx
zv = h5read(filename, "/Model/zg_coord")   # length nvz

# 3. Read a centre field → reshape from 1D to 2D
T_flat = h5read(filename, "/Centers/T")         # length ncx*ncz
T      = reshape(Float64.(T_flat), ncx, ncz)    # 2D array

# 4. Read a vertex field
sxz_flat = h5read(filename, "/Vertices/sxz")
sxz      = reshape(Float64.(sxz_flat), nvx, nvz)

# 5. Read velocity (staggered — note +1 on one dimension)
Vx = reshape(Float64.(h5read(filename, "/VxNodes/Vx")), nvx, nvz + 1)
Vz = reshape(Float64.(h5read(filename, "/VzNodes/Vz")), nvx + 1, nvz)
```

### Interpolating Between Grids

Vertex → centre (for 2nd invariant computation):
```julia
# Average sxz from 4 surrounding vertices to cell centre
sxz_c = 0.25 .* (sxz[1:end-1, 1:end-1] .+ sxz[2:end, 1:end-1] .+
                  sxz[1:end-1, 2:end]   .+ sxz[2:end, 2:end])
```

Velocity → centre:
```julia
Vxc = 0.5 .* (Vx[1:end-1, 2:end-1] .+ Vx[2:end, 2:end-1])
Vzc = 0.5 .* (Vz[2:end-1, 1:end-1] .+ Vz[2:end-1, 2:end])
```

### Computing Invariants at Centres

```julia
τxx = reshape(Float64.(h5read(filename, "/Centers/sxxd")), ncx, ncz)
τzz = reshape(Float64.(h5read(filename, "/Centers/szzd")), ncx, ncz)
τyy = -(τxx .+ τzz)   # plane strain: τ_yy = -(τ_xx + τ_zz)

sxz_c = 0.25 .* (sxz[1:end-1,1:end-1] .+ sxz[2:end,1:end-1] .+
                  sxz[1:end-1,2:end]   .+ sxz[2:end,2:end])

τII = sqrt.(0.5 .* (τxx.^2 .+ τyy.^2 .+ τzz.^2) .+ sxz_c.^2)
```

### Masking Air Cells

Phase = -1 marks air (free surface). Mask before plotting:
```julia
ph   = reshape(Float64.(h5read(filename, "/VizGrid/compo")), ncx, ncz)
mask = ph .== -1.0
T[mask] .= NaN    # NaN values are invisible in plots
```

## Scaling and Units

All values in HDF5 are **already in SI units** (the C code multiplies by the scaling factors before writing). No further de-scaling is needed when reading.

- Coordinates: metres (m)
- Stress: Pascals (Pa)
- Strain rate: s⁻¹
- Viscosity: Pa·s
- Temperature: Kelvin (K) — subtract 273.15 for °C
- Velocity: m/s
- Time (`Params[0]`): seconds — divide by `3.1558e7` for Myr
- Grain size: metres (m)

## Julia/Makie Project Setup

### Dependencies

The visualisation toolkit lives in `JuliaVisualisation/` and uses:
- **HDF5** — read output files
- **CairoMakie** — static PNG/PDF figures
- **GLMakie** — interactive OpenGL windows (optional)
- **FFMPEG** — video encoding
- **Colors**, **ColorSchemes** — colour maps
- **MathTeXEngine** — LaTeX labels
- **CSV**, **DataFramesMeta** — marker history tables

### Installation

```julia
cd("JuliaVisualisation")
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

### Module Structure

`JuliaVisualisation/src/` provides:

| File | Exports |
|------|---------|
| `LoadDataFiles.jl` | `ExtractData(file, path)`, `ExtractField(file, field, size, mask_air, mask)`, `ReadFile(path, step, scales, options, PlotOnTop)`, `Print2Disk(f, path, field, istep, Mak; res=4)` |
| `ContoursQuivers.jl` | `AddCountourQuivers!(...)`, `PrincipalStress(τxx, τzz, τxz, P)` |

## Available Plot Fields

The main script `Main_Visualisation_Makie_MD7.jl` supports 22+ fields:

| Field | Source dataset | Physical meaning |
|-------|---------------|-----------------|
| Phases | `VizGrid/compo` | Material phase map (custom palette) |
| Density | `Centers/rho_n` or `Vertices/rho_s` | ρ [kg/m³] |
| Viscosity | `Centers/eta_n` or `Vertices/eta_s` | η [Pa·s], log₁₀ scale |
| Stress II | computed from `sxxd`, `szzd`, `sxz` | τ_II [Pa] |
| Normal stress xx | `Centers/sxxd` | τ'_xx [Pa] |
| Normal stress zz | `Centers/szzd` | τ'_zz [Pa] |
| Strain rate II | computed from `exxd`, `ezzd`, `exz` | ε̇_II [s⁻¹], log₁₀ |
| Plastic strain rate | `Centers/eII_pl` | ε̇_pl [s⁻¹] |
| Pressure | `Centers/P` | P [Pa] |
| Divergence | `Centers/divu` | ∇·v [s⁻¹] |
| Temperature | `Centers/T` | T [K or °C] |
| Grain size | `Centers/d` | d [m], log₁₀ |
| Melt fraction | `Centers/X` | X [0–1] |
| Anisotropy factor | `Centers/ani_fac` | δ |
| Velocity Vx | `VxNodes/Vx` | m/s |
| Velocity Vz | `VzNodes/Vz` | m/s |
| Velocity magnitude | computed from Vx, Vz | |v| |
| Topography | `Topo/z_grid` | Surface elevation [m] |
| Cohesion | `Centers/cohesion` | C [Pa] |
| Friction | `Centers/friction` | φ [rad] |
| Accumulated strain | `Centers/strain` | ε_total |
| Time series | `TimeSeries/*` | Evolution plots |

## Overlay Options

Configure in `PlotOnTop` struct:

| Overlay | What it shows |
|---------|---------------|
| Phase contours | Material boundaries |
| Temperature contours | Isotherms |
| Fabric / director | Anisotropy orientation ticks (nx, nz) |
| Principal stress | Eigenvectors of stress tensor (σ₁, σ₃) |
| Principal strain rate | Eigenvectors of strain rate tensor |
| Velocity vectors | Quiver plot of (Vx, Vz) |
| Topography | Surface line overlay on field plots |
| Melt fraction contours | Melt boundary lines |

## Multi-Panel Layout

`Main_Visualisation_Makie_MD7_3x3.jl` creates a 3×3 grid of subplots showing the same field at 9 different time steps for rapid temporal comparison. Default field: ε̇_II.

## Marker Trajectory Analysis

### Selecting Markers

`Select_Markers_MD7.jl`:
1. Define spatial bounds: `xmin`, `xmax`, `zmin`, `zmax`
2. Set starting file index: `file_select_start`
3. Script reads `Particles<NNNNN>.gzip.h5` and filters markers by position
4. Tracks matched marker IDs across subsequent time steps
5. Produces arrays: Time × Marker ID for position, temperature, pressure, phase

### Reading Marker History

`Read_Marker_History_MD7.jl`:
- Reads `MarkerData.csv`
- Fields: T [K], P [Pa], x [m], z [m], phase
- Useful for P–T path reconstruction of individual particles

## Batch Export and Video

### PNG Export

```julia
Print2Disk(f, path, field_name, istep, Mak; res=4)
```
- `res` controls DPI multiplier (default 4 → high resolution)
- Saves to `path/field_name_NNNNN.png`

### Video Generation

Use FFMPEG.jl to stitch PNGs into MP4:
```julia
using FFMPEG
FFMPEG.ffmpeg_exe(`-framerate 15 -i output_%05d.png -c:v libx264 -pix_fmt yuv420p movie.mp4`)
```

## Workflow

1. Run MDOODZ — produces `Output*.gzip.h5` files in run directory
2. Open `Main_Visualisation_Makie_MD7.jl`
3. Set `path` to output directory
4. Choose field and overlays
5. Set time step range (`file_start`, `file_step`, `file_end`)
6. Run script — produces figures or interactive window
7. Optionally export PNGs and stitch into video
