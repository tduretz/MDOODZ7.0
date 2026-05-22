---
name: skill-model-setup
description: Create and configure MDOODZ simulations — .c callback files, .txt parameter files, SetPhase, SetTemperature, SetGrainSize, SetSurfaceZCoord, boundary conditions, MdoodzSetup, Coordinates, POSITION enum, geometry helpers, and scenario authoring workflow.
---

# Model Setup and Configuration

## The Paired-File Pattern

Every MDOODZ simulation is defined by two files in `SETS/`:

1. **C callback file** (`<Name>.c`): Contains callback functions that define initial conditions and boundary conditions. Has a `main()` that populates a `MdoodzSetup` struct and calls `RunMDOODZ`.
2. **Parameter file** (`<Name>.txt`): Plain-text file with all runtime parameters — scaling, resolution, switches, material properties.

### Creating a New Scenario

```bash
# 1. Copy an existing pair as a starting point
cp SETS/RiftingChenin.c SETS/MyModel.c
cp SETS/RiftingChenin.txt SETS/MyModel.txt

# 2. Edit callbacks in MyModel.c (phase geometry, BCs, initial T, etc.)
# 3. Tune parameters in MyModel.txt (resolution, materials, switches)

# 4. Build and run
make build SET=MyModel
make run SET=MyModel
```

## Entry Point Pattern

Every `.c` file follows this canonical structure:

```c
#include "mdoodz.h"

// ... callback function implementations ...

int main(int nargs, char *args[]) {
    char *input = GetSetupFileName(nargs, args);
    MdoodzSetup setup = {
        .BuildInitialTopography = &(BuildInitialTopography_ff){
            .SetSurfaceZCoord = SetSurfaceZCoord,
        },
        .SetParticles = &(SetParticles_ff){
            .SetPhase       = SetPhase,
            .SetTemperature = SetTemperature,
            .SetGrainSize   = SetGrainSize,
            // ... other particle callbacks ...
        },
        .SetBCs = &(SetBCs_ff){
            .SetBCVx    = SetPureShearBCVx,   // or custom
            .SetBCVz    = SetPureShearBCVz,
            .SetBCPType = SetBCPType,
            .SetBCT     = SetBCT,
        },
        .MutateInput = NULL,  // optional
    };
    RunMDOODZ(input, &setup);
    free(input);
}
```

## The MdoodzSetup Struct

Groups all callback function pointers:

| Field | Type | Purpose |
|-------|------|---------|
| `BuildInitialTopography` | `BuildInitialTopography_ff*` | Surface topography callbacks |
| `SetParticles` | `SetParticles_ff*` | Particle initialisation callbacks |
| `SetBCs` | `SetBCs_ff*` | Boundary condition callbacks |
| `MutateInput` | `MutateInput_f*` | Optional: modify params programmatically after `.txt` parsing |

## Callback Function Signatures

### Particle Initialisation Callbacks (SetParticles_ff)

| Callback | Signature | Returns | Purpose |
|----------|-----------|---------|---------|
| `SetPhase` | `(MdoodzInput*, Coordinates)` | `int` | Phase ID at this position |
| `SetTemperature` | `(MdoodzInput*, Coordinates)` | `double` | Temperature (scaled) |
| `SetGrainSize` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Grain size (scaled) |
| `SetPorosity` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Porosity |
| `SetDensity` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Initial density |
| `SetXComponent` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Chemical composition |
| `SetPressure` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Initial pressure |
| `SetNoise` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Random perturbation |
| `SetDualPhase` | `(MdoodzInput*, Coordinates, int phase)` | `int` | Secondary phase for mixing |
| `SetAnisoAngle` | `(MdoodzInput*, Coordinates, int phase, double predefined)` | `double` | Anisotropy angle |
| `SetDefGrad` | `(MdoodzInput*, Coordinates, int phase)` | `Tensor2D` | Initial deformation gradient |
| `SetTxx` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Initial deviatoric stress xx |
| `SetTzz` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Initial deviatoric stress zz |
| `SetTxz` | `(MdoodzInput*, Coordinates, int phase)` | `double` | Initial deviatoric stress xz |
| `AdjustPhaseToTemperature` | `(MdoodzInput*, Coordinates, double T, int phase)` | `int` | Phase adjustment after T init |

### Topography Callbacks (BuildInitialTopography_ff)

| Callback | Signature | Returns | Purpose |
|----------|-----------|---------|---------|
| `SetSurfaceZCoord` | `(MdoodzInput*, double x_coord)` | `double` | Surface elevation at x |
| `SetSurfacePhase` | `(MdoodzInput*, double x_coord)` | `int` | Phase of surface material |

### Boundary Condition Callbacks (SetBCs_ff)

| Callback | Signature | Returns | Purpose |
|----------|-----------|---------|---------|
| `SetBCVx` | `(MdoodzInput*, POSITION, Coordinates)` | `SetBC` | Horizontal velocity BC |
| `SetBCVz` | `(MdoodzInput*, POSITION, Coordinates)` | `SetBC` | Vertical velocity BC |
| `SetBCT` | `(MdoodzInput*, POSITION, Coordinates, double gridT)` | `SetBC` | Temperature BC |
| `SetBCC` | `(MdoodzInput*, POSITION, Coordinates, double gridX)` | `SetBC` | Chemical BC |
| `SetBCPType` | `(MdoodzInput*, POSITION)` | `char` | Pressure BC type |
| `FixTemperature` | `(MdoodzInput*, double pressure)` | `double` | Fixed T at given P |

### The SetBC Return Struct

```c
typedef struct {
    double value;  // BC value (velocity, temperature, etc.)
    char   type;   // BC type code (see table below)
} SetBC;
```

### MutateInput Callback

```c
void MutateInput(MdoodzInput *input);
```

Called after `.txt` file parsing, before simulation starts. Use to programmatically override parameters.

## The Coordinates Struct

```c
typedef struct {
    double l;  // distance parameter
    double k;  // index parameter
    double x;  // horizontal position (scaled, non-dimensional)
    double z;  // vertical position (scaled, non-dimensional)
} Coordinates;
```

**Important**: Coordinates inside callbacks are in **scaled (non-dimensional) units**. To convert to physical units: `x_physical = coordinates.x * input->scaling.L`.

## The POSITION Enum

Used in boundary condition callbacks to identify where on the boundary the point is:

| Value | Meaning |
|-------|---------|
| `NE` | North-East corner |
| `NW` | North-West corner |
| `SE` | South-East corner |
| `SW` | South-West corner |
| `N` | North edge |
| `S` | South edge |
| `W` | West edge |
| `E` | East edge |
| `INTERNAL` | Interior point |
| `free_surface` | Free surface point |

## Boundary Condition Type Codes

| Code | Meaning |
|------|---------|
| `0` | Dirichlet — matches the physical boundary (Vx: left/right, Vz: bottom/top) |
| `11` | Dirichlet — does not match the physical boundary (Vx: bottom/top, Vz: left/right) |
| `2` | Neumann — does not match the physical boundary |
| `13` | Neumann — matches the physical boundary |
| `-2` | Periodic in x direction |
| `-1` | Not a BC point (interior) |
| `30` | Not calculated (air above free surface) |

## Built-in BC Templates

Pre-built boundary conditions you can use directly:

| Function | Description |
|----------|-------------|
| `SetPureShearBCVx` / `SetPureShearBCVz` | Pure shear (extension/compression) |
| `SetSimpleShearBCVx` / `SetSimpleShearBCVz` | Simple shear |
| `SetPureOrSimpleShearBCVx` / `SetPureOrSimpleShearBCVz` | Combined pure + simple shear |

Usage: `.SetBCVx = SetPureShearBCVx` in the `SetBCs_ff` struct.

### free_surface and Top Boundary Behaviour

| `free_surface` | Top boundary condition | Use for |
|----------------|----------------------|---------|
| `0` | Free-slip (v·n = 0) — impenetrable wall | Closed-box benchmarks (Blankenbach, Rayleigh-Bénard) |
| `1` | Zero traction (σ·n = 0) — material can flow through | Geodynamic models with topography evolution |

**Warning**: Using `free_surface=1` for closed-box convection benchmarks produces ~35% error in Vrms because the top boundary is not impenetrable. When `free_surface=1` is enabled, extra air cells are added above the domain, and `BuildInitialTopography` must be provided in the setup struct.

## Geometry Helpers

### Ellipse

```c
typedef struct {
    double centreX, centreZ;
    double radiusX, radiusZ;
    double angle;
} Ellipse;

bool IsEllipseCoordinates(Coordinates coords, Ellipse ellipse, double scalingL);
```

### Rectangle

```c
typedef struct {
    double centreX, centreZ;
    double sizeX, sizeZ;
    double angle;
} Rectangle;

bool IsRectangleCoordinates(Coordinates coords, Rectangle rect, double scalingL);
```

**Example** — placing an elliptical inclusion in `SetPhase`:
```c
int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    Ellipse inclusion = {
        .centreX = 0.0 / input->scaling.L,
        .centreZ = -30e3 / input->scaling.L,
        .radiusX = 5e3 / input->scaling.L,
        .radiusZ = 5e3 / input->scaling.L,
        .angle   = 0.0,
    };
    if (IsEllipseCoordinates(coordinates, inclusion, input->scaling.L)) {
        return 1;  // inclusion phase
    }
    return 0;  // background phase
}
```

## The .txt Parameter File

### Section Reference

#### RESTART
| Parameter | Description |
|-----------|-------------|
| `istep` | Step number for restart |
| `irestart` | 0 = fresh start, 1 = restart from checkpoint |

#### OUTPUT FILES
| Parameter | Description |
|-----------|-------------|
| `writer` | Enable output writing (0/1) |
| `writer_step` | Write output every N steps |
| `writer_markers` | Write marker data (0/1) |
| `writer_debug` | Write debug fields (0/1) |

#### SCALES
| Parameter | Description | Typical values |
|-----------|-------------|----------------|
| `eta` | Viscosity scale [Pa·s] | 1e20–1e23 |
| `L` | Length scale [m] | 1e3–1e6 |
| `V` | Velocity scale [m/s] | 1e-10–1e-8 |
| `T` | Temperature scale [K] | 1e2 |

#### SPACE-TIME
| Parameter | Description |
|-----------|-------------|
| `Nx`, `Nz` | Grid resolution (horizontal, vertical) |
| `Nt` | Number of time steps |
| `xmin`, `xmax`, `zmin`, `zmax` | Domain bounds [m] |
| `dt` | Initial time step [s] |
| `Courant` | Courant number (typically 0.3–0.5) |
| `penalty` | Pressure penalty parameter |
| `eta_average` | Viscosity averaging: 0=arithmetic, 1=harmonic, 2=geometric |
| `interp_stencil` | Interpolation stencil width |

#### SWITCHES
| Parameter | Values | Description |
|-----------|--------|-------------|
| `mechanical` | 0/1 | Enable mechanical solver |
| `thermal` | 0/1 | Enable thermal solver |
| `elastic` | 0/1 | Enable elasticity |
| `free_surface` | 0/1 | Enable free surface |
| `free_surface_stab` | 0/1 | Free surface stabilisation (FSSA) |
| `pure_shear_ALE` | -1/0/1 | ALE mode: 1=stretch, -1=fixed box |
| `shear_style` | 0/1 | 0=pure shear, 1=periodic simple shear (auto-sets `periodic_x=1`) |
| `periodic_x` | 0/1 | Periodic boundaries in x (set automatically by `shear_style`) |
| `shear_heating` | 0/1 | Shear heating source term |
| `adiab_heating` | 0/1 | Adiabatic heating |
| `subgrid_diffusion` | 0/1/2 | Particle subgrid thermal diffusion mode |
| `surface_processes` | 0/1/2/3 | 0=none, 1=diffusion, 2=fill, 3=diffusion+source |
| `anisotropy` | 0/1 | Enable mechanical anisotropy |
| `constant_dt` | 0/1 | 0=adaptive, 1=fixed time step |
| `RK` | 1/2/4 | Runge-Kutta order for advection |
| `stress_rotation` | 0/1/2 | Stress objectivity: 0=none, 1=Jaumann, 2=analytical |
| `conserv_interp` | 0/1 | Conservative velocity interpolation (Taras scheme) |

#### SETUP DEPENDANT
| Parameter | Description |
|-----------|-------------|
| `bkg_strain_rate` | Background strain rate [s⁻¹] (negative = extension) |
| `user0`–`user9` | User-defined parameters accessible in callbacks |

#### GRAVITY
| Parameter | Description |
|-----------|-------------|
| `gx` | Horizontal gravity [m/s²] (usually 0) |
| `gz` | Vertical gravity [m/s²] (usually -9.81) |

#### PHASE PROPERTIES
```
Nb_phases = 3   // Total number of phases

ID   = 0        // Phase index (0-based)
rho  = 2700.00  // Reference density [kg/m³]
G    = 1e10     // Shear modulus [Pa]
Cp   = 1050     // Heat capacity [J/kg/K]
k    = 2.5      // Thermal conductivity [W/m/K]
Qr   = 1e-6     // Radiogenic heat production [W/m³]
C    = 1e7      // Cohesion [Pa]
phi  = 30       // Friction angle [degrees]
alp  = 32.0e-6  // Thermal expansivity [K⁻¹]
bet  = 1e-11    // Compressibility [Pa⁻¹]
cstv = 0        // Constant viscosity flag (0=off, >0=database ID)
pwlv = 40       // Power law viscosity flag (0=off, >0=database ID)
linv = 40       // Linear viscosity flag (0=off, >0=database ID)
gbsv = 0        // Grain boundary sliding flag
expv = 40       // Peierls/exponential creep flag
aniso_factor = 1.0  // Anisotropy factor (1.0 = isotropic)
aniso_angle  = 45   // Initial anisotropy angle [degrees]
```

## Key Source Files

- Public API and types: `MDLIB/include/mdoodz.h`
- Problem initialisation: `MDLIB/Setup.c`
- Example scenario: `SETS/RiftingChenin.c` + `SETS/RiftingChenin.txt`
- Configuration parsing: `MDLIB/InputOutput.c`

## Verification Workflow

Before running a production simulation, create coarse-resolution verification `.txt` files to catch divergence or instability early.

### Two-Tier Verification Pattern

1. **Copy the production `.txt` to create verification variants**:
   ```bash
   cp SETS/MyModel.txt SETS/MyModel_quick.txt
   cp SETS/MyModel.txt SETS/MyModel_medium.txt
   ```

2. **Edit the verification files** — reduce resolution and time steps, keep all physics identical:

   | Parameter | Production | Quick | Medium |
   |-----------|-----------|-------|--------|
   | `Nx` | 200+ | 101 | 101 |
   | `Nz` | 160+ | 81 | 81 |
   | `Nt` | 1000+ | 10 | 100 |

3. **Run quick verification** (10 steps) — catches immediate crashes, NaN, bad parameters:
   ```bash
   make build SET=MyModel TXT=MyModel_quick.txt
   make run SET=MyModel > logs_quick.txt 2>&1
   ```

4. **Check the log** for problems:
   ```bash
   grep -i 'nan\|diverge\|error\|exit' logs_quick.txt
   grep 'Ending' logs_quick.txt          # verify it reached the end
   grep 'T (°C)' logs_quick.txt           # check temperature range is physical
   ```

5. **Run medium verification** (100 steps) — catches slower instabilities (e.g., `eta_vp` too low with melting):
   ```bash
   make build SET=MyModel TXT=MyModel_medium.txt
   make run SET=MyModel > logs_medium.txt 2>&1
   ```

6. **If both pass**, proceed to production resolution.

### What To Watch For

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| NaN in output | Unstable rheology or strain rate | Check `eta_vp`, `min_eta`/`max_eta`, flow law indices |
| Negative viscosity | `eta_vp` too low with melting | Increase `eta_vp` (try 2.5e20 for lithospheric scale) |
| Values ~1e60 | Numerical blow-up | Same as above; also check `dt`, `Courant` |
| Temperature out of range | Bad thermal BCs or initial T | Verify `SetBCT` callback and `SetTemperature` |
| Solver does not converge | Viscosity contrast too high | Reduce `max_eta`/`min_eta` ratio, try `Picard2Newton` |
