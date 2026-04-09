## Why

The SolVi benchmark (Schmid & Podladchikov 2003) reveals that MDOODZ achieves ~1.0 convergence order for velocity and only ~0.75 for pressure near viscosity discontinuities. The dominant error source is the FD stencil at cells straddling the inclusion boundary — arithmetic averaging of viscosity across cell faces smears the discontinuity.

Experimental testing (`TESTS/SolViMarkerComparison.cpp`) has ruled out two other approaches:
- **Increasing markers per cell** (4×4 → 8×8 → 16×16): gives 0% Vx improvement and only ~12% P reduction at fixed resolution, while the P convergence order actually *worsens* (0.75 → 0.59). Cost: up to 6.8× runtime. The error is dominated by the FD truncation at the discontinuity, not by marker interpolation noise.
- **Changing `eta_average`** (arithmetic/harmonic/geometric): controls marker-to-node phase mixing *within* each cell, not cell-face interpolation in the stencil. Geometric gives 18% better absolute P but no convergence order gain.

Literature (Deubelbeiss & Kaus 2008, Schmid & Podladchikov 2003) shows harmonic averaging at cell faces can improve pressure convergence order by +0.3–0.5 with negligible runtime cost. This is the remaining lever that hasn't been tried.

## What Changes

- Add a new model parameter (`cell_avg`) to select the cell-face viscosity averaging method used when the Stokes stencil reads viscosity from neighbouring cells: `0` = direct (current behaviour, no change), `1` = harmonic mean at cell faces
- Modify the Stokes assembly (`StokesAssemblyDecoupled.c`) to compute harmonic averages of the constitutive tensor components (`D11_n`, `D33_s`, etc.) at cell faces before inserting them into the stencil, when `cell_avg = 1`
- Modify the Jacobian assembly (`FD_Jacobian.c`) to use the same averaging, ensuring Newton convergence is consistent
- Add SolVi L2 benchmark tests comparing `cell_avg = 0` vs `cell_avg = 1` at multiple resolutions
- Document the parameter and findings in skills

## Capabilities

### New Capabilities
- `cell-face-averaging`: The cell-face harmonic averaging feature — parameter, stencil modification, Jacobian consistency, and L2 benchmark verification

### Modified Capabilities
- `skill-solvers`: Add documentation of cell-face averaging parameter and its effect on convergence orders
- `skill-testing-guide`: Add experimental comparison data (L2 tables for cell_avg=0 vs cell_avg=1)

## Impact

- **Core solver code**: `MDLIB/StokesAssemblyDecoupled.c` (~2390 lines) — the stencil reads `D11E`, `D33N` etc. from `mesh->D_xx_n[index]`. The averaging inserts a harmonic mean step between the mesh array lookup and the stencil coefficient computation
- **Jacobian**: `MDLIB/FD_Jacobian.c` — must use the same averaged values for consistent Newton linearisation
- **Constitutive tensor**: `MDLIB/StokesRoutines.c` — computes `D11_n`, `D33_s` from `eta_n`, anisotropy factor, director. These arrays are the input to the averaging
- **Parameter reading**: `MDLIB/InputOutput.c` — add `cell_avg` parameter
- **Model struct**: `MDLIB/mdoodz-private.h` — add `cell_avg` field
- **Tests**: `TESTS/SolViBenchmarkTests.cpp` or new test file — L2 comparison at 41/51/81 resolutions
- **No breaking changes**: `cell_avg` defaults to `0`, preserving all existing behaviour
