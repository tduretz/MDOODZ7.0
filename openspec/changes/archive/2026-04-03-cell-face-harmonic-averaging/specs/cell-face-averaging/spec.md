## ADDED Requirements

### Requirement: cell_avg parameter controls cell-face viscosity averaging
The system SHALL read an integer parameter `cell_avg` from the `.txt` parameter file with default value `0`. When `cell_avg = 0`, the Stokes stencil SHALL use direct (unaveraged) D-tensor values from neighbouring grid points, preserving current behaviour. When `cell_avg = 1`, the stencil SHALL use harmonic-averaged D-tensor values at cell faces.

#### Scenario: Default behaviour unchanged
- **WHEN** `cell_avg` is not specified in the `.txt` file
- **THEN** the solver SHALL behave identically to the current code (direct D-tensor reads, no averaging)

#### Scenario: Parameter value 0 explicit
- **WHEN** `cell_avg = 0` is set in the `.txt` file
- **THEN** the solver SHALL use direct D-tensor reads (same as default)

#### Scenario: Parameter value 1 enables harmonic averaging
- **WHEN** `cell_avg = 1` is set in the `.txt` file
- **THEN** the Stokes stencil SHALL compute harmonic means of D11, D22, D33 at each cell face before inserting them into the stencil coefficients

### Requirement: Harmonic averaging applied to x-momentum stencil
In `Xmomentum_InnerNodesDecoupled`, when `cell_avg = 1`, the stencil SHALL replace the direct reads `D11E = mesh->D11_n[iPrE]` and `D11W = mesh->D11_n[iPrW]` with a single harmonic mean value `D11E = D11W = 2·D11_E·D11_W / (D11_E + D11_W)`. The same SHALL apply to `D33N` and `D33S` using their respective mesh values.

#### Scenario: X-momentum with uniform viscosity
- **WHEN** `cell_avg = 1` and all cells have the same viscosity η
- **THEN** D11E = D11W = 2η (identical to the non-averaged case) and L2 errors SHALL be within 1% of the `cell_avg = 0` result

#### Scenario: X-momentum with viscosity discontinuity
- **WHEN** `cell_avg = 1` and the SolVi inclusion benchmark is run (η_c/η_m = 1000)
- **THEN** the pressure L2 error SHALL be lower than with `cell_avg = 0` at the same resolution

### Requirement: Harmonic averaging applied to z-momentum stencil
In `Zmomentum_InnerNodesDecoupled`, when `cell_avg = 1`, the stencil SHALL replace `D22S`/`D22N` with their harmonic mean and `D33E`/`D33W` with their harmonic mean, using the same formula as the x-momentum stencil.

#### Scenario: Z-momentum symmetry preserved
- **WHEN** `cell_avg = 1` and the SolVi benchmark is run with a centred inclusion
- **THEN** the Vx and Vz L2 errors SHALL be of comparable magnitude (symmetry not broken by the averaging)

### Requirement: Division-by-zero protection
When computing the harmonic mean, if both D values sum to less than 1e-30, the result SHALL be 0.0 instead of performing the division.

#### Scenario: Air cell neighbour
- **WHEN** `cell_avg = 1` and one neighbour is an air cell (D = 0)
- **THEN** the harmonic mean SHALL be 0.0 and the existing `inE`/`inW`/`inN`/`inS` flags SHALL zero out the contribution in the stencil (no crash, no NaN)

### Requirement: Jacobian consistency
The Jacobian assembly SHALL use the same cell-face averaging as the Stokes assembly. Since `FD_Jacobian.c` calls the same `Xmomentum_InnerNodesDecoupled` and `Zmomentum_InnerNodesDecoupled` functions, this SHALL be automatic with no separate modification.

#### Scenario: Newton convergence with cell_avg=1
- **WHEN** `cell_avg = 1` and the SolVi benchmark is run with Newton iteration
- **THEN** the solver SHALL converge (final residual < 1e-8) and the number of Newton iterations SHALL not increase by more than 50% compared to `cell_avg = 0`

### Requirement: SolVi L2 benchmark comparison
A test SHALL compare L2 errors for `cell_avg = 0` and `cell_avg = 1` on the SolVi benchmark at resolutions 41×41, 51×51, and 81×81. The test SHALL report L2(Vx), L2(P), and convergence orders for both modes.

#### Scenario: Pressure convergence order improvement
- **WHEN** the SolVi benchmark is run at 41×41 and 81×81 with `cell_avg = 1`
- **THEN** the pressure convergence order SHALL be at least 0.8 (improved from ~0.75 with `cell_avg = 0`)

#### Scenario: Velocity accuracy not degraded
- **WHEN** `cell_avg = 1` on SolVi at 51×51
- **THEN** the Vx L2 error SHALL not be more than 20% worse than with `cell_avg = 0`

### Requirement: Parameter field added to model struct
A field `int cell_avg` SHALL be added to the `params` struct in `mdoodz-private.h`. It SHALL be read from the `.txt` file in `InputOutput.c` with default value 0.

#### Scenario: Struct field accessible in stencil
- **WHEN** the stencil assembly functions receive the `model` parameter
- **THEN** `model.cell_avg` SHALL be accessible and determine the averaging mode

### Requirement: Only diagonal D-tensor terms averaged
Harmonic averaging SHALL apply only to D11_n, D22_n, and D33_s. Off-diagonal terms (D12_n, D13_n, D14_n, D23_n, D24_n, D31_s, D32_s, D34_s) SHALL remain as direct reads.

#### Scenario: Isotropic case cross-terms zero
- **WHEN** `anisotropy = 0` and `cell_avg = 1`
- **THEN** cross-terms are zero regardless of averaging mode, so behaviour SHALL be identical to averaging all D components
