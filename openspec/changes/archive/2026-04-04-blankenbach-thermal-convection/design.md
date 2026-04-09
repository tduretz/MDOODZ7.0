## Context

The test suite has 47 GTests but zero thermo-mechanical coupling tests. The Blankenbach et al. (1989) Case 1a benchmark is the community standard for validating isoviscous Rayleigh-Bénard convection codes. Published diagnostics (Nu = 4.884, Vrms = 42.865, T_max = 0.4266 at mid-depth) exist from multiple independent codes.

MDOODZ already has all the physics: Boussinesq buoyancy (`ρ = ρ₀(1 - αΔT)` in RheologyDensity.c), Stokes solver, energy equation with advection (Main_DOODZ.c marker advection of T), and Dirichlet thermal BCs (SetBCT callback with type=1). No new physics code is needed.

HDF5 output stores: T on centres `(Nx-1)×(Nz-1)` as `Centers/T`, Vx on `Nx×(Nz+1)` as `VxNodes/Vx`, Vz on `(Nx+1)×Nz` as `VzNodes/Vz`. All in SI units (already re-scaled).

## Goals / Non-Goals

**Goals:**
- Validate Boussinesq thermal-buoyancy coupling, temperature advection, and long-integration stability
- Assert Nusselt number, RMS velocity, and mid-depth max temperature against published Blankenbach Case 1a values
- Produce a presentable gnuplot temperature-field visualization with velocity vectors
- Keep CI runtime under ~60 seconds

**Non-Goals:**
- No temperature-dependent viscosity (Case 1b/2a) — isoviscous only
- No grid convergence study (single resolution is sufficient for the benchmark)
- No time-series tracking of transient evolution — only steady-state diagnostics
- No new TestHelpers utilities beyond what's needed for Nu/Vrms computation

## Decisions

### 1. Non-dimensional Blankenbach parameters via MDOODZ scaling

The Blankenbach benchmark is defined non-dimensionally: unit square, Ra = αρgΔTL³/(κη) = 10⁴. MDOODZ uses its own scaling system (η₀, L₀, V₀, T₀). We choose SI parameters that naturally give Ra = 10⁴:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| L = 1e6 m | Domain size | Convenient mantle scale |
| η = 1e22 Pa·s | Viscosity | Typical mantle |
| ρ = 3300 kg/m³ | Density | Standard mantle |
| α = 1e-5 K⁻¹ | Thermal expansion | Standard olivine |
| ΔT = 1000 K | T_bot - T_top | user1 - user0 |
| Cp = 1250 J/kg/K | Heat capacity | Chosen to set κ |
| k = 3.3 W/m/K | Conductivity | κ = k/(ρCp) ≈ 8e-7 m²/s |
| g = 10 m/s² | Gravity | Standard |

Check: Ra = αρgΔTL³/(κη) = 1e-5 × 3300 × 10 × 1000 × (1e6)³ / (8e-7 × 1e22) = 3.3e23 / 8e15 ≈ 4.1e7. That's too high.

**Revised approach:** Use smaller domain or adjust parameters. The simplest way: set L = 1e5 m (100 km) with appropriately adjusted parameters, or directly use near-unity scaling: η=1, L=1, g=Ra, α=1, ρ=1, ΔT=1, κ=1. MDOODZ's scaling allows this.

Actually, the cleanest approach: choose physical parameters so Ra = 10⁴ exactly. With L=1e6 m, η=1e22, ρ=3300, α=1e-5, ΔT=1000, g=10, k=3.3, Cp=1250: κ = k/(ρCp) = 3.3/(3300×1250) = 8e-7. Ra = 1e-5 × 3300 × 10 × 1000 × 1e18 / (8e-7 × 1e22) = 3.3e17 / 8e15 = 41.25. Too low with L=1e6.

Let me fix this properly. Ra = αρgΔTL³/(κη). Target Ra = 1e4.

Choose: α=3e-5, ρ=3300, g=10, ΔT=2500 K, L=1e6 m, η=1e21 Pa·s, k=3.564 W/m/K, Cp=1000.
κ = k/(ρCp) = 3.564/(3300×1000) = 1.08e-6.
Ra = 3e-5 × 3300 × 10 × 2500 × 1e18 / (1.08e-6 × 1e21) = 2.475e21 / 1.08e15 = 2.29e6. Still wrong scale.

The fundamental issue is that mantle-realistic parameters give either Ra ≈ 10⁶–10⁸. For Ra = 10⁴, we need lab-scale or artificial parameters. The simplest: use MDOODZ scaling to make internal values near-unity.

**Final approach:** Set SCALES to make the problem near-unity internally. Use:
- `eta = 1`, `L = 1`, `V = 1`, `T = 1` — unit scaling
- Domain: xmin=-0.5, xmax=0.5, zmin=-0.5, zmax=0.5 (unit square, L=1 in physical units)
- ρ = 1, α = 1, g = Ra = 1e4, ΔT = 1 (T_top=0, T_bot=1)
- η = 1 (constant), k = 1, Cp = 1 → κ = k/(ρCp) = 1
- Ra = αρgΔTL³/(κη) = 1 × 1 × 1e4 × 1 × 1 / (1 × 1) = 1e4 ✓

This is the standard non-dimensional Blankenbach formulation. The "gravity" parameter absorbs Ra.

**Alternative considered:** Physical SI parameters. Rejected — the algebra to match Ra exactly is error-prone and the non-dimensional formulation is the standard approach used by Blankenbach et al. themselves.

### 2. Resolution: 41×41 for CI, with option for 81×81

At 41×41, published results show Nu accurate to ~1–2% for Case 1a. This balances CI runtime (<60s) with accuracy. The test threshold will be 3% for Nu and Vrms.

**Alternative considered:** 21×21 (faster but ~5% error), 81×81 (more accurate but ~4× slower). 41×41 is the sweet spot.

### 3. Steady-state detection: fixed number of time steps

Run for a large number of steps (Nt=500) with a moderately large dt. The simulation will reach thermal steady state when the convective cell stabilises. We validate diagnostics from the final output only.

**Alternative considered:** Adaptive steady-state detection (compare Nu at consecutive outputs). Adds complexity and the fixed-step approach is simpler and deterministic for CI.

### 4. Nusselt number from finite-difference dT/dz at top

Nu = -(1/ΔT) × mean(dT/dz|_{z=top}) × L. Since T is stored at cell centres, dT/dz at the top boundary is computed using the top row of centres and the imposed T_top Dirichlet BC:

dT/dz|_top ≈ (T_top - T_centre_top_row) / (dz/2)

Integrate over x and normalise by L and ΔT. This is a standard one-sided finite-difference approach.

**Alternative considered:** Reading heat flux from solver internals. Not available in HDF5 output; the FD approach from stored T is standard.

### 5. RMS velocity from Vx and Vz grids

Vrms = sqrt(mean(Vx² + Vz²)) over the domain. Since Vx and Vz live on staggered grids (Vx: Nx×(Nz+1), Vz: (Nx+1)×Nz), we compute RMS on their respective grids and combine: Vrms = sqrt(mean(Vx²) + mean(Vz²)). For a uniform grid this is equivalent to interpolating to centres.

**Alternative considered:** Interpolate both to centres first. Unnecessary for the RMS integral on a uniform grid.

### 6. Gnuplot visualization: temperature field with velocity arrows

The gnuplot script reads the final HDF5 output (via h5dump to ASCII or direct binary), plots the temperature field as a colour map using `pm3d` with a red-yellow-white heat palette, and overlays velocity vectors using `with vectors`. Output as PNG.

The script will be standalone (runnable manually, not part of CI assertions) and will live in `TESTS/BlankenBench/plot_convection.gp`.

**Alternative considered:** Julia/Makie visualisation. While more powerful, gnuplot has no dependencies beyond gnuplot itself and matches the user's request for a presentable plot with a heat colourmap.

### 7. Initial condition: conductive profile with perturbation

Start with a linear temperature profile T(z) = T_bot + (T_top - T_bot) × (z - z_bot)/H plus a small cos(πx)×sin(πz) perturbation to seed the single convection cell. Without the perturbation, the symmetric initial condition would take very long to break symmetry.

### 8. Free-slip BCs on all walls

All four walls get free-slip velocity BCs (standard for Blankenbach). Top and bottom get Dirichlet T BCs (T=0, T=1). Left and right get zero heat flux (Neumann, type=0 in SetBCT).

## Risks / Trade-offs

- **[Risk] CI runtime ~30–60s** → Acceptable for a coupled benchmark. The 41×41 grid with 500 steps is the minimum viable resolution. If too slow, reduce to Nt=300 and accept slightly looser thresholds.
- **[Risk] Steady state not reached in 500 steps** → The non-dimensional diffusion time is L²/κ = 1. With dt ≈ 2e-4 (Courant-limited), 500 steps = 0.1 diffusion times. May need dt ≈ 5e-3 (larger, implicit thermal solver handles it) for 500 steps = 2.5 diffusion times, sufficient for steady state.
- **[Risk] Nusselt number accuracy at 41×41** → Published community results at 41×41 typically achieve ~1% for Case 1a. Our 3% threshold is conservative.
- **[Risk] Gnuplot not installed in CI** → The plot script is manual-only, not part of GTest. No CI dependency.
- **[Risk] HDF5 field names for T/Vx/Vz assumed** → Verified from HDF5Output.c: `Centers/T`, `VxNodes/Vx`, `VzNodes/Vz`.
