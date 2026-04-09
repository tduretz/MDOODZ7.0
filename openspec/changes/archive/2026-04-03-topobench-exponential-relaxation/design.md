## Context

The CI suite lacks quantitative transient free-surface verification. The existing `SinkingBlock` test checks only that the surface deflects (Vz ≠ 0, P > 0) — a qualitative smoke test. TopoBench Case 1 (Crameri et al., 2012) is the standard benchmark: a sinusoidal cosine perturbation h₀·cos(kx) on a viscous half-space decays exponentially under gravity. MDOODZ already ships a production scenario (`SETS/TopoBenchCase1.c/.txt`) with the required `SetSurfaceZCoord` callback, BCs, and two-phase setup.

## Goals / Non-Goals

**Goals:**
- Add `TopoBenchRelaxation` test: run multi-step simulation, extract amplitude from `Topo/z_grid` at each output step, compare against analytical exp(−t/τᵣ) decay
- Add `TopoBenchConvergence` test: run at 2–3 resolutions, verify monotonic error decrease
- Use finite-depth corrected relaxation time τᵣ = 4πη/(ρgλ) · coth(2πH/λ)
- Document analytical solution in AnalyticalSolutions.md

**Non-Goals:**
- Multi-layer (Case 2/3) benchmarks
- Viscosity contrast between lithosphere and mantle in the analytical comparison (use uniform η in the test)
- Free-surface algorithm improvements
- Sub-percent accuracy targets (free surface is inherently noisy from marker chain)

## Decisions

### 1. Uniform viscosity for the test (single phase η)

TopoBenchCase1 has two phases: mantle (η=1e21) and lithosphere (η=1e23). The analytical relaxation formula assumes a **uniform viscosity** half-space. The 100× viscosity contrast in the lithosphere makes the effective relaxation time deviate from the simple formula.

**Decision:** Use uniform viscosity (both phases η=1e21, drop the lithosphere layer) for the CI test. This gives a clean comparison with the analytical solution.

**Alternative considered:** Keep two phases and use an effective viscosity — rejected because there's no simple closed-form for layered viscosity relaxation, and the point is verifying the free-surface numerics, not the rheology.

### 2. Extract amplitude via max(z_grid)

The `Topo/z_grid` array (Nx floats) stores the surface elevation at each column. For a cosine perturbation, the amplitude at time t is simply max(z_grid). We compare the ratio h(t)/h₀ against exp(−t/τᵣ).

**Alternative considered:** Fit a cosine to the z_grid profile — more robust but unnecessary for a single-wavelength perturbation. max(z_grid) is exact for a pure cosine.

### 3. Finite-depth correction

The standard half-space formula τᵣ = 4πη/(ρgλ) assumes infinite depth. The TopoBench domain has finite depth H. The correction factor is coth(2πH/λ) ≈ coth(kH).

For H=700 km, λ=2800 km: kH = 2π·700/2800 = π/2 ≈ 1.57, coth(π/2) ≈ 1.09. This is a 9% correction — significant enough to matter at the 10–20% accuracy level.

**Decision:** Use τᵣ = (4πη)/(ρgλ) · coth(2πH/λ) as the analytical reference.

### 4. Domain and resolution choices

**Base test (TopoBenchRelaxation):**
- Domain: 2800 km × 700 km (one full wavelength, same as TopoBenchCase1 but symmetric about x=0 for cleaner cosine)
- Resolution: 51×51 (minimum stable for free surface, fast CI)
- 20 time steps, dt=5e10 s (~1.6 kyr), writing every step
- Total simulated time: ~1e12 s ≈ 32 kyr (roughly τᵣ/40, so ~2.5% decay)

Actually, with τᵣ ≈ (4π·1e21)/(3300·10·2.8e6)·1.09 ≈ 1.48e11 s, 20 steps of dt=5e10 s covers t=1e12 s ≈ 6.8 τᵣ, so amplitude decays to ~0.1%. That's too much — the surface will be flat.

**Revised:** Use dt=1e10 s, 20 steps → t_final = 2e11 s ≈ 1.35 τᵣ → amplitude decays to ~26%. Good dynamic range for testing.

**Convergence resolutions:** 31×16, 51×26, 101×51 (keeping aspect ratio ~2:1 matching domain 2800×700).

### 5. Threshold: relative error < 15%

Free-surface marker chain introduces noise from interpolation and remeshing. Combined with time-stepping error, expect ~5–10% error at moderate resolution. A 15% threshold gives margin while catching real regressions.

**Convergence:** Require monotonic error decrease across resolutions. Order ≥ 0.5 (free surface convergence is slower than interior Stokes due to marker chain).

### 6. New test file (TopoBenchTests.cpp)

Create a dedicated file rather than adding to `FreeSurfaceTests.cpp`. The TopoBench setup needs different callbacks (sinusoidal topo, free-slip BCs, uniform phase) — a separate test fixture keeps things clean.

## Risks / Trade-offs

- **[Risk] Free-surface marker noise** → The marker chain discretizes the surface with ~23 markers/cell. At low resolution, amplitude extraction from max(z_grid) may be noisy. Mitigation: 51×26 base resolution gives adequate marker density.
- **[Risk] CI time — multi-step at 3 resolutions** → ~15–30s total. Mitigation: keep resolution modest (max 101×51), keep time steps to 20.
- **[Risk] Analytical formula accuracy at small depth/wavelength ratio** → coth correction helps but the analytical solution still assumes Stokes flow with flat bottom BC. The numerical solution has a rigid bottom boundary. Mitigation: 15% threshold absorbs this.
- **[Tradeoff] Uniform viscosity simplifies test but doesn't exercise layered rheology** → Acceptable because the purpose is verifying free-surface numerics, not rheological coupling.
