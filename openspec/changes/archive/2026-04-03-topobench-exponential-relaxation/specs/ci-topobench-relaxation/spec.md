## ADDED Requirements

### Requirement: Topographic relaxation amplitude decay
The test suite SHALL verify that a sinusoidal free-surface perturbation decays exponentially under gravity, matching the analytical relaxation time τᵣ = 4πη/(ρgλ) · coth(2πH/λ).

#### Scenario: Amplitude decay follows analytical exponential
- **WHEN** a simulation runs with `free_surface=1`, uniform viscosity η=1e21 Pa·s, ρ=3300 kg/m³, g=10 m/s², initial perturbation h₀·cos(2πx/λ) with h₀=7 km, λ=2800 km, domain depth H=700 km, for 20 time steps
- **THEN** the relative error between the measured amplitude ratio h(t)/h₀ and the analytical exp(−t/τᵣ) SHALL be less than 15% at each output step after the first

#### Scenario: Topography remains finite and smooth
- **WHEN** the TopoBench simulation completes all time steps
- **THEN** all values in `Topo/z_grid` SHALL be finite (no NaN or Inf) at every output step

#### Scenario: Amplitude decreases monotonically
- **WHEN** the TopoBench simulation runs for 20 steps
- **THEN** the peak surface elevation max(z_grid) SHALL decrease monotonically from step to step

### Requirement: Topographic relaxation grid convergence
The test suite SHALL verify that the relaxation error decreases with grid refinement, demonstrating numerical convergence of the free-surface implementation.

#### Scenario: Monotonic error decrease with resolution
- **WHEN** the TopoBench simulation runs at resolutions 31×16, 51×26, and 101×51 (aspect ratio ~2:1 matching the 2800×700 km domain)
- **THEN** the mean relative error of h(t)/h₀ vs analytical SHALL decrease monotonically from coarsest to finest grid

#### Scenario: Convergence order
- **WHEN** pair-wise convergence orders are computed as log(error_coarse / error_fine) / log(h_coarse / h_fine)
- **THEN** the convergence order SHALL be at least 0.5 for each resolution pair
