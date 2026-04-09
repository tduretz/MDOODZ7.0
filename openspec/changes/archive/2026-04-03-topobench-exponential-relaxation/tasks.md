## 1. Parameter Files

- [x] 1.1 Create `TESTS/TopoBench/TopoBenchRelaxation.txt` — base 51×51 grid, uniform η=1e21, ρ=3300, g=−10, free_surface=1, free_surface_stab=1, 20 steps dt=1e10 s, writer_step=1, domain 2800×750 km
- [x] 1.2 Create `TESTS/TopoBench/TopoBenchConvergence31.txt` — 31×16 variant
- [x] 1.3 Create `TESTS/TopoBench/TopoBenchConvergence101.txt` — 101×51 variant

## 2. Test Implementation

- [x] 2.1 Create `TESTS/TopoBenchTests.cpp` with `TopoBench` fixture: SetSurfaceZCoord (h₀·cos(2πx/λ)), SetPhase (single phase), free-slip BCs, SetBCPType
- [x] 2.2 Add `TopoBenchRelaxation` test: run simulation, extract max(z_grid) at each step, verify exponential decay against τᵣ = 4ηk/(ρg)·coth(kH), assert relative error < 15%, monotonic decay, finite values
- [x] 2.3 Add `TopoBenchConvergence` test: run at 3 resolutions, compute mean relative error at each, assert monotonic decrease and overall convergence order ≥ 0.3

## 3. Build Integration

- [x] 3.1 Register `TopoBenchTests` executable in `TESTS/CMakeLists.txt`

## 4. Build and Validate

- [x] 4.1 Build and run: both tests pass. Relaxation mean error ~1.9%, convergence order 0.34

## 5. Documentation

- [x] 5.1 Add §4 to `TESTS/AnalyticalSolutions.md` documenting the TopoBench relaxation analytical solution, equal-viscosity interface correction (2× factor), and measured accuracy
- [x] 5.2 Update `skill-testing-guide` with new test entries
