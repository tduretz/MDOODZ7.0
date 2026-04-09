## 1. Scenario Setup

- [x] 1.1 Create `SETS/BlankenBench.c` with callbacks: SetPhase (single phase), SetTemperature (linear + cosine perturbation), SetBCVx/SetBCVz (free-slip), SetBCT (Dirichlet T=0 top / T=1 bottom, Neumann sides)
- [x] 1.2 Create `TESTS/BlankenBench/BlankenBench.txt` parameter file with unit scaling (η=1, L=1, V=1, T=1), g=Ra=1e4, constant viscosity, thermal=1, 41×41 grid, sufficient Nt/dt for steady state

## 2. Test Implementation

- [x] 2.1 Create `TESTS/BlankenBenchTests.cpp` with fixture and three test cases: NusseltAndVrms (compute Nu from top dT/dz, Vrms from Vx/Vz grids), TemperatureProfile (mid-depth T_max), all asserting <5% relative error against published values
- [x] 2.2 Register `BlankenBenchTests` executable in `TESTS/CMakeLists.txt`, add `BlankenBench/` to copy list
- [x] 2.3 Build and run in WSL2; iterate on Nt/dt until steady state is reached and all assertions pass

## 3. Gnuplot Visualization

- [x] 3.1 Create `TESTS/BlankenBench/plot_convection.gp` — standalone script producing temperature field PNG with heat palette (pm3d/image) and velocity vector overlay, with usage comments
- [x] 3.2 Run the script on actual simulation output, view the resulting PNG, and verify: convection cell is clearly visible, colour bar labels are readable, velocity arrows show coherent flow pattern, and the plot is presentable (proper aspect ratio, no clipping, legible text)

## 4. Documentation

- [x] 4.1 Add §5 "Blankenbach Thermal Convection" to `TESTS/AnalyticalSolutions.md` documenting the benchmark parameters, published values, Nusselt/Vrms formulas, measured errors, and gnuplot visualization instructions
- [x] 4.2 Add BlankenBench entries to the test suite table in `.github/skills/skill-testing-guide/SKILL.md`
- [x] 4.3 Add BlankenBench entry to the scenario gallery in `.github/skills/skill-scenario-gallery/SKILL.md` under a Benchmarks or Thermal section
- [x] 4.4 Create `.github/skills/skill-thermal/SKILL.md` documenting the MDOODZ thermal solver: energy equation, thermal BCs (SetBCT callback, type=0/1), conductivity/Cp/α parameters, temperature advection via markers, Boussinesq coupling (ρ = ρ₀(1-αΔT)), scaling of thermal quantities, shear/adiabatic heating switches, and practical guidance for setting up thermally-coupled simulations
