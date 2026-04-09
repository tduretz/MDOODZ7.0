# MDOODZ 7.0 — Copilot Instructions

## Project Overview

MDOODZ is a 2D Visco-Elasto-Plastic Thermo-Mechanical (VEPTM) geodynamic modelling code written in C. It solves the Stokes equations coupled with the energy equation using a staggered-grid finite-difference method with marker-in-cell advection.

## Repository Layout

- `MDLIB/` — Core C library (36 modules): solver, rheology, advection, I/O, free surface, anisotropy
- `SETS/` — 79 predefined simulation scenarios (paired `.c` + `.txt` files)
- `JuliaVisualisation/` — Julia/Makie scripts for post-processing HDF5 output
- `.github/skills/` — Domain-specific Copilot skills (see below)

## Key Conventions

- **Internal units are non-dimensional.** The code scales all quantities using a reference viscosity (η), length (L), velocity (V), and temperature (T). HDF5 output is written in SI units (already re-scaled).
- **Staggered grid naming**: `_n` = cell centres, `_s` = cell vertices, `_0` = previous time step.
- **Creep abbreviations**: `pwl` = power-law/dislocation, `lin` = linear/diffusion, `exp` = Peierls, `gbs` = grain-boundary sliding, `cst`/`cstv` = constant viscosity.
- **Phase = -1** means air (free surface). Mask these cells with NaN before plotting.

## When to Use Skills

Before answering questions about MDOODZ, check whether a relevant skill exists. The skills contain verified, codebase-specific knowledge that is more accurate than general reasoning.

| User is asking about... | Use skill |
|------------------------|-----------|
| Building, compiling, dependencies, CMake | skill-build-and-run |
| Creating a new simulation, .c/.txt files, callbacks, boundary conditions | skill-model-setup |
| Viscosity, creep laws, flow laws, plasticity, elasticity, material properties | skill-rheology |
| Partial melting, melt fraction, solidus/liquidus, melt weakening | skill-melting |
| Anisotropy, fabric, director field, foliation | skill-anisotropy |
| Solver convergence, Newton, Picard, time stepping, penalty | skill-solvers |
| Governing equations, Stokes, energy, scaling, advection, free surface | skill-physics-theory |
| Thermal solver, energy equation, thermal BCs, conductivity, heating, Boussinesq | skill-thermal |
| Variable names, struct fields, what a code symbol means | skill-code-glossary |
| Which scenario to use, example setups, benchmarks | skill-scenario-gallery |
| Plotting, HDF5 output, Julia visualisation, extracting data | skill-visualisation |
| Validating .txt parameters, checking physical ranges | skill-parameter-validation |
| Writing tests, running tests, CI tests, visual regression tests, GTest | skill-testing-guide |

## Code Style

- C11 with OpenMP (optional). No C++.
- Build system: CMake 3.16+. Convenience targets via top-level `makefile`.
- External dependencies: SuiteSparse (CHOLMOD/UMFPACK), HDF5, BLAS/LAPACK.
- Simulation setup uses callback functions registered through the `MdoodzSetup` struct in `MDLIB/include/mdoodz.h`.
