## Why

MDOODZ7.0 is a research-grade Visco-Elasto-Plastic Thermo-Mechanical solver with 36 core modules, 100+ preconfigured scenarios, and a Julia visualization toolkit — but virtually no structured documentation exists. New users (students, collaborators, reviewers) must reverse-engineer C source code to understand model setup, rheology, solver behaviour, and output interpretation. A set of Copilot skills would provide contextual, queryable documentation that lives alongside the code and can guide users through common workflows without requiring deep familiarity with the codebase.

## What Changes

- Add a new `.github/skills/` skill for each major documentation domain (model setup, rheology, anisotropy, solvers, visualization, build/run, physics theory, scenario gallery)
- Each skill contains a `SKILL.md` with structured domain knowledge, code pointers, and workflow guidance
- Skills reference the attached anisotropic rifting paper as domain context where relevant (anisotropy, rheology, rifting scenarios)
- No changes to existing source code, build system, or runtime behaviour

## Capabilities

### New Capabilities

- `skill-model-setup`: Skill explaining how to create and configure simulations — the paired `.c` callback file and `.txt` parameter file pattern, all callback signatures (`SetPhase`, `SetTemperature`, `SetGrainSize`, `SetSurfaceZCoord`, `SetHorizontalVelocity`, `SetVerticalVelocity`, `SetBCPType`, `SetBCT`), key data structures, and step-by-step guidance for authoring new scenarios
- `skill-rheology`: Skill documenting the rheological framework — flow laws (dislocation creep, diffusion creep, Peierls, exponential), viscosity computation pipeline, elastic and plastic branches, phase-dependent material parameters, and how rheological properties are assigned through `.txt` configuration
- `skill-anisotropy`: Skill covering the anisotropy module — director-based fabric tracking, viscous anisotropy formulation, anisotropy factor definitions (δ_n, δ_s), coupling with the Stokes solver, relation to the anisotropic rifting methodology in the attached paper, and relevant scenario examples (RiftingAnisotropy, SimpleShearAniso*, AnisoHomo*)
- `skill-solvers`: Skill explaining the numerical solver architecture — staggered-grid finite difference discretisation, Stokes assembly (decoupled velocity–pressure), Newton–Raphson nonlinear iteration, thermal solver, sparse linear algebra (SuiteSparse/UMFPACK), convergence criteria, and the main time-stepping loop in `Main_DOODZ.c`
- `skill-visualisation`: Skill for the Julia/Makie visualisation toolkit — reading HDF5 output, available field types (phases, temperature, velocity, stress, strain rate, grain size, fabric), overlay options, multi-panel layouts, marker trajectory analysis, and batch export workflows
- `skill-build-and-run`: Skill documenting build and execution — CMake configuration, dependency installation (SuiteSparse, HDF5, BLAS/LAPACK), OpenMP parallelisation, make targets (`build`, `build-dev`, `run`), environment setup (`env.cmake`), and common troubleshooting
- `skill-physics-theory`: Skill summarising the governing equations and numerical methods — conservation of momentum/energy/mass, constitutive relations, non-dimensionalisation (scaling parameters η, L, V, T), marker-in-cell advection, free surface algorithm, and references to the theoretical background from the anisotropic rifting paper
- `skill-scenario-gallery`: Skill cataloguing the 100+ preconfigured scenarios in `SETS/` — organised by geodynamic application (rifting, subduction, collision, shear, thermal, chemical, benchmarks), with brief descriptions of purpose, key parameters, and expected outputs for each category
- `skill-parameter-validation`: Skill encoding physically reasonable bounds and sanity checks for `.txt` configuration parameters — covers scaling values (η typically 1e18–1e25 Pa·s, L ~1e3–1e6 m, V ~1e-12–1e-8 m/s), material properties (density 2700–3300 kg/m³, thermal conductivity 1–5 W/m/K, heat capacity 800–1200 J/kg/K), domain geometry (Nx/Nz resolution vs. domain size), boundary condition consistency, switch combinations that are incompatible, and common pitfalls (e.g., Courant number violations from dt vs. grid spacing, unrealistic initial temperatures). Provides a checklist-style validation workflow users can invoke when preparing or reviewing a `.txt` file.
- `skill-code-glossary`: Skill providing a comprehensive code-to-physics glossary — maps code variable names (e.g., `sxxd`, `eta_n`, `alp`, `bet`, `Eii`, `Tii`), struct fields, function names, and solver processes to their physical meaning, equations, and SI units. Organised as a quick-reference lookup covering: tensor components and invariants, material property symbols, grid array naming conventions (_n = cell centres, _s = vertices), solver process stages (what is physically happening at each step of the time loop), and common abbreviations (pwl = power law, lin = linear/diffusion, exp = exponential/Peierls, gbs = grain boundary sliding, cstv = constant viscosity). Enables users to read the C source and immediately understand the physical interpretation.

### Modified Capabilities

_(none — no existing specs to modify)_

## Impact

- **New files**: 10 new skill directories under `.github/skills/`, each containing a `SKILL.md`
- **No source code changes**: All additions are documentation/tooling only
- **No build changes**: Skills are metadata files, not compiled
- **Dependencies**: None — skills are standalone markdown consumed by Copilot
- **Users affected**: Anyone interacting with the codebase through Copilot or reading skill files directly
