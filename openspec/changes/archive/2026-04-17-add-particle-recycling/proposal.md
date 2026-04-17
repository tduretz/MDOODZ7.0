## Why

In closed-box simulations (e.g. BlankenBench convection), particles never leave the domain so no dead particles (`phase == -1`) accumulate. The active reseeding paths (`CountPartCell` mode 1, `CountPartCell_OLD` mode 0) only recycle dead particles — they never **deactivate** excess particles in over-populated cells. This means `Nb_part` grows monotonically until it hits `Nb_part_max` (4.1× initial), at which point the code calls `exit(190)` and crashes. The only deactivation logic exists in `CountPartCell2`, which is dead code never called from `Main_DOODZ.c`.

## What Changes

- **Add excess-particle deactivation to active reseeding paths**: Port the deactivation logic from dead `CountPartCell2` (lines ~650–665 of `ParticleReseeding.c`) into both `CountPartCell` (mode 1) and `CountPartCell_OLD` (mode 0) in `ParticleRoutines.c`. When a cell has more particles than `max_part_cell`, mark the extras as `phase = -1` so they become available for recycling.
- **Clean up dead code in `ParticleReseeding.c`**: `CountPartCell2` is never called. `ind_part_reuse` / `ind_part_reuse2` are allocated and freed without being consumed. `inc_reuse` is declared but never incremented. `AddPartCell` and `AddPartVert` (with recycling) are also never called. Remove or clearly mark this dead code.
- **Validate particle count stability in BlankenBench**: BlankenBench requires ~100k steps to reach thermal steady state, which is too expensive for CI. The existing CI smoke test (500 steps) may not be enough to trigger the particle overflow — we may need more steps, or can accelerate convergence by tuning settings (larger `dt`, higher `Courant`, coarser grid). Once particle recycling is working, the manual steady-state test (`NusseltAndVrms`, gated by `BLANKENBACH_STEADY=1`) should run to completion without crashing, and its results should match the published Blankenbach et al. (1989) values (Nu = 4.884, Vrms = 42.865). Add a particle count stability check to the CI test to confirm `Nb_part` is not growing monotonically — adjusting step count or time stepping as needed to reliably trigger the reseeding path.
- **Create `skill-particle-reseeding`**: A Copilot skill documenting how particles are managed — reseeding modes, recycling, deactivation, the markers struct, `Nb_part` / `Nb_part_max`, and common failure modes.

## Capabilities

### New Capabilities
- `particle-excess-deactivation`: Logic to mark excess particles in over-populated cells as `phase = -1`, enabling recycling in closed-box simulations
- `skill-particle-reseeding`: Copilot skill documenting the particle reseeding system — modes 0/1, recycling logic, deactivation, markers struct, failure modes

### Modified Capabilities
- `ci-blankenbach-convection`: Add `Nb_part` stability assertion to the 500-step CI smoke test. The full steady-state test (100k steps, manual only) should reach completion and match Blankenbach et al. (1989) published values once recycling is fixed.

## Impact

- `MDLIB/ParticleRoutines.c` — add deactivation of excess particles in `CountPartCell` and `CountPartCell_OLD`
- `MDLIB/ParticleReseeding.c` — remove or annotate dead code (`CountPartCell2`, `AddPartCell`, `AddPartVert`)
- `TESTS/BlankenBenchTests.cpp` — add `Nb_part` stability check to 500-step CI test; full steady-state validation (Nu, Vrms vs Blankenbach 1989) remains manual-only
- `.github/skills/skill-particle-reseeding/SKILL.md` — new skill file
- No breaking changes. Existing simulations with outflow (rifting, subduction) are unaffected since they already produce dead particles naturally.
