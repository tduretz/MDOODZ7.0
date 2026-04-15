## 1. Deactivation in Active Reseeding Paths

- [x] 1.1 Add deactivation pass to `CountPartCell` (mode 1) in `ParticleRoutines.c`: after the reseeding loop, scan all cells using `part_cell[][]` and mark particles beyond `min_part_cell + 4` as `phase = -1`. Skip cells with `BCt.type == 30`.
- [x] 1.2 Add deactivation pass to `CountPartCell_OLD` (mode 0) in `ParticleRoutines.c`: after the parallel reseeding section, add a serial loop using `ind_per_cell[][]` with the same `min_part_cell + 4` threshold and `BCt.type == 30` guard.
- [x] 1.3 Add LOG_INFO for deactivated particle count in both functions (matching `CountPartCell2` pattern).

## 2. Dead Code Cleanup

- [x] 2.1 Delete `CountPartCell2` function from `ParticleReseeding.c`. Removed file from CMakeLists.txt (entire file was just CountPartCell2).
- [x] 2.2 Delete unused `inc_reuse` variable and any `ind_part_reuse` / `ind_part_reuse2` allocations that are never consumed in `ParticleReseeding.c`. (Covered by 2.1 — whole file removed from build.)
- [x] 2.3 Delete `AddPartCell` and `AddPartVert` functions from `ParticleRoutines.c` (kept `AddPartCell2` which is actively used).
- [x] 2.4 Remove `CountPartCell2` declaration from `mdoodz-private.h`. (`AddPartCell`/`AddPartVert` had no header declarations.)

## 3. Nb_part in HDF5 Output

- [x] 3.1 Added `Nb_part` as integer dataset in "Model" group in `WriteOutputHDF5` in `HDF5Output.c`.
- [x] 3.2 Added `readScalarInt` helper in `TestHelpers.h` to read integer datasets from HDF5.

## 4. BlankenBench Test Updates

- [x] 4.1 Added `Nb_part` stability assertion to `ConvectionDevelops`: reads `Nb_part` from final output, asserts ≤ 2× initial count.
- [ ] 4.2 Tune `BlankenBenchSteady.txt` to reach thermal steady state in ~1000 steps (increase `Courant`, `dt_max`, adjust grid resolution).
- [ ] 4.3 Verify the `NusseltAndVrms` steady-state test runs to completion without `exit(190)` and matches Blankenbach values (Nu = 4.884, Vrms = 42.865) within 5%.

## 5. Skill: particle-reseeding

- [x] 5.1 Created `.github/skills/skill-particle-reseeding/SKILL.md`.
- [x] 5.2 Added `skill-particle-reseeding` entry to the skill table in `.github/copilot-instructions.md`.

## 6. Validation

- [x] 6.1 Full build passes (mdoodz + BlankenBenchTests, zero errors).
- [x] 6.2 `ConvectionDevelops` passes with Nb_part=50951 (ratio 1.99×, under 2× limit). Deactivation active (11 particles culled on final step).
