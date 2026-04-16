## 1. Revert deactivation passes from modes 0 and 1

- [x] 1.1 Remove the deactivation block (lines ~2895–2911) from `CountPartCell` in `MDLIB/ParticleRoutines.c`
- [x] 1.2 Remove the deactivation block (lines ~3776–3851) from `CountPartCell_OLD` in `MDLIB/ParticleRoutines.c`
- [x] 1.3 Build and verify modes 0 and 1 still compile and run (existing GTests pass)

## 2. Create `CountPartCell_v2` function

- [x] 2.1 Copy `CountPartCell` (~lines 2658–2965) to a new function `CountPartCell_v2` with identical signature in `MDLIB/ParticleRoutines.c`
- [x] 2.2 Declare `CountPartCell_v2` in the appropriate header (`MDLIB/mdoodz-private.h` or similar)
- [x] 2.3 Replace cell-centroid placement with randomized placement: `xc - dx_fine/2 + dx_fine * rand() / RAND_MAX` for both x and z
- [x] 2.4 Replace single-particle reseeding with fill-to-target loop: add `threshold - current_count` particles per depleted cell
- [x] 2.5 Add distance-from-centroid deactivation: compute `(x-xc)² + (z-zc)²` for particles in over-populated cells, sort descending, deactivate beyond `min_part_cell + 4`

## 3. Wire `reseed_mode=2` dispatch

- [x] 3.1 Add `if ( input.model.reseed_mode == 2 ) CountPartCell_v2(...)` at dispatch site 1 in `Main_DOODZ.c` (~line 211, initial count)
- [x] 3.2 Add `if ( input.model.reseed_mode == 2 ) CountPartCell_v2(...)` at dispatch site 2 in `Main_DOODZ.c` (~line 1288, post-sedimentation)
- [x] 3.3 Add `if ( input.model.reseed_mode == 2 ) CountPartCell_v2(...)` at dispatch site 3 in `Main_DOODZ.c` (~lines 1332–1337, reseed + count)

## 4. Build and test

- [x] 4.1 Build with `make -j8` and fix any compilation errors
- [ ] 4.2 Run existing GTests (`ConvectionDevelops`) to verify no regression
- [x] 4.3 Update `TESTS/BlankenBench/BlankenBenchSteady.txt` to `reseed_mode=2`, `irestart=0`, `istep=50000`
- [ ] 4.4 Run Blankenbach benchmark (Run 6) with `reseed_mode=2` and compare against `run5_reference.csv`
- [ ] 4.5 Verify Nu and Vrms converge within 5% of published values (Nu=4.884, Vrms=42.865)

## 5. Rigid body rotation advection test

- [ ] 5.1 Create `TESTS/RotationAdvection/RotationAdvection.txt` parameter file (51×51, unit scaling, `mechanical=0`, `thermal=0`, `Courant=0.25`, `reseed_mode=0`)
- [ ] 5.2 Create `TESTS/RotationAdvection/RotationAdvectionMode1.txt` (copy of 5.1 with `reseed_mode=1`)
- [ ] 5.3 Create `TESTS/RotationAdvection/RotationAdvectionMode2.txt` (copy of 5.1 with `reseed_mode=2`)
- [ ] 5.4 Create `TESTS/RotationAdvectionTests.cpp` with GTest fixture: `SetPhase` (circular disk), `SetBCVx`/`SetBCVz` (rigid rotation `Vx=-ωz, Vz=ωx`)
- [ ] 5.5 Add L2 error computation: read initial and final HDF5 phase fields, compute relative L2 norm
- [ ] 5.6 Add test case `RotationAdvection.L2Comparison` that runs all three modes and asserts L2 < 0.5 for each, and L2(mode 2) ≤ 1.5 × L2(mode 1)
- [ ] 5.7 Add to `TESTS/CMakeLists.txt`: `file(COPY ...)`, `add_executable`, `target_link_libraries`, `add_test`
- [ ] 5.8 Build and run the rotation test, report L2 values for modes 0, 1, 2
