## Why

The `CountPartCell` reseeding logic (mode 1) introduces sustained numerical noise that prevents the Blankenbach benchmark from converging to steady state. Run 5 (101×101, 25k steps) shows Nu_top oscillating between 1.9–4.6 with no damping, while the target is 4.88. The reseeding algorithm has three structural problems: (1) deactivation selects particles by index order, biased toward removing the most recently advected particles that carry accurate thermal information; (2) only one particle is added per depleted cell per step, making recovery from depletion slow; (3) particles are always placed at cell centers, creating systematic spatial concentration bias. Reference data from Run 5 is saved in `TESTS/BlankenBench/run5_reference.csv` (257 outputs).

Note: `AssignMarkerProperties` already supports grid interpolation for T, P, stress, etc. via `Centers2Particle` when `direct_neighbour=0` (the default). The property assignment itself is not the problem.

## What Changes

- **New `reseed_mode=2`**: Add a new `CountPartCell_v2` function with improved reseeding, selectable via `reseed_mode=2` in the `.txt` parameter file. Modes 0 and 1 remain untouched for backward compatibility.
- **Deactivation selection** (mode 2): Distance-from-centroid ordering — deactivate particles farthest from the cell center first, preserving spatially central particles that best represent the cell
- **Multi-particle reseeding** (mode 2): When a cell is below the reseeding threshold, add enough particles to reach the threshold in a single pass instead of adding only one per step
- **Randomized placement** (mode 2): Place new particles at random positions within the cell instead of always at the cell centroid, improving spatial coverage
- **Revert deactivation passes**: Remove the deactivation blocks previously added to `CountPartCell` (mode 1) and `CountPartCell_OLD` (mode 0) to restore them to their original behaviour

## Capabilities

### New Capabilities
- `improved-reseeding`: New `reseed_mode=2` with distance-from-centroid deactivation, multi-particle fill-to-target, and randomized placement within cells
- `rotation-advection-test`: Rigid body rotation benchmark — prescribe solid-body rotation velocity field, advect a compositional anomaly for one full revolution, measure L2 error of composition field vs initial condition for reseed_mode 0, 1, and 2

### Modified Capabilities

(none — no existing spec-level behavior changes; modes 0 and 1 restored to original)

## Impact

- **Code**: `MDLIB/ParticleRoutines.c` — new `CountPartCell_v2` function; revert deactivation blocks in `CountPartCell` and `CountPartCell_OLD`
- **Code**: `MDLIB/Main_DOODZ.c` — add `reseed_mode == 2` dispatch at three call sites
- **Code**: `MDLIB/InputOutput.c` — document `reseed_mode=2` option (no parse changes needed, already reads as int)
- **Testing**: Blankenbach steady-state test (`TESTS/BlankenBenchTests.cpp`) is the primary validation — Run 6 with `reseed_mode=2` should show reduced oscillations compared to `run5_reference.csv`
- **Testing**: New rigid body rotation GTest (`TESTS/RotationAdvectionTests.cpp`) — quantitative L2 comparison of composition field after one revolution across all three reseed_modes
- **Risk**: Low — modes 0 and 1 are restored to original, mode 2 is a new code path with no side effects on existing simulations
- **New parameter value**: `reseed_mode=2` in `.txt` files activates the improved reseeding. Default remains 0.
