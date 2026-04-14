## Why

P2Mastah (particle-to-grid interpolation) accounts for 29% of wall time at 16 threads (2.56s/step) and is bandwidth-saturated — flat from 8 to 16 threads. The root cause is redundant memory traffic: ~35 separate P2Mastah calls per timestep each loop over all ~12M particles, recomputing identical grid weights and re-reading the same particle coordinates. Fusing calls that share a centroid type into a single pass would reduce total particle memory traffic by ~8×, breaking through the bandwidth plateau.

## What Changes

- New `P2Mastah_Fused()` function that scatters multiple fields per particle in a single pass over the particle array, computing grid weights once and writing to N destination grids per particle visit
- New `interp_mode = 3` value: fused atomic scatter (extends the existing 0/1/2 mode selection)
- Call-site grouping: consecutive P2Mastah calls sharing a centroid type are replaced by a single `P2Mastah_Fused()` call with a field descriptor array
- `interp_mode = 0/1/2` remain fully functional and unchanged — mode 3 is opt-in
- Existing `P2Mastah()` is not modified — `P2Mastah_Fused()` is a new function alongside it

## Capabilities

### New Capabilities
- `fused-p2g-interp`: Fused multi-field particle-to-grid interpolation — single-pass scatter per centroid type, field descriptor API, backward-compatible mode selection, EC2 benchmark validation

### Modified Capabilities
- `interp-buffer-pool`: Adding `interp_mode = 3` to the mode selection; pool must pre-allocate destination buffers for fused scatter batches

## Impact

- **Code**: `MDLIB/ParticleRoutines.c` (new `P2Mastah_Fused` function), `MDLIB/Main_DOODZ.c` (call-site refactoring — ~35 individual calls grouped into ~6 fused calls), `MDLIB/include/mdoodz.h` (new structs/prototypes)
- **Parameters**: `interp_mode = 3` added to `.txt` config (default remains 0)
- **CI**: Existing interp mode 1/2 tests unaffected; new mode 3 test fixtures required for all 19 base scenarios
- **Performance target**: `interp_s` reduction from ~2.56s to <1.5s at 16 threads on c5ad.4xlarge (1001×801), validated via EC2 highres benchmark sweep
- **Memory**: Mode 3 should use comparable or less memory than mode 2 (no additional thread-local buffers)
- **Numerical**: Mode 3 atomic scatter produces machine-epsilon differences vs mode 0 (same as mode 2), not bit-identical
