## Context

MDOODZ uses marker-in-cell (MIC) advection with periodic reseeding to maintain particle coverage. The reseeding logic lives in `MDLIB/ParticleRoutines.c` in two functions: `CountPartCell` (mode 1, lines 2658–2965) and `CountPartCell_OLD` (mode 0, lines 2970–3988). Both operate on a fine mesh (2× resolution) and share three structural problems: index-order deactivation, single-particle reseeding, and cell-centroid placement.

Rather than modifying these legacy functions, we will create a new `CountPartCell_v2` function activated by `reseed_mode=2`. This preserves backward compatibility — modes 0 and 1 remain untouched. The deactivation passes previously added to modes 0 and 1 will be reverted to restore their original behaviour.

Property assignment (`AssignMarkerProperties`) already supports grid interpolation via `Centers2Particle` when `direct_neighbour=0` (the default), so T/P/stress fields are correctly set for new particles. The noise is from the deactivation and placement, not from property interpolation.

Blankenbach Run 5 (101×101, 25k steps, `reseed_mode=1`) shows sustained Nu_top oscillations (1.9–4.6) that never damp. The particle count saturates at 459k, meaning reseeding and deactivation are in continuous equilibrium — a steady noise pump.

## Goals / Non-Goals

**Goals:**
- Reduce reseeding-induced noise to allow the Blankenbach benchmark to converge within 5% of published values (Nu=4.884, Vrms=42.865)
- Implement improved reseeding as a new `reseed_mode=2` (`CountPartCell_v2`) to allow A/B comparison
- Revert deactivation passes in modes 0 and 1 to restore backward compatibility
- Deactivate spatially redundant (farthest-from-centroid) particles instead of index-ordered ones
- Fill depleted cells to target count in a single pass
- Randomize new particle positions within cells

**Non-Goals:**
- Modifying `CountPartCell` (mode 1) or `CountPartCell_OLD` (mode 0) beyond reverting the deactivation pass
- Changing the deactivation threshold (`min_part_cell + 4`)
- Modifying `AssignMarkerProperties` or property interpolation
- Thread-safety changes

## Decisions

### 1. Deactivation ordering: distance-from-centroid (farthest first)

**Choice:** Sort particles in each over-populated cell by distance from the fine-cell centroid; deactivate the farthest ones beyond the threshold.

**Rationale:** Particles near the cell center contribute most to cell-averaged quantities during particle-to-grid (P2G) interpolation. Particles near cell edges are partially redundant with neighbours in adjacent cells. Farthest-first removal preserves the most representative particles.

**Alternatives considered:**
- *Random selection*: Simpler but still destroys well-positioned particles. Distance-based is deterministic and spatially motivated.
- *Oldest-first (by creation time)*: Would require tracking particle age. More complex with unclear benefit.
- *Lowest-weight (by interpolation weight)*: Requires computing bilinear weights. Expensive and tightly couples reseeding to the interpolation scheme.

**Implementation:** For each over-populated cell, compute `(x - xc)² + (z - zc)²` for each particle in the cell, sort descending, and deactivate entries beyond `min_part_cell + 4`. Use a simple insertion sort — cell populations are small (typically 16–24), so O(n²) is negligible.

### 2. Multi-particle fill-to-target

**Choice:** When a cell is below the reseeding threshold, add `(threshold - current_count)` particles in one pass instead of 1.

**Rationale:** With threshold=2 and target recovery to ~4 particles, at most 2–3 particles are added per cell per step. This eliminates the multi-step lag where a depleted cell runs with 1 particle for several steps, degrading interpolation accuracy in that cell each step.

**Alternatives considered:**
- *Fill to `min_part_cell`*: Would add up to 16 particles at once, which is wasteful — most cells only need 1–3 to get above threshold. Stick with filling to the threshold value.
- *Gradual increase (e.g., double each step)*: Adds complexity without clear benefit over direct fill.

**Implementation:** Replace the single `AddPartCell` / recycling block with a loop over `needed = threshold - current_count`. Each iteration reuses the existing nearest-phase search and `AssignMarkerProperties` call, but with a different random position.

### 3. Randomized placement within cell

**Choice:** Place new particles at random positions within the fine cell bounds using `rand() / RAND_MAX` jitter, matching the existing `PutPartInBox` pattern (line 1133).

**Rationale:** Cell-centroid placement creates spatial clustering. Random placement gives better coverage, especially when multiple particles are added in a single pass (Decision 2). The existing codebase already uses `rand()` for initial particle placement — no new RNG infrastructure needed.

**Implementation:**
```
new_x = xc - dx_fine/2 + (dx_fine * rand() / RAND_MAX)
new_z = zc - dz_fine/2 + (dz_fine * rand() / RAND_MAX)
```
Where `dx_fine = mesh->dx/2`, `dz_fine = mesh->dz/2` (fine cell dimensions). This places particles uniformly within the fine cell bounds.

### 4. New function `CountPartCell_v2` activated by `reseed_mode=2`

**Choice:** Create a new standalone function `CountPartCell_v2` with the same signature as `CountPartCell`, based on a copy of its structure (~300 lines). Wire it via `reseed_mode == 2` at all three dispatch sites in `Main_DOODZ.c`. Revert the deactivation blocks previously added to `CountPartCell` (mode 1) and `CountPartCell_OLD` (mode 0).

**Rationale:** A new mode preserves backward compatibility — existing simulations using modes 0 or 1 are unaffected. It enables direct A/B testing (same setup, just change `reseed_mode` from 1 to 2). The deactivation passes we added to modes 0/1 were a stopgap; removing them restores those functions to their original tested state.

**Alternatives considered:**
- *Modify mode 1 in place*: Simpler but breaks backward compatibility and makes it impossible to A/B compare.
- *Runtime flag within mode 1*: Adds conditionals to an already complex function; harder to reason about.

**Implementation:**
- Copy `CountPartCell` (~lines 2658–2965) to a new function `CountPartCell_v2` with identical signature
- Apply all three improvements (distance deactivation, fill-to-target, random placement) in the new function
- Add `if ( input.model.reseed_mode == 2 ) CountPartCell_v2(...)` at lines ~211, ~1288, ~1332–1337 in `Main_DOODZ.c`
- Remove the deactivation block (lines 2895–2911) from `CountPartCell`
- Remove the deactivation block (lines 3776–3851) from `CountPartCell_OLD`

## Risks / Trade-offs

**[Sorting overhead in deactivation]** → Mitigated by small cell populations (16–24 particles). Insertion sort on 24 elements is ~300 comparisons — negligible compared to the Stokes solve. No allocation needed if we sort in-place on the `part_cell[kc][]` array.

**[Random placement may put particle near cell boundary]** → Acceptable. Particles near boundaries are valid — they'll be counted in the correct cell and contribute to interpolation normally. The deactivation ordering (farthest-first) will preferentially remove boundary particles in over-populated cells anyway.

**[`rand()` is not thread-safe]** → Not a concern. The reseeding loop in `CountPartCell_v2` is serial (no OpenMP pragma). If parallelism is added later, switch to `rand_r()` with per-thread seeds.

**[`rand()` has no explicit seed]** → Acceptable for reseeding. Reproducibility is not required for particle placement noise — only the statistical properties matter. The existing `PutPartInBox` also uses unseeded `rand()`.

**[Fill-to-target may add particles that are immediately deactivated]** → Not possible. Reseeding triggers at `< 2`, filling to 2. Deactivation triggers at `> min_part_cell + 4 = 20`. The gap (2→20) prevents immediate removal.

**[Code duplication from copying CountPartCell]** → Acceptable tradeoff for backward compatibility. The new function will diverge intentionally. If mode 2 proves superior long-term, modes 0/1 can be deprecated in a separate change.

### 5. Rigid body rotation test for quantitative reseeding comparison

**Choice:** Create a new GTest (`RotationAdvectionTests.cpp`) that prescribes rigid body rotation `Vx = -ω·z, Vz = ω·x`, advects a circular compositional disk for one full revolution, and compares the L2 error of the cell-centre phase field against the initial condition. Run the same scenario three times with `reseed_mode` 0, 1, and 2.

**Rationale:** The Blankenbach benchmark tests the full coupled system — reseeding noise is convolved with thermal diffusion and the Stokes solve. A prescribed-velocity advection test isolates the marker advection + reseeding loop, and the analytical solution is trivially the initial condition after one revolution. Any L2 error is directly attributable to the particle method (advection discretization + reseeding artifacts). This is a standard benchmark in computational geodynamics (Gerya 2010, Thielmann et al. 2014).

**Design:**
- Domain: [-0.5, 0.5]² with unit scaling (`eta=1, L=1, V=1, T=1`)
- Grid: 51×51, 4×4 markers/cell
- Anomaly: circular disk of phase 1 at offset (x0, z0) with radius R, surrounded by phase 0
- Velocity: `Vx = -ω·z`, `Vz = ω·x` via `SetBCVx`/`SetBCVz` callbacks returning the prescribed field (type=0 Dirichlet everywhere)
- `mechanical=0` to skip the Stokes solve — pure advection test
- Total time: `T = 2π/ω` for one revolution; ω chosen so Courant ≤ 0.25
- L2 metric: relative L2 of cell-centre phase proportion field, comparing final HDF5 output to initial
- Three runs per test: the test function loops over `reseed_mode` 0, 1, 2, modifying the `.txt` file or using three separate `.txt` files

**Alternatives considered:**
- *Blankenbach-only validation*: The Blankenbach benchmark already tests reseeding indirectly via Nu/Vrms, but it's expensive (hours) and the noise signal is mixed with solver/thermal effects. The rotation test provides a clean, fast (minutes) diagnostic.
- *Zalesak slotted disk*: Sharper interface test but requires discontinuous composition tracking that MDOODZ doesn't support well. A smooth circular disk is sufficient to measure reseeding diffusion.
- *van Keken RT benchmark*: Tests advection + Stokes coupling. Good complement but doesn't isolate reseeding quality. Can be added later.

**Implementation notes:**
- MDOODZ needs `mechanical=0` to skip the Stokes solve. If this flag doesn't exist, use a constant-viscosity setup with `pure_shear_ALE=0`, `bkg_strain_rate=0`, and set velocity via BC callbacks. The velocity field will be "incorrect" from the Stokes perspective but the markers will still be advected with the prescribed field.
- The `SetBCVx`/`SetBCVz` callbacks must return `type=0` (Dirichlet) at all nodes with the analytical rotation velocity. This overwrites the Stokes solution at every step.
- Phase proportions on cell centres are computed by `P2Mastah` and stored in HDF5. The L2 comparison reads these from the initial and final output files.
