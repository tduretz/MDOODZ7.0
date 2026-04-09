# Proposal: blankenbach-steady-state

## Problem

The existing BlankenBench test (Suite 18) runs only 500 steps (~78 s) as a CI smoke test. It verifies that convection starts but cannot validate the code against the published Blankenbach et al. (1989) Case 1a reference values because the simulation is far from steady state (Nu = 0.32 at step 500 vs. published 4.884).

There is no quantitative benchmark that tests the full thermo-mechanical coupling loop against community-agreed reference values.

## Solution

Add a **long-running** Blankenbach steady-state test that:

1. Uses a separate `.txt` parameter file with high step count (Nt = 100,000) and writes only the final output (`writer_step = 100000`)
2. A C++ GTest reads the final HDF5 output and computes:
   - **Nusselt number** (Nu): surface heat flux integral
   - **Vrms**: root-mean-square velocity over the domain
3. Compares against published reference values with resolution-appropriate tolerances (~5% at 41×41)

This test is **not** added to CI/ctest — it's built alongside the other tests but excluded from `add_test()`. The user runs it manually, leaves it overnight, and checks the result.

## Capabilities

1. **steady-state-parameter-file**: A `BlankenBench/BlankenBenchSteady.txt` with Nt=100000, writer_step=100000
2. **nusselt-vrms-test**: A new GTest (`BlankenBench.NusseltNumber`, `BlankenBench.VrmsMatch`) that computes Nu and Vrms from the final HDF5 output and asserts against the Blankenbach reference values

## Scope

**In scope:**
- New parameter file with high Nt for steady-state run
- New GTest cases computing Nu, Vrms from the final output
- Tolerances calibrated for 41×41 resolution (~5%)
- Documentation in README.md Suite 18 section

**Out of scope:**
- Grid convergence study (running at multiple resolutions)
- Corner heat flux comparison (q1, q2)
- Any modification to the MDLIB C code
- Adding this to CI (too slow)

## Non-goals

- We are NOT modifying the existing 500-step smoke test
- We are NOT adding early-termination / convergence-detection logic to MDOODZ
- We are NOT running higher-Ra cases (Case 1b, 1c, 2a)

## Impact

- Adds 2 new test cases (not in CI) to the BlankenBench suite
- Adds 1 new parameter file
- Updates README.md documentation
- No changes to production code (MDLIB)
