## 1. Add persistent thermal DirectSolver to Main_DOODZ.c

- [x] 1.1 Declare `DirectSolver ThermalSolver;` alongside the existing `CholmodSolver` in `Main_DOODZ.c`
- [x] 1.2 Initialize `cholmod_start(&ThermalSolver.c)` and set `ThermalSolver.Analyze = 1` before the timestep loop
- [x] 1.3 Set `ThermalSolver.c.nthreads_max` using the same `cholmod_nthreads` logic as the Stokes solver
- [x] 1.4 Add `cholmod_free_factor(&ThermalSolver.Lfact, &ThermalSolver.c)` and `cholmod_finish(&ThermalSolver.c)` after the timestep loop exits

## 2. Update EnergyDirectSolve signature and internals

- [x] 2.1 Change `EnergyDirectSolve` signature to accept `DirectSolver *thermal_solver` parameter
- [x] 2.2 Update the declaration in `mdoodz-private.h`
- [x] 2.3 Remove local `cholmod_common c; cholmod_start(&c);` and `cholmod_finish(&c)` — use `thermal_solver->c` instead
- [x] 2.4 Replace local `Afact` with `thermal_solver->Lfact`: on first call (`Analyze == 1`) store the factor from `FactorEnergyCHOLMOD`, on subsequent calls reuse it
- [x] 2.5 Remove `cholmod_free_factor(&Afact, &c)` at end of function — factor lifetime is now managed by the caller

## 3. Split FactorEnergyCHOLMOD into analyze+factorize paths

- [x] 3.1 Add `int analyze` parameter to `FactorEnergyCHOLMOD` (or accept the existing `cholmod_factor*` when `Analyze == 0`)
- [x] 3.2 When `analyze == 1`: run `cholmod_analyze` then `cholmod_factorize`, return new factor (current behavior)
- [x] 3.3 When `analyze == 0`: accept existing `cholmod_factor*`, run only `cholmod_factorize` on it, skip `cholmod_analyze`
- [x] 3.4 Update the declaration in `mdoodz-private.h` if needed

## 4. Update ThermalSteps signature and call site

- [x] 4.1 Change `ThermalSteps` signature to accept `DirectSolver *thermal_solver`
- [x] 4.2 Update the declaration in `mdoodz-private.h`
- [x] 4.3 Pass `thermal_solver` through to its internal `EnergyDirectSolve` call
- [x] 4.4 Update the `ThermalSteps` call site in `Main_DOODZ.c` (line ~254, initial_cooling) to pass `&ThermalSolver`

## 5. Update call sites in Main_DOODZ.c

- [x] 5.1 Update the `EnergyDirectSolve` call in the timestep loop (line ~1037) to pass `&ThermalSolver`
- [x] 5.2 After the first successful thermal solve, set `ThermalSolver.Analyze = 0`

## 6. Verify and benchmark

- [x] 6.1 Build and run locally at low resolution (201×161) — confirm no crashes, no warnings, correct temperature output
- [x] 6.2 Compare output against a baseline run to verify bit-for-bit identical results
- [x] 6.3 Run EC2 benchmark at 1001×801 with threads 1/2/4/6/8/12/16 and compare `thermal_s` against baseline
