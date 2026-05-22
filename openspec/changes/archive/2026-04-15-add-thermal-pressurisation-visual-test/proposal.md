## Why

Visual verification plots for the thermal validation tests (ThermoElastic / thermal pressurisation, OceanicCooling) should be auto-generated and embedded in the VISUAL_TESTS README. The plots already exist as PNGs in `VISUAL_TESTS/img/` but are not referenced in the README, so there is no visual audit trail. Additionally, the ThermoElastic plot currently shows P(t) and T(t) separately — it should be a **P vs T** plot with an error panel (expected ~0.37% error over 10 steps), matching the Julia visualisation in `Main_Visualisation_Makie_MD7.jl` (field `:PTtime`). The analytical relationship is $P = P_0 + \frac{\alpha}{\beta}(T - T_0)$ with $\alpha = 10^{-5}$ K$^{-1}$, $\beta = 10^{-10}$ Pa$^{-1}$.

## What Changes

- **Improve ThermoElastic gnuplot**: Replace the current T(t) + P(t) two-panel layout with a P-vs-T plot plus a pressure error panel (matching the reference screenshot and the Julia `:PTtime` field). Update `Validation.cpp` to write the data columns needed for P-vs-T plotting.
- **Increase ThermoElastic steps to 10**: The CI test currently runs 5 steps. The expected error (0.37%) is calibrated for 10 steps. Add a 10-step visual-test variant to show the error behaviour over more steps.
- **Add CI test comparing thermal_solver 0 vs 1**: The existing `ThermoElastic` and `ThermoElasticPCG` CI tests run independently. Add an explicit CI assertion that CHOLMOD and PCG produce equivalent results for thermal pressurisation (L2 error between solvers < tolerance), similar to the existing cross-comparison pattern in OceanicCooling tests.
- **Embed thermal plots in README**: Add sections to `VISUAL_TESTS/readme.md` for ThermoElastic (thermal pressurisation), OceanicCooling, and ThermalAccuracy — all three already generate PNGs but none are in the README.

## Capabilities

### New Capabilities
- `thermal-pressurisation-plot`: P-vs-T gnuplot with error panel for the ThermoElastic benchmark, auto-generated during `make run-vis`, embedded in VISUAL_TESTS README
- `thermal-solver-comparison-ci`: CI test asserting CHOLMOD and PCG thermal solvers produce equivalent ThermoElastic results

### Modified Capabilities

## Impact

- `VISUAL_TESTS/Validation.cpp` — modify `PlotThermoElastic()` data output format
- `VISUAL_TESTS/ThermoElastic.gnu` — replace with P-vs-T + error layout
- `VISUAL_TESTS/readme.md` — add three new sections with embedded images
- `TESTS/ThermalVerificationTests.cpp` — potentially add solver cross-comparison assertion
- `TESTS/Thermal/ThermoElastic.txt` — may need 10-step variant for visual test
- No breaking changes. Existing CI tests unaffected.
