## OUTCOME: DESIGN INVALIDATED — ALL TASKS CANCELLED

Implementation was attempted (April 2026). Tasks 1.1–6.1 were completed,
but at task 6.3 the `cell_avg=1` mode caused CHOLMOD failure:
"matrix not positive definite". All code changes were reverted.

**Root cause**: Setting `D11E = D11W = harmonic_mean` makes the operator
locally constant across each pair of cells, destroying the viscosity
contrast and breaking SPD structure. In this staggered grid the D values
at cell centres already represent the natural face values for the FD
differencing scheme — there is no separate "face viscosity" to average.

**Actions taken after failure:**
- [x] Reverted all code changes (mdoodz.h, InputOutput.c,
      StokesAssemblyDecoupled.c, SolViMarkerComparison.cpp)
- [x] Deleted test parameter files (SolViRes*_cellavg1.txt)
- [x] Documented negative result in `TESTS/SolViMarkerComparison.cpp` header
- [x] Updated `skill-solvers/SKILL.md` with cell-face averaging failure
- [x] Updated `skill-testing-guide/SKILL.md` with "ruled out" summary

**Conclusion**: All three approaches to improving P convergence order have
been tested and ruled out: (1) more markers, (2) eta_average, (3) D-tensor
harmonic averaging. Improving beyond ~0.75 would require fundamentally
different stencil methods (immersed-boundary, ghost-fluid).

---

## Original Tasks (all cancelled)

## 1. Parameter Infrastructure

- [ ] ~~1.1 Add `int cell_avg` field to the `params` struct in `MDLIB/mdoodz-private.h`~~
- [ ] ~~1.2 Add `ReadInt2(fin, "cell_avg", 0)` in `MDLIB/InputOutput.c` alongside `eta_average`~~
- [ ] ~~1.3 Verify the parameter is printed during startup (add to parameter summary printf block)~~

## 2. X-Momentum Stencil Modification

- [ ] ~~2.1 In `Xmomentum_InnerNodesDecoupled` (`MDLIB/StokesAssemblyDecoupled.c`), add harmonic averaging of `D11E`/`D11W` gated on `model.cell_avg == 1`, with div-by-zero guard (sum < 1e-30 → 0.0)~~
- [ ] ~~2.2 In the same function, add harmonic averaging of `D33N`/`D33S` gated on `model.cell_avg == 1`~~

## 3. Z-Momentum Stencil Modification

- [ ] ~~3.1 In `Zmomentum_InnerNodesDecoupled` (`MDLIB/StokesAssemblyDecoupled.c`), add harmonic averaging of `D22S`/`D22N` gated on `model.cell_avg == 1`, with div-by-zero guard~~
- [ ] ~~3.2 In the same function, add harmonic averaging of `D33E`/`D33W` gated on `model.cell_avg == 1`~~

## 4. Test Parameter Files

- [ ] ~~4.1 Create `TESTS/SolViBenchmark/SolViRes41_cellavg1.txt` (copy of SolViRes41.txt with `cell_avg = 1`)~~
- [ ] ~~4.2 Create `TESTS/SolViBenchmark/SolViRes51_cellavg1.txt` (copy of SolViRes51.txt with `cell_avg = 1`)~~
- [ ] ~~4.3 Create `TESTS/SolViBenchmark/SolViRes81_cellavg1.txt` (copy of SolViRes81.txt with `cell_avg = 1`)~~

## 5. L2 Comparison Test

- [ ] ~~5.1 Add `CellFaceAveraging` test in `TESTS/SolViMarkerComparison.cpp`~~
- [ ] ~~5.2 Assert P convergence order ≥ 0.8 with `cell_avg = 1`~~
- [ ] ~~5.3 Assert Vx L2 not more than 20% worse than `cell_avg = 0`~~

## 6. Build & Verify

- [ ] ~~6.1 Build with `cmake -DTEST=ON` and verify no compilation errors~~
- [ ] ~~6.2 Run existing SolVi tests (`SolViBenchmarkTests`) to confirm `cell_avg = 0` (default) produces identical results (no regression)~~
- [ ] ~~6.3 Run the new `CellFaceAveraging` test and record measured L2 values and convergence orders~~

## 7. Documentation

- [ ] ~~7.1 Update `skill-solvers/SKILL.md`~~
- [ ] ~~7.2 Update `skill-testing-guide/SKILL.md`~~
- [ ] ~~7.3 Update `skill-model-setup/SKILL.md`~~
- [ ] ~~7.4 Update `skill-code-glossary/SKILL.md`~~
- [ ] ~~7.5 Add findings to `TESTS/SolViMarkerComparison.cpp` header comment block~~
