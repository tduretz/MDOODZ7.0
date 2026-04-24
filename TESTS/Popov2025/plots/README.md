# Popov et al. (2025) — paper-figure reproduction

A gnuplot script that renders a qualitative counterpart of Fig. 5 of

> Popov, A. A., Berlie, N., Kaus, B. J. P. (2025). *A dilatant visco-elasto-
> viscoplasticity model with globally continuous tensile cap: stable
> two-field mixed formulation.* Geosci. Model Dev., 18, 7035–7058.
> doi:[10.5194/gmd-18-7035-2025](https://doi.org/10.5194/gmd-18-7035-2025)

The script reads the HDF5 output of the 0D integration tests in
`TESTS/Popov2025Tests.cpp` via the `extract_popov2025` helper. It is a
**developer aid**: NOT invoked by CI and NOT wired into the
`VISUAL_TESTS/` regression pipeline.

## Invocation order

1. **Build** the test binary and the extractor:
   ```bash
   cmake -B ./cmake-build -DTEST=ON
   cmake --build ./cmake-build
   ```

2. **Run** the 0D integration tests:
   ```bash
   make -C cmake-build run-tests          # includes Popov2025Tests
   # or directly:
   (cd cmake-build/TESTS && ./Popov2025Tests)
   ```

3. **Extract** ASCII `.dat` files from the HDF5 output directories:
   ```bash
   cd cmake-build/TESTS
   ./extract_popov2025 Popov0D_VolumetricExtension
   ./extract_popov2025 Popov0D_DeviatoricShear
   ./extract_popov2025 Popov0D_MixedStrain
   ```
   The extractor emits a `trajectory.dat` with columns
   `(step, t[s], sII[Pa], P[Pa], eII_pl[1/s], divu_pl[1/s])` for the
   centre cell of each run.

4. **Plot** with gnuplot from `TESTS/Popov2025/plots/`:
   ```bash
   cd TESTS/Popov2025/plots
   gnuplot fig5_0d_stress.gp
   ```
   The PNG lands in `out/fig5_0d_stress.png` (gitignored).  A committed
   reference copy lives at `reference/fig5_0d_stress.png` and is used by
   `TESTS/AnalyticalSolutions.md` §6.1b.

## Overrides

`fig5_0d_stress.gp` defaults `BUILD_DIR` to `../../../cmake-build/TESTS`.
Override on the command line if your build tree lives elsewhere:

```bash
gnuplot -e 'BUILD_DIR="/path/to/build/TESTS"' fig5_0d_stress.gp
```

## What the plot shows

The four panels correspond to paper Fig. 5 a/b/c/d:

- **(a) Volumetric extension** — `p → p_T = −0.5 MPa`, `τII ≈ 0` throughout.
- **(b) Deviatoric shear** — `τII` yields on the Drucker-Prager envelope
  (≈ 1 MPa) at ≈ 20 yr; `p` grows post-yield from dilatant ψ = 10° coupling.
- **(c) Mixed strain** — `p` dips to ≈ −0.35 MPa while `τII` grows through
  yield; the trajectory sweeps through the cap ↔ DP delimiter.
- **(d) Meridional `P–τII` trajectories** — three characteristic shapes:
  extension slides onto the cap tip, shear climbs the `τII` axis and bends
  right, mixed forms the classic S-curve.

The pressure plotted is the integration-point value
`p_local = p* + K·θ̇_vp·dt` (paper Eq. 31), reconstructed in gnuplot from
the HDF5 `P` (global trial) and `divu_pl` fields. This is why the Fig. 5a
saturation sits exactly at `p_T`, not at the one-step-lagged trial value.
