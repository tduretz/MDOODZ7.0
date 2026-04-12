## 1. Create the C scenario file

- [x] 1.1 Copy `SETS/MLPS_Ellipses.c` to `SETS/RiftingCombinedYield.c`
- [x] 1.2 Update the `RunMDOODZ` call to reference `"RiftingCombinedYield.txt"`
- [x] 1.3 Verify all callbacks (`SetPhase`, `SetTemperature`, `SetDensity`, `SetBCT`) and ellipse geometry are preserved unchanged

## 2. Create the production .txt parameter file

- [x] 2.1 Copy `SETS/MLPS_Ellipses.txt` to `SETS/RiftingCombinedYield.txt`
- [x] 2.2 Add `plast = 2`, `T_st = -10e6`, and `eta_vp = 1.0e19` to all 9 phase blocks (IDs 0–8)
- [x] 2.3 Add melting parameters: `melting = 1`, `melt_weak = 30.0`, `force_melt_weak = 1` in the global switches section
- [x] 2.4 Add per-phase `melt` values: `melt = 11` for crustal phases (0, 4, 6, 7, 8) and `melt = 41` for mantle phases (1, 2, 3, 5)
- [x] 2.5 Update solver settings: `Newton = 0`, `Picard2Newton = 1`, `Picard2Newton_tol = 1e-1`, `nit_max = 20`, `line_search = 1`
- [x] 2.6 Confirm domain, resolution (`Nx=200`, `Nz=160`), elastic/thermal/shear_heating switches are correct
- [x] 2.7 Confirm no anisotropy parameters are present

## 3. Create verification .txt files

- [x] 3.1 Copy `RiftingCombinedYield.txt` to `SETS/RiftingCombinedYield_quick.txt`, set `Nx=101`, `Nz=81`, `Nt=10`
- [x] 3.2 Copy `RiftingCombinedYield.txt` to `SETS/RiftingCombinedYield_medium.txt`, set `Nx=101`, `Nz=81`, `Nt=100`

## 4. Register in CMake

- [x] 4.1 Add `add_set(RiftingCombinedYield)` to `SETS/CMakeLists.txt` (or wherever sets are registered)

## 5. Build and quick verification

- [x] 5.1 Run `make build` and confirm `RiftingCombinedYield` compiles without errors
- [x] 5.2 Run with `RiftingCombinedYield_quick.txt` (10 steps) — verify no NaN or divergence, HDF5 output produced
- [x] 5.3 Run with `RiftingCombinedYield_medium.txt` (100 steps) — verify stable over longer evolution, no thermal runaway or melt instabilities
