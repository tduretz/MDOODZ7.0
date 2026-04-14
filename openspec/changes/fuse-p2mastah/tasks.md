## 1. P2MastahField Struct and Declarations

- [x] 1.1 Define `P2MastahField` struct in `mdoodz-private.h` with fields: `double *src`, `double *dst`, `char *BCtype`, `int flag`, `int avg`, `int prop`, `int stencil` (per design D3)
- [x] 1.2 Declare `P2Mastah_Fused` function prototype in `mdoodz-private.h`: `void P2Mastah_Fused(params *, markers, grid *, int centroid, int interp_stencil, P2MastahField *, int nfields, InterpBufPool *)`

## 2. Extend InterpBufPool for Mode 3

- [x] 2.1 Add mode 3 fields to `InterpBufPool` struct in `mdoodz-private.h`: `int max_fused_fields` (16), `int max_centroid_size`, `double *BMWM_fused[16]`
- [x] 2.2 Extend `InterpBufPoolInit` in `ParticleRoutines.c`: when `interp_mode == 3`, compute `max_centroid_size = max(sizes[0..3])`, allocate 16 `BMWM_fused` arrays of `max_centroid_size` doubles each, also allocate `WM[4]` and `BMWM[4]` (for singleton fallback to mode 2 path)
- [x] 2.3 Extend `InterpBufPoolFree` in `ParticleRoutines.c`: when `interp_mode == 3`, free all 16 `BMWM_fused` arrays

## 3. P2Mastah_Fused Implementation

- [x] 3.1 Implement `P2Mastah_Fused` setup phase in `ParticleRoutines.c`: derive `Nx`, `Nz`, `X_vect`, `Z_vect` from centroid; compute `dx_itp`, `dz_itp` from `interp_stencil`; validate all `fields[f].stencil == interp_stencil`; resolve `pool_idx` from centroid
- [x] 3.2 Implement buffer init: memset shared `WM[pool_idx]` to zero; memset `BMWM_fused[0..nfields-1]` to zero (only up to `Nx*Nz` bytes per each, not `max_centroid_size`)
- [x] 3.3 Implement overflow guard: if `nfields > pool->max_fused_fields`, log error via `LOG_ERR` and fall back — call individual `P2Mastah` for each field descriptor using mode 2 path, then return
- [x] 3.4 Implement scatter loop (`#pragma omp parallel for`): for each active particle (`phase[k] != -1`), compute `ip`, `jp`, `dxm`, `dzm` once; for centre stencil node compute `w_` and atomically accumulate `w_` into `WM[kp]` once, then for each field `f` compute `mark_val` from flag/src and atomically accumulate `mark_val * w_` into `BMWM_fused[f][kp]`
- [x] 3.5 Implement 9-cell stencil neighbors (N, S, E, W, NE, NW, SE, SW) in scatter loop: for each neighbor, compute `kp`, `dxm`, `dzm`, `w_`; atomically accumulate `w_` into `WM[kp]` once; for each field `f` atomically accumulate `mark_val_f * w_` into `BMWM_fused[f][kp]` (handle `periodic_x` for E/W/NE/NW/SE/SW)
- [x] 3.6 Implement `mark_val` computation per field: handle `flag == 0` (read `fields[f].src[phase[k]]`), `flag == 1` (read `fields[f].src[k]`), `flag == -1` (compute `2*nx²*nz²`), `flag == -2` (compute `nx*nz*(-nx²+nz²)`), `flag == -3` (compute `acos(nx)`); apply harmonic (`1/val`) and geometric (`log(val)`) averaging transforms per `fields[f].avg`
- [x] 3.7 Implement normalization loop (`#pragma omp parallel for`): for each node `i` in `0..Nx*Nz-1`, skip if `WM[i] < 1e-30` or `NodeType==30/31` (use `fields[0].BCtype` since all share centroid), otherwise for each field `f`: `fields[f].dst[i] = BMWM_fused[f][i] / WM[i]`; apply harmonic inverse (`1/val`) if `avg==1`, geometric inverse (`exp(val)`) if `avg==2`
- [x] 3.8 Implement periodic-x boundary fixup: if `centroid == 0` and `model->periodic_x == 1`, for each field `f` with `prop == 0`, average the left and right boundary columns of `fields[f].dst`

## 4. P2Mastah Mode 2 Fallback for Mode 3 Singletons

- [x] 4.1 Verify that existing `P2Mastah` handles `imode >= 2` identically (i.e., `pool->interp_mode == 3` triggers the `else` branch same as mode 2) — no code changes expected, just confirm

## 5. Call-Site Refactoring — Pre-Solve Hot Path (Main_DOODZ.c lines 539–594)

- [x] 5.1 Refactor Batch 1 (cent, stencil=1, lines 539–547): wrap `Cp`, `Qr`, `T0_n`, `divth0_n` in `if (pool && pool->interp_mode == 3)` block building a 4-field descriptor array and calling `P2Mastah_Fused(..., cent, 1, ...)`, with `else` preserving original 4 `P2Mastah` calls
- [x] 5.2 Refactor Batch 2a (cent, stencil=interp, lines 554–565): build descriptor array with `p0_n`, `phi0_n`, and conditionally `sxxd0`, `szzd0` (if `elastic == 1`); call `P2Mastah_Fused(..., cent, interp_stencil, ...)`; preserve original calls in `else`
- [x] 5.3 Refactor Batch 2b (vert, stencil=interp, lines 555–566): build descriptor array with `p0_s`, and conditionally `sxz0` (if `elastic == 1`); call `P2Mastah_Fused(..., vert, interp_stencil, ...)`; preserve original calls in `else`
- [x] 5.4 Keep the `InterpCentroidsToVerticesDouble` / `InterpVerticesToCentroidsDouble` / `ShearModCompExpGrid` barrier between Batch 2 and Batch 3 untouched
- [x] 5.5 Refactor Batch 3a (cent, stencil=interp, lines 578–580): if `anisotropy == 1`, build 3-field descriptor array for `d1_n`, `d2_n`, `angle_n` (flags -1, -2, -3); call `P2Mastah_Fused(..., cent, interp_stencil, ...)`; preserve original calls in `else`
- [x] 5.6 Refactor Batch 3b (vert, stencil=interp, lines 581–583): if `anisotropy == 1`, build 3-field descriptor array for `d1_s`, `d2_s`, `angle_s` (flags -1, -2, -3); call `P2Mastah_Fused(..., vert, interp_stencil, ...)`; preserve original calls in `else`
- [x] 5.7 Keep the `FiniteStrainAspectRatio` barrier between Batch 3 and Batch 4 untouched
- [x] 5.8 Refactor Batch 4a (cent, stencil=interp, lines 587–594): build 3-field descriptor array for `X0_n`, `noise_n`, `d0_n`; call `P2Mastah_Fused(..., cent, interp_stencil, ...)`; preserve original calls in `else`
- [x] 5.9 Refactor Batch 4b (vert, stencil=interp, lines 588–590): build 2-field descriptor array for `X0_s`, `noise_s`; call `P2Mastah_Fused(..., vert, interp_stencil, ...)`; preserve original calls in `else`
- [x] 5.10 Keep kx/kz singleton calls (lines 542–543, vxnodes/vznodes) as existing `P2Mastah` — no changes needed

## 6. Local Validation

- [x] 6.1 Build with `make build` — verify zero compilation errors and zero warnings in new code
- [x] 6.2 Run full CI test suite locally with `interp_mode = 0` (default): `cd cmake-build && ctest --output-on-failure` — all tests pass (no regressions)
- [x] 6.3 Run BlankenBench locally with `interp_mode = 3` and `thermal_solver = 1`: verify completion and results within tolerance
- [x] 6.4 Spot-check numerical equivalence: BlankenBench 500 steps mode 0 vs mode 3 — T range identical, velocity/Nu within expected transient variance

## 7. CI Test Fixture

- [x] 7.1 Create BlankenBench optimized `.txt` variant with `interp_mode = 3` and `thermal_solver = 1`
- [x] 7.2 Add `FusedP2Mastah` test case in `BlankenBenchTests.cpp` using `BlankenBenchOptimized.txt` — reuses existing BlankenBench fixture and CMake target
- [x] 7.3 Verify CI test passes: 20/20 tests pass including FusedP2Mastah (177.64s total)

## 8. EC2 Benchmark (Experiment 7)

- [ ] 8.1 Push branch, start EC2 instance, pull and build on EC2
- [ ] 8.2 Run benchmark sweep: `--scenario RiftingComprehensive --resolutions "highres" --threads "1 2 4 6 8 12 16" --steps 10 --sed 's/^thermal_solver.*/thermal_solver = 1/' --sed 's/^interp_mode.*/interp_mode = 3/'`
- [ ] 8.3 Download results, verify `avg_interp_s < 1.5s` at 16 threads and no >10% regression at 1 thread vs mode 2 baseline
- [ ] 8.4 Document results as Experiment 7 in `skill-benchmarking/SKILL.md` with hypothesis, implementation summary, results table, and conclusion
