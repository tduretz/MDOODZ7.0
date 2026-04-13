## 1. Parameter Reading

- [x] 1.1 Add `interp_mode` (int, default 0) field to the `params` struct in `mdoodz.h` (near `interp_stencil`)
- [x] 1.2 Read `interp_mode` in `InputOutput.c` using `ReadInt2` alongside existing interpolation parameters
- [x] 1.3 Clamp invalid `interp_mode` values to 0 after reading (anything other than 0, 1, 2 → treat as 0)

## 2. InterpBufPool Struct and Lifecycle

- [x] 2.1 Define `InterpBufPool` struct in `mdoodz-private.h` with `WM[4]`, `BMWM[4]`, `sizes[4]`, `Wm[4]`, `BmWm[4]`, `nthreads`, `interp_mode` fields
- [x] 2.2 Implement `InterpBufPoolInit(grid *mesh, int interp_mode, int nthreads)` in `ParticleRoutines.c` — allocate `WM`/`BMWM` for all 4 centroid sizes; for mode 1 also allocate thread-local `Wm[c][nthreads][size]` and `BmWm[c][nthreads][size]`; for mode 2 set `Wm`/`BmWm` to NULL
- [x] 2.3 Implement `InterpBufPoolFree(InterpBufPool *pool)` in `ParticleRoutines.c` — free all allocated arrays and the pool struct itself; handle `pool == NULL` gracefully
- [x] 2.4 Declare `InterpBufPoolInit` and `InterpBufPoolFree` in `mdoodz-private.h`

## 3. P2Mastah Signature and Call Sites

- [x] 3.1 Add `InterpBufPool *pool` as the last parameter to `P2Mastah` declaration in `mdoodz-private.h`
- [x] 3.2 Update `P2Mastah` definition in `ParticleRoutines.c` to accept the new `pool` parameter
- [x] 3.3 Update all ~30 `P2Mastah` call sites in `Main_DOODZ.c` to pass `pool`
- [x] 3.4 Update `CountPartCell` signature to accept `InterpBufPool *pool` and pass it to its `P2Mastah` call
- [x] 3.5 Update `CountPartCell` call sites in `Main_DOODZ.c` to pass `pool`

## 4. Pool Lifecycle in Main_DOODZ.c

- [x] 4.1 Add `InterpBufPool *pool = NULL;` declaration and `InterpBufPoolInit` call after grid initialization, before first `P2Mastah` call
- [x] 4.2 Add `InterpBufPoolFree(pool)` at simulation cleanup

## 5. Mode 1 — Persistent Buffers with Thread-Local Scatter

- [x] 5.1 Add centroid-to-pool-index mapping (`pool_idx` switch) at top of `P2Mastah`
- [x] 5.2 Add mode dispatch: when `pool != NULL && pool->interp_mode == 1`, use `pool->WM[pool_idx]`/`pool->BMWM[pool_idx]` (memset-zeroed) instead of `DoodzCalloc` for `WM`/`BMWM`; use `pool->Wm[pool_idx][k]`/`pool->BmWm[pool_idx][k]` (memset-zeroed) instead of per-call thread-local allocations
- [x] 5.3 For mode 1: keep `Wm_ph` per-call alloc/free unchanged when `prop==1`
- [x] 5.4 For mode 1: skip `DoodzFree` for `WM`, `BMWM`, `Wm`, `BmWm` at function end (pool owns them)
- [x] 5.5 Build and run locally with `interp_mode = 0` — verify bit-identical output (backward compat)
- [x] 5.6 Build and run locally with `interp_mode = 1` — verify bit-identical output to mode 0

## 6. Mode 2 — Atomic Scatter

- [x] 6.1 Add mode 2 branch: when `pool != NULL && pool->interp_mode == 2`, use `pool->WM[pool_idx]`/`pool->BMWM[pool_idx]` (memset-zeroed) as direct scatter targets
- [x] 6.2 In mode 2 scatter loop: replace `Wm[thread_num][kp] += ...` and `BmWm[thread_num][kp] += ...` with `#pragma omp atomic` writes to `WM[kp]` and `BMWM[kp]`
- [x] 6.3 In mode 2 scatter loop for `prop==1`: replace `Wm_ph[thread_num][p][kp] += ...` with `#pragma omp atomic` writes to `mesh->phase_perc_s[p][kp]` or `mesh->phase_perc_n[p][kp]`
- [x] 6.4 In mode 2: skip the reduction loop (lines ~2258–2271) since `WM`/`BMWM` already contain accumulated values
- [x] 6.5 In mode 2: skip `Wm_ph` allocation entirely (no thread-local phase arrays needed)
- [x] 6.6 Build and run locally with `interp_mode = 2` — verify results match mode 0 within machine epsilon

## 7. SETS Parameter Files

- [x] 7.1 Add `interp_mode = 0` to all 4 RiftingComprehensive `.txt` files (after `interp_stencil` line) so `--sed` substitution works for benchmarking

## 8. CI Tests

- [x] 8.1 Create `interp_mode = 1` variant `.txt` files for all 19 test suites (except BlankenBench) — copy existing `.txt` and append `interp_mode = 1`
- [x] 8.2 Create `interp_mode = 2` variant `.txt` files for all 19 test suites (except BlankenBench) — copy existing `.txt` and append `interp_mode = 2`
- [x] 8.3 Add GTest cases `<TestName>InterpMode1` for each of the 19 suites — same assertions as baseline (bit-identical for mode 1)
- [x] 8.4 Add GTest cases `<TestName>InterpMode2` for each of the 19 suites — relaxed tolerance assertions (1e-10 relative)
- [x] 8.5 Build and run full CI test suite locally — all mode 0, 1, 2 variants pass

## 9. EC2 Benchmarking

- [ ] 9.1 Run EC2 benchmark: RiftingComprehensive 1001×801, 10 steps, threads 1/2/4/6/8/12/16 with `interp_mode = 2` (atomic scatter — expect largest gains)
- [ ] 9.2 Run EC2 benchmark: same config with `interp_mode = 1` (persistent buffers only — isolate allocation benefit)
- [ ] 9.3 Download results, compare against baseline (Experiment 4 CHOLMOD+PCG data), document as Experiment 5 in skill-benchmarking SKILL.md
