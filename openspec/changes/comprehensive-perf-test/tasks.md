## 1. Core C Changes — Time-based Termination and perf.csv Expansion

- [x] 1.1 Add `t_end` parameter: read `t_end` in `InputOutput.c` via `ReadDou2(fin, "t_end", 0.0)`, add field to model struct in `mdoodz.h`
- [x] 1.2 Add `time_unit` parameter: read `time_unit` in `InputOutput.c` via `ReadInt2(fin, "time_unit", 0)` (0=Ma, 1=Ka, 2=yr), add field to model struct
- [x] 1.3 Modify main loop in `Main_DOODZ.c`: change termination condition to check `model.t_end > 0 ? model.time < model.t_end / scaling.t : model.step <= model.Nt`
- [x] 1.4 Add per-subsystem timing variables in `Main_DOODZ.c`: `dt_thermal`, `dt_advection`, `dt_free_surface`, `dt_reseeding`, `dt_melting`, `dt_anisotropy`, `dt_gse`, `dt_output`
- [x] 1.5 Instrument thermal solver block with `omp_get_wtime()` → `dt_thermal`
- [x] 1.6 Instrument advection block (full block including sub-steps) → `dt_advection`
- [x] 1.7 Instrument free surface operations within advection → `dt_free_surface`
- [x] 1.8 Instrument reseeding (`CountPartCell`) → `dt_reseeding`
- [x] 1.9 Instrument melting blocks (`MeltFractionGrid` + `UpdateAlphaCp`) → `dt_melting`
- [x] 1.10 Instrument anisotropy blocks (`UpdateAnisoFactor` + grid viscosity outside Stokes loop) → `dt_anisotropy`
- [x] 1.11 Instrument GSE blocks (both `UpdateParticleGrainSize` calls) → `dt_gse`
- [x] 1.12 Instrument output blocks (HDF5 write + breakpoint) → `dt_output`
- [x] 1.13 Update perf.csv header and fprintf to include: `time_ma`, `thermal_s`, `advection_s`, `free_surface_s`, `reseeding_s`, `melting_s`, `anisotropy_s`, `gse_s`, `output_s`
- [x] 1.14 Compute `time_ma` from `model.time * scaling.t / (1e6 * 365.25 * 24 * 3600)` (adjust column name for `time_unit`)
- [x] 1.15 Extend `mdoodz_log_set_step()` to also accept model time; update metadata prefix to include `T<value><unit>` (e.g. `T0.42Ma`)
- [x] 1.16 Build and verify: run RiftingBasic 3 steps, confirm new perf.csv columns exist with zeros for disabled subsystems

## 2. Scenario — RiftingComprehensive

- [ ] 2.1 Merge `origin/add-meltin-vp-model` into current branch (or rebase onto it)
- [ ] 2.2 Copy `SETS/RiftingCombinedYield.c` → `SETS/RiftingComprehensive.c`, update `RunMDOODZ()` call to use `"RiftingComprehensive.txt"`
- [ ] 2.3 Copy `SETS/RiftingCombinedYield.txt` → `SETS/RiftingComprehensive.txt` (200×160 default), set `t_end = 4.73e14`, `Nt = 99999`
- [ ] 2.4 Add anisotropy parameters to mantle phases (IDs 1, 2, 3, 5): `aniso = 1`, `aniso_fstrain = 1`, per-phase anisotropy factors
- [ ] 2.5 Add GSE to phase 1 (dry olivine): `gsel = 1` with grain-size parameters
- [ ] 2.6 Create `RiftingComprehensive_lowres.txt` — 100×80, same physics, tuned dt
- [ ] 2.7 Create `RiftingComprehensive_medres.txt` — 501×401, same physics, tuned dt
- [ ] 2.8 Create `RiftingComprehensive_highres.txt` — 1000×800, same physics, tuned dt
- [ ] 2.9 Register `RiftingComprehensive` in `SETS/CMakeLists.txt` via `add_set(RiftingComprehensive)`
- [ ] 2.10 Build the scenario: `cmake -DSET=RiftingComprehensive -DOPT=ON -DOMP=ON`

## 3. Model Stability Validation

- [ ] 3.1 Run lowres (100×80) for 3 timesteps — verify exit code 0, no NaN
- [ ] 3.2 Run default (200×160) for 3 timesteps — verify exit code 0, no NaN
- [ ] 3.3 Run medres (501×401) for 3 timesteps — verify exit code 0, no NaN
- [ ] 3.4 Run highres (1000×800) for 3 timesteps — verify exit code 0, no NaN
- [ ] 3.5 If any resolution crashes: diagnose, adjust parameters (eta_vp, min_eta, nit_max, anisotropy factors, dt), and re-validate all resolutions
- [ ] 3.6 Verify perf.csv has non-zero values for thermal_s, advection_s, melting_s, anisotropy_s, gse_s columns

## 4. Benchmark Pipeline — --scenario, --resolutions, --validate

- [ ] 4.1 Add `--scenario` flag to `misc/benchmark.sh` (default: `RiftingBasic`), build with `-DSET=<scenario>`
- [ ] 4.2 Add `--resolutions` flag: iterate over `<Scenario>_<res>.txt` files instead of `--grids` patching
- [ ] 4.3 Implement `--validate` mode: run 3 steps per resolution (low→high), 1 thread, bail on first failure
- [ ] 4.4 Update `make_bench_txt` to handle resolution-based runs (copy .txt file directly instead of patching Nx/Nz)
- [ ] 4.5 Update summary.csv awk to handle expanded perf.csv columns
- [ ] 4.6 Update `misc/benchmark-report.sh` to display new subsystem timing columns and model time
- [ ] 4.7 Test: `./misc/benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default" --threads "1 4" --validate`

## 5. AWS Preflight and Environment

- [ ] 5.1 Create `misc/benchmark-preflight.sh`: check AWS CLI, BENCH_EC2_INSTANCE, BENCH_SSH_KEY, BENCH_S3_BUCKET, SSH connectivity
- [ ] 5.2 Run preflight — if it fails, stop and guide user through setup
- [ ] 5.3 Verify: AWS CLI configured, EC2 instance accessible, SSH works, S3 bucket writable

## 6. AWS Automation Scripts

- [ ] 6.1 Create `misc/benchmark-aws.sh`: preflight → start instance → wait SSH → rsync repo → remote exec → download summary → stop instance
- [ ] 6.2 Add trap handler for EXIT/INT/TERM: always stop instance on exit
- [ ] 6.3 Add `--instance-id`, `--ssh-key`, `--s3-bucket`, `--scenario`, `--resolutions` flags
- [ ] 6.4 Create `misc/benchmark-ec2-setup.sh`: install deps (apt), build, run `benchmark.sh --validate --scenario <name> --resolutions <list> --threads "1 2 4 6 8 12 16"`
- [ ] 6.5 Wrap remote execution in `timeout 6h`
- [ ] 6.6 Add incremental S3 upload after each resolution completes
- [ ] 6.7 Add S3 upload on timeout (partial results)
- [ ] 6.8 Test locally: run `benchmark-ec2-setup.sh` in WSL to verify build + benchmark logic works

## 7. Skill Documentation

- [ ] 7.1 Update `.github/skills/skill-benchmarking/SKILL.md` with AWS Benchmarking section: IAM permissions, SSH key setup, env vars, instance recommendations, quick-start workflow, cost guidance
- [ ] 7.2 Add `*.pem` to `.gitignore` if not already present
- [ ] 7.3 Update `.github/copilot-instructions.md` skill table if needed

## 8. End-to-End Verification

- [ ] 8.1 Run full local benchmark: `./misc/benchmark.sh --scenario RiftingComprehensive --resolutions "lowres default" --threads "1 2 4" --validate`
- [ ] 8.2 Generate report: `./misc/benchmark-report.sh benchmark-results/<latest>/`
- [ ] 8.3 Verify report includes: model time column, all subsystem timing breakdowns, thread scaling, grid scaling
- [ ] 8.4 Run AWS benchmark (if preflight passed): `./misc/benchmark-aws.sh`
- [ ] 8.5 Verify results appear in S3 and instance is stopped after completion
