## Context

MDOODZ currently has `perf.csv` output per timestep with columns: `step, wall_s, rheology_s, assembly_s, solve_s, nit, n_particles, neq_mom, neq_cont, peak_rss_mb, user_cpu_s, sys_cpu_s`. Only three subsystems (rheology, Stokes assembly, Stokes solve) are individually timed — the remaining wall time (thermal solve, advection, reseeding, melting, anisotropy, GSE, free surface, I/O) is unaccounted for.

The benchmark pipeline (`misc/benchmark.sh`) is hardcoded to `RiftingBasic`, which exercises only basic creep + elasticity. The `RiftingCombinedYield` model on `origin/add-meltin-vp-model` already has 9 phases with melting + combined yield + visco-plasticity but no anisotropy or GSE.

There is no AWS automation — benchmarks must be run manually.

## Goals / Non-Goals

**Goals:**
- Create a `RiftingComprehensive` scenario that exercises every major subsystem simultaneously
- Expand `perf.csv` to capture per-subsystem timing for: thermal, advection, free surface, reseeding, melting, anisotropy, GSE, I/O — in addition to existing rheology/assembly/solve
- Allow `benchmark.sh` to run any scenario (not just RiftingBasic)
- Automate the full lifecycle on AWS EC2: provision → build → benchmark → upload to S3 → self-stop
- Run sweeps across 4 resolutions × 7 thread counts with a 6-hour timeout
- Terminate simulations by model time (t_end) not timestep count, so all resolutions reach the same physical state
- Include model time in perf.csv and log metadata for cross-resolution event correlation

**Non-Goals:**
- CI integration — this is an on-demand performance analysis tool, not part of CI
- Correctness validation — no reference solution comparison (separate concern)
- Multi-node / MPI benchmarking — single-machine OpenMP only
- Automated instance type selection — user chooses the EC2 instance

## Decisions

### 1. Scenario based on RiftingCombinedYield + anisotropy + GSE

**Decision**: Fork `RiftingCombinedYield.c/.txt` into `RiftingComprehensive.c/.txt`. Add `aniso = 1, aniso_fstrain = 1` on mantle phases (IDs 1, 2, 3, 5) and `gsel = 1` on dry olivine mantle phase (ID 1). Keep crust phases isotropic (realistic: only mantle develops fabric).

**Alternative considered**: Enabling anisotropy on all phases. Rejected — unrealistic and may cause convergence issues in crustal phases where anisotropy parameters are untested.

**Rationale**: Mantle-only anisotropy is physically motivated and exercises the full anisotropy pipeline (director tracking, ViscosityConciseAniso, finite strain) without destabilising the model.

### 2. Four resolution tiers as separate .txt files

**Decision**: Provide `RiftingComprehensive.txt` (200×160 default), `_lowres.txt` (100×80), `_medres.txt` (501×401), `_highres.txt` (1000×800). All share the same `.c` file — only grid size and writer frequency differ.

**Alternative considered**: A single .txt with runtime patching by benchmark.sh. Rejected — resolution changes may need adjusted `dt` or Courant number, and having explicit files lets users run them standalone outside the benchmark pipeline.

**Rationale**: Each resolution tier can be individually tuned (dt, Courant) for stability. The benchmark script selects which ones to run.

### 3. Expand perf.csv with per-subsystem columns

**Decision**: Add timing columns to `perf.csv`:

| New column | What it measures |
|-----------|-----------------|
| `thermal_s` | `EnergyDirectSolve` |
| `advection_s` | Full advection block (including sub-stepping) |
| `free_surface_s` | Free surface operations within advection |
| `reseeding_s` | `CountPartCell` reseeding |
| `melting_s` | `MeltFractionGrid` + `UpdateAlphaCp` |
| `anisotropy_s` | `UpdateAnisoFactor` + `NonNewtonianViscosityGridAniso` (outside Stokes loop) |
| `gse_s` | `UpdateParticleGrainSize` (both calls) |
| `output_s` | HDF5 write + breakpoint |

The header becomes: `step,wall_s,time_ma,rheology_s,assembly_s,solve_s,thermal_s,advection_s,free_surface_s,reseeding_s,melting_s,anisotropy_s,gse_s,output_s,nit,n_particles,neq_mom,neq_cont,peak_rss_mb,user_cpu_s,sys_cpu_s`

**Alternative considered**: Separate timing log file. Rejected — keeping everything in one CSV makes analysis straightforward.

**Rationale**: Adds 9 columns (8 timing + model time) but all are doubles. The CSV remains easily parseable. The `benchmark-report.sh` will need updating to display the new columns.

### 3a. Model time in perf.csv and log metadata

**Decision**: Add a `time_ma` column to `perf.csv` (model time in Ma). Computed as `model.time * scaling.t / (1e6 * 365.25 * 24 * 3600)`. This enables:
- Cross-resolution comparison at the same physical time
- Correlating events of note (e.g. onset of melting, strain localisation) by model time

Also extend `mdoodz_log_set_step()` to accept model time, and include it in the log metadata prefix:
```
[    6.565|S0001|T0.42Ma|N00] INFO  | ...
```

The time unit is configurable via a new `.txt` parameter:
```
time_unit = 0   / 0=Ma (default), 1=Ka, 2=yr
```
The perf.csv column name changes accordingly (`time_ma`, `time_ka`, `time_yr`). Default is Ma.

### 3b. Time-based simulation termination

**Decision**: Add a `t_end` parameter to `.txt` files (in seconds, like `dt`):
```
t_end = 4.73e14   / 15 Ma in seconds
```

The main loop condition changes from:
```c
for (; input.model.step <= input.model.Nt; input.model.step++)
```
to:
```c
for (; (input.model.t_end > 0 ? input.model.time < input.model.t_end / scaling.t 
        : input.model.step <= input.model.Nt); input.model.step++)
```

When `t_end > 0`, the simulation runs until model time reaches `t_end`, ignoring `Nt`. When `t_end = 0` (default), the existing `Nt`-based loop is preserved — no breaking change.

The RiftingComprehensive `.txt` files will use `t_end = 4.73e14` (15 Ma) and set `Nt = 99999` as a safety cap.

**Alternative considered**: Using Ma directly in the .txt file. Rejected — the code already uses seconds for `dt` and other time parameters; consistency is more important.

**Rationale**: Different resolutions with adaptive dt will take different numbers of timesteps to reach 15 Ma. Comparing by step number is meaningless for benchmarking — comparing by model time is the only way to correlate physical events across resolutions.

### 4. `benchmark.sh --scenario` flag and validation mode

**Decision**: Add `--scenario <name>` flag (default: `RiftingBasic`). The script builds with `-DSET=<name>` and looks for `<name>.txt` in the SETS directory. The `make_bench_txt` function patches Nx/Nz/Nt as before.

For the comprehensive benchmark, the user specifies `--scenario RiftingComprehensive` and the script iterates over the resolution-specific .txt files instead of patching a single file.

Add `--resolutions` flag: `--resolutions "lowres default medres highres"`. Each maps to `RiftingComprehensive_<res>.txt`. The grid size comes from the .txt file itself (not from `--grids`).

Add `--validate` flag: runs a quick validation pass before the full benchmark. The validation:
1. Runs each resolution (100×80 → 200×160 → 501×401 → 1000×800) for **3 timesteps** with 1 thread
2. Proceeds in order from lowest to highest — if a resolution crashes, stops immediately and reports which resolution failed and the error
3. Only after all resolutions pass does the full benchmark sweep (resolutions × thread counts × 500 steps) begin
4. This catches convergence issues (bad `eta_vp`, unstable anisotropy params, etc.) in seconds rather than hours

The AWS script (`benchmark-aws.sh`) always runs with `--validate` by default.

**Important**: The validation pass is also part of the implementation workflow. Before the scenario is considered ready for benchmarking, the agent must:
1. Build the scenario
2. Run validation (3 steps at each resolution, lowest to highest)
3. If any resolution crashes or diverges — diagnose and fix (adjust `eta_vp`, `min_eta`, `nit_max`, anisotropy factors, `dt`, etc.)
4. Re-run validation until all resolutions pass
5. Only then proceed to commit the scenario files and move to the benchmark automation tasks
6. After model validation passes, run `misc/benchmark-preflight.sh` to verify AWS connectivity (CLI, EC2 instance, SSH key, S3 bucket)
7. If preflight fails, stop and guide the user through setup — do not proceed to AWS automation tasks

This ensures the model is stable and the infrastructure is ready before any long-running benchmark is ever launched.

**Alternative considered**: Auto-detecting resolution from .txt file name. Rejected — too fragile.

### 5. AWS preflight checks — gate before implementation proceeds

**Decision**: Before any AWS-related tasks are implemented (and before the benchmark automation is applied), the agent must verify that the full AWS environment is operational:

1. **AWS CLI installed and configured**: `aws sts get-caller-identity` succeeds
2. **EC2 instance exists and is accessible**: `aws ec2 describe-instances --instance-ids $BENCH_EC2_INSTANCE` returns a valid instance
3. **SSH key exists and has correct permissions**: `BENCH_SSH_KEY` file exists, is readable, and `ssh -o ConnectTimeout=10` to the instance succeeds (instance must be started temporarily for this check)
4. **S3 bucket exists and is writable**: `aws s3 ls s3://$BENCH_S3_BUCKET` succeeds, and a test upload/delete works
5. **Required env vars are set**: `BENCH_EC2_INSTANCE`, `BENCH_S3_BUCKET`, `BENCH_SSH_KEY` are all non-empty

If any check fails, the agent must stop and tell the user exactly what is missing — with instructions on how to set it up (referencing the skill doc). Implementation must not proceed until all preflight checks pass.

This is codified as a `misc/benchmark-preflight.sh` script that can be run standalone. The `benchmark-aws.sh` script calls it automatically at startup. The agent runs it during the implementation workflow.

### 6. AWS EC2 automation as a self-contained script

**Decision**: Create `misc/benchmark-aws.sh` that:
1. **Starts** a stopped EC2 instance by instance ID (user provides via `--instance-id` or env var `BENCH_EC2_INSTANCE`)
2. **Waits** for SSH to become available
3. **Rsyncs** the repo to the instance
4. **Runs** a remote setup script (`misc/benchmark-ec2-setup.sh`) that: installs deps (apt), builds MDOODZ, runs `benchmark.sh --scenario RiftingComprehensive`, uploads results to S3
5. **Downloads** the summary back to local machine
6. **Stops** the instance (self-pause — billing stops)

A companion `misc/benchmark-ec2-setup.sh` runs on the instance itself. It includes a 6-hour `timeout` wrapper and always uploads partial results to S3 before exiting.

**Alternative considered**: Using AWS CDK/CloudFormation for instance management. Rejected — over-engineered for a single benchmark machine. Plain `aws ec2 start-instances` / `stop-instances` is sufficient.

**Alternative considered**: Terminating the instance after runs. Rejected — stopping preserves the disk so deps don't need reinstalling each time.

### 7. S3 result storage

**Decision**: Upload results to `s3://<bucket>/mdoodz-benchmarks/<hostname>/<timestamp>/`. Structure mirrors the local `benchmark-results/` layout. Bucket name provided via `--s3-bucket` or env var `BENCH_S3_BUCKET`.

Upload happens:
- After each resolution completes (incremental)
- On timeout (partial results)
- On completion (final)

### 8. Skill documentation for AWS prerequisites

**Decision**: Update `.github/skills/skill-benchmarking/SKILL.md` with a new "AWS Benchmarking" section covering:
- Required IAM permissions (`ec2:StartInstances`, `ec2:StopInstances`, `ec2:DescribeInstances`, `s3:PutObject`, `s3:GetObject`)
- SSH key pair setup (key name, `.pem` file location)
- S3 bucket creation
- Instance preparation (first-time setup vs. subsequent runs)
- Environment variables (`BENCH_EC2_INSTANCE`, `BENCH_S3_BUCKET`, `BENCH_SSH_KEY`)

## Risks / Trade-offs

**[Anisotropy + melting convergence]** → The combination of anisotropy, melting, and combined yield has not been tested together. Mitigation: start with the default resolution (200×160) and validate convergence before running the full sweep. Use `nit_max = 20` with Picard-to-Newton switching.

**[6-hour timeout may truncate high-res runs]** → At 1000×800 with 500 timesteps, high-res may not finish. Mitigation: the script uploads partial results incrementally and the report generator handles incomplete perf.csv gracefully.

**[AWS cost if instance left running]** → User forgets to stop. Mitigation: the script always stops the instance in a trap handler (even on error/Ctrl-C). Document this in the skill.

**[perf.csv column expansion breaks existing report scripts]** → Mitigation: `benchmark-report.sh` will be updated simultaneously. Old perf.csv files with fewer columns will still parse (awk handles missing fields).

**[SSH key security]** → The `.pem` key file must not be committed to the repo. Mitigation: `.gitignore` already covers `*.pem`. The skill documents this explicitly.
