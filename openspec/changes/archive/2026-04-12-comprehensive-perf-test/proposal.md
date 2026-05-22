## Why

MDOODZ has no single integration-level scenario that exercises all major subsystems (anisotropy, melting, combined-yield plasticity, visco-plasticity, GSE, elasticity, thermal, free surface) together. Performance benchmarking currently uses RiftingBasic — a minimal setup that barely stresses the solver. We need a realistic, feature-comprehensive scenario that can serve as a performance baseline across machines and resolutions, showing how MDOODZ performs and how time is distributed across different parts of the code. The RiftingCombinedYield model on `origin/add-meltin-vp-model` is nearly there but lacks anisotropy and GSE, and runs only 1 timestep.

## What Changes

- Create a new scenario `RiftingComprehensive` (`.c` + `.txt`) based on `RiftingCombinedYield` from `origin/add-meltin-vp-model`:
  - Enable **anisotropy** (`aniso = 1`, `aniso_fstrain = 1`) on mantle phases (dry olivine, wet olivine) to activate the director fabric tracking subsystem
  - Enable **grain-size evolution** (`gsel = 1`) on at least one mantle phase to exercise GSE
  - Set **t_end = 15 Ma** (model time termination) instead of a fixed timestep count — this ensures all resolutions evolve to the same physical state regardless of adaptive dt, enabling meaningful cross-resolution comparison
  - Provide **four resolution variants**:
    - `_lowres.txt` — 100×80 (quick sanity check)
    - `_default.txt` — 200×160 (baseline benchmark)
    - `_medres.txt` — 501×401 (scaling study)
    - `_highres.txt` — 1000×800 (stress test)
- Create **AWS EC2 automation** (`misc/benchmark-aws.sh`):
  - Start an EC2 instance (e.g. c6i.4xlarge, 16 vCPUs)
  - Upload code via rsync/scp, install deps, build
  - Run benchmark suite across resolutions with thread counts: 1, 2, 4, 6, 8, 12, 16
  - Upload results (logs, perf.csv, REPORT.md, HDF5 outputs) to an S3 bucket
  - **Self-pause** (stop) the instance when finished — billing stops immediately
  - Enforce a **time limit of ~6 hours** per run; abort and upload partial results if exceeded
- Ensure **comprehensive per-subsystem timing** in `perf.csv`: every major subsystem (rheology, Stokes assembly, Stokes solve, thermal solve, advection, free surface, anisotropy update, melting, GSE, particle reseeding) must be individually timed and logged so the benchmark captures where time is actually spent
- Add **model time** to `perf.csv` and log metadata so that events can be correlated by physical time across resolutions. Default unit is Ma; configurable via `.txt` parameter (Ma, Ka, yr)
- Integrate with `misc/benchmark.sh`: add `--scenario` flag so the comprehensive scenario can be used instead of RiftingBasic
- Create a **Copilot skill** (`skill-benchmarking`) update documenting AWS prerequisites: AWS CLI configuration, EC2 SSH key pair, S3 bucket name, IAM permissions, instance type recommendations

This is **not CI** — it is a standalone benchmarking workflow run on-demand to understand MDOODZ performance characteristics (wall time breakdown per subsystem, thread scaling efficiency, memory usage vs. resolution).

## Capabilities

### New Capabilities
- `scenario-rift-comprehensive`: The new RiftingComprehensive scenario definition — phase layout, material properties (with anisotropy + GSE), boundary conditions, thermal setup, and four resolution parameter files
- `perf-scenario-integration`: Integration of the comprehensive scenario into the benchmark pipeline — `--scenario` flag for `benchmark.sh`, adapted `make_bench_txt` patching, summary.csv compatibility
- `aws-benchmark-automation`: AWS EC2 lifecycle script — instance start, code upload, build, benchmark execution, S3 result upload, instance self-pause, timeout enforcement
- `skill-aws-benchmarking`: Copilot skill documenting AWS dependencies (CLI keys, SSH key pair, S3 bucket, IAM policy, instance types) and usage workflow

### Modified Capabilities
- `scenario-rift-combined-yield`: The existing spec covers the base RiftingCombinedYield model. The new scenario extends it with anisotropy and GSE phases, so requirements for phase definitions and material property ranges are expanding.

## Impact

- **New files**: `SETS/RiftingComprehensive.c`, `SETS/RiftingComprehensive.txt` (default), `SETS/RiftingComprehensive_lowres.txt`, `SETS/RiftingComprehensive_medres.txt`, `SETS/RiftingComprehensive_highres.txt`, `misc/benchmark-aws.sh`
- **Modified files**: `SETS/CMakeLists.txt` (register new scenario), `misc/benchmark.sh` (add `--scenario` flag), `.github/skills/skill-benchmarking/SKILL.md` (AWS docs)
- **Dependencies**: Requires the `add-meltin-vp-model` branch changes to be merged first (combined yield, melting infrastructure). Requires AWS CLI configured with appropriate IAM credentials and an S3 bucket.
- **Runtime**: Low-res ~minutes; default ~15–30 min; medium ~1–2 hours; high-res potentially several hours. Full sweep (4 resolutions × 7 thread counts) runs sequentially with 6-hour timeout enforcement.
