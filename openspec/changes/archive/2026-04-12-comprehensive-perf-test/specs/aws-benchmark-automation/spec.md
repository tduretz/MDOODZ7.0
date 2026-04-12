## ADDED Requirements

### Requirement: Preflight check script
`misc/benchmark-preflight.sh` SHALL verify the full AWS environment is operational before any benchmark runs. It SHALL check:
1. AWS CLI is installed and configured (`aws sts get-caller-identity` succeeds)
2. `BENCH_EC2_INSTANCE` env var is set and the instance exists (`aws ec2 describe-instances`)
3. `BENCH_SSH_KEY` env var points to a readable `.pem` file
4. `BENCH_S3_BUCKET` env var is set and the bucket is accessible (`aws s3 ls`)
5. SSH connectivity to the instance (with a 10-second timeout)

Each check SHALL print a pass/fail status. If any check fails, the script SHALL exit with code 1 and print which prerequisite is missing with setup instructions.

#### Scenario: All checks pass
- **WHEN** AWS CLI is configured, instance exists, SSH key is valid, and S3 bucket is accessible
- **THEN** the script exits with code 0 and prints "All preflight checks passed"

#### Scenario: Missing AWS CLI
- **WHEN** `aws` command is not found
- **THEN** the script prints "FAIL: AWS CLI not installed. Install with: pip install awscli" and exits with code 1

#### Scenario: Missing environment variable
- **WHEN** `BENCH_EC2_INSTANCE` is not set
- **THEN** the script prints "FAIL: BENCH_EC2_INSTANCE not set. Export your EC2 instance ID." and exits with code 1

#### Scenario: SSH key not found
- **WHEN** `BENCH_SSH_KEY` points to a non-existent file
- **THEN** the script prints "FAIL: SSH key not found at <path>" and exits with code 1

### Requirement: EC2 lifecycle orchestration script
`misc/benchmark-aws.sh` SHALL orchestrate the full benchmark lifecycle on a remote EC2 instance:
1. Run `benchmark-preflight.sh` — abort if it fails
2. Start the EC2 instance (`aws ec2 start-instances`)
3. Wait for the instance to be running and SSH-accessible
4. Rsync the repository to the instance
5. Execute the remote setup/run script via SSH
6. Download summary results to local machine
7. Stop the instance (`aws ec2 stop-instances`) — billing stops

The script SHALL accept `--instance-id`, `--ssh-key`, `--s3-bucket` flags that override the corresponding environment variables. It SHALL also accept `--scenario` (default: `RiftingComprehensive`) and `--resolutions` flags to pass to the remote benchmark.

#### Scenario: Full lifecycle completes
- **WHEN** `./misc/benchmark-aws.sh` is run with valid credentials
- **THEN** the instance starts, code is uploaded, benchmark runs, results are uploaded to S3, and the instance is stopped

#### Scenario: Preflight failure aborts
- **WHEN** `benchmark-preflight.sh` fails
- **THEN** the script prints the preflight errors and exits without starting the instance

#### Scenario: Instance stops on error
- **WHEN** the benchmark crashes mid-run
- **THEN** the instance is still stopped (via trap handler) and partial results are uploaded to S3

### Requirement: Trap handler ensures instance stops
`benchmark-aws.sh` SHALL register a trap handler for `EXIT`, `INT`, and `TERM` signals. The handler SHALL always stop the EC2 instance, even if the script is interrupted with Ctrl-C or encounters an error.

#### Scenario: Ctrl-C stops instance
- **WHEN** the user presses Ctrl-C during a benchmark run
- **THEN** the trap handler stops the EC2 instance before the script exits

### Requirement: Remote execution script
`misc/benchmark-ec2-setup.sh` SHALL run on the EC2 instance and:
1. Install dependencies if not already present (apt: build-essential, cmake, libsuitesparse-dev, libhdf5-dev, libblas-dev, liblapack-dev)
2. Build the scenario with `cmake -DOPT=ON -DOMP=ON -DSET=<scenario>`
3. Run `benchmark.sh --validate --scenario <scenario> --resolutions <resolutions> --threads "1 2 4 6 8 12 16"`
4. Upload results to S3 incrementally (after each resolution completes)
5. Generate the benchmark report

The script SHALL be wrapped in `timeout 6h` to enforce the time limit. On timeout, it SHALL upload partial results before exiting.

#### Scenario: Dependencies already installed
- **WHEN** the script runs on an instance that was previously set up
- **THEN** it skips dependency installation and proceeds to build

#### Scenario: 6-hour timeout
- **WHEN** the benchmark exceeds 6 hours
- **THEN** the `timeout` command kills the run, partial results are uploaded to S3, and the exit code indicates timeout

### Requirement: S3 result upload
Results SHALL be uploaded to `s3://<bucket>/mdoodz-benchmarks/<hostname>/<timestamp>/`. The upload SHALL include: `summary.csv`, `system.txt`, `REPORT-*.md`, all `perf.csv` files, and `stdout.log` files. HDF5 output files SHALL be uploaded only if present.

Upload SHALL happen:
- After each resolution completes (incremental)
- On timeout (partial results)
- On completion (final)

#### Scenario: Incremental upload after resolution
- **WHEN** the lowres benchmark completes
- **THEN** its `perf.csv` and `stdout.log` are uploaded to S3 before the default resolution starts

#### Scenario: Partial upload on timeout
- **WHEN** the benchmark times out during medres
- **THEN** results for lowres and default (completed) plus partial medres are uploaded to S3

### Requirement: Instance self-pause after completion
After the remote script finishes (success, failure, or timeout), the EC2 instance SHALL be stopped either by `benchmark-aws.sh` (from local) or by the remote script itself as a fallback. The instance SHALL NOT be terminated — stopping preserves the disk for future runs.

#### Scenario: Instance stopped after successful run
- **WHEN** `benchmark-aws.sh` completes successfully
- **THEN** `aws ec2 describe-instances` shows the instance in `stopped` state
