## ADDED Requirements

### Requirement: AWS prerequisites section in benchmarking skill
`.github/skills/skill-benchmarking/SKILL.md` SHALL include an "AWS Benchmarking" section that documents all prerequisites for running benchmarks on EC2.

#### Scenario: Skill documents IAM permissions
- **WHEN** a user reads the AWS Benchmarking section
- **THEN** they find the required IAM permissions listed: `ec2:StartInstances`, `ec2:StopInstances`, `ec2:DescribeInstances`, `s3:PutObject`, `s3:GetObject`, `s3:ListBucket`, `sts:GetCallerIdentity`

### Requirement: SSH key setup documentation
The skill SHALL document how to create and configure an EC2 SSH key pair, including key generation, `chmod 400`, and setting the `BENCH_SSH_KEY` environment variable.

#### Scenario: SSH key instructions complete
- **WHEN** a user follows the SSH key setup instructions
- **THEN** they can successfully SSH to the EC2 instance using the documented key

### Requirement: Environment variable documentation
The skill SHALL document all required environment variables:
- `BENCH_EC2_INSTANCE`: EC2 instance ID (e.g. `i-0abc123def456`)
- `BENCH_S3_BUCKET`: S3 bucket name (e.g. `mdoodz-benchmarks`)
- `BENCH_SSH_KEY`: Path to the `.pem` SSH key file
- `BENCH_EC2_USER`: SSH username (default: `ubuntu`)

#### Scenario: Env vars documented with examples
- **WHEN** a user reads the environment variable section
- **THEN** each variable has a description, example value, and how to export it

### Requirement: Instance type recommendations
The skill SHALL recommend EC2 instance types suitable for MDOODZ benchmarking, with guidance on core count, memory, and cost trade-offs.

#### Scenario: Instance recommendations present
- **WHEN** a user reads the instance type section
- **THEN** they find at least 3 recommended instance types with vCPU count, RAM, and approximate hourly cost

### Requirement: Quick start workflow
The skill SHALL include a quick-start workflow showing the complete sequence from first-time setup to running a benchmark and retrieving results.

#### Scenario: Quick start covers full lifecycle
- **WHEN** a user follows the quick-start steps
- **THEN** they can: set up env vars, run `benchmark-preflight.sh`, run `benchmark-aws.sh`, and find results in S3

### Requirement: Cost and billing guidance
The skill SHALL document that stopped instances only incur EBS storage costs, and SHALL recommend stopping (not terminating) instances to preserve installed dependencies.

#### Scenario: Billing guidance present
- **WHEN** a user reads the cost section
- **THEN** they understand that compute billing stops on instance stop, and only storage costs (~$0.08/GB/month) continue
