## Why

The `add-gmg-stokes-solver` change shipped a working GMG-preconditioned FGMRES path for Stokes (26/26 tests green, Newton and Picard SolVi equivalence at ≤ 1e-8 vs CHOLMOD), but deliberately deferred three classes of work: the performance and memory regression harness (§12), the unfinished integration fixtures for `DISABLED_*` benchmarks (§15.16, §11.2–11.8), and the public-facing defence of the effort. Without those, the change is merged but not *justified* — the memory-scaling and wall-time numbers that made the whole project worthwhile are not measured, the remaining benchmarks that would lock down correctness under anisotropy / free surface / Newton have no live fixtures, and there is no single document a reviewer or collaborator can read to understand why GMG matters for MDOODZ.

This change closes that loop. The headline deliverable is a single printable **defence** document that ties theory → technology → implementation → measurements into a coherent narrative (primer-series style, with a V-cycle animation and clean process diagrams). The supporting engineering makes the numbers in that defence real: an EC2-backed perf-and-memory harness that measures α ≤ 1.2 and ≥ 3× speedup at 201², and the sweep through the unfinished tasks from the previous change.

## What Changes

- **NEW — defence document**: `defence.md` + `defence.pdf` alongside the existing reading library (primer / codes / stokes_internals / journey / anisotropy_seismology / phenomena). Narrative walkthrough of the GMG change covering motivation, theory, implementation choices (D1–D11), measurement results, test receipts, and limitations. Embeds the V-cycle animation and 3–5 clean process diagrams.
- **NEW — V-cycle animation**: the `gmg_dump_vcycle = 1` instrumented-run plus a Python post-processor that produces an animated GIF showing residual and error at every V-cycle step (pre-smooth, restrict, coarse solve, prolongate, post-smooth), with the spatial power spectrum of the residual on a companion panel.
- **NEW — reusable EC2 performance harness**: a scripted benchmark runner that provisions on AWS via the existing MacBook AWS CLI setup, runs MDOODZ at a controlled set of grid sizes under `lin_solver = 0` and `lin_solver = 3`, captures wall time and peak RSS, and emits a reproducible results artefact (CSV + regression plots). Structured for re-use on later changes (not one-off).
- **NEW — grid-sweep performance tests**: 41², 81², 161², 321², 801² memory-scaling points for the α ≤ 1.2 regression; 201² wall-time comparison for the ≥ 3× bound; OpenMP strong-scaling sweep 1 → 32 threads at 201².
- **NEW — integration fixtures for the remaining `DISABLED_*` benchmarks**: `.txt` twins for each (free-surface relaxation, anisotropic shear band, high-contrast SolKz, power-law shear zone, conjugate shear bands), dual-solver comparison fixtures asserting L2 equivalence vs CHOLMOD within 1e-8.
- **IN SCOPE — engine fixes surfaced by the harness or fixtures**: if enabling a benchmark or running the perf sweep exposes a latent bug in MDOODZ core (assembly, rheology, BC handling, marker machinery, etc.), fixing it is part of this change rather than deferred. The defence document calls out every such fix explicitly so the reader sees the full cost/benefit picture of the GMG work, including which latent bugs its test discipline surfaced.
- **MODIFIED — skill-solvers spec**: entry documenting `lin_solver = 3`, the `gmg_*` parameters, the dual-solver fixture pattern, and the convergence-factor diagnostic. Carried from the deferred task 15.3–15.4 / 16.3–16.4 in the previous change.
- **NEW — process diagrams for the defence**: 3–5 figures covering the multigrid V-cycle, the decoupled + bridge architecture (MDOODZ assembler ↔ `StokesAssemblyGMG` matvec ↔ V-cycle), the Picard→Newton phasing, the active-mask handling above the marker-chain free surface, and the FGMRES-with-GMG-preconditioner control flow. Intended to be printable and intuitive, not overloaded.

This change is **behaviourally non-breaking by default**: the perf harness is additive, the defence document is a reading artefact, the fixtures exercise existing code, and the skill-solvers spec entry is documentation. The one situation where production behaviour *may* change is if the new fixtures surface a latent bug in an engine path whose fix corrects previously-wrong output; any such fix is called out individually in the defence with before/after comparisons so the change is not silent.

## Capabilities

### New Capabilities

- `gmg-performance-harness`: reusable EC2-backed benchmark runner for MDOODZ perf + memory measurements. Covers host provisioning via AWS CLI, deterministic grid-sweep runs under both `lin_solver` values, wall-time and peak-RSS capture, CSV-plus-plot output artefact format, and a stable interface that future changes can re-invoke for their own regression checks. Not MDOODZ-internal — it lives in a new top-level `benchmarks/ec2/` directory with its own scripts and documentation.
- `gmg-defence-materials`: the defence document, V-cycle animation, and process diagrams as first-class spec-level deliverables. Distinct from the existing reading library primers in that it carries normative content (measurements, test receipts) that other work will cite as evidence.

### Modified Capabilities

- `skill-solvers`: the documentation skill's behaviour contract gains a section on `lin_solver = 3` and the `gmg_*` parameter family, carried over from the deferred docs task in the previous change.

## Impact

**New files**

- `defence.md` + `defence.pdf` + `figs/defence_*.png` — the defence document and its diagrams.
- `figs/vcycle_*.gif` (and per-frame PNGs) — the V-cycle animation.
- `benchmarks/ec2/run_perf_sweep.sh`, `benchmarks/ec2/provision.sh`, `benchmarks/ec2/analyse.py`, `benchmarks/ec2/README.md` — the perf harness.
- `TESTS/PerfHarness/*.cpp` — small test fixtures the harness drives (one per grid size, parameterised over `lin_solver`).
- `SETS/SolViBenchmark/SolViRes51_anisotropy_gmg.txt`, `SETS/TopoBench/TopoRelax_gmg.txt`, `SETS/ShearBand/SingleBand_gmg.txt`, `SETS/ShearBand/ConjugateBands_gmg.txt`, and their CHOLMOD twins — the `.txt` twins for the remaining benchmarks.
- `TESTS/DisabledBenchmarks/*.cpp` — the integration fixtures that enable the previously-disabled tests.

**Modified files**

- `MDLIB/include/mdoodz.h` + `MDLIB/InputOutput.c` — add the `gmg_dump_vcycle` debug parameter (default `0`, no production impact).
- `MDLIB/MultigridStokes.c` — dump residual + error snapshots at every V-cycle operator when `gmg_dump_vcycle = 1` is active.
- `openspec/specs/skill-solvers/spec.md` — a MODIFIED delta adding the `lin_solver = 3` section.
- `TESTS/CMakeLists.txt` — register the new fixtures.
- **Any MDLIB file needed for an engine fix** — if a benchmark or harness run surfaces a latent bug in core MDOODZ (e.g. `StokesAssemblyDecoupled.c`, `RheologyParticles.c`, `FreeSurface.c`, `AnisotropyRoutines.c`), the fix lands here. Each fix is explicitly listed in the defence document with the originating test, root cause, and before/after receipt.

**APIs / user-facing**

- One new debug `.txt` parameter (`gmg_dump_vcycle = 1`, default off) enabling V-cycle instrumentation.
- No change to HDF5 output schema, rheology interfaces, or solver dispatch logic.

**External dependencies**

- **AWS CLI** already present on the MacBook agent workstation; the perf harness uses it but the MDOODZ build itself gains no AWS dependency. The local build and CI remain AWS-free.
- No new scientific-computing dependencies.

**Risks**

- **EC2 runs produce numbers that vary with instance type and AWS noise**. Mitigated by (a) pinning a specific instance type in the harness config, (b) running each measurement N=5 times and reporting median + IQR, (c) rejecting runs where thermal throttling is detected.
- **Defence document scope creep**. Mitigated by writing the proposal to a strict page budget (~25 pages printable) and scoping the V-cycle animation to a single small problem.
- **Disabled-benchmark fixtures expose latent bugs**. This is a feature, not a risk — the whole point of enabling them is to find any remaining correctness gaps. Engine fixes are in scope for this change (see "What Changes"); each is called out individually in the defence document. If a surfaced bug is large enough that fixing it would blow the change's page/time budget, it is deferred to a targeted follow-up; the fixture stays `DISABLED_` and the deferral is logged in STATUS.md with the specific scope boundary.
- **Performance numbers disappoint**. Also in scope: if α > 1.2 or the 201² speedup is < 3×, the defence reports that honestly and scopes the gap as future work. The defence is about honest accounting, not marketing.
