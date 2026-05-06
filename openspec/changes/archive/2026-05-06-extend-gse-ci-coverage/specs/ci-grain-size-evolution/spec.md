## ADDED Requirements

### Requirement: Analytical L2 benchmark covers the coupled dislocation+diffusion regime

The CI test suite SHALL include a 5-strain-rate sweep test (`RheologyCreep.GrainSizeSweepCoupled`) that runs the calcite paleowattmeter scenario with **both** dislocation creep (`pwlv = 15`, Renner et al. 2002) and diffusion creep (`linv = 15`, Calcite Herwegh 2003) active simultaneously, in addition to the wattmeter (`gs = 10`). The test SHALL L2-compare the resulting grid grain-size field at each of five strain rates against an **analytical coupled steady state** computed inside the test by solving the 1D fixed-point on $\tau_{II}$ that satisfies $\dot\varepsilon_{II}^{tot} = \dot\varepsilon_{II}^{pwl}(\tau_{II}) + \dot\varepsilon_{II}^{lin}(\tau_{II}, d_{ss})$ with $d_{ss} = (B_g\dot\varepsilon_{II}^{pwl}\tau_{II}p/A_g)^{-1/(p+1)}$. The test's solver SHALL use a different iteration scheme than MDOODZ's (e.g. secant on $\tau_{II}$ rather than bisection-then-Newton on $\eta_{ve}$) so that a self-consistent-but-wrong MDOODZ result is exposed by L2 disagreement rather than masked.

#### Scenario: Coupled sweep converges to analytical d_ss across 4 decades

- **WHEN** five single-step (`Nt = 1`) homogeneous pure-shear simulations run with `pwlv = 15`, `linv = 15`, `gs = 10`, `elastic = 0`, `thermal = 0`, fixed temperature, single phase (`user1 = 0`), and `bkg_strain_rate ∈ {10⁻¹⁶, 10⁻¹⁵, 10⁻¹⁴, 10⁻¹³, 10⁻¹²}` s⁻¹
- **THEN** for each strain rate, the relative L2 error between `mesh.d_n` (from `Output00001.gzip.h5`) and a constant-vector at the analytical $d_{ss}$ (solved from the 1D coupled fixed-point at that Eii_total) SHALL be less than `5e-3`
- **AND** the mean of `mesh.d_n` SHALL be within 1% of the analytical $d_{ss}$
- **AND** `min(mesh.d_n) > 0` and finite
- **AND** `max(mesh.d_n) / min(mesh.d_n) − 1 < 1e-3` (homogeneity)

#### Scenario: Test's analytical solver converges

- **WHEN** the test's coupled fixed-point solver runs at any of the five strain rates
- **THEN** it SHALL converge to a $\tau_{II}$ satisfying $|Eii_{total} - Eii_{pwl}(\tau_{II}) - Eii_{lin}(\tau_{II}, d_{ss})| / Eii_{total} < 10^{-12}$ within 50 iterations
- **AND** the resulting $d_{ss}$ SHALL be larger than the dislocation-only $d_{ss}$ at the same Eii_total (because diffusion takes a fraction of the strain rate, leaving less for dislocation, which the wattmeter rewards with a larger steady-state grain size)

#### Scenario: Single base fixture file is used for the coupled sweep

- **WHEN** the test source tree is inspected
- **THEN** exactly one fixture SHALL exist for the coupled sweep, [TESTS/RheologyCreep/GrainSizeSweepCoupledBase.txt](TESTS/RheologyCreep/GrainSizeSweepCoupledBase.txt), with `pwlv = 15`, `linv = 15`, `gs = 10`
- **AND** the test SHALL run this base file 5 times with `bkg_strain_rate` and `writer_subfolder` injected per iteration via MDLIB's `MutateInput` callback hook
- **AND** the dislocation-only `RheologyCreep.GrainSizeSweep` test SHALL similarly share its base ([TESTS/RheologyCreep/GrainSizeSteadyState.txt](TESTS/RheologyCreep/GrainSizeSteadyState.txt)) — no per-strain-rate fixture files SHALL exist

#### Scenario: Coupled-sweep test failure points to the coupled-Jacobian path

- **WHEN** `RheologyCreep.GrainSizeSweepCoupled` fails while `RheologyCreep.GrainSizeSweep` (dislocation only) still passes
- **THEN** the failure SHALL be diagnostic of a regression in the coupled-mechanism Jacobian path of `LocalIterationViscoElasticGrainSize` (specifically the `dfdeta += params.m_lin * Eii_lin * dddeta / d_ve` term and the diffusion strain-rate evaluation), not in the dislocation-only path

### Requirement: 2D smoke test gates the integrated PinchSwellGSE code path

The CI test suite SHALL include a 2D smoke test (`RheologyCreep.PinchSwellGSESmoke`) that runs a downsized version of the Schmalholz & Duretz (2017) pinch-and-swell scenario through MDOODZ's full integrated pipeline (Stokes solve, advection, P2G interpolation, harmonic averaging, particle reseeding) and asserts that grain size remains physically valid and develops the expected spatial heterogeneity.

#### Scenario: Downsized PinchSwellGSE completes with healthy d_n field

- **WHEN** [TESTS/RheologyCreep/PinchSwellGSESmoke.txt](TESTS/RheologyCreep/PinchSwellGSESmoke.txt) — a copy of [SETS/PinchSwellGSE.txt](SETS/PinchSwellGSE.txt) with `Nx = Nz = 51` and `Nt = 10` — is run by `RheologyCreep.PinchSwellGSESmoke`
- **THEN** the simulation SHALL exit with status 0 and `Output00010.gzip.h5` SHALL exist
- **AND** all values in `mesh.d_n` SHALL be strictly positive and finite (catches NaN propagation in the integrated code path)

#### Scenario: Heterogeneity emerges in the layer

- **WHEN** the same simulation output is read
- **THEN** `max(mesh.d_n) / min(mesh.d_n)` SHALL exceed 5 (catches "scenario silently produced uniform field" regressions — the wattmeter MUST drive heterogeneity given the cosine-perturbed pinch-swell stress field)

#### Scenario: Grain reduction occurs in the layer phase

- **WHEN** the same simulation output is read
- **THEN** `min(mesh.d_n)` SHALL be less than the layer phase's `gs_ref` (1×10⁻³ m) — confirming the wattmeter actively reduced grain size where stress is highest, not bypassed silently

#### Scenario: Smoke test runtime is bounded

- **WHEN** `RheologyCreep.PinchSwellGSESmoke` runs on a typical CI machine
- **THEN** total wall time SHALL be under 30 seconds (51×51 × 10 steps; expected ~5 s)

### Requirement: Coverage extensions documented in `TESTS/AnalyticalSolutions.md`

The new coupled sweep and 2D smoke tests SHALL be documented in [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) as a "Coverage extensions" subsection appended to §8 (rather than as a new top-level section), since they validate the same paleowattmeter physics on different code paths.

#### Scenario: §8 contains a Coverage extensions subsection

- **WHEN** [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** §8 SHALL contain a subsection titled "Coverage extensions" (or equivalent) describing both `RheologyCreep.GrainSizeSweepCoupled` and `RheologyCreep.PinchSwellGSESmoke`, what each catches beyond the headline sweep, the parameter file paths, the analytical or invariant references, and the measured numerical results

#### Scenario: Summary table includes the new tests

- **WHEN** the summary table at the bottom of [TESTS/AnalyticalSolutions.md](TESTS/AnalyticalSolutions.md) is read
- **THEN** it SHALL contain at least two new rows: one for the coupled sweep (analytical reference: paleowattmeter steady state), and one for the 2D smoke (reference: invariants on `d_n` finite/positive/heterogeneous)
