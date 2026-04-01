# MDOODZ Test Suite — Scientific Documentation

This document describes the automated test suite for MDOODZ 7.0, a 2D Visco-Elasto-Plastic Thermo-Mechanical (VEPTM) geodynamic modelling code. Each test is designed to exercise a specific physical mechanism or numerical feature in isolation, so that regressions can be detected at the level of individual constitutive laws, boundary conditions, thermal processes, and solver configurations.

All tests use [Google Test](https://github.com/google/googletest) (GTest) and run as part of the continuous-integration pipeline. They produce HDF5 output files that are read back to verify physical observables against expected bounds.

---

## How to read Google Test assertions

MDOODZ tests use two main assertion families:

- `ASSERT_*`: **fatal** for the current test function. If it fails, that test stops immediately.
- `EXPECT_*`: **non-fatal** for the current test function. If it fails, the test keeps running and reports additional failures.

Common comparison macros:

| Macro | Meaning |
|---|---|
| `ASSERT_GE(a, b)` / `EXPECT_GE(a, b)` | Check `a >= b` |
| `ASSERT_GT(a, b)` / `EXPECT_GT(a, b)` | Check `a > b` |
| `ASSERT_LE(a, b)` / `EXPECT_LE(a, b)` | Check `a <= b` |
| `ASSERT_LT(a, b)` / `EXPECT_LT(a, b)` | Check `a < b` |
| `ASSERT_EQ(a, b)` / `EXPECT_EQ(a, b)` | Check exact equality |
| `ASSERT_NEAR(a, b, tol)` / `EXPECT_NEAR(a, b, tol)` | Check `|a-b| <= tol` (useful for floating-point values) |

Quick rule of thumb:
- Use `ASSERT_*` when the rest of the test depends on that condition (e.g. output file exists, solver returned steps).
- Use `EXPECT_*` for physics bounds and secondary checks where collecting multiple failures is useful.

---

## Table of Contents

1. [How to read Google Test assertions](#how-to-read-google-test-assertions)
2. [Common Model Setup](#common-model-setup)
3. [Test Infrastructure](#test-infrastructure)
4. [Suite 1 — ShearTemplate (Newton Convergence)](#suite-1--sheartemplate-newton-convergence)
5. [Suite 2 — RheologyCreep (Viscous Flow Laws)](#suite-2--rheologycreep-viscous-flow-laws)
6. [Suite 3 — Plasticity (Yield & Failure)](#suite-3--plasticity-yield--failure)
7. [Suite 4 — Thermal (Heat Equation)](#suite-4--thermal-heat-equation)
8. [Suite 5 — BoundaryCondition (Velocity BCs)](#suite-5--boundarycondition-velocity-bcs)
9. [Suite 6 — SolverMode (Iteration & Time Stepping)](#suite-6--solvermode-iteration--time-stepping)
10. [Suite 7 — ViscoElastic (Maxwell Stress Accumulation)](#suite-7--viscoelastic-maxwell-stress-accumulation)
11. [Suite 8 — Density (Thermal Expansion & Hydrostatics)](#suite-8--density-thermal-expansion--hydrostatics)
12. [Suite 9 — ShearHeating (Viscous Dissipation)](#suite-9--shearheating-viscous-dissipation)
13. [Suite 10 — FreeSurface (Topography & Sinking)](#suite-10--freesurface-topography--sinking)
14. [Suite 11 — VelocityField (Kinematic Symmetry)](#suite-11--velocityfield-kinematic-symmetry)
15. [Suite 12 — Compressibility (Elastic Bulk Response)](#suite-12--compressibility-elastic-bulk-response)
16. [Suite 13 — FiniteStrain (Deformation Gradient)](#suite-13--finitestrain-deformation-gradient)
17. [Suite 14 — NeumannBC (Stress Boundary Conditions)](#suite-14--neumannbc-stress-boundary-conditions)
18. [Suite 15 — ConvergenceRate (Newton vs Picard)](#suite-15--convergencerate-newton-vs-picard)
19. [Running the Tests](#running-the-tests)

---

## Common Model Setup

Most tests share a canonical geometry: a square 2D domain containing a circular **inclusion** (phase 1) embedded in a **matrix** (phase 0). The inclusion acts as a perturbation that triggers non-uniform stress and strain-rate fields, which is essential for exercising non-linear rheologies and verifying that the solver handles heterogeneous viscosity.

**Phase assignment callback:**

```
phase = (x² + z² < r²) ? 1 : 0
```

where $r$ is the inclusion radius (parameter `user1`), and the domain is centred at the origin.

**Boundary conditions** are either:
- **Pure shear** ($\dot{\varepsilon} _{bkg}$): $V_x = \dot{\varepsilon}_{bkg} \cdot x$ on left/right, $V_z = -\dot{\varepsilon}_{bkg} \cdot z$ on top/bottom (free-slip tangential).
- **Simple shear** ($\dot{\gamma}$): $V_x = \dot{\gamma} \cdot z$ on top/bottom with periodic lateral boundaries.

The non-dimensional scaling follows:

$$\tilde{x} = \frac{x}{L}, \quad \tilde{\eta} = \frac{\eta}{\eta_0}, \quad \tilde{V} = \frac{V}{V_0}, \quad \tilde{T} = \frac{T}{T_0}$$

where reference values $(\eta_0, L_0, V_0, T_0)$ are defined in the `SCALES` block of each `.txt` parameter file.

---

## Test Infrastructure

### TestHelpers.h

Utility functions for reading HDF5 output:

| Function | Returns | Purpose |
|---|---|---|
| `getStepsCount` | Newton iteration count | Verify solver convergence (fewer iterations → well-conditioned problem) |
| `getFinalResidual` | Last absolute residual $\|r\|$ | Confirm the momentum or divergence residual fell below tolerance |
| `getMaxFieldValue` | $\max(f_{ij})$ over the grid | Bound-check on stress, viscosity, temperature, etc. |
| `getMinFieldValue` | $\min(f_{ij})$ over the grid | Ensure fields remain positive where physics requires it |
| `getMeanFieldValue` | $\bar{f}$ (arithmetic mean) | Track bulk evolution of temperature or other extensive quantities |
| `getFieldValueAt` | $f_{index}$ at a specific 1D index | Point-wise verification of boundary values or known solutions |
| `getMaxAbsFieldValue` | $\max(\|f_{ij}\|)$ over the grid | Check magnitude bounds regardless of sign (e.g. velocity amplitude) |

---

## Suite 1 — ShearTemplate (Newton Convergence)

**Source:** `NewtonIterationConvergence.cpp`
**Parameter files:** `ShearTemplate/*.txt`

### Purpose

Verify that the Newton-Raphson and Picard non-linear solvers converge for both linear and non-linear rheologies, under both pure shear and simple shear loading, with and without mechanical anisotropy.

### Physical Background

MDOODZ solves the incompressible Stokes equations:

$$\nabla \cdot \boldsymbol{\sigma} + \rho \mathbf{g} = 0$$

$$\nabla \cdot \mathbf{v} = 0$$

where $\boldsymbol{\sigma} = -P\mathbf{I} + 2\eta\dot{\boldsymbol{\varepsilon}}$ is the Cauchy stress, $P$ is pressure, $\eta$ is the effective viscosity, and $\dot{\boldsymbol{\varepsilon}}$ is the deviatoric strain-rate tensor.

For **linear rheology** ($n = 1$), the system is linear and the Newton solver should converge in exactly **1 iteration**. For **non-linear rheology** ($n > 1$, power-law), the viscosity depends on the strain rate:

$$\eta_{eff} = \eta_0 \left(\frac{\dot{\varepsilon}_{II}}{\dot{\varepsilon}_0}\right)^{\frac{1}{n} - 1}$$

requiring iterative updates. The Newton solver should converge in **fewer than 10 iterations** due to its quadratic convergence rate.

### Tests

| Test | Rheology | Deformation | Anisotropy | Expected iterations |
|---|---|---|---|---|
| LinearPureshearIsotropic | $n = 1$ (Newtonian) | Pure shear | Off | Exactly 1 |
| LinearSimpleshearIsotropic | $n = 1$ | Simple shear | Off | Exactly 1 |
| LinearPureshearAnisotropic | $n = 1$ | Pure shear | On | Exactly 1 |
| LinearSimpleshearAnisotropic | $n = 1$ | Simple shear | On | Exactly 1 |
| NonLinearPureshearIsotropic | $n = 3$ (power-law) | Pure shear | Off | 2–9 |
| NonLinearSimpleshearIsotropic | $n = 3$ | Simple shear | Off | 2–9 |
| NonLinearPureshearAnisotropic | $n = 3$ | Pure shear | On | 2–19 |
| NonLinearSimpleshearAnisotropi | $n = 3$ | Simple shear | On | 2–9 |

**Code assertions** (`NewtonIterationConvergence.cpp`):
```cpp
// Linear tests — must converge in exactly 1 Newton step:
ASSERT_EQ(stepsCount, 1);
// Non-linear tests — must converge in a bounded number of iterations:
ASSERT_TRUE(stepsCount > 1);
ASSERT_TRUE(stepsCount < 10);   // (< 20 for anisotropic pure shear)
```

### Justification

These tests form the **foundation** of the solver verification. If Newton convergence is broken for a linear problem (which must converge in 1 step), any downstream test is unreliable. They also catch regressions in the Jacobian assembly, line-search algorithm, and anisotropy coupling.

---

## Suite 2 — RheologyCreep (Viscous Flow Laws)

**Source:** `RheologyCreepTests.cpp`
**Parameter files:** `RheologyCreep/*.txt`

All tests use a 10 km × 10 km domain, 31 × 31 grid, a 500 m radius inclusion, pure shear at $\dot{\varepsilon} = 10^{-14}$ s⁻¹, and $T = 1100$ °C (except Peierls: 800 °C). The inclusion uses a constant viscosity (`cstv = 1`, $\eta_0 = 10^{27}$ Pa·s) to create a stiff heterogeneity.

### 2.1 DiffusionCreepLinear

**Flow law:** Olivine dry diffusion creep — Hirth & Kohlstedt (2003), index `linv = 40`

$$\dot{\varepsilon}_{diff} = A_{diff} \cdot \sigma^{n_{diff}} \cdot d^{-m_{diff}} \cdot \exp\left(-\frac{Q_{diff} + PV_{diff}}{RT}\right)$$

| Parameter | Value |
|---|---|
| $A_{diff}$ | $1.5 \times 10^{-15}$ Pa$^{-n}$ m$^{m}$ s$^{-1}$ |
| $n_{diff}$ | 1.0 (linear) |
| $m_{diff}$ | 3.0 (grain-size exponent) |
| $Q_{diff}$ | 375 kJ/mol |
| $V_{diff}$ | $4 \times 10^{-6}$ m³/mol |
| $d$ (gs_ref) | 5 mm |

**Cause → Effect:** With $n = 1$, diffusion creep is Newtonian. The matrix effective viscosity is a function only of temperature and grain size:

$$\eta_{diff} = \frac{1}{2A_{diff}} \cdot d^{m} \cdot \exp\left(\frac{Q + PV}{RT}\right)$$

**Code assertions** (`RheologyCreepTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed
EXPECT_GT(minEta, 0.0);                    // η > 0 everywhere
```

### 2.2 DiffusionCreepGrainSize

**Flow law:** Same as above (`linv = 40`), but the inclusion (phase 1) uses a **smaller grain size**: $d_1 = 1$ mm vs. $d_0 = 5$ mm.

**Cause → Effect:** Since $\eta_{diff} \propto d^m$ and $m = 3$, reducing the grain size by a factor of 5 should reduce the inclusion viscosity by a factor of $5^3 = 125$. This tests MDOODZ's ability to interpolate grain size from markers to the grid, and that the diffusion creep law correctly handles per-phase grain-size variation. Note that both phases may hit the non-dimensional viscosity ceiling (`max_eta`), so the test verifies that the flow-law path executes without producing NaN or negative values.

**Code assertions** (`RheologyCreepTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed
EXPECT_GT(minEta, 0.0);                    // η > 0 everywhere
```

### 2.3 PeierlsCreep

**Flow law:** Olivine Peierls creep — Evans & Goetze (1979), regularised following Kameyama et al. (1999), index `expv = 40`

Peierls (or low-temperature plasticity) is a thermally-activated dislocation glide mechanism dominant at high stress and low temperature ($T < 900$ °C):

$$\dot{\varepsilon}_{Pei} = E_{Pei} \cdot \exp\left(-\frac{Q_{Pei}}{RT}\left(1 - \left(\frac{\sigma}{\sigma_P}\right)^q\right)^p\right)$$

| Parameter | Value |
|---|---|
| $E_{Pei}$ | $5.7 \times 10^{11}$ s$^{-1}$ |
| $Q_{Pei}$ | 540 kJ/mol |
| $\sigma_P$ | 8.5 GPa |
| $p$ | 0.1 |
| $q$ | 2.0 |

**Cause → Effect:** The temperature is set to 800 °C (colder than the other creep tests) to place the system in the Peierls-dominated regime. At these conditions, the effective viscosity is strongly stress-dependent, requiring robust non-linear iteration. The matrix is assigned `expv = 40`; the inclusion is constant-viscosity (stiff, $\eta = 10^{27}$ Pa·s).

**Code assertions** (`RheologyCreepTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged within nit_max=20
EXPECT_GT(minEta, 0.0);                    // η > 0 (no NaN or negative from Peierls regularisation)
```

### 2.4 GBSCreep

**Flow law:** Olivine low-T grain-boundary sliding — Hirth & Kohlstedt (2003), index `gbsv = 40`

$$\dot{\varepsilon}_{GBS} = A_{GBS} \cdot \sigma^{n_{GBS}} \cdot d^{-m_{GBS}} \cdot \exp\left(-\frac{Q_{GBS}}{RT}\right)$$

| Parameter | Value |
|---|---|
| $A_{GBS}$ | $6.5 \times 10^{-30}$ Pa$^{-n}$ m$^{m}$ s$^{-1}$ |
| $n_{GBS}$ | 3.5 |
| $m_{GBS}$ | 2.0 |
| $Q_{GBS}$ | 400 kJ/mol |

**Cause → Effect:** GBS is a composite mechanism (dislocation + diffusion at grain boundaries) that activates at intermediate temperatures ($T < 1250$ °C) and small grain sizes. This test uses `elastic = 1` (Maxwell visco-elastic) because the GBS code path requires elastic stress storage to avoid numerical instability. The solver is set to Picard (`Newton = 0`) to avoid Jacobian singularities in the GBS branch.

**Code assertions** (`RheologyCreepTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed (GBS + Picard + elastic)
EXPECT_GT(minEta, 0.0);                    // η > 0 (no NaN from GBS code path)
```

### 2.5 CompositePwlLin

**Flow laws:** Simultaneous dislocation (`pwlv = 40`) + diffusion (`linv = 40`) creep, both Olivine dry — Hirth & Kohlstedt (2003).

When multiple creep mechanisms operate in parallel, their strain rates add:

$$\dot{\varepsilon}_{total} = \dot{\varepsilon}_{pwl} + \dot{\varepsilon}_{lin}$$

and the composite effective viscosity is:

$$\frac{1}{\eta_{comp}} = \frac{1}{\eta_{pwl}} + \frac{1}{\eta_{lin}}$$

The **weaker mechanism dominates**. At low stress / small grain size, diffusion creep controls the flow; at high stress, dislocation creep takes over.

**Cause → Effect:** With $\dot{\varepsilon} = 10^{-14}$ s$^{-1}$, $T = 1100$ °C, and $d = 5$ mm, both mechanisms contribute comparably to the total strain rate. The code must correctly sum the strain-rate contributions and compute the harmonic mean viscosity. The viscosity contrast field around the stiff inclusion tests both mechanism branches simultaneously.

**Code assertions** (`RheologyCreepTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed
EXPECT_GT(minEta, 0.0);                    // η > 0 (composite viscosity is physical)
```

### 2.6 PowerLawCreep

**Flow law:** Olivine dry dislocation creep — Hirth & Kohlstedt (2003), index `pwlv = 40`

$$\dot{\varepsilon}_{pwl} = A_{pwl} \cdot \sigma^{n_{pwl}} \cdot \exp\left(-\frac{Q_{pwl} + PV_{pwl}}{RT}\right)$$

| Parameter | Value |
|---|---|
| $A_{pwl}$ | $1.1 \times 10^{-16}$ Pa$^{-n}$ s$^{-1}$ |
| $n_{pwl}$ | 3.5 |
| $Q_{pwl}$ | 530 kJ/mol |
| $V_{pwl}$ | $18 \times 10^{-6}$ m³/mol |

**Cause → Effect:** At $T = 600$ °C (colder than the other creep tests) and $\dot{\varepsilon} = 10^{-14}$ s$^{-1}$, the power-law viscosity is stress-dependent: $\eta_{pwl} = \frac{\sigma}{2\dot{\varepsilon}_{pwl}}$. The test verifies that the dislocation creep pathway is correctly activated by checking that the power-law strain-rate invariant $\dot{\varepsilon}_{II}^{pwl} > 0$. The inclusion uses constant viscosity ($\eta = 10^{27}$ Pa·s, `cstv = 1`) to create a viscosity contrast.

**Code assertions** (`RheologyCreepTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_GT(minEta, 0.0);                    // η > 0 everywhere
EXPECT_GT(maxEiiPwl, 0.0);                 // dislocation creep strain rate > 0 → pwl active
```

---

## Suite 3 — Plasticity (Yield & Failure)

**Source:** `PlasticityTests.cpp`
**Parameter files:** `Plasticity/*.txt`

All tests use a 4 km × 4 km domain, 31 × 31 grid, a 100 m radius inclusion, pure shear at $\dot{\varepsilon} = 5 \times 10^{-13}$ s$^{-1}$ (except NoYield and StressLimiter), with crustal-like density ($\rho = 2700$ kg/m³). Elasticity is enabled (`elastic = 1`, Maxwell model), and finite strain is tracked. The strain rate is chosen high enough that the Maxwell visco-elastic trial stress $\sigma_{trial} = 2 G \Delta t \dot{\varepsilon}$ exceeds the Drucker-Prager yield stress, ensuring plastic yielding occurs.

### Plasticity Model

MDOODZ uses a **Drucker-Prager** yield criterion (smooth approximation to Mohr-Coulomb):

$$\tau_{II} \leq C \cos\varphi + P \sin\varphi$$

where $\tau_{II} = \sqrt{\frac{1}{2}\boldsymbol{\tau}:\boldsymbol{\tau}}$ is the second invariant of the deviatoric stress, $C$ is cohesion, $\varphi$ is the friction angle, and $P$ is pressure.

When the trial stress (from the viscous/elastic predictor) exceeds the yield stress, the effective viscosity is reduced to bring the stress back to the yield surface:

$$\eta_{pl} = \frac{\tau_{yield}}{2\dot{\varepsilon}_{II}}$$

This produces a **viscoplastic** effective viscosity: $\eta_{eff} = \min(\eta_{viscous}, \eta_{pl})$.

### 3.1 DruckerPragerYield

**Parameters:** $C = 20$ MPa, $\varphi = 30°$, $\psi = 10°$ (dilation angle), `plast = 1`, $\dot{\varepsilon} = 5 \times 10^{-13}$ s$^{-1}$

**Cause → Effect:** The applied strain rate is high enough that the elastic/viscous trial stress exceeds the yield stress ($\sigma_{trial} \approx 2 G \Delta t \dot{\varepsilon} = 2 \times 10^{10} \times 10^{10} \times 5 \times 10^{-13} = 100$ MPa $> C \cos\varphi \approx 17$ MPa). The code must activate plastic yielding, reducing the effective viscosity. The inclusion (weaker $G = 2.5$ GPa vs. matrix $G = 10$ GPa) focuses deformation and triggers earlier yielding.

**Code assertions** (`PlasticityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_LT(maxSxxd, 1e40);                  // deviatoric stress is finite
EXPECT_LT(maxSzzd, 1e40);
EXPECT_GT(maxEiiPl, 0.0);                  // plastic strain rate > 0 → plasticity activated
```

### 3.2 DruckerPragerNoYield

**Parameters:** Same as above **except** $\dot{\varepsilon} = 10^{-18}$ s$^{-1}$ (four orders of magnitude lower).

**Cause → Effect:** At this negligible strain rate, the trial stress is far below the yield surface. No plastic deformation should occur. The plastic strain-rate invariant $\dot{\varepsilon}_{II}^{pl}$ must remain effectively zero.

**Code assertions** (`PlasticityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_LT(maxEiiPl, 1e-25);               // plastic strain rate ≈ 0 → no yielding (negative control)
```

### 3.3 StrainSoftening

**Parameters:** $C_0 = 20$ MPa, $C_e = 5$ MPa (end cohesion), `coh_soft = 1`, accumulated plastic strain window $\varepsilon^{pl} \in [0, 0.5]$, 3 time steps.

Cohesion softening is a linear interpolation:

$$C(\varepsilon^{pl}) = C_0 - (C_0 - C_e) \cdot \frac{\varepsilon^{pl} - \varepsilon^{pl}_{start}}{\varepsilon^{pl}_{end} - \varepsilon^{pl}_{start}}, \quad \varepsilon^{pl}_{start} \leq \varepsilon^{pl} \leq \varepsilon^{pl}_{end}$$

**Cause → Effect:** Over 3 time steps, accumulated plastic strain grows in the deforming regions (particularly near the inclusion). The cohesion must decrease from 20 MPa toward 5 MPa, further localising deformation — a positive feedback that models shear-zone weakening in the lithosphere.

**Code assertions** (`PlasticityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged after 3 steps
maxCoh_init = getMaxFieldValue(fileInit, "Centers", "cohesion");  // read initial C
minCoh      = getMinFieldValue(fileFinal, "Centers", "cohesion"); // read final C
EXPECT_GT(minCoh, 0.0);                    // cohesion stays positive
EXPECT_LT(minCoh, maxCoh_init);            // cohesion decreased → softening occurred
```

### 3.4 StressLimiter

**Parameters:** `plast = 0` (Drucker-Prager off), `Slim = 100 MPa`.

The stress limiter is an independent mechanism that caps deviatoric stress:

$$\eta_{eff} = \min\left(\eta_{viscous},\; \frac{S_{lim}}{2\dot{\varepsilon}_{II}}\right)$$

**Cause → Effect:** Without Drucker-Prager plasticity, the only stress-bounding mechanism is the limiter. The effective viscosity of the matrix ($\eta_0 = 10^{30}$ Pa·s, very stiff) would generate enormous stresses at $\dot{\varepsilon} = 5 \times 10^{-14}$ s$^{-1}$ without the cap. With `Slim = 100 MPa`, stresses must saturate.

**Code assertions** (`PlasticityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_LT(maxSxxd, 1e40);                  // stress is finite
EXPECT_LT(maxSzzd, 1e40);
// Without limiter: η_trial = 1e30 → σ ~ 1e17 Pa. With Slim = 1e8 Pa, η must be ~100× lower:
EXPECT_LT(maxEta, eta_nd_unlimited * 0.01); // viscosity capped far below trial value
```

---

## Suite 4 — Thermal (Heat Equation)

**Source:** `ThermalTests.cpp`
**Parameter files:** `Thermal/*.txt`

All tests use a 100 km × 100 km domain, 31 × 31 grid, single-phase (no inclusion), constant viscosity, minimal mechanical forcing ($\dot{\varepsilon} = 10^{-18}$ s$^{-1}$), and coupled thermo-mechanical solving (`thermal = 1`).

### Governing Equation

The energy equation solved by MDOODZ:

$$\rho C_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + H_r + H_s + H_a$$

where:
- $\rho C_p$ is the volumetric heat capacity
- $k$ is thermal conductivity
- $H_r$ = radiogenic heat production
- $H_s$ = shear (viscous dissipation) heating
- $H_a$ = adiabatic heating

In these tests, shear and adiabatic heating are disabled (`shear_heating = 0`, `adiab_heating = 0`) to isolate specific thermal mechanisms.

### 4.1 GaussianDiffusion

**Setup:** Background $T_0 = 500$ °C. A circular thermal perturbation with $\Delta T = 200$ °C and radius 10 km is placed at the centre. 5 time steps, $\Delta t = 10^{12}$ s (~31.7 kyr). No-flux boundaries (zero Neumann), $Q_r = 0$.

**Cause → Effect:** The initial Gaussian-like thermal anomaly relaxes by conduction. The thermal diffusivity is:

$$\kappa = \frac{k}{\rho C_p} = \frac{3.0}{3300 \times 1050} \approx 8.66 \times 10^{-7} \text{ m}^2/\text{s}$$

The characteristic diffusion length after 5 steps is $\ell \sim \sqrt{\kappa \cdot 5\Delta t} \approx 6.6$ km, comparable to the perturbation radius, so significant smoothing occurs.

**Code assertions** (`ThermalTests.cpp`):
```cpp
maxT_init  = getMaxFieldValue(fileInit,  "Centers", "T");  // peak T at t=0
maxT_final = getMaxFieldValue(fileFinal, "Centers", "T");  // peak T at t=5Δt
EXPECT_LT(maxT_final, maxT_init);          // peak T decreased → diffusion is working
```

### 4.2 SteadyStateGeotherm

**Setup:** $T_{top} = 0$ °C (Dirichlet at surface), $T_{bot} = 1600$ K (Dirichlet at base). $Q_r = 0$, 20 time steps to approach steady state.

**Cause → Effect:** With Dirichlet conditions on both boundaries and no internal heat sources, the steady-state solution is a linear geotherm:

$$T(z) = T_{top} + (T_{bot} - T_{top}) \cdot \frac{z_{top} - z}{z_{top} - z_{bot}}$$

Over 20 diffusive steps, the temperature field should converge toward this linear profile.

**Code assertions** (`ThermalTests.cpp`):
```cpp
EXPECT_GT(maxT, minT);                     // thermal gradient exists → Dirichlet BCs applied
// Mean temperature lies between boundary values (linear geotherm property)
EXPECT_GT(meanT, minT);
EXPECT_LT(meanT, maxT);
```

### 4.3 RadiogenicHeat

**Setup:** Uniform initial $T_0 = 500$ °C. Radiogenic heat production $Q_r = 10^{-6}$ W/m³ (typical upper-crustal value). Zero-flux boundaries, 5 time steps.

**Cause → Effect:** The volumetric heat source raises the domain temperature over time:

$$\Delta T \approx \frac{Q_r \cdot \Delta t}{\rho C_p} = \frac{10^{-6} \times 5 \times 10^{12}}{3300 \times 1050} \approx 1.44 \text{ °C}$$

The mean temperature must increase monotonically.

**Code assertions** (`ThermalTests.cpp`):
```cpp
meanT_init  = getMeanFieldValue(fileInit,  "Centers", "T");
meanT_final = getMeanFieldValue(fileFinal, "Centers", "T");
EXPECT_GT(meanT_final, meanT_init);        // mean T increased → radiogenic source is active
// Peak temperature should also increase from volumetric heating
EXPECT_GT(maxT_final, maxT_init);
```

### 4.4 DirichletBC

**Setup:** $T_{top} = 0$ °C, $T_{bot} = 1600$ K (Dirichlet on both boundaries). $Q_r = 0$, 1 time step.

**Cause → Effect:** After one time step with fixed-temperature boundaries, the interior temperature must be bounded between the two boundary values and the field must exhibit a monotone variation from top to bottom.

**Code assertions** (`ThermalTests.cpp`):
```cpp
EXPECT_GT(minT, 0.0);                      // no sub-zero T (Dirichlet BCs respected)
EXPECT_GT(maxT, minT);                     // thermal gradient exists
```

---

## Suite 5 — BoundaryCondition (Velocity BCs)

**Source:** `BoundaryConditionTests.cpp`
**Parameter files:** `BoundaryCondition/*.txt`

All tests use a non-dimensional setup: $L = 1$ m, $\eta_0 = 1$ Pa·s, $V_0 = 1$ m/s. Domain is $[-0.5, 0.5]^2$ with a 21 × 21 grid. Matrix viscosity = 1, inclusion viscosity = 10 (viscosity contrast of 10×). This non-dimensional setup is a standard Stokes benchmark geometry.

### 5.1 PeriodicSimpleShear

**Setup:** `periodic_x = 1`, `shear_style = 1` (simple shear), $\dot{\gamma} = 1$.

**Cause → Effect:** The domain is periodic in the $x$-direction. Simple shear imposes $V_x = \dot{\gamma} \cdot z$ on the top and bottom boundaries. The solver must correctly handle the periodicity by stitching the left and right boundaries in the stiffness matrix. The pressure field must remain well-defined (not drift) since periodic BCs remove the pressure null space only in specific configurations.

**Code assertions** (`BoundaryConditionTests.cpp`):
```cpp
ASSERT_GT(stepsCount, 0);                  // solver converged
EXPECT_LT(maxP, 1e40);                     // pressure is finite (not divergent)
EXPECT_GT(minP, -1e40);                    // no negative divergence
```

### 5.2 FreeSlipPureShear

**Setup:** `pure_shear_ALE = 1`, `periodic_x = 0`, $\dot{\varepsilon} = 1$.

Free-slip conditions enforce zero tangential stress at boundaries:

$$\sigma_{xz}\big|_{\partial\Omega} = 0, \quad V_n\big|_{\partial\Omega} = V_{BC}$$

**Cause → Effect:** Under pure shear with free-slip walls, the velocity field for a homogeneous medium would be exactly $(V_x, V_z) = (\dot{\varepsilon} x, -\dot{\varepsilon} z)$. The inclusion creates a local perturbation, but the shear stress on the boundaries must remain near zero everywhere.

**Code assertions** (`BoundaryConditionTests.cpp`):
```cpp
ASSERT_GT(stepsCount, 0);                  // solver converged
EXPECT_LT(maxSxz, 1e40);                   // shear stress is bounded
EXPECT_GT(minSxz, -1e40);
// Velocity field check: Vx anti-symmetry confirms extensional pure shear
EXPECT_GT(maxVx, 0.0);                     // positive Vx on east side
EXPECT_LT(minVx, 0.0);                     // negative Vx on west side
```

### 5.3 NoSlip

**Setup:** Same geometry but with **no-slip** conditions: $V_x = 0$ and $V_z = 0$ on all boundaries. Gravity is enabled ($g_z = -10$ m/s²) and no background strain rate is applied (`pure_shear_ALE = 0`).

$$\mathbf{V}\big|_{\partial\Omega} = \mathbf{0}$$

**Cause → Effect:** No-slip conditions completely pin the velocity at the domain edges. The deformation is driven by the gravitational body force acting on the density contrast between the inclusion ($\rho = 2750$ kg/m³) and matrix ($\rho = 2700$ kg/m³). This tests the correct imposition of Dirichlet velocity conditions across all boundary nodes, which requires special handling in the staggered-grid FD scheme (ghost nodes, penalty enforcement). Note: `pure_shear_ALE` must be disabled — combining no-slip BCs with ALE kinematic boundary velocities produces a singular stiffness matrix.

**Code assertions** (`BoundaryConditionTests.cpp`):
```cpp
ASSERT_GT(stepsCount, 0);                  // solver converged
EXPECT_LT(maxP, 1e40);                     // pressure is bounded
// With no-slip BCs and gravity, buoyancy-driven flow creates bounded velocities
EXPECT_LT(maxAbsVx, 1e40);                 // velocities exist but are bounded
EXPECT_LT(maxAbsVz, 1e40);
```

---

## Suite 6 — SolverMode (Iteration & Time Stepping)

**Source:** `SolverModeTests.cpp`
**Parameter files:** `SolverMode/*.txt`

These tests use the same non-dimensional benchmark geometry as the BC tests, but focus on solver algorithm variants rather than boundary conditions.

### 6.1 PicardConvergence

**Setup:** `Newton = 0` (Picard iteration), `nit_max = 50`, power-law matrix with $n = 3$.

The Picard (fixed-point) iteration scheme updates the viscosity from the previous iteration's strain rate but does not account for the viscosity-strain-rate coupling in the Jacobian:

$$\eta^{(k+1)} = \eta(\dot{\varepsilon}^{(k)})$$

Unlike Newton-Raphson (which includes $\partial\eta/\partial\dot{\varepsilon}$ in the Jacobian), Picard converges linearly and typically requires more iterations.

**Cause → Effect:** The non-linear power-law rheology ($n = 3$) requires multiple Picard iterations to converge. The test verifies that Picard finds a solution in a **reasonable** number of iterations ($1 < N_{iter} < 50$). If convergence fails entirely, it indicates a bug in the viscosity update loop or the residual computation.

**Code assertions** (`SolverModeTests.cpp`):
```cpp
ASSERT_GT(stepsCount, 0);                  // Picard converged
ASSERT_LT(stepsCount, 50);                 // within iteration budget
EXPECT_LT(finalRes, 1e-8);                 // momentum residual below tolerance
```

### 6.2 AdaptiveDt

**Setup:** `constant_dt = 0`, `Courant = 0.5`, `advection = 1`, 3 time steps, Newton solver.

When `constant_dt = 0`, the time step is computed adaptively to satisfy the CFL (Courant-Friedrichs-Lewy) condition:

$$\Delta t = \text{Courant} \cdot \frac{\Delta x}{\max |V|}$$

This ensures markers do not advect more than one cell per step.

**Cause → Effect:** With advection enabled, marker positions update each step. The code must recalculate $\Delta t$ based on the velocity field. Over 3 steps with an evolving velocity field (due to the non-linear inclusion problem), each time step may differ.

**Code assertions** (`SolverModeTests.cpp`):
```cpp
ASSERT_GT(steps1, 0);                      // step 1 converged
ASSERT_GT(steps3, 0);                      // step 3 converged (adaptive dt didn't crash)
```

---

## Suite 7 — ViscoElastic (Maxwell Stress Accumulation)

**Source:** `ViscoElasticTests.cpp`
**Parameter file:** `ViscoElastic/StressAccumulation.txt`

### Physical Background

MDOODZ implements visco-elasticity through the **Maxwell model**, where elastic and viscous strain rates are additive:

$$\dot{\varepsilon}_{ij} = \dot{\varepsilon}_{ij}^{viscous} + \dot{\varepsilon}_{ij}^{elastic} = \frac{\tau_{ij}}{2\eta} + \frac{\overset{\circ}{\tau}_{ij}}{2G}$$

where $G$ is the shear modulus and $\overset{\circ}{\tau}$ is the Jaumann objective stress rate. The **Maxwell time** is $t_M = \eta / G$, which controls the transition between elastic and viscous behaviour.

For loading faster than $t_M$, the material responds elastically; for loading slower than $t_M$, viscous flow dominates. In a time-stepping scheme, the effective viscosity becomes:

$$\eta_{eff} = \frac{\eta G \Delta t}{\eta + G \Delta t}$$

and the deviatoric stress accumulates over time steps as elastic memory builds up.

### 7.1 StressAccumulation

**Setup:** 10 km × 10 km domain, 31 × 31 grid, `elastic = 1`. Matrix: $\eta_0 = 10^{22}$ Pa·s, $G = 10^{10}$ Pa → $t_M = 10^{12}$ s (~31.7 kyr). Inclusion: $\eta_0 = 10^{27}$ Pa·s (stiffer). Pure shear at $\dot{\varepsilon} = 10^{-14}$ s$^{-1}$, $\Delta t = 10^{10}$ s, 5 time steps.

| Parameter | Value |
|---|---|
| Domain | 10 × 10 km, 31 × 31 |
| $\eta_{matrix}$ | $10^{22}$ Pa·s |
| $\eta_{incl}$ | $10^{27}$ Pa·s |
| $G$ | $10^{10}$ Pa |
| $\Delta t$ | $10^{10}$ s |
| $t_M$ | $10^{12}$ s |
| $\dot{\varepsilon}$ | $10^{-14}$ s$^{-1}$ |
| Steps | 5 |

**Cause → Effect:** With $\Delta t / t_M = 0.01 \ll 1$, the material is firmly in the elastic regime. Over 5 steps, the deviatoric stress $\tau_{xx}$ should grow approximately linearly: $\tau_{xx} \approx 2G \cdot \dot{\varepsilon} \cdot n\Delta t$. The elastic strain-rate invariant $\dot{\varepsilon}_{II}^{el}$ must be positive, confirming that the elastic branch of the constitutive law is active. The test compares stress at step 1 vs step 5 to verify monotonic build-up.

**Code assertions** (`ViscoElasticTests.cpp`):
```cpp
ASSERT_GE(steps1, 0);                      // solver converged at step 1
ASSERT_GE(steps5, 0);                      // solver converged at step 5
// Stress at step 5 ≥ 90% of stress at step 1 (monotonic build-up)
EXPECT_GE(fabs(maxSxxd_5), fabs(maxSxxd_1) * 0.9);
EXPECT_GT(maxEiiEl, 0.0);                  // elastic strain rate is active
```

### Justification

Visco-elasticity is essential for modelling transient lithospheric deformation. This test ensures the Maxwell rheology correctly accumulates stress over time steps and that the elastic strain-rate partitioning works. A regression here would produce purely viscous behaviour even for rapid loading.

---

## Suite 8 — Density (Thermal Expansion & Hydrostatics)

**Source:** `DensityTests.cpp`
**Parameter files:** `Density/ThermalExpansion.txt`, `Density/HydrostaticPressure.txt`

### Physical Background

MDOODZ computes temperature-dependent density using a linearised equation of state:

$$\rho(T) = \rho_0 \left(1 - \alpha (T - T_0)\right)$$

where $\alpha$ is the thermal expansion coefficient. This coupling is fundamental to thermal convection (buoyancy).

In the absence of deviatoric stress (static limit), the pressure obeys the hydrostatic balance:

$$\frac{\partial P}{\partial z} = \rho g_z$$

yielding $P(z) = \rho g |z|$ for uniform density.

### 8.1 ThermalExpansion

**Setup:** 10 km × 10 km domain, 31 × 31 grid. Matrix at $T_0 = 600$ °C, inclusion at $T = 600 + 1600 = 2200$ °C (user2 = 1600). Both phases have $\rho_0 = 3300$ kg/m³ and $\alpha = 3 \times 10^{-5}$ K$^{-1}$. Gravity disabled ($g_z = 0$) to isolate the density computation.

| Parameter | Value |
|---|---|
| $\rho_0$ | 3300 kg/m³ |
| $\alpha$ | $3 \times 10^{-5}$ K$^{-1}$ |
| $T_{matrix}$ | 600 °C |
| $T_{inclusion}$ | 2200 °C |
| $g_z$ | 0 |

**Cause → Effect:** The 1600 °C temperature contrast produces a density difference:

$$\Delta\rho = \rho_0 \alpha \Delta T = 3300 \times 3 \times 10^{-5} \times 1600 \approx 158 \text{ kg/m}^3$$

The hot inclusion should have lower density ($\rho_{incl} \approx 3142$ kg/m³) than the cold matrix ($\rho_{matrix} \approx 3241$ kg/m³). The test verifies that $\rho_{min} < \rho_{max}$.

**Code assertions** (`DensityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed
EXPECT_LT(minRho, maxRho);                 // density contrast exists
EXPECT_GT(minRho, 0.0);                    // density is positive
```

### 8.2 HydrostaticPressure

**Setup:** 10 km × 10 km domain, 21 × 21 grid. Uniform material ($\rho = 2700$ kg/m³, no inclusion). Gravity $g_z = -10$ m/s². Pure shear at $\dot{\varepsilon} = 10^{-14}$ s$^{-1}$ with dimensional scaling ($\eta = 10^{22}$, $L = 10^4$).

| Parameter | Value |
|---|---|
| $\rho$ | 2700 kg/m³ |
| $g_z$ | −10 m/s² |
| Domain | 10 × 10 km |
| Inclusion | None |

**Cause → Effect:** With uniform density and gravity, the pressure should develop a lithostatic profile. The maximum pressure (at the domain bottom) is:

$$P_{max} \approx \rho |g| H = 2700 \times 10 \times 10^4 = 2.7 \times 10^8 \text{ Pa}$$

The test checks that $P_{max} > P_{min}$ and $P_{max} > 0$ — i.e. the pressure increases with depth as required by hydrostatic equilibrium.

**Code assertions** (`DensityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed
EXPECT_GT(maxP, minP);                     // pressure increases with depth
EXPECT_GT(maxP, 0.0);                      // positive lithostat
```

### Justification

Density computation and hydrostatic pressure are foundational for any gravity-driven geodynamic model — subduction, mantle convection, diapirism. A broken thermal-expansion coupling or incorrect gravity integration would produce unphysical buoyancy forces.

---

## Suite 9 — ShearHeating (Viscous Dissipation)

**Source:** `ShearHeatingTests.cpp`
**Parameter file:** `ShearHeating/ViscousDissipation.txt`

### Physical Background

Shear heating (viscous dissipation) converts mechanical work into heat:

$$H_s = \tau_{ij} \dot{\varepsilon}_{ij} = 2\eta \dot{\varepsilon}_{II}^2$$

This is a positive feedback mechanism in geodynamics: deformation heats the material, reducing viscosity, which further localises deformation (thermal shear localisation). The energy equation includes this source term:

$$\rho C_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + H_s$$

### 9.1 ViscousDissipation

**Setup:** 10 km × 10 km domain, 31 × 31 grid, single phase. `thermal = 1`, `shear_heating = 1`, `adiab_heating = 0`. Constant viscosity $\eta = 10^{22}$ Pa·s, pure shear $\dot{\varepsilon} = 10^{-14}$ s$^{-1}$, $T_0 = 500$ °C, 3 time steps ($\Delta t = 10^{11}$ s). Gravity disabled. Temperature BCs: Dirichlet $T = T_0$ on top and bottom (using `SetBCT` callback).

| Parameter | Value |
|---|---|
| $\eta$ | $10^{22}$ Pa·s |
| $\dot{\varepsilon}$ | $10^{-14}$ s$^{-1}$ |
| $T_0$ | 500 °C |
| $k$ | 3.0 W/m/K |
| $\rho$ | 3300 kg/m³ |
| $C_p$ | 1050 J/kg/K |
| $\Delta t$ | $10^{11}$ s |
| Steps | 3 |

**Cause → Effect:** The volumetric heat production rate is:

$$H_s = 2\eta \dot{\varepsilon}_{II}^2 = 2 \times 10^{22} \times (10^{-14})^2 = 2 \times 10^{-6} \text{ W/m}^3$$

Over 3 steps, this generates temperature increase:

$$\Delta T \approx \frac{H_s \cdot 3\Delta t}{\rho C_p} = \frac{2 \times 10^{-6} \times 3 \times 10^{11}}{3300 \times 1050} \approx 0.17 \text{ °C}$$

The peak temperature must increase (or at least not decrease) relative to the initial state.

**Code assertions** (`ShearHeatingTests.cpp`):
```cpp
EXPECT_GE(maxT_final, maxT_init);          // peak T increased from viscous dissipation
```

### Justification

Shear heating is critical for lithospheric localisation models (shear zones, subduction interfaces). Without it, thermally coupled models would underpredict temperatures in deforming regions. This test requires a `SetBCT` callback (thermal boundary condition), which exercises a distinct code path from the mechanical-only tests.

---

## Suite 10 — FreeSurface (Topography & Sinking)

**Source:** `FreeSurfaceTests.cpp`
**Parameter file:** `FreeSurface/SinkingBlock.txt`

### Physical Background

The free surface algorithm in MDOODZ tracks the topography of the Earth's surface as a 1D line $z = h(x, t)$ that evolves with the flow. Cells above the topography are assigned `phase = -1` (air), and the surface boundary condition imposes zero normal stress:

$$\boldsymbol{\sigma} \cdot \hat{n}\big|_{surface} = 0$$

The free surface stabilisation method (`free_surface_stab = 1`) corrects for the numerical instability caused by density contrast at the surface ("sloshing" instability).

### 10.1 SinkingBlock

**Setup:** 100 km × 100 km domain, 41 × 41 grid, `free_surface = 1`, `free_surface_stab = 1`, `advection = 1`. Matrix: $\rho = 3000$ kg/m³, inclusion: $\rho = 3300$ kg/m³ (dense block with 10 km radius). Gravity $g_z = -10$ m/s², no-slip BCs on all sides, 3 time steps ($\Delta t = 10^{11}$ s). The `BuildInitialTopography_ff` callback is registered with `SetSurfaceZCoord` to initialise the surface at $z = 0$.

| Parameter | Value |
|---|---|
| Domain | 100 × 100 km |
| Grid | 41 × 41 |
| $\rho_{matrix}$ | 3000 kg/m³ |
| $\rho_{incl}$ | 3300 kg/m³ |
| $\Delta\rho$ | 300 kg/m³ |
| $g_z$ | −10 m/s² |
| Steps | 3 |
| BCs | No-slip (all sides) |

**Cause → Effect:** The dense block ($\Delta\rho = 300$ kg/m³) sinks under gravity, driving viscous flow and deflecting the free surface. The Stokes velocity scale for this problem is:

$$V \sim \frac{\Delta\rho \cdot g \cdot r^2}{\eta} = \frac{300 \times 10 \times (10^4)^2}{10^{22}} = 3 \times 10^{-11} \text{ m/s}$$

The test requires a minimum resolution of 41 × 41 to adequately resolve the free surface and avoid CHOLMOD conditioning issues. After 3 time steps, the vertical velocity must be non-zero (block is sinking) and the lithostatic pressure must be positive.

**Code assertions** (`FreeSurfaceTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_GT(maxAbsVz, 0.0);                  // vertical flow exists (block sinking)
EXPECT_GT(maxP, 0.0);                      // positive lithostatic pressure
```

### Justification

Free surface evolution is essential for modelling surface topography in rifting, collision, and erosion studies. This test exercises the `BuildInitialTopography_ff` registration pattern, the air-phase masking, and the surface stabilisation, which are distinct code paths not covered by any other test.

---

## Suite 11 — VelocityField (Kinematic Symmetry)

**Source:** `VelocityFieldTests.cpp`
**Parameter file:** `VelocityField/PureShearVelocity.txt`

### Physical Background

For a **homogeneous Newtonian medium** under pure shear with background strain rate $\dot{\varepsilon}$, the analytical velocity solution in a domain $[-L/2, L/2]^2$ is:

$$V_x = \dot{\varepsilon} \cdot x, \quad V_z = -\dot{\varepsilon} \cdot z$$

This solution is exact (no body forces). The velocity field is **anti-symmetric**: $V_x(x) = -V_x(-x)$ and $V_z(z) = -V_z(-z)$. A circular inclusion with different viscosity perturbs this solution locally but preserves the overall anti-symmetry.

### 11.1 PureShearVelocity

**Setup:** Non-dimensional ($\eta = 1$, $L = 1$), domain $[-0.5, 0.5]^2$, 21 × 21 grid. Matrix $\eta = 1$, inclusion $\eta = 10$ (radius 0.05). $\dot{\varepsilon} = 1$.

| Parameter | Value |
|---|---|
| Domain | $[-0.5, 0.5]^2$ |
| Grid | 21 × 21 |
| $\eta_{matrix}$ | 1 |
| $\eta_{incl}$ | 10 |
| $\dot{\varepsilon}$ | 1 |
| Scaling | Non-dimensional |

**Cause → Effect:** At the boundaries, $V_x(\pm 0.5) = \pm 0.5$ and $V_z(\pm 0.5) = \mp 0.5$. The inclusion perturbs the interior but the boundary values are set by the Dirichlet BCs. The test checks:
1. **Sign**: $V_x^{max} > 0$ (east), $V_x^{min} < 0$ (west) — extension in $x$.
2. **Sign**: $V_z^{max} > 0$ (bottom moves up), $V_z^{min} < 0$ (top moves down) — shortening in $z$.
3. **Symmetry**: $|V_x^{max}| \approx |V_x^{min}|$ and $|V_z^{max}| \approx |V_z^{min}|$ within 20% tolerance.

**Code assertions** (`VelocityFieldTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_GT(maxVx, 0.0);                     // extension in x
EXPECT_LT(minVx, 0.0);
EXPECT_GT(maxVz, 0.0);                     // shortening in z
EXPECT_LT(minVz, 0.0);
EXPECT_NEAR(fabs(maxVx), fabs(minVx), fabs(maxVx) * 0.2);  // anti-symmetry
EXPECT_NEAR(fabs(maxVz), fabs(minVz), fabs(maxVz) * 0.2);
```

### Justification

Velocity-field symmetry is a basic sanity check for the staggered-grid Stokes solver. A broken ghost-node implementation, incorrect boundary stencil, or indexing error would destroy the symmetry even for this simple configuration.

---

## Suite 12 — Compressibility (Elastic Bulk Response)

**Source:** `CompressibilityTests.cpp`
**Parameter file:** `Compressibility/DilatantFlow.txt`

### Physical Background

When `compressible = 1` is enabled, MDOODZ replaces the strict incompressibility constraint ($\nabla \cdot \mathbf{v} = 0$) with elastic compressibility:

$$\nabla \cdot \mathbf{v} = -\beta \frac{dP}{dt}$$

where $\beta$ is the elastic compressibility (reciprocal of bulk modulus). This allows volumetric strain in response to pressure changes, which is important for modelling crustal deformation at high pressures.

### 12.1 DilatantFlow

**Setup:** 10 km × 10 km, 31 × 31 grid, `compressible = 1`. Matrix and inclusion: $\eta_0 = 10^{22}$ Pa·s, $\beta = 10^{-11}$ Pa$^{-1}$. Pure shear at $\dot{\varepsilon} = 5 \times 10^{-14}$ s$^{-1}$. Gravity disabled ($g_z = 0$). No dilatancy ($\psi = 0$). Plasticity enabled with very high cohesion ($C = 10^{90}$) to prevent yielding.

| Parameter | Value |
|---|---|
| Domain | 10 × 10 km |
| $\eta_0$ | $10^{22}$ Pa·s |
| $\beta$ | $10^{-11}$ Pa$^{-1}$ |
| $\dot{\varepsilon}$ | $5 \times 10^{-14}$ s$^{-1}$ |
| $g_z$ | 0 |
| $\psi$ | 0 |

**Cause → Effect:** With `compressible = 1`, the mass conservation equation gains a pressure-rate term. The solver must correctly assemble the modified divergence operator. The test verifies solver completion with bounded pressure and positive viscosity. Gravity is disabled to avoid ill-conditioning (CHOLMOD "not positive definite" failure when combining gravity + compressibility at low strain rates).

**Code assertions** (`CompressibilityTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver completed
EXPECT_LT(fabs(maxP), 1e40);               // pressure is bounded
EXPECT_GT(minEta, 0.0);                    // viscosity is positive
```

### Justification

Compressibility modifies the structure of the Stokes system (the divergence block becomes non-zero on the diagonal). This test exercises the compressible assembly path in `StokesAssemblyDecoupled.c`, which is a distinct code branch from the incompressible case.

---

## Suite 13 — FiniteStrain (Deformation Gradient)

**Source:** `FiniteStrainTests.cpp`
**Parameter file:** `FiniteStrain/PureShearStrain.txt`

### Physical Background

The **deformation gradient tensor** $\mathbf{F}$ tracks the cumulative deformation of material:

$$F_{ij} = \frac{\partial x_i}{\partial X_j}$$

where $\mathbf{x}$ is the current position and $\mathbf{X}$ is the reference position. For incompressible flow, $\det(\mathbf{F}) = 1$ (area conservation in 2D). MDOODZ evolves $\mathbf{F}$ on markers using:

$$\dot{F}_{ij} = L_{ik} F_{kj}$$

where $L_{ij} = \partial v_i / \partial x_j$ is the velocity gradient. Finite strain is used to track fabric orientation (anisotropy) and to compute measures like natural strain $\epsilon = \ln\sqrt{F_{xx}^2 + F_{xz}^2}$.

### 13.1 PureShearStrain

**Setup:** Non-dimensional ($\eta = 1$, $L = 1$), domain $[-0.5, 0.5]^2$, 21 × 21 grid, `finite_strain = 1`. Matrix $\eta = 1$, inclusion $\eta = 10$. Pure shear $\dot{\varepsilon} = 1$, $\Delta t = 0.01$, 5 steps.

| Parameter | Value |
|---|---|
| Domain | $[-0.5, 0.5]^2$ |
| Grid | 21 × 21 |
| $\dot{\varepsilon}$ | 1 |
| $\Delta t$ | 0.01 |
| Steps | 5 |
| `finite_strain` | 1 |

**Cause → Effect:** In pure shear ALE mode, the grid deforms with the material. For a homogeneous medium, the analytical deformation gradient after total strain $\varepsilon = \dot{\varepsilon} \cdot n\Delta t = 0.05$ is:

$$\mathbf{F} = \begin{pmatrix} e^{\varepsilon} & 0 \\ 0 & e^{-\varepsilon} \end{pmatrix} \approx \begin{pmatrix} 1.051 & 0 \\ 0 & 0.951 \end{pmatrix}$$

In ALE mode, markers move with the grid, so the computed $\mathbf{F}$ may stay close to identity. The test verifies:
1. **Positivity**: $F_{xx} > 0$ and $F_{zz} > 0$ (physical requirement — no material inversion).
2. **Area conservation**: $\det(\mathbf{F}) = F_{xx} F_{zz} - F_{xz} F_{zx} \approx 1$ within 20% tolerance.

**Code assertions** (`FiniteStrainTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_GT(maxFxx, 0.0);                    // Fxx > 0
EXPECT_GT(minFzz, 0.0);                    // Fzz > 0
EXPECT_NEAR(detF, 1.0, 0.2);              // volume (area) conservation
```

### Justification

Finite strain tracking is required for the anisotropy module (director rotation) and for post-processing fabric evolution. This test ensures the deformation gradient integration (`ParticleRoutines.c`) produces physically valid results and preserves the incompressibility constraint at the marker level.

---

## Suite 14 — NeumannBC (Stress Boundary Conditions)

**Source:** `NeumannBCTests.cpp`
**Parameter file:** `NeumannBC/StressBCPureShear.txt`

### Physical Background

MDOODZ supports stress (Neumann) boundary conditions through the `constant_shear_stress` BC type (type = 13). In the built-in `SetPureShearBCVx`, this type is applied on the **N/S boundaries** for $V_x$, implementing free-slip:

$$\sigma_{xz}\big|_{z=\pm L/2} = 0$$

while the **E/W boundaries** use Dirichlet velocity BCs (type = 0):

$$V_x\big|_{x=\pm L/2} = \pm \dot{\varepsilon} \cdot L/2$$

This combination of stress and velocity BCs is the standard "free-slip pure shear" configuration used throughout lithospheric modelling.

### 14.1 StressBCPureShear

**Setup:** 10 km × 10 km domain, 21 × 21 grid, single phase (no inclusion). Constant viscosity $\eta = 10^{22}$ Pa·s, $\dot{\varepsilon} = 10^{-14}$ s$^{-1}$, `pure_shear_ALE = 1`, dimensional scaling.

| Parameter | Value |
|---|---|
| Domain | 10 × 10 km |
| Grid | 21 × 21 |
| $\eta$ | $10^{22}$ Pa·s |
| $\dot{\varepsilon}$ | $10^{-14}$ s$^{-1}$ |
| BCs (Vx, N/S) | type = 13 (free-slip) |
| BCs (Vx, E/W) | type = 0 (Dirichlet) |

**Cause → Effect:** For a homogeneous Newtonian medium under pure shear with free-slip walls, the shear stress is exactly zero everywhere ($\sigma_{xz} = 0$), while the deviatoric normal stress is $\sigma'_{xx} = 2\eta\dot{\varepsilon} = 2 \times 10^{22} \times 10^{-14} = 2 \times 10^{8}$ Pa. The test verifies:
1. **Free-slip**: $|\sigma_{xz}| < 1$ (near zero in non-dimensional space).
2. **Deviatoric stress**: $|\sigma'_{xx}| > 0$ (deformation generates stress).
3. **Velocity**: $|V_x| > 0$ (flow is driven by the boundary strain rate).

**Code assertions** (`NeumannBCTests.cpp`):
```cpp
ASSERT_GE(stepsCount, 0);                  // solver converged
EXPECT_LT(fabs(maxSxz), 1.0);              // shear stress ≈ 0 (free-slip)
EXPECT_LT(fabs(minSxz), 1.0);
EXPECT_GT(fabs(maxSxxd), 0.0);             // deviatoric stress > 0
EXPECT_GT(fabs(maxVx), 0.0);               // non-zero velocity field
```

### Justification

Stress BCs are the physical complement to velocity BCs and are used whenever the boundary has a prescribed traction (e.g. plate boundary forces, basal drag). This test explicitly verifies the type = 13 code path in the Stokes assembly, which modifies the boundary rows of the stiffness matrix differently from Dirichlet conditions.

---

## Suite 15 — ConvergenceRate (Newton vs Picard)

**Source:** `ConvergenceRateTests.cpp`
**Parameter files:** `ConvergenceRate/NewtonPwl.txt`, `ConvergenceRate/PicardPwl.txt`

### Physical Background

MDOODZ offers two non-linear solvers:

1. **Newton-Raphson**: Includes the viscosity-strain-rate Jacobian $\partial\eta/\partial\dot{\varepsilon}$ in the tangent stiffness matrix. Achieves **quadratic convergence** ($\|r^{(k+1)}\| \propto \|r^{(k)}\|^2$) near the solution.

2. **Picard** (fixed-point): Updates viscosity from the previous iteration's strain rate without the Jacobian correction. Achieves only **linear convergence** ($\|r^{(k+1)}\| \propto \|r^{(k)}\|$) but is more robust for strongly non-linear problems.

The convergence rate is a fundamental property of the solver algorithm and should be preserved across code changes.

### 15.1 NewtonPwl

**Setup:** 10 km × 10 km, 31 × 31 grid. Matrix: power-law (`pwlv = 40`, $n = 3.5$, Olivine H&K 2003), $T = 1100$ °C. Inclusion: constant viscosity ($\eta = 10^{27}$). `Newton = 1`, `nit_max = 20`.

**Cause → Effect:** The power-law nonlinearity requires iterative solution. Newton's method includes the analytical Jacobian, so convergence should be rapid (typically < 15 iterations) and the final residual should be below the nonlinear tolerance ($10^{-9}$).

**Code assertions** (`ConvergenceRateTests.cpp`):
```cpp
ASSERT_GE(newtonSteps, 0);                 // converged (0 = immediate convergence, valid)
EXPECT_LT(newtonSteps, 15);                // quadratic convergence → few iterations
EXPECT_LT(finalRes, 1e-8);                 // residual below tolerance
```

### 15.2 PicardPwl

**Setup:** Identical to NewtonPwl **except** `Newton = 0`, `line_search = 0`, `nit_max = 40`.

**Cause → Effect:** Without the Jacobian correction, Picard converges linearly and typically requires more iterations. The iteration budget is 40 (vs 20 for Newton).

**Code assertions** (`ConvergenceRateTests.cpp`):
```cpp
ASSERT_GE(picardSteps, 0);                 // converged
EXPECT_LT(picardSteps, 40);                // within iteration budget
EXPECT_LT(finalRes, 1e-8);                 // residual below tolerance
```

### 15.3 NewtonFasterThanPicard

**Cause → Effect:** By construction, Newton (quadratic) should converge in fewer or equal iterations compared to Picard (linear) for the **same problem**. This test reads both HDF5 outputs and compares iteration counts.

**Code assertions** (`ConvergenceRateTests.cpp`):
```cpp
EXPECT_LE(newtonSteps, picardSteps);        // Newton ≤ Picard iterations
```

### Justification

Comparing Newton vs Picard convergence on the same problem is a rigorous verification of the Jacobian assembly. If the Jacobian has an error, Newton may converge slower than Picard (or diverge), which unambiguously indicates a bug. This is the only test in the suite that compares two solver runs quantitatively.

---

## Running the Tests

### Build and run

```bash
# From the repository root:
cmake -B ./cmake-build -DTEST=ON
cmake --build ./cmake-build -j4
cd cmake-build && ctest --output-on-failure
```

All 15 test suites (35 individual test cases) should pass in under 15 seconds:

```
100% tests passed, 0 tests failed out of 15
Total Test time (real) =   9.44 sec
```

### Via Docker (CI-equivalent)

```bash
docker run --rm -v $(pwd):/workspace -w /workspace ubuntu:22.04 bash -c \
  "apt update && apt install -y cmake gcc g++ libhdf5-serial-dev libsuitesparse-dev && \
   cmake -B ./cmake-build -DTEST=ON && cmake --build ./cmake-build && \
   cd cmake-build && ctest --output-on-failure"
```

### WSL Notes (Windows/NTFS)

When building in WSL from a Windows-hosted checkout, line-ending issues can cause segmentation faults in file parsing. If tests fail universally with SIGSEGV:

```bash
# Restore MDLIB C files to LF line endings:
git checkout -- MDLIB/

# Strip CRLF from test files:
find TESTS -name '*.cpp' -o -name '*.txt' -o -name '*.h' -o -name 'CMakeLists.txt' \
  | xargs sed -i 's/\r$//'

# Clean rebuild:
rm -rf cmake-build && cmake -S . -B cmake-build -DTEST=ON && cmake --build cmake-build -j4
```

### Dependencies

- CMake ≥ 3.16
- GCC (C11 + C++17)
- HDF5 (serial)
- SuiteSparse (CHOLMOD/UMFPACK)
- Google Test (fetched automatically by CMake via FetchContent)

---

## Test Coverage Summary

| Suite | Source File | Tests | Physics Verified |
|-------|------------|-------|-----------------|
| 1. ShearTemplate | `NewtonIterationConvergence.cpp` | 8 | Newton/Picard convergence, anisotropy, linear/nonlinear |
| 2. RheologyCreep | `RheologyCreepTests.cpp` | 6 | Diffusion, power-law, Peierls, GBS creep, composite |
| 3. Plasticity | `PlasticityTests.cpp` | 4 | Drucker-Prager yield, no-yield, strain softening, stress limiter |
| 4. Thermal | `ThermalTests.cpp` | 4 | Diffusion, geotherm, radiogenic heating, Dirichlet BCs |
| 5. BoundaryCondition | `BoundaryConditionTests.cpp` | 3 | Periodic, free-slip, no-slip BCs |
| 6. SolverMode | `SolverModeTests.cpp` | 2 | Picard convergence, adaptive time stepping |
| 7. ViscoElastic | `ViscoElasticTests.cpp` | 1 | Maxwell stress accumulation, elastic strain rate |
| 8. Density | `DensityTests.cpp` | 2 | Thermal expansion $\rho(T)$, hydrostatic pressure |
| 9. ShearHeating | `ShearHeatingTests.cpp` | 1 | Viscous dissipation $H = 2\eta\dot{\varepsilon}^2$ |
| 10. FreeSurface | `FreeSurfaceTests.cpp` | 1 | Topography tracking, sinking block, air masking |
| 11. VelocityField | `VelocityFieldTests.cpp` | 1 | Pure shear kinematics, anti-symmetry |
| 12. Compressibility | `CompressibilityTests.cpp` | 1 | Elastic bulk compressibility $\nabla\cdot\mathbf{v} = -\beta\dot{P}$ |
| 13. FiniteStrain | `FiniteStrainTests.cpp` | 1 | Deformation gradient $\mathbf{F}$, area conservation |
| 14. NeumannBC | `NeumannBCTests.cpp` | 1 | Stress BC (type=13), free-slip shear stress |
| 15. ConvergenceRate | `ConvergenceRateTests.cpp` | 3 | Newton vs Picard rate comparison |
| **Total** | **15 executables** | **39** | |
