# Blankenbach Benchmark (Blankenbach et al., 1989)

## Modifications

### Input file

The following sections describe the essential modifications relative to the input file for the Blankenbach benchmark provided by Roman. Additional minor changes are not discussed here but can be inspected in the input file in the `/TESTS/BlankenBench` directory.

#### Scaling parameters

The model is defined using dimensional material and domain parameters. The following characteristic scales are used:

```text
/***** SCALES *****/
eta = 1.0e23
L   = 1.0e6
V   = 1.0e-12
T   = 1.0e3
```

These values are mutually consistent because the thermal diffusivity is

$\kappa = \frac{k}{\rho C_p} = \frac{5}{4000\cdot1250} = 10^{-6}\ \mathrm{m^2\,s^{-1}},
$

and therefore

$
V_{\mathrm{scale}} = \frac{\kappa}{L} = 10^{-12}\ \mathrm{m\,s^{-1}}.
$

It would also be possible to formulate the benchmark using unit scaling parameters and nondimensional input quantities. In that case, particular care would be required when specifying the buoyancy term because the Rayleigh number would have to enter the nondimensional momentum equation explicitly.

The dimensional formulation was retained because the material and model parameters are easier to interpret physically.

#### Linear solver

The penalty parameter is set to

```text
penalty = 1e6
```

#### Nonlinear solver

The constant-viscosity Stokes problem does not require nonlinear rheological iterations. These iterations are therefore effectively disabled using

```text
/***** NON-LINEAR SOLVER *****/
Newton      = 0
line_search = 0
nit_max     = 1
```

#### Marker handling

The marker-reseeding mode is set to

```text
reseed_mode = 2
```

Tests with different marker densities and with reseeding disabled showed that the marker configuration was not the primary source of the remaining oscillations in the Nusselt number.

#### Material phase

The material phase with `ID = 0` is assigned the benchmark properties:

```text
ID   = 0
rho  = 4000.0
Cp   = 1250.0
k    = 5.0
Qr   = 0.0
alp  = 2.5e-5
bet  = 0.0
drho = 0.0
C    = 1.0e90
```

Plastic yielding is suppressed by assigning a very large cohesion. Elasticity is disabled separately using the corresponding model switch.

For these parameters, the Rayleigh number is

$ Ra = \frac{\rho g\alpha\Delta T L^3}{\kappa\eta} = \frac{\rho^2 C_p g\alpha\Delta T L^3}{k\eta}.
$

Consequently, changing the constant phase viscosity gives

```text
eta0 = 1e23 Pa s  -> Ra = 1e4
eta0 = 1e22 Pa s  -> Ra = 1e5
eta0 = 1e21 Pa s  -> Ra = 1e6
```

The characteristic viscosity used for scaling remains $10^{23}\ \mathrm{Pa\,s}$, while `eta0` defines the physical viscosity of the material phase.

#### Additional settings

The following additional settings are used:

```text
interp_mode    = 3
thermal_solver = 1
```

### MDOODZ source files

#### `HDF5Output.c`

Datasets for the root-mean-square velocity and the upper and lower Nusselt numbers were added to the `TimeSeries` group written by `WriteOutputHDF5()`:

```text
Vrms
Nu_top
Nu_bottom
```

#### `InputOutput.c`

The three new diagnostic histories were added to the breakpoint writer and reader:

```text
Vrms_time
Nu_top_time
Nu_bottom_time
```

The serialization of all time-series fields in the breakpoint files was also changed. Previously, `Nt + 1` entries were written and loaded. The number of stored entries is now determined from the current model step:

```c
nTime = (size_t)model.step + 1;
```

when writing a breakpoint, and

```c
nTime = (size_t)model->step + 1;
```

when loading a breakpoint.

A breakpoint therefore contains the complete time-series history from step zero through the current step. When restarting a model, these historical values are restored and the continued simulation appends new values to the same arrays. The final value of `Nt` can consequently be increased during a restart without erasing or misaligning the existing time-series data.

#### `Main_DOODZ.c`

Calls to `ComputeNusseltNumber()` were added for the initial state and immediately after the thermal solve.

This ensures that the Nusselt numbers are evaluated from the newly calculated Eulerian temperature field before the subsequent update of the marker temperatures.

#### `mdoodz-private.h`

Scalar and time-series fields for the root-mean-square velocity and the upper and lower Nusselt numbers were added to the grid structure. The declaration of `ComputeNusseltNumber()` was also added.

#### `MemoryAllocFree.c`

Memory allocation and deallocation were added for the new time-series fields:

```text
Vrms_time
Nu_top_time
Nu_bottom_time
```

#### `RheologyDensity.c`

The root-mean-square velocity is calculated from the velocity components interpolated from the staggered velocity nodes to the centres of the active cells.

A separate `ComputeNusseltNumber()` function calculates the upper and lower Nusselt numbers from the Eulerian cell-centred temperature field. Because the physical boundaries lie half a grid spacing from the nearest temperature centroids, the boundary-normal temperature gradients are calculated using the corresponding half-cell distance.

`LogTimeSeries()` stores the calculated values of $V_{\mathrm{RMS}}$, $\mathrm{Nu}_{\mathrm{top}}$, and $\mathrm{Nu}_{\mathrm{bottom}}$ in their respective time-series arrays.

Tests showed that changing the marker density, disabling reseeding, using a single OpenMP thread, and reducing the Courant number had little effect on the remaining Nusselt-number oscillations. Increasing the spatial grid resolution produced the clearest reduction in their amplitude. The remaining oscillations are therefore most consistent with the limited spatial resolution of the thermal boundary gradients.

### Julia visualisation

The file `Main_Visualisation_Makie_MD7_Blanken.jl` was added to visualise the complete time series of

- the nondimensional root-mean-square velocity $V_{\mathrm{RMS}}$;
- the upper Nusselt number $\mathrm{Nu}_{\mathrm{top}}$;
- the lower Nusselt number $\mathrm{Nu}_{\mathrm{bottom}}$.

# Results

## Benchmark values

The reference values for the constant-viscosity benchmark cases of Blankenbach et al. (1989) are:

| $Ra$ | Case | $V_{\mathrm{RMS}}$ | $\mathrm{Nu}$ |
|---:|:---:|---:|---:|
| $10^4$ | 1a | 42.864947 | 4.884409 |
| $10^5$ | 1b | 193.21454 | 10.534095 |
| $10^6$ | 1c | 833.98977 | 21.972465 |

The benchmark values of $V_{\mathrm{RMS}}$ are nondimensional. The dimensional value stored by MDOODZ is converted according to

$
V_{\mathrm{RMS}}^* = \frac{V_{\mathrm{RMS}}^{\mathrm{dim}}}{V_{\mathrm{scale}}}.
$

## Statistical quantities

At every time step, the mean Nusselt number is calculated as

$
\mathrm{Nu}_{\mathrm{mean}}(t) = \frac{\mathrm{Nu}_{\mathrm{top}}(t) + \mathrm{Nu}_{\mathrm{bottom}}(t)}{2}.
$

The reported Nusselt number is its temporal average over the selected steady-state interval:

$
\overline{\mathrm{Nu}} = \operatorname{mean}\left(\mathrm{Nu}_{\mathrm{mean}}\right).
$

The temporal standard deviation is

$
\sigma_{\mathrm{Nu}} = \operatorname{std}\left(\mathrm{Nu}_{\mathrm{mean}}\right).
$

The signed imbalance is defined as

$
I_{\mathrm{signed}} = \operatorname{mean}\left(\mathrm{Nu}_{\mathrm{top}} - \mathrm{Nu}_{\mathrm{bottom}}\right).
$

It measures the systematic difference between the upper and lower heat fluxes:

- $I_{\mathrm{signed}}>0$: the upper Nusselt number is larger on average;
- $I_{\mathrm{signed}}<0$: the lower Nusselt number is larger on average;
- $I_{\mathrm{signed}}\approx0$: there is no substantial systematic bias.

Positive and negative instantaneous deviations can cancel in the signed imbalance.

The absolute imbalance is defined as

$
I_{\mathrm{abs}} = \operatorname{mean}\left(\left|\mathrm{Nu}_{\mathrm{top}} - \mathrm{Nu}_{\mathrm{bottom}}\right|\right).
$

It measures the typical instantaneous disagreement between the upper and lower Nusselt numbers, irrespective of its sign. Consequently,

$
\left|I_{\mathrm{signed}}\right|\leq I_{\mathrm{abs}}.
$

All averages and standard deviations reported below are evaluated over the selected steady-state part of the respective time series.

## Low-Rayleigh-number case: $Ra=10^4$

| Cells | $V_{\mathrm{RMS}}$ | $\mathrm{Nu}$ | $\sigma_{\mathrm{Nu}}$ | $I_{\mathrm{signed}}$ | $I_{\mathrm{abs}}$ |
|---:|---:|---:|---:|---:|---:|
| $40^2$  | 42.886 | 4.8364 | 0.0298 | −0.0006172 | 0.04799169 |
| $80^2$  | 43.014 | 4.8578 | 0.0161 | −0.0012847 | 0.02293156 |
| $160^2$ | 43.038 | 4.8616 | 0.0104 | −0.0035791 | 0.01028848 |

The calculated Nusselt number approaches the benchmark value as the grid resolution increases. Both its temporal variability and the absolute imbalance between the upper and lower boundaries decrease.

The root-mean-square velocity does not converge monotonically towards the benchmark value. Its relative difference increases from approximately $0.05\%$ for the $40^2$ model to approximately $0.40\%$ for the $160^2$ model.

## High-Rayleigh-number case: $Ra=10^6$

| Cells | $V_{\mathrm{RMS}}$ | $\mathrm{Nu}$ | $\sigma_{\mathrm{Nu}}$ | $I_{\mathrm{signed}}$ | $I_{\mathrm{abs}}$ |
|---:|---:|---:|---:|---:|---:|
| $80^2$  | 827.891 | 21.362 | 0.286 | −0.2878 | 0.463 |
| $160^2$ | 839.724 | 21.838 | 0.214 | −0.3378 | 0.366 |

Increasing the resolution substantially improves the calculated Nusselt number. Its relative difference from the benchmark decreases from approximately $2.78\%$ for the $80^2$ model to approximately $0.61\%$ for the $160^2$ model. The temporal variability and the absolute imbalance also decrease.

The root-mean-square velocity crosses the benchmark value. It is approximately $0.73\%$ below the benchmark for the $80^2$ model and approximately $0.69\%$ above it for the $160^2$ model. At least one additional, higher-resolution calculation is therefore required to determine its convergence trend reliably.
