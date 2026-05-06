# Grain Size Steady-State Benchmark — Calcite Paleowattmeter (Austin & Evans 2002)
#
# Plots the analytical d_ss(tau_II) curve over a stress range, with the MDOODZ
# measurement overlaid as a single point. Generated from grain_size_benchmark.dat
# which the GTest case `RheologyCreep.GrainSizeSteadyState` writes after running.
#
# Run from cmake-build-test/TESTS:
#     gnuplot ../../TESTS/RheologyCreep/plot_grain_size.gp
#
# Output: grain_size_benchmark.png (sibling to director_benchmark.png, etc.)

set terminal pngcairo size 900,650 enhanced font "Arial,12"
set output "grain_size_benchmark.png"

set title "Calcite paleowattmeter — Austin & Evans (2002)"
set xlabel "{/Symbol t}_{II} (Pa)"
set ylabel "Steady-state grain size d_{ss} (m)"

set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set xrange [1e4:1e9]
set yrange [1e-6:1e-2]
set grid
set key bottom left

# Constants (must match GrainSizeSteadyState.txt and the analytical computed in the test)
T_K   = 623.15       # 350 °C + 273.15
R     = 8.314        # J/mol/K
# Calcite paleowattmeter — gs=10 — Austin & Evans (2002)
p     = 3.0
Kg    = 2.5e9 * 10**(-6.0*p)
Qg    = 175e3
gam   = 1.0
lam   = 0.1
cg    = pi
Bg    = lam / (cg * gam)
Ag    = Kg * exp(-Qg / (R * T_K))
# Calcite power-law — pwlv=15 — Renner et al. (2002)
n_pwl = 4.7
Q_pwl = 297e3
A_pwl = 1.5849e-25
F_pwl = (1.0/6.0) * 2**(1.0/n_pwl) * 3**((n_pwl-1.0)/(2.0*n_pwl))
B_pwl = F_pwl * A_pwl**(-1.0/n_pwl) * exp(Q_pwl / (n_pwl * R * T_K))

# Closed-form d_ss along the single-mechanism dislocation steady state. Each test
# scenario imposes a different bkg_strain_rate Eii, which fixes tau via Renner+02:
#     Eii_pwl = (tau / (2*B_pwl))^n
# Substituting into the Austin & Evans wattmeter d_ss = (Bg * Eii_pwl * tau * p / Ag)^(-1/(p+1))
# gives a *single* curve in (tau, d) space that all 5 MDOODZ sweep points must land on:
d_ss(tau) = ( Bg * (tau/(2*B_pwl))**n_pwl * tau * p / Ag )**(-1.0/(p+1.0))

plot d_ss(x) with lines linewidth 2.5 linecolor rgb "#1f77b4" \
        title "Analytical d_{ss}({/Symbol t}_{II}) (Austin&Evans 2002 + Renner+02 dislocation, T=623 K)", \
     "grain_size_benchmark.dat" using 3:5 with points pointtype 7 pointsize 1.8 \
        linecolor rgb "#d62728" \
        title "MDOODZ sweep — 5 strain rates, {/Symbol e}'_{II}={10^{-16}, 10^{-15}, 10^{-14}, 10^{-13}, 10^{-12}} s^{-1}"

set output
