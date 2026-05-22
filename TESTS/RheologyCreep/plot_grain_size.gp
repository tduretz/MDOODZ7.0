# Grain Size Steady-State Benchmark — calcite paleowattmeter (Austin & Evans 2007/2009)
#
# Plots d_ss vs the IMPOSED total strain rate Eii_total. In this view the
# dislocation-only and coupled scenarios trace VISIBLY DISTINCT curves:
# adding diffusion creep means a fraction of Eii_total is absorbed by diffusion
# rather than dislocation, so dislocation creep "sees" a smaller effective Eii
# at any given Eii_total → larger steady-state grain size at any given Eii_total.
#
# (Note: in the alternative (tau_II, d_ss) view, both scenarios trace the SAME
# curve because the wattmeter formula d_ss(tau, Eii_pwl(tau)) is mechanism-
# independent given tau. The (Eii_total, d_ss) view exposes the coupled-physics
# difference where it actually lives — in the partition between mechanisms.)
#
# Run from cmake-build-test/TESTS:
#     gnuplot ../../TESTS/RheologyCreep/plot_grain_size.gp
#
# Output: grain_size_benchmark.png

set terminal pngcairo size 900,650 enhanced font "Arial,12"
set output "grain_size_benchmark.png"

set title "Calcite paleowattmeter — Austin and Evans (2007, 2009) at T = 623 K"
set xlabel "Imposed total strain rate {/Symbol e}'_{II} (s^{-1})"
set ylabel "Steady-state grain size d_{ss} (m)"

set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set xrange [1e-17:1e-11]
set yrange [1e-5:1e-3]

# Paper-style ticks + dense grid network — minor tics inward (at 2, 3, ..., 9 per
# decade), mirrored on top/right. Matches the look of plot_grain_size_paper_view.gp.
set tics in scale 2.0, 1.0
set xtics mirror
set ytics mirror
set mxtics 10
set mytics 10
set style line 81 linetype 1 linecolor rgb "#bbbbbb" linewidth 0.6
set style line 82 linetype 1 linecolor rgb "#e0e0e0" linewidth 0.4
set grid xtics ytics back linestyle 81
set grid mxtics mytics back linestyle 82

set key bottom left

# Constants — must match GrainSizeSteadyState.txt and the test's analytical
T_K   = 623.15
R     = 8.314
# Calcite paleowattmeter — gs=10
p     = 3.0
Kg    = 2.5e9 * 10**(-6.0*p)
Qg    = 175e3
gam   = 1.0
lam   = 0.1
cg    = pi
Bg    = lam / (cg * gam)
Ag    = Kg * exp(-Qg / (R * T_K))
# Calcite power-law — pwlv=15
n_pwl = 4.7
Q_pwl = 297e3
A_pwl = 1.5849e-25
F_pwl = (1.0/6.0) * 2**(1.0/n_pwl) * 3**((n_pwl-1.0)/(2.0*n_pwl))
B_pwl = F_pwl * A_pwl**(-1.0/n_pwl) * exp(Q_pwl / (n_pwl * R * T_K))

# Dislocation-only analytical: closed form. tau = 2*B_pwl*Eii^(1/n), then wattmeter.
tau_disl(Eii) = 2.0 * B_pwl * Eii**(1.0/n_pwl)
d_disl(Eii)   = ( Bg * Eii * tau_disl(Eii) * p / Ag )**(-1.0/(p+1.0))

plot d_disl(x) with lines linewidth 2.5 linecolor rgb "#1f77b4" \
        title "Analytical d_{ss} — dislocation only (Renner 2002 + Austin and Evans 2007/2009)", \
     "grain_size_benchmark.dat" using 1:5 with points pointtype 7 pointsize 1.8 \
        linecolor rgb "#1f77b4" \
        title "MDOODZ — dislocation-only sweep (linv = 0)", \
     "grain_size_benchmark_coupled.dat" index 0 using 1:3 with lines linewidth 2.5 linecolor rgb "#d62728" \
        title "Analytical d_{ss} — coupled (+ Herwegh 2003 diffusion creep)", \
     "grain_size_benchmark_coupled.dat" index 1 using 1:3 with points pointtype 9 pointsize 1.8 \
        linecolor rgb "#d62728" \
        title "MDOODZ — coupled sweep (linv = 15)"

set output
