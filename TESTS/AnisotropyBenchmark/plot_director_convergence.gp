#!/usr/bin/env gnuplot
#
# Director evolution benchmark: analytical vs MDOODZ results at multiple dt
#
# Shows θ(t) trajectory — analytical curve vs real MDOODZ HDF5 output
# at 5 time step sizes, demonstrating first-order convergence visually.
#
# Prerequisites:
#   cd build/TESTS
#   ./AnisotropyBenchmarkTests --gtest_filter="*DtConvergence"
#   (generates director_trajectory_dt*.dat)
#
# Usage:
#   gnuplot AnisotropyBenchmark/plot_director_convergence.gp
#
# Output: director_benchmark.png

reset
set terminal pngcairo size 800,600 enhanced font "Arial,13"
set output "director_benchmark.png"

# --- Physical parameters ---
theta0 = pi/4.0      # initial director angle (rad) = 45°
gdot   = 1.0         # shear rate: dudz = 2 * bkg_strain_rate = 1.0

# Analytical solution: θ(t) = atan(tan(θ₀) - γ̇·t)
analytical(t) = atan(tan(theta0) - gdot*t) * 180.0/pi

set xlabel "Time t" font ",13"
set ylabel "Director angle θ (°)" font ",13"
set title "Director Evolution: Analytical vs MDOODZ" font ",15"
set key left bottom font ",11" spacing 1.3
set xrange [0:0.55]
set yrange [22:48]
set grid

plot analytical(x) with lines lw 3 lc rgb "#000000" title "Analytical", \
     "director_trajectory_dt250.dat"  using 1:2 with linespoints pt 7 ps 1.4 lw 1.5 lc rgb "#E41A1C" title "Δt = 0.250 (2 steps)", \
     "director_trajectory_dt100.dat"  using 1:2 with linespoints pt 5 ps 1.1 lw 1.3 lc rgb "#377EB8" title "Δt = 0.100 (5 steps)", \
     "director_trajectory_dt050.dat"  using 1:2 with linespoints pt 9 ps 0.9 lw 1.1 lc rgb "#4DAF4A" title "Δt = 0.050 (10 steps)", \
     "director_trajectory_dt025.dat"  using 1:2 with linespoints pt 13 ps 0.7 lw 0.9 lc rgb "#984EA3" title "Δt = 0.025 (20 steps)", \
     "director_trajectory_dt0125.dat" using 1:2 with lines lw 1.8 lc rgb "#FF7F00" title "Δt = 0.0125 (40 steps)"
set output

print "Generated: director_benchmark.png"
