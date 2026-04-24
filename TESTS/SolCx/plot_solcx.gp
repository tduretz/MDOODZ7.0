#!/usr/bin/env gnuplot
#
# SolCx grid-convergence plot.
#
# Reads solcx_convergence.dat (written by SolCxBenchmarkTests.GridConvergence)
# and produces solcx_convergence.png: log-log L2(Vx), L2(Vz), L1(P) vs h at
# three resolutions (21, 41, 81), with dashed 1st- and 2nd-order reference lines.
#
# Prerequisites:
#   cd cmake-build/TESTS
#   ./SolCxBenchmarkTests --gtest_filter=SolCxBenchmark.GridConvergence
#   (writes solcx_convergence.dat to the current directory)
#
# Usage:
#   gnuplot SolCx/plot_solcx.gp
#
# Output: SolCx/solcx_convergence.png
#
# Reference: Duretz, May, Gerya, Tackley (2011) "Discretization errors and
# free surface stabilization in the finite difference and marker-in-cell
# method for applied geodynamics", G-Cubed 12, Q07004.

reset
set terminal pngcairo size 800,600 enhanced font "Arial,13"
set output "SolCx/solcx_convergence.png"

set title "SolCx Grid Convergence (step viscosity 1 to 10^6, staggered FD)" font ",14"
set xlabel "h = 1/N" font ",13"
set ylabel "L2 / L1 error" font ",13"
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set key bottom right spacing 1.2
set grid xtics ytics mxtics mytics ls 0 lc rgb "#bbbbbb"
set xrange [0.008 : 0.08]
set yrange [1e-4 : 3e-1]

# Reference slopes: anchor at (h=5e-2, y=2e-1) for slope-1; y=1e-1 for slope-2
ref1(x) = 2e-1 * (x / 5e-2)**1
ref2(x) = 1e-1 * (x / 5e-2)**2

plot \
    "solcx_convergence.dat" using 1:2 with linespoints lw 2 pt 7 ps 1.2 lc rgb "#d62728" title "L_2(V_x)", \
    "solcx_convergence.dat" using 1:3 with linespoints lw 2 pt 5 ps 1.2 lc rgb "#1f77b4" title "L_2(V_z)", \
    "solcx_convergence.dat" using 1:4 with linespoints lw 2 pt 9 ps 1.2 lc rgb "#2ca02c" title "L_1(P)", \
    ref1(x) with lines dt 2 lw 1.5 lc rgb "#888888" title "slope 1", \
    ref2(x) with lines dt 3 lw 1.5 lc rgb "#888888" title "slope 2"
