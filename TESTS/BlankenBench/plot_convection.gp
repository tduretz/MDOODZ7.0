#!/usr/bin/env gnuplot
#
# Blankenbach Case 1a — Thermal convection benchmark visualisation
#
# Plots the temperature field as a colour map with velocity vector overlay,
# showing the convection cell pattern for the Blankenbach Ra=10^4 benchmark.
#
# Prerequisites:
#   cd build/TESTS
#   ./BlankenBenchTests                           # produces BlankenBench/Output00500.gzip.h5
#   ./extract_blankenbach BlankenBench/Output00500.gzip.h5
#   (generates BlankenBench/temperature.dat and BlankenBench/velocity.dat)
#
# Usage:
#   gnuplot BlankenBench/plot_convection.gp
#
# Output: blankenbach_convection.png

reset
set terminal pngcairo size 900,800 enhanced font "Arial,13"
set output "blankenbach_convection.png"

# --- Layout ---
set xlabel "x (km)" font ",13"
set ylabel "z (km)" font ",13"
set title "Blankenbach Case 1a: Temperature and Velocity (Ra = 10^4)" font ",15"
set size ratio 1
set xrange [-50:50]
set yrange [-50:50]

# --- Temperature colour map ---
set pm3d map interpolate 2,2
set palette defined ( \
    0.0 "#313695", 0.1 "#4575B4", 0.2 "#74ADD1", \
    0.3 "#ABD9E9", 0.4 "#E0F3F8", 0.5 "#FFFFBF", \
    0.6 "#FEE090", 0.7 "#FDAE61", 0.8 "#F46D43", \
    0.9 "#D73027", 1.0 "#A50026" )
set cbrange [0:1]
set cblabel "T / {/Symbol D}T" font ",13" offset 1,0
set colorbox

# --- Velocity vectors ---
# velocity.dat has: x(km) z(km) vx_norm vz_norm
# every 2 = plot every 2nd vector for readability

splot "BlankenBench/temperature.dat" using 1:2:3 with pm3d notitle, \
      "BlankenBench/velocity.dat" every 2:2 using 1:2:(0):($3*6.0):($4*6.0):(0) \
      with vectors head size 0.8,20 filled lc rgb "#000000" lw 1.0 notitle

set output
print "Generated: blankenbach_convection.png"
