#!/usr/bin/env gnuplot
#
# Two-panel pattern check for Blankenbach benchmark
#

reset
set terminal pngcairo size 1400,600 enhanced font "Arial,13"
set output "BlankenBench/pattern_check.png"

set multiplot layout 1,2

# --- Panel 1: Temperature + velocity vectors ---
set xlabel "x (km)"
set ylabel "z (km)"
set title "Step 26000: Temperature + Velocity (0.5 {/Symbol t}_{diff})" font ",14"
set size ratio 1
set xrange [-50:50]
set yrange [-50:50]
set pm3d map interpolate 2,2
set palette defined ( \
    0.0 "#313695", 0.1 "#4575B4", 0.2 "#74ADD1", \
    0.3 "#ABD9E9", 0.4 "#E0F3F8", 0.5 "#FFFFBF", \
    0.6 "#FEE090", 0.7 "#FDAE61", 0.8 "#F46D43", \
    0.9 "#D73027", 1.0 "#A50026" )
set cbrange [0:1]
set cblabel "T / {/Symbol D}T" offset 1,0
splot "BlankenBench/temperature.dat" using 1:2:3 with pm3d notitle, \
      "BlankenBench/velocity.dat" every 3:3 using 1:2:(0):($3*5.0):($4*5.0):(0) \
      with vectors head size 0.6,20 filled lc rgb "#000000" lw 0.8 notitle

# --- Panel 2: Isotherms (Blankenbach paper style) ---
set title "Step 26000: Isotherms (T/DT = 0.1, 0.2, ..., 0.9)" font ",14"
unset pm3d
unset colorbox
set view map
set contour base
set cntrparam levels incremental 0.1, 0.1, 0.9
unset surface
set style data lines
set key off
splot "BlankenBench/temperature.dat" using 1:2:3 with lines lc rgb "black" lw 1.2 notitle

unset multiplot
set output
print "Generated: BlankenBench/pattern_check.png"
