#!/usr/bin/env gnuplot
#
# 4-panel pattern evolution for Blankenbach benchmark
#

reset
set terminal pngcairo size 1600,1400 enhanced font "Arial,12"
set output "BlankenBench/pattern_evolution.png"

set multiplot layout 2,2 title "Blankenbach Case 1a (Ra=10^4) — Pattern Evolution" font ",16"

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
set cblabel "T/{/Symbol D}T" offset 1,0

# Panel 1: step 500
set xlabel "x (km)"
set ylabel "z (km)"
set title "Step 500 (0.01 {/Symbol t}_{diff})" font ",13"
splot "BlankenBench/temperature_00500.dat" using 1:2:3 with pm3d notitle, \
      "BlankenBench/velocity_00500.dat" every 3:3 using 1:2:(0):($3*5.0):($4*5.0):(0) \
      with vectors head size 0.5,20 filled lc rgb "#000000" lw 0.6 notitle

# Panel 2: step 5000
set title "Step 5000 (0.10 {/Symbol t}_{diff})" font ",13"
splot "BlankenBench/temperature_05000.dat" using 1:2:3 with pm3d notitle, \
      "BlankenBench/velocity_05000.dat" every 3:3 using 1:2:(0):($3*5.0):($4*5.0):(0) \
      with vectors head size 0.5,20 filled lc rgb "#000000" lw 0.6 notitle

# Panel 3: step 13000
set title "Step 13000 (0.25 {/Symbol t}_{diff})" font ",13"
splot "BlankenBench/temperature_13000.dat" using 1:2:3 with pm3d notitle, \
      "BlankenBench/velocity_13000.dat" every 3:3 using 1:2:(0):($3*5.0):($4*5.0):(0) \
      with vectors head size 0.5,20 filled lc rgb "#000000" lw 0.6 notitle

# Panel 4: step 26000
set title "Step 26000 (0.50 {/Symbol t}_{diff})" font ",13"
splot "BlankenBench/temperature_26000.dat" using 1:2:3 with pm3d notitle, \
      "BlankenBench/velocity_26000.dat" every 3:3 using 1:2:(0):($3*5.0):($4*5.0):(0) \
      with vectors head size 0.5,20 filled lc rgb "#000000" lw 0.6 notitle

unset multiplot
set output
print "Generated: BlankenBench/pattern_evolution.png"
