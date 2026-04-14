# InterpEquivalence.gnu — 2D heatmap of |T_mode0 - T_mode3|
set terminal pngcairo enhanced size 700,600 font "Helvetica,14"
set output filename

set xlabel "x [m]"
set ylabel "z [m]"
set cblabel "|T_{mode0} - T_{mode3}| [K]"

set title "Interpolation Equivalence: Mode 0 vs Mode 3"

set pm3d map
set palette defined (0 "white", 1 "yellow", 2 "red", 3 "dark-red")
set size ratio -1

splot data using 1:2:3 notitle
