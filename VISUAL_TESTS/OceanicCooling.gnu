# OceanicCooling.gnu — Depth vs temperature: analytical, CHOLMOD, PCG
set terminal pngcairo enhanced size 700,600 font "Helvetica,13"
if (!exists("filename")) filename='../VISUAL_TESTS/img/oceanic_cooling.png'
set output filename

set xlabel "Temperature [K]"
set ylabel "Depth [km]"
set title "Oceanic Cooling — Half-Space erf Verification"
set grid
set key top right

if (!exists("data")) data='oceanic_cooling.dat'

plot data using 2:($1/1e3) with lines lw 2 dt 2 lc rgb "#000000" title "Analytical (erf)", \
     data using 3:($1/1e3) with lines lw 2 lc rgb "#2166ac" title "CHOLMOD", \
     data using 4:($1/1e3) with linespoints lw 1.5 pt 7 ps 0.6 dt 3 lc rgb "#b2182b" title "PCG"
