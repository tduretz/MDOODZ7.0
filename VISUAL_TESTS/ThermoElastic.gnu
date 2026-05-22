# ThermoElastic.gnu — Two-panel: (1) T vs time, (2) mean P vs time
set terminal pngcairo enhanced size 900,500 font "Helvetica,13"
if (!exists("filename")) filename='../VISUAL_TESTS/img/thermoelastic.png'
set output filename

if (!exists("data")) data='thermoelastic.dat'

set multiplot layout 1,2 title "Thermo-Elastic Coupling Verification"

set format x "%.1e"
set xtics rotate by -30

# Panel 1: Temperature vs time
set xlabel "Time [s]"
set ylabel "Mean Temperature [K]"
set title "Temperature evolution"
set grid
set key top left

plot data using 1:2 with linespoints lw 2 pt 7 ps 0.8 lc rgb "#2166ac" title "CHOLMOD", \
     data using 1:3 with linespoints lw 2 pt 5 ps 0.8 lc rgb "#b2182b" title "PCG", \
     data using 1:4 with lines lw 2 dt 2 lc rgb "#000000" title "Analytical"

# Panel 2: Pressure vs time
set xlabel "Time [s]"
set ylabel "Mean Pressure [Pa]"
set title "Pressure response"
set key bottom right

plot data using 1:5 with linespoints lw 2 pt 7 ps 0.8 lc rgb "#2166ac" title "CHOLMOD", \
     data using 1:6 with linespoints lw 2 pt 5 ps 0.8 lc rgb "#b2182b" title "PCG", \
     data using 1:7 with lines lw 2 dt 2 lc rgb "#000000" title "Analytical"

unset multiplot
