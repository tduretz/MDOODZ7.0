# nbpart_timeseries.gnu — Nb_part and Nb_active over time for all 3 modes

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1400,500 enhanced font "Arial,12"
set output sprintf("%s/nbpart_timeseries.png", outdir)

set multiplot layout 1,2 title "Particle Counts Over Time (Normal Density: 4×4 / cell)" font "Arial,14"

set grid
set key top left font "Arial,10"

# --- Left panel: Total Nb_part ---
set title "Total Array Size (Nb\\_part)" font "Arial,13"
set xlabel "Time Step"
set ylabel "Nb\\_part"

plot "reseed_nbpart_Mode0.dat" using 1:2 with linespoints lw 2 pt 5 ps 0.6 lc rgb "#E31A1C" title "Mode 0", \
     "reseed_nbpart_Mode1.dat" using 1:2 with linespoints lw 2 pt 7 ps 0.6 lc rgb "#FF7F00" title "Mode 1", \
     "reseed_nbpart_Mode2.dat" using 1:2 with linespoints lw 2 pt 9 ps 0.6 lc rgb "#1F78B4" title "Mode 2"

# --- Right panel: Active particles ---
set title "Active Particles (phase ≠ -1)" font "Arial,13"
set xlabel "Time Step"
set ylabel "Nb\\_active"

plot "reseed_nbpart_Mode0.dat" using 1:3 with linespoints lw 2 pt 5 ps 0.6 lc rgb "#E31A1C" title "Mode 0", \
     "reseed_nbpart_Mode1.dat" using 1:3 with linespoints lw 2 pt 7 ps 0.6 lc rgb "#FF7F00" title "Mode 1", \
     "reseed_nbpart_Mode2.dat" using 1:3 with linespoints lw 2 pt 9 ps 0.6 lc rgb "#1F78B4" title "Mode 2"

unset multiplot
