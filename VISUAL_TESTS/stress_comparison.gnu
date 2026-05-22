# stress_comparison.gnu — The case for Mode 2
# Left: Normal density — Mode 0 Nb_part grows; modes 1/2 stable
# Right: High density (6×6) — Mode 2 actively deactivates excess particles

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1400,500 enhanced font "Arial,12"
set output sprintf("%s/stress_comparison.png", outdir)

set multiplot layout 1,2 title "The Case for Mode 2: Controlled Particle Population" font "Arial,14"

set grid
set key top left font "Arial,11"

# --- Left panel: Normal density Nb_part (mode 0 grows) ---
set title "Normal Density (4×4): Nb\\_part Growth" font "Arial,13"
set xlabel "Time Step"
set ylabel "Nb\\_part"

plot "reseed_nbpart_Mode0.dat" using 1:2 with linespoints lw 2.5 pt 5 ps 0.7 lc rgb "#E31A1C" title "Mode 0 — unbounded growth", \
     "reseed_nbpart_Mode1.dat" using 1:2 with linespoints lw 2.5 pt 7 ps 0.7 lc rgb "#FF7F00" title "Mode 1 — stable", \
     "reseed_nbpart_Mode2.dat" using 1:2 with linespoints lw 2.5 pt 9 ps 0.7 lc rgb "#1F78B4" title "Mode 2 — stable"

# --- Right panel: High density Nb_active (mode 2 deactivates) ---
set title "High Density (6×6): Active Particles" font "Arial,13"
set xlabel "Time Step"
set ylabel "Nb\\_active (phase ≠ -1)"

plot "reseed_nbpart_HighDensMode1.dat" using 1:3 with linespoints lw 2.5 pt 7 ps 0.7 lc rgb "#FF7F00" title "Mode 1 — no deactivation", \
     "reseed_nbpart_HighDensMode2.dat" using 1:3 with linespoints lw 2.5 pt 9 ps 0.7 lc rgb "#1F78B4" title "Mode 2 — active cleanup"

unset multiplot
