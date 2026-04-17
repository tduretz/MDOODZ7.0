# cell_uniformity.gnu — Per-fine-cell particle count distribution
# Shows that mode 2 deactivation produces a tighter, more uniform distribution

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1400,500 enhanced font "Arial,12"
set output sprintf("%s/cell_uniformity.png", outdir)

set multiplot layout 1,2 title "Per-Cell Particle Count Distribution (After 50 Steps)" font "Arial,14"

set grid
set style data boxes
set style fill solid 0.6 border -1
set boxwidth 0.35

# --- Left panel: Normal density (4×4 particles/cell) ---
set title "Normal Density (4×4 / cell)" font "Arial,13"
set xlabel "Active Particles per Fine Cell"
set ylabel "Number of Fine Cells"
set key top right font "Arial,10"

plot "reseed_cellhist_Mode0.dat" using ($1-0.35):2 with boxes lw 1 lc rgb "#E31A1C" title "Mode 0", \
     "reseed_cellhist_Mode1.dat" using 1:2 with boxes lw 1 lc rgb "#FF7F00" title "Mode 1", \
     "reseed_cellhist_Mode2.dat" using ($1+0.35):2 with boxes lw 1 lc rgb "#1F78B4" title "Mode 2"

# --- Right panel: High density (6×6 particles/cell), Mode 1 vs Mode 2 ---
set title "High Density (6×6 / cell)" font "Arial,13"
set xlabel "Active Particles per Fine Cell"
set ylabel "Number of Fine Cells"
set key top left font "Arial,10"
set xrange [0.5:11.5]
set yrange [0:900]

# Threshold annotation
set arrow 1 from 8,0 to 8,900 nohead lw 2 dt 2 lc rgb "#666666"
set label 1 "deactivation\nthreshold" at 5.5,800 font "Arial,10" tc rgb "#666666"
set arrow 2 from 7.5,780 to 7.95,780 head lw 1.5 lc rgb "#666666"

# Value labels for clipped bars
set label 2 "1400" at 7.5,870 center font "Arial,11,Bold" tc rgb "#1F78B4"
set arrow 3 from 7.7,860 to 8.15,860 head lw 1 lc rgb "#1F78B4"

plot "reseed_cellhist_HighDensMode1.dat" using ($1-0.175):2 with boxes lw 1 lc rgb "#FF7F00" title "Mode 1", \
     "reseed_cellhist_HighDensMode2.dat" using ($1+0.175):2 with boxes lw 1 lc rgb "#1F78B4" title "Mode 2"

unset multiplot
