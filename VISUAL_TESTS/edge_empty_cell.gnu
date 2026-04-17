# edge_empty_cell.gnu — Bulk reseeding comparison: Mode 1 vs Mode 2 (low density)

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1800,500 enhanced font "Arial,11"
set output sprintf("%s/edge_empty_cell.png", outdir)

set multiplot layout 1,3 title "Edge Case: Low Particle Density (min\\_part\\_cell=4)" font "Arial,13"

Ncx = 20
dx = 1.0 / (2.0 * Ncx)
xmin = -0.5
zmin = -0.5

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set size ratio 1

set style line 10 lc rgb "#CCCCCC" lw 0.3
do for [i=0:2*Ncx] {
  set arrow from (xmin + i*dx), -0.5 to (xmin + i*dx), 0.5 nohead ls 10 front
  set arrow from -0.5, (zmin + i*dx) to 0.5, (zmin + i*dx) nohead ls 10 front
}

set palette defined (0 "grey70", 1 "red")
set cbrange [0:1]
unset colorbox

set title "Step 0 — Initial (Sparse)"
set xlabel "x"
set ylabel "z"
plot "reseed_edge_init_EdgeLowPart.dat" using 1:2:4 with points pt 7 ps 0.5 lc palette notitle

set title "Step 50 — After Reseeding (Mode 1)"
plot "reseed_edge_final_EdgeLowPart_Mode1.dat" using 1:2:4 with points pt 7 ps 0.5 lc palette notitle

set title "Step 50 — After Reseeding (Mode 2)"
plot "reseed_edge_final_EdgeLowPart.dat" using 1:2:4 with points pt 7 ps 0.5 lc palette notitle

unset multiplot
