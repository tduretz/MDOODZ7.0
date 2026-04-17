# edge_multiphase.gnu — 3-phase multiphase: dominant-phase selection in interface cells

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1200,500 enhanced font "Arial,11"
set output sprintf("%s/edge_multiphase.png", outdir)

set multiplot layout 1,2 title "Edge Case: 3-Phase Multiphase Reseeding" font "Arial,13"

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

# Phase colors: 0=blue, 1=orange, 2=green
set palette defined (0 "steelblue", 1 "dark-orange", 2 "#228B22")
set cbrange [0:2]
set cbtics ("phase 0" 0, "phase 1" 1, "phase 2" 2)

set title "Step 0 — Initial (3 Phases)"
set xlabel "x"
set ylabel "z"
plot "reseed_edge_init_EdgeMultiphase.dat" using 1:2:3 with points pt 7 ps 0.4 lc palette notitle

set title "Step 50 — After Reseeding"
plot "reseed_edge_final_EdgeMultiphase.dat" using 1:2:3 with points pt 7 ps 0.4 lc palette notitle

unset multiplot
