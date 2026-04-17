# reseed_placement.gnu — 3-panel: reseeded particles (generation==1) for modes 0/1/2
# Highlights the different placement strategies

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1600,500 enhanced font "Arial,11"
set output sprintf("%s/reseed_placement.png", outdir)

set multiplot layout 1,3 title "Reseeding Placement Strategy Comparison" font "Arial,13"

# Fine-mesh grid overlay
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

# Color: generation==0 grey, generation==1 red
set palette defined (0 "grey70", 1 "red")
set cbrange [0:1]
unset colorbox

set title "Mode 0 (CountPartCell\\_OLD)"
set xlabel "x"
set ylabel "z"
plot "reseed_all_final_Mode0.dat" using 1:2:4 with points pt 7 ps 0.3 lc palette notitle

set title "Mode 1 (CountPartCell) — Centroid"
plot "reseed_all_final_Mode1.dat" using 1:2:4 with points pt 7 ps 0.3 lc palette notitle

set title "Mode 2 (CountPartCell\\_v2) — Random"
plot "reseed_all_final_Mode2.dat" using 1:2:4 with points pt 7 ps 0.3 lc palette notitle

unset multiplot
