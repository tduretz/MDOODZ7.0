# advection_primer.gnu — 2-panel: before/after advection showing cell depletion
# Shows why reseeding is needed: particles drift, some cells empty out

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1200,500 enhanced font "Arial,11"
set output sprintf("%s/advection_primer.png", outdir)

set multiplot layout 1,2 title "Why Reseeding? Particle Drift Under Rigid-Body Rotation" font "Arial,13"

# Fine-mesh grid lines: 2*Ncx cells, Ncx=20, dx_fine = 1.0/40 = 0.025
Ncx = 20
dx = 1.0 / (2.0 * Ncx)
xmin = -0.5
zmin = -0.5

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set size ratio 1
set key off

# Draw fine-mesh grid
set style line 10 lc rgb "#CCCCCC" lw 0.3
do for [i=0:2*Ncx] {
  set arrow from (xmin + i*dx), -0.5 to (xmin + i*dx), 0.5 nohead ls 10 front
  set arrow from -0.5, (zmin + i*dx) to 0.5, (zmin + i*dx) nohead ls 10 front
}

set title "Step 0 — Initial Distribution"
set xlabel "x"
set ylabel "z"
plot "reseed_all_mid_Mode0.dat" using 1:2:(($3==0) ? 1 : 2) with points pt 7 ps 0.3 lc variable notitle

set title "Step 25 — After Advection (No Reseeding Shown)"
plot "reseed_all_mid_Mode1.dat" using 1:2:(($3==0) ? 1 : 2) with points pt 7 ps 0.3 lc variable notitle

unset multiplot
