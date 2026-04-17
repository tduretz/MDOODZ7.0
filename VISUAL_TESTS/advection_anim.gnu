# advection_anim.gnu — Animated GIF: particle drift under rigid-body rotation
# Heatmap of per-cell particle count shows depletion (red) and pileup (blue)
# Uses Mode 0 data where Nb_part grows from 6400 to 8716 (no deactivation)

if (!exists("outdir")) outdir = "."

set terminal gif animate delay 30 loop 0 size 600,550 enhanced font "Arial,12"
set output sprintf("%s/advection_anim.gif", outdir)

# Diverging palette centred on 4 (expected count): blue=depleted, red=pileup
set palette defined ( 0 "#2166AC", 1 "#4393C3", 2 "#92C5DE", \
                      3 "#D1E5F0", 4 "#FFFFBF", \
                      5 "#FEE08B", 6 "#FDAE61", \
                      7 "#F46D43", 8 "#D73027", 9 "#A50026" )
set cbrange [0:9]
set cblabel "Particles / fine cell" font "Arial,10"

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set size ratio 1
set xlabel "x"
set ylabel "z"
unset key

# 51 frames: steps 0, 1, 2, ..., 50
do for [step=0:50:1] {
  set title sprintf("Particle Drift (Mode 0) — Step %d / 50\nBlue = depleted cells, Red = overpopulated cells", step) font "Arial,12"
  gridname = sprintf("reseed_grid_Mode0_%03d.dat", step)
  plot gridname using 1:2:3 with image notitle
}

unset output
