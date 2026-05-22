# mode_comparison_anim.gnu — 3-panel animated GIF comparing Mode 0 vs 1 vs 2
# Heatmap of per-cell particle count: shows how each mode manages density

if (!exists("outdir")) outdir = "."

set terminal gif animate delay 30 loop 0 size 1600,500 enhanced font "Arial,11"
set output sprintf("%s/mode_comparison_anim.gif", outdir)

# Diverging palette centred on 4: blue=depleted, red=pileup
set palette defined ( 0 "#2166AC", 1 "#4393C3", 2 "#92C5DE", \
                      3 "#D1E5F0", 4 "#FFFFBF", \
                      5 "#FEE08B", 6 "#FDAE61", \
                      7 "#F46D43", 8 "#D73027", 9 "#A50026" )
set cbrange [0:9]

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set size ratio 1
unset key
set xlabel "x"
set ylabel "z"
unset colorbox

# Panel width 0.28, gap between panels, colorbox on far right
pw = 0.28
x0 = 0.05
x1 = x0 + pw + 0.02
x2 = x1 + pw + 0.02
py = 0.12
ph = 0.76

do for [step=0:50:1] {
  set multiplot title sprintf("Per-Cell Particle Count — Step %d / 50", step) font "Arial,13"

  unset colorbox
  set lmargin at screen x0
  set rmargin at screen x0 + pw
  set bmargin at screen py
  set tmargin at screen py + ph
  set title "Mode 0 (no deactivation)" font "Arial,11"
  plot sprintf("reseed_grid_Mode0_%03d.dat", step) using 1:2:3 with image notitle

  set lmargin at screen x1
  set rmargin at screen x1 + pw
  set title "Mode 1 (centroid, no deactivation)" font "Arial,11"
  plot sprintf("reseed_grid_Mode1_%03d.dat", step) using 1:2:3 with image notitle

  set colorbox vertical user origin (x2 + pw + 0.015), py size 0.015, ph
  set cblabel "Particles / fine cell" font "Arial,10"
  set lmargin at screen x2
  set rmargin at screen x2 + pw
  set title "Mode 2 (random + deactivation)" font "Arial,11"
  plot sprintf("reseed_grid_Mode2_%03d.dat", step) using 1:2:3 with image notitle

  unset multiplot
}

unset output
