# closed_box_crash_anim.gnu — Animated GIF showing closed-box convection
# Mode 0 and 1 grow unbounded and crash (exit 190);
# Mode 2 stays stable thanks to deactivation.
# Panel titles show Nb_part / Nb_part_max; crashed modes show "CRASHED".

if (!exists("outdir")) outdir = "."

set terminal gif animate delay 15 loop 0 size 1600,550 enhanced font "Arial,11"
set output sprintf("%s/closed_box_crash_anim.gif", outdir)

# Diverging palette: blue=depleted, yellow=nominal(4), red=overpopulated
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

# Panel layout
pw = 0.28
x0 = 0.05
x1 = x0 + pw + 0.02
x2 = x1 + pw + 0.02
py = 0.12
ph = 0.72

# Initial Nb_part = 10*10*16 = 1600; Nb_part_max = ceil(4.1 * 1600) = 6560
nb_max = 6560
total_steps = 3000
wstep = 30

# Detect crash: read the last step available for each mode
last0 = system("tail -1 reseed_nbpart_ClosedBoxMode0.dat | awk '{print $1}'") + 0
last1 = system("tail -1 reseed_nbpart_ClosedBoxMode1.dat | awk '{print $1}'") + 0
last2 = system("tail -1 reseed_nbpart_ClosedBoxMode2.dat | awk '{print $1}'") + 0

# Read the final Nb_part at crash point for each mode
crashnb0 = system(sprintf("awk '$1==%d {print $2}' reseed_nbpart_ClosedBoxMode0.dat", last0)) + 0
crashnb1 = system(sprintf("awk '$1==%d {print $2}' reseed_nbpart_ClosedBoxMode1.dat", last1)) + 0

do for [frame=0:100] {
  step = frame * wstep

  # --- Read Nb_part for each mode ---
  nb0_str = system(sprintf("awk '$1==%d {print $2}' reseed_nbpart_ClosedBoxMode0.dat", step))
  nb1_str = system(sprintf("awk '$1==%d {print $2}' reseed_nbpart_ClosedBoxMode1.dat", step))
  nb2_str = system(sprintf("awk '$1==%d {print $2}' reseed_nbpart_ClosedBoxMode2.dat", step))

  nb0 = (strlen(nb0_str) > 0) ? nb0_str + 0 : 1600
  nb1 = (strlen(nb1_str) > 0) ? nb1_str + 0 : 1600
  nb2 = (strlen(nb2_str) > 0) ? nb2_str + 0 : 1600

  pct0 = 100.0 * nb0 / nb_max
  pct1 = 100.0 * nb1 / nb_max
  pct2 = 100.0 * nb2 / nb_max

  # Determine if mode has crashed by this step
  crashed0 = (step > last0) ? 1 : 0
  crashed1 = (step > last1) ? 1 : 0
  crashed2 = (step > last2) ? 1 : 0

  set multiplot title sprintf("Closed-Box Convection — Step %d / %d", step, total_steps) font "Arial,13"

  # --- Mode 0 panel ---
  unset colorbox
  set lmargin at screen x0
  set rmargin at screen x0 + pw
  set bmargin at screen py
  set tmargin at screen py + ph
  if (crashed0) {
    set title sprintf("Mode 0 — CRASHED at step %d (%d / %d)", last0, int(crashnb0), nb_max) font "Arial,10" tc rgb "red"
    # Show last available frame frozen
    f0 = sprintf("reseed_grid_ClosedBoxMode0_%03d.dat", last0)
  } else {
    set title sprintf("Mode 0 — %d / %d (%.0f%%)", int(nb0), nb_max, pct0) font "Arial,10" tc rgb "black"
    f0 = sprintf("reseed_grid_ClosedBoxMode0_%03d.dat", step)
  }
  if (system(sprintf("test -f %s && echo 1 || echo 0", f0)) + 0) {
    plot f0 using 1:2:3 with image notitle
  }

  # --- Mode 1 panel ---
  set lmargin at screen x1
  set rmargin at screen x1 + pw
  if (crashed1) {
    set title sprintf("Mode 1 — CRASHED at step %d (%d / %d)", last1, int(crashnb1), nb_max) font "Arial,10" tc rgb "red"
    f1 = sprintf("reseed_grid_ClosedBoxMode1_%03d.dat", last1)
  } else {
    set title sprintf("Mode 1 — %d / %d (%.0f%%)", int(nb1), nb_max, pct1) font "Arial,10" tc rgb "black"
    f1 = sprintf("reseed_grid_ClosedBoxMode1_%03d.dat", step)
  }
  if (system(sprintf("test -f %s && echo 1 || echo 0", f1)) + 0) {
    plot f1 using 1:2:3 with image notitle
  }

  # --- Mode 2 panel (with colorbox) ---
  set colorbox vertical user origin (x2 + pw + 0.015), py size 0.015, ph
  set cblabel "Particles / fine cell" font "Arial,10"
  set lmargin at screen x2
  set rmargin at screen x2 + pw
  if (crashed2) {
    set title sprintf("Mode 2 — CRASHED at step %d", last2) font "Arial,10" tc rgb "red"
    f2 = sprintf("reseed_grid_ClosedBoxMode2_%03d.dat", last2)
  } else {
    set title sprintf("Mode 2 — %d / %d (%.0f%%)", int(nb2), nb_max, pct2) font "Arial,10" tc rgb "black"
    f2 = sprintf("reseed_grid_ClosedBoxMode2_%03d.dat", step)
  }
  if (system(sprintf("test -f %s && echo 1 || echo 0", f2)) + 0) {
    plot f2 using 1:2:3 with image notitle
  }

  unset multiplot
}

unset output
