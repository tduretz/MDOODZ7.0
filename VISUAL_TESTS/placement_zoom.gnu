# placement_zoom.gnu — Close-up of ONLY reseeded particles (generation==1)
# Makes the centroid-grid pattern (mode 1) vs random scatter (mode 2) obvious

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1600,700 enhanced font "Arial,12"
set output sprintf("%s/placement_zoom.png", outdir)

set multiplot layout 1,3 title "Reseeded Particle Placement: Centroid Grid vs Random" font "Arial,14"

Ncx = 20
dx = 1.0 / (2.0 * Ncx)
xmin_g = -0.5
zmin_g = -0.5

# Zoom to bottom-right corner where rotation depletes cells
xlo = 0.2;  xhi = 0.5
zlo = -0.5; zhi = -0.2

set xrange [xlo:xhi]
set yrange [zlo:zhi]
set size ratio 1

# Draw fine-mesh grid
set style line 10 lc rgb "#E0E0E0" lw 0.5
do for [i=0:2*Ncx] {
  xv = xmin_g + i*dx
  if (xv >= xlo && xv <= xhi) {
    set arrow from xv, zlo to xv, zhi nohead ls 10 front
  }
}
do for [j=0:2*Ncx] {
  zv = zmin_g + j*dx
  if (zv >= zlo && zv <= zhi) {
    set arrow from xlo, zv to xhi, zv nohead ls 10 front
  }
}

set key top left font "Arial,10"

# --- Mode 0: closest-neighbour ---
set title "Mode 0 — Closest-Neighbour" font "Arial,12"
set xlabel "x"
set ylabel "z"
plot "reseed_all_final_Mode0.dat" \
       using ($4==0 ? $1 : 1/0):($4==0 ? $2 : 1/0) \
       with points pt 7 ps 0.5 lc rgb "#CCCCCC" title "original", \
     "reseed_all_final_Mode0.dat" \
       using ($4==1 ? $1 : 1/0):($4==1 ? $2 : 1/0) \
       with points pt 7 ps 2.0 lc rgb "#E31A1C" title "reseeded"

# --- Mode 1: centroid placement ---
set title "Mode 1 — Centroid (Grid-Aligned)" font "Arial,12"
plot "reseed_all_final_Mode1.dat" \
       using ($4==0 ? $1 : 1/0):($4==0 ? $2 : 1/0) \
       with points pt 7 ps 0.5 lc rgb "#CCCCCC" title "original", \
     "reseed_all_final_Mode1.dat" \
       using ($4==1 ? $1 : 1/0):($4==1 ? $2 : 1/0) \
       with points pt 7 ps 2.0 lc rgb "#E31A1C" title "reseeded"

# --- Mode 2: random placement ---
set title "Mode 2 — Random Within Cell" font "Arial,12"
plot "reseed_all_final_Mode2.dat" \
       using ($4==0 ? $1 : 1/0):($4==0 ? $2 : 1/0) \
       with points pt 7 ps 0.5 lc rgb "#CCCCCC" title "original", \
     "reseed_all_final_Mode2.dat" \
       using ($4==1 ? $1 : 1/0):($4==1 ? $2 : 1/0) \
       with points pt 7 ps 2.0 lc rgb "#E31A1C" title "reseeded"

unset multiplot
