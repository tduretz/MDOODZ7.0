# deactivation_strategy.gnu — 2-panel: mode 1 vs mode 2 deactivation
# Uses high-density (6×6) data where mode 2 deactivation is dramatic

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1200,500 enhanced font "Arial,11"
set output sprintf("%s/deactivation_strategy.png", outdir)

set multiplot layout 1,2 title "Deactivation at High Density (6×6 / cell): Blue = Active, Red = Deactivated" font "Arial,13"

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

set key top right

set title "Mode 1 — No Deactivation"
set xlabel "x"
set ylabel "z"
plot "reseed_hd_final_HighDensMode1.dat" using ($3!=-1 ? $1 : 1/0):($3!=-1 ? $2 : 1/0) with points pt 7 ps 0.2 lc rgb "steelblue" title "active", \
     "reseed_hd_final_HighDensMode1.dat" using ($3==-1 ? $1 : 1/0):($3==-1 ? $2 : 1/0) with points pt 7 ps 0.5 lc rgb "red" title "deactivated"

set title "Mode 2 — Distance-Based Deactivation"
plot "reseed_hd_final_HighDensMode2.dat" using ($3!=-1 ? $1 : 1/0):($3!=-1 ? $2 : 1/0) with points pt 7 ps 0.2 lc rgb "steelblue" title "active", \
     "reseed_hd_final_HighDensMode2.dat" using ($3==-1 ? $1 : 1/0):($3==-1 ? $2 : 1/0) with points pt 7 ps 0.5 lc rgb "red" title "deactivated"

unset multiplot
