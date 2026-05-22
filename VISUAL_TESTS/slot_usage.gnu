# slot_usage.gnu — Array index usage: mode 1 vs mode 2 at high density
# Plots particle index vs position, showing how mode 2 recycles dead slots

if (!exists("outdir")) outdir = "."

set terminal pngcairo size 1200,500 enhanced font "Arial,11"
set output sprintf("%s/slot_usage.png", outdir)

set multiplot layout 1,2 title "Array Slot Usage (6×6 / cell): Mode 1 vs Mode 2" font "Arial,13"

set yrange [-0.5:0.5]
set key top right

# Three colors: active-original (blue), active-reseeded (orange), dead (light red)

set title "Mode 1 — Slot Reuse (No Deactivation)"
set xlabel "Array Index"
set ylabel "z position"
plot "reseed_hd_final_HighDensMode1.dat" using ($3!=-1 && $4==0 ? $5 : 1/0):($3!=-1 && $4==0 ? $2 : 1/0) with points pt 7 ps 0.15 lc rgb "steelblue" title "original", \
     "reseed_hd_final_HighDensMode1.dat" using ($3!=-1 && $4==1 ? $5 : 1/0):($3!=-1 && $4==1 ? $2 : 1/0) with points pt 7 ps 0.4 lc rgb "#FF7F00" title "reseeded", \
     "reseed_hd_final_HighDensMode1.dat" using ($3==-1 ? $5 : 1/0):($3==-1 ? $2 : 1/0) with points pt 7 ps 0.15 lc rgb "#FFCCCC" title "dead (phase=-1)"

set title "Mode 2 — Slot Recycling + Active Deactivation"
plot "reseed_hd_final_HighDensMode2.dat" using ($3!=-1 && $4==0 ? $5 : 1/0):($3!=-1 && $4==0 ? $2 : 1/0) with points pt 7 ps 0.15 lc rgb "steelblue" title "original", \
     "reseed_hd_final_HighDensMode2.dat" using ($3!=-1 && $4==1 ? $5 : 1/0):($3!=-1 && $4==1 ? $2 : 1/0) with points pt 7 ps 0.4 lc rgb "#FF7F00" title "reseeded", \
     "reseed_hd_final_HighDensMode2.dat" using ($3==-1 ? $5 : 1/0):($3==-1 ? $2 : 1/0) with points pt 7 ps 0.15 lc rgb "#FFCCCC" title "dead (phase=-1)"

unset multiplot
