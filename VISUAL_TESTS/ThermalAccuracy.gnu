# ThermalAccuracy.gnu — Grouped bar chart of L2 errors (CHOLMOD vs PCG)
set terminal pngcairo enhanced size 800,500 font "Helvetica,14"
set output filename

set style data histogram
set style histogram clustered gap 1
set style fill solid 0.7 border -1
set boxwidth 0.8

set ylabel "Relative L2 error"
set xlabel ""
set logscale y
set format y "10^{%L}"
set yrange [1e-8:1e-1]
set grid ytics

set key top right

# data format: "Test" "Solver" L2_error
# We plot CHOLMOD (rows 1,3) and PCG (rows 2,4) as grouped bars

set title "Thermal Solver Accuracy: CHOLMOD vs PCG"

plot data using 3:xtic(1) every 2::0 title "CHOLMOD" lc rgb "#2166ac", \
     data using 3:xtic(1) every 2::1 title "PCG" lc rgb "#b2182b"
