# PCGConvergence.gnu — Log-scale convergence plot for PCG thermal solver
set terminal pngcairo enhanced size 800,500 font "Helvetica,14"
set output filename

set xlabel "Iteration"
set ylabel "Relative residual"
set logscale y
set format y "10^{%L}"
set grid

set title "PCG Thermal Solver Convergence"
set key top right

plot data using 1:2 with linespoints title "Step 0 (cold start)" lw 2 pt 7 ps 0.8 lc rgb "#d73027", \
     data using 1:3 with linespoints title "Step 4 (warm start)" lw 2 pt 5 ps 0.8 lc rgb "#4575b4"
