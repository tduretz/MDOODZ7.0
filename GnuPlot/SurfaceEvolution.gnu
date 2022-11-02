set terminal png size 600,600;
set output 'SurfaceEvolution.png';
set ylabel "height [m]"
set xlabel "{/:Italic t} [kyr]"
plot 'SurfaceEvolution.dat'