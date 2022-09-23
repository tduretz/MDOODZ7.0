reset
set contour
unset surface
set cntrparam levels incr 0,1,4
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"
set cblabel "({/:Italic P}) [Pa]"
set title 'P map + Tii contour'
set dgrid3d 100,100,1
set table $COUNTOUR
splot data u 1:2:4
unset table
unset key
if (!exists("filename")) filename='../VISUAL_TESTS/img/ShearTemplate.png'
set terminal png size 800,600;
set output filename;
set style textbox noborder
load 'inferno.pal'
plot data with image, data using 1:2:($5/10):($6/10) with vectors lc -1 filled, $COUNTOUR w l lc rgb "white", '' u 1:2:3 every 50 w labels boxed tc rgb "white"