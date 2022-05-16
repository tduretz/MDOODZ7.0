reset
set contour
unset surface
set cntrparam levels incr 0,1,4
set xrange [-0.5:0.5]
set yrange [-0.5:0.5]
set dgrid3d 100,100,1
set table $COUNTOUR
splot 'ShearTemplate.dat' u 1:2:4
unset table
unset key
set style textbox noborder
load 'inferno.pal'
plot 'ShearTemplate.dat' with image, 'ShearTemplate.dat' using 1:2:($5/10):($6/10) with vectors lc -1 filled, $COUNTOUR w l lc rgb "white", '' u 1:2:3 every 50 w labels boxed tc rgb "white"