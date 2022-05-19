reset
unset table
unset key
set xrange [-33000:33000]
set yrange [0:31000]
set ylabel "z (m)"
set xlabel "y (m)"
set size ratio -1

if (!exists("filename")) filename='../VISUAL_TESTS/img/ShearHeatingDuretz14.png'
set terminal png size 1000,1000;
set output filename;
set multiplot layout 2,1
set title 'ShearHeatingDuretz14 step 5: eII'
plot 'ShearHeatingDuretz14.dat' u 1:2:3 with image
set title 'ShearHeatingDuretz14 step 5: sII'
plot 'ShearHeatingDuretz14.dat' u 1:2:4 with image
unset multiplot