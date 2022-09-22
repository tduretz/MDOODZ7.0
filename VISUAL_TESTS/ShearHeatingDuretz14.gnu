reset
unset table
unset key
set xrange [-33000:33000]
set yrange [0:31000]
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"
set size ratio -1
load 'turbo.pal'

if (!exists("filename")) filename='../VISUAL_TESTS/img/ShearHeatingDuretz14.png'
set terminal png size 800,800;
set output filename;
set multiplot layout 2,1
set cblabel "log_{10} {/:Italic ε}II"
plot 'ShearHeatingDuretz14.dat' u 1:2:3 with image
set cblabel "log_{10} {/:Italic σ}II"
plot 'ShearHeatingDuretz14.dat' u 1:2:4 with image
unset multiplot