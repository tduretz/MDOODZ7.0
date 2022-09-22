reset
unset table
unset key
set xrange [-150000:150000]
set yrange [-180000:10000]
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"
set size ratio -1
load 'turbo.pal'

if (!exists("filename")) filename='../VISUAL_TESTS/img/RiftingChenin.png'
set terminal png size 1000,1000;
set output filename;
set multiplot layout 2,1
set cblabel "log_{10} {/:Italic ε}II"
set cbrange [-17:-14]
plot 'RiftingChenin.dat' u 1:2:3 with image
set cblabel "log_{10} {/:Italic σ}II"
set cbrange [5:9]
plot 'RiftingChenin.dat' u 1:2:4 with image
unset multiplot