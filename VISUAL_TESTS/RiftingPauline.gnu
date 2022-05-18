reset
unset table
unset key
set xrange [-150000:150000]
set yrange [-180000:10000]
set size ratio -1

if (!exists("filename")) filename='../VISUAL_TESTS/img/RiftingPauline.png'
set terminal png size 1000,1000;
set output filename;
set multiplot layout 2,1
set title 'RiftingPauline step 50: eII'
plot 'RiftingPauline.dat' u 1:2:3 with image
set title 'RiftingPauline step 50: sII'
plot 'RiftingPauline.dat' u 1:2:4 with image
unset multiplot