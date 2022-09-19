set terminal gif animate delay 40
set output filename
set cbrange [-12:-8.5]
set xrange [-0.4:0.4]
set yrange [-0.17:0.17]
set ylabel "z (m)"
set xlabel "x (m)"
set cblabel "grain size (log)"
set rmargin at screen 0.8

set key center tmargin

do for [i=1:10] {
    plot data index i with image title columnheader(1)
}
