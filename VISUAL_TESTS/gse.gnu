set terminal gif animate delay 40
set output filename
set xrange [-0.4:0.4]
set yrange [-0.17:0.17]
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"
set cblabel "log_{10} ({/:Italic d}) [μm]"
set cbrange [1:3]
set rmargin at screen 0.8
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')

set key center tmargin

do for [i=0:9] {
    plot data index i with image title columnheader(1)
}
