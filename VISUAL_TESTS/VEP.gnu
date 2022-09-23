set terminal gif animate delay 40 size 800,400
set output filename
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"
set cblabel "{/:Italic P} [MPa]" offset 1
set cbrange [-6:4]
set xrange [-2000:2000]
set yrange [-1000:1000]
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

do for [i=0:10] {
    plot data index i with image title columnheader(1)
}
