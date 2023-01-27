reset
set terminal gif animate delay 40
set contour
unset surface
set cntrparam levels incr -3,1,0
set output filename
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"
set cblabel "{/:Italic x} [%]"
set cbrange [0:0.01]
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
    set dgrid3d 50,50,1
    set table $COUNTOUR
    splot data u 1:2:4
    unset table
    unset key
    set style textbox noborder
    plot data index i with image title columnheader(1), $COUNTOUR w l lc rgb "white"
}
