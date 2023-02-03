reset
set terminal gif animate delay 40
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

do for [i=0:9] {
    set contour
    unset surface
    set cntrparam levels incr -3,0.5,1
    set key center tmargin
    set style textbox noborder
    set dgrid3d 50,50,1
    set table $COUNTOUR
    splot data u 1:2:4 index i
    unset table
    unset key
    plot data index i with image title columnheader(1), $COUNTOUR w l lc rgb "white", '' u 1:2:3 every 50 w labels boxed tc rgb "red" font "Arial,12,bold"
}
