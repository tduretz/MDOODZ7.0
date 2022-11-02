set terminal gif animate delay 40 size 1300,1000
set output filename
set ylabel "{/:Italic z} [m]"
set xlabel "{/:Italic x} [m]"


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

do for [i=1:frames-1] {
    print i
    set multiplot layout 2,1
        set cblabel "log_{10} {/:Italic ε}II"
        set cbrange [-18:-14.5]
        plot data index i u 1:2:3 with image title columnheader(1)
        set cblabel "log_{10} {/:Italic σ}II"
        set cbrange [5:9]
        plot data index i u 1:2:4 with image
    unset multiplot
}
