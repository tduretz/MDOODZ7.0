# set title "Dev. stress components"
# set xlabel "Angle"
# set ylabel "Stress magnitude"
# set key autotitle columnheader
# set style data lines
# set style line 1 lw 1 lc rgb "blue"
# set style line 2 lw 1 lc rgb "green"
# set style line 3 pt 1 ps 1.0 lc rgb "red"
# set style line 4 pt 2 ps 1.0 lc rgb "black"
# set grid x y mx my
# set grid lt -1 lc "black" lw 1 , lt -1 lc bgnd lw 16
# p for [col=2:5] "data4plot.txt" u 1:col w lp

### lines and points in the same plot for-loop
reset session

LINECOLORS = "black blue  green red  black"
LINEWIDTHS = '0     2.0   2.0   2.0   2.0    '
POINTTYPES = '0     0     0     0     0      '
POINTSIZES = '0     1     1     1.0   1.0    '
TITLES     = 'angle Txx   Tzz   Txz   Tii   '

myLinecolor(i) = word(LINECOLORS,i)
myLinewidth(i) = real(word(LINEWIDTHS,i))
myPointtype(i) = int(word(POINTTYPES,i))
myPointsize(i) = real(word(POINTSIZES,i))
myLinetype(i) = myLinewidth(i) == 0 ? -2 : 1
myTitle(i) = word(TITLES,i)

set samples 31
set key out
set grid lt -1 lc "gray" lw 1 , lt -1 lc bgnd lw 16

# p for [col=2:5] "data4plot.txt" u 1:col w lp
plot for [i=2:words(TITLES)]  "data4plot.txt" u 1:i w lp pt myPointtype(i) ps myPointsize(i) \
    lt myLinetype(i) lw myLinewidth(i) lc rgb myLinecolor(i) title myTitle(i)
### end of code