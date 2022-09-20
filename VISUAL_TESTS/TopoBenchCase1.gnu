set terminal png size 600,600;
set output '../VISUAL_TESTS/img/TopoBenchCase1.png';
set ylabel "height [m]"
set xlabel "{/:Italic t} [Kyr]"
plot 7e3*exp(-0.2139e-11*(x*31556952000000)) title "Analytical - Crameri et al. (2012)", 'TopoBenchCase1.dat' title "MDOODZ TopoBenchCase1"