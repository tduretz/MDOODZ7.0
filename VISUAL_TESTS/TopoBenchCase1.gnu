set terminal png size 1000,1000;
set output '../VISUAL_TESTS/img/TopoBenchCase1.png';
set ylabel "height (m)"
set xlabel "time (s)"
plot 7e3*exp(-0.2139e-11*x) title "Analytical - Crameri et al. (2012)", 'TopoBenchCase1.dat' title "MDOODZ TopoBenchCase1"