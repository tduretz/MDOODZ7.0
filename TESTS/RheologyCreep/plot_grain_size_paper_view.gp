# Grain Size Steady-State Benchmark — Schmalholz & Duretz 2017 Fig. 2 layout
#
# Same physics as plot_grain_size.gp but in the paper's chosen axes:
#   x = grain size d (µm), y = differential stress 2*tau_II (MPa), both log
#
# Reproduces the deformation-map view from S&D 2017 Fig. 2 for visual paper-
# comparison. Renner 2002 dislocation + Herwegh 2003 diffusion + Austin & Evans 2007/2009
# wattmeter, T = 350 °C — all matching the paper's "Renner 2002 & Herwegh 2003"
# panel.
#
# Run from cmake-build-test/TESTS:
#     gnuplot ../../TESTS/RheologyCreep/plot_grain_size_paper_view.gp
#
# Output: grain_size_paper_view.png

set terminal pngcairo size 900,700 enhanced font "Arial,12"
set output "grain_size_paper_view.png"

set encoding utf8
set title "Calcite deformation map at T = 350 C - Schmalholz & Duretz 2017 Fig. 2 layout" noenhanced

# Paper has stress on y, grain size on x — flipped vs our other plot
set xlabel "Grain size d (µm)" noenhanced
set ylabel "Differential stress 2·{/Symbol t}_{II} (MPa)"

set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"
set xrange [1:3000]
set yrange [1:1000]

# Paper-style ticks: minor tics on log axes (at 2, 3, ..., 9 per decade), all
# pointing inward, mirrored on top and right axes. Matches Schmalholz & Duretz
# 2017 Fig. 2. `set tics scale` makes major/minor tics more visible than gnuplot's
# pngcairo default (which renders them very small).
set tics in scale 2.0, 1.0
set xtics mirror
set ytics mirror
set mxtics 10
set mytics 10

# Full grid at both major and minor tick positions — paper Fig. 2 style.
# Major gridlines slightly darker, minor gridlines faint. Both are drawn behind
# the data (`back`).
set style line 81 linetype 1 linecolor rgb "#bbbbbb" linewidth 0.6
set style line 82 linetype 1 linecolor rgb "#e0e0e0" linewidth 0.4
set grid xtics ytics back linestyle 81
set grid mxtics mytics back linestyle 82
set key bottom right

# Regime labels — same positions as Schmalholz & Duretz 2017 Fig. 2
set label "Diffusion creep" at 4, 4 rotate by 70 font ",13" textcolor rgb "#222222" front
set label "Dislocation creep" at 400, 70 font ",13" textcolor rgb "#222222" front

# Iso-contour labels placed directly on the rising legs (paper style — no legend entries)
set label "10^{-10} s^{-1}" at 1.6, 2.5  font ",10" textcolor rgb "#444444" front
set label "10^{-12} s^{-1}" at 6,   2.5  font ",10" textcolor rgb "#444444" front
set label "10^{-14} s^{-1}" at 27,  2.5  font ",10" textcolor rgb "#444444" front

# Layout signature
set label "Layout: Schmalholz & Duretz (2017), J. Struct. Geol. 103, Fig. 2" \
    at screen 0.99, screen 0.012 right font ",8" textcolor rgb "#666666" noenhanced front

# .dat is in SI (m, Pa). Convert to (µm, MPa) inline:
#   d in µm  =  d_SI * 1e6
#   sigma in MPa = 2 * tau_SI / 1e6     (paper plots differential stress sigma = 2*tau_II)

# Block 0: paleowattmeter samples — coupled (Renner 2002 + Herwegh 2003 + wattmeter)
# Block 1: 5 MDOODZ measurement points
# Block 2: iso-contour Eii = 1e-10 s^-1
# Block 3: iso-contour Eii = 1e-12 s^-1
# Block 4: iso-contour Eii = 1e-14 s^-1

plot "grain_size_benchmark_coupled.dat" index 2 using ($1*1e6):($2*2/1e6) \
        with lines linewidth 1.0 linecolor rgb "#444444" notitle, \
     "grain_size_benchmark_coupled.dat" index 3 using ($1*1e6):($2*2/1e6) \
        with lines linewidth 1.0 linecolor rgb "#444444" notitle, \
     "grain_size_benchmark_coupled.dat" index 4 using ($1*1e6):($2*2/1e6) \
        with lines linewidth 1.0 linecolor rgb "#444444" \
        title "Coupled creep iso-contours ({/Symbol e}'_{II} labelled inline)", \
     "grain_size_benchmark_coupled.dat" index 0 using ($3*1e6):($2*2/1e6) \
        with lines linewidth 3.0 linecolor rgb "#000000" \
        title "Paleowattmeter: Renner 2002 and Herwegh 2003", \
     "grain_size_benchmark_coupled.dat" index 1 using ($3*1e6):($2*2/1e6) \
        with points pointtype 9 pointsize 2.0 linecolor rgb "#d62728" \
        title "MDOODZ — coupled sweep (linv = 15)"

set output
