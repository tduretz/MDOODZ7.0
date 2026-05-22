#!/usr/bin/env gnuplot
#
# Finite-strain anisotropy factor evolution: MDOODZ vs analytical reference.
#
# δ = min(FS_AR, ani_fac_max) where FS_AR is the strain-ellipse aspect ratio
# from the polar decomposition of the deformation gradient. Under simple
# shear with shear rate γ̇:
#
#   FS_AR(γ) = 1 + γ²/2 + γ·sqrt(γ²/4 + 1)            (closed form, derived
#                                                      in TESTS/AnalyticalSolutions.md)
#
# Three curves are drawn:
#   - solid  black: unbounded FS_AR(γ) curve (no clamp)
#   - dashed black: clamped δ(γ) = min(FS_AR, ani_fac_max)  ← MDOODZ form
#   - red triangles: MDOODZ measurements (mean of Centers/ani_fac per step)
#
# Reference markers:
#   - grey dashed horizontal: y = ani_fac_max (saturation level)
#   - grey dashed vertical:   x = γ_sat (saturation strain, where FS_AR=ani_fac_max)
#
# The MDOODZ markers should sit on the dashed black line at every step. The
# kink at (γ_sat, ani_fac_max) is the visual signature of the min() clamp.
#
# Prerequisites:
#   cd build/TESTS
#   ./AnisotropyBenchmarkTests --gtest_filter="*AniFstrainSimpleShear"
#   (generates aniso_factor_evolution.dat)
#
# Usage:
#   gnuplot ../../TESTS/AnisotropyBenchmark/plot_aniso_factor.gp
#
# Output: aniso_factor_evolution.png

reset
set terminal pngcairo size 900,600 enhanced font "Arial,13"
set output "aniso_factor_evolution.png"

# --- Reference parameters (must match the test) ---
ani_fac_max_test = 4.0   # not the value used in the .dat, just for the
                         # illustrative clamp line

# Closed-form FS_AR for simple shear
fsar(g) = 1.0 + 0.5*g*g + g*sqrt(0.25*g*g + 1.0)
# Inverse: γ at which FS_AR(γ) = D, derived from FS_AR² - 1 = γ²·FS_AR
gamma_sat(D) = sqrt(D*D - 1.0) / sqrt(D)

# Clamped form δ(γ) = min(FS_AR, ani_fac_max)
delta_clamped(g, D) = (fsar(g) < D) ? fsar(g) : D

g_sat = gamma_sat(ani_fac_max_test)

set xlabel "Shear strain γ" font ",13"
set ylabel "Anisotropy factor δ" font ",13"
set title "Finite-strain anisotropy: δ = min(FS\\_AR(γ), ani\\_fac\\_max)" font ",15"
set key left top font ",11" spacing 1.3
set xrange [-0.2:5.2]
set yrange [0.8:30]
set logscale y
set grid
set format y "%g"

set arrow from 0,ani_fac_max_test to 5,ani_fac_max_test \
    nohead linewidth 1 dashtype 3 linecolor rgb "#888888" front
set arrow from g_sat,0.8 to g_sat,30 \
    nohead linewidth 1 dashtype 3 linecolor rgb "#888888" front
set label "ani\\_fac\\_max = 4" at 0.1,4.6 textcolor rgb "#888888" font ",11" front
set label sprintf("γ_{sat} ≈ %.2f", g_sat) at g_sat+0.08,0.95 \
    textcolor rgb "#888888" font ",11" front

plot fsar(x) with lines linewidth 2.0 dashtype 1 linecolor rgb "#000000" \
        title "FS\\_AR(γ) = 1 + γ²/2 + γ√(γ²/4+1)  (unbounded)", \
     delta_clamped(x, ani_fac_max_test) with lines linewidth 2.5 dashtype 2 \
        linecolor rgb "#000000" \
        title "δ = min(FS\\_AR, 4)  (MDOODZ form)", \
     "aniso_factor_evolution.dat" using 3:6 with points pointtype 9 \
        pointsize 1.6 linecolor rgb "#d62728" \
        title "MDOODZ — mean Centers/ani\\_fac"
set output

print "Generated: aniso_factor_evolution.png"
