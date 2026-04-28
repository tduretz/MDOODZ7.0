# Reproduces Popov et al. (2025, GMD) Fig. 5 — 0D stress integration tests.
# Reads trajectory.dat files produced by extract_popov2025 from the three 0D
# test output directories.
#
# Usage (from TESTS/Popov2025/plots/):
#   gnuplot -e 'BUILD_DIR="../../../cmake-build/TESTS"' fig5_0d_stress.gp
#
# All data comes from MDOODZ 7.0 runs with parameters matching paper
# Table 1 "0D Fig. 5 a, b" exactly; see TESTS/Popov2025/Popov0D_*.txt.

if (!exists("BUILD_DIR")) BUILD_DIR = "../../../cmake-build/TESTS"

VOL = BUILD_DIR."/Popov0D_VolumetricExtension/trajectory.dat"
DEV = BUILD_DIR."/Popov0D_DeviatoricShear/trajectory.dat"
MIX = BUILD_DIR."/Popov0D_MixedStrain/trajectory.dat"

set terminal pngcairo size 1200,940 enhanced font "Helvetica,11"
set output "out/fig5_0d_stress.png"

# ----- helpers --------------------------------------------------------------

yr  = 365.25*86400.0
# MDOODZ step 0 contains an initialization artifact for the stress invariant
# (sII ≈ 6e14 Pa before any Stokes solve has run); replace it with 0 so the
# plot starts cleanly from the origin, matching the paper.  P at step 0 is
# zero already, so no correction needed.
clean_sII(sII_val, step) = (step == 0 ? 0.0 : sII_val)

# Panel (d) plots the meridional (P, τ_II) trajectories for each test, picking
# the pressure variable that LANDS ON YIELD for that loading geometry:
#   * Shear  parks on the DP wing → P* (Centers/P, column 4) is on yield.
#   * Mixed  parks on the cap arc / DP wing → P* (column 4) is on yield.
#   * Volumetric Extension parks at the cap APEX (vertical tangent), where P*
#     overshoots the apex by K·θ̇_vp·dt ≈ 90 kPa.  For that one trajectory we
#     plot the reconstructed p_local = P* + K·θ̇_vp·dt (column 7 of
#     trajectory.dat, written by extract_popov2025) so the trajectory lands on
#     the apex as the paper's Fig 5(d) does.  Paper Eq. 31's two-field
#     formulation explicitly defines this offset between p* and p_local at the
#     vertical-tangent singularity.

# Yield-surface geometry (paper Eqs. 13–17), paper Table 1 "0D" values:
phi_deg = 30.0
psi_deg = 10.0
C_mc    = 1.0e6
p_T     = -5.0e5
kk   = sin(phi_deg*pi/180.0)
kkq  = sin(psi_deg*pi/180.0)
cc   = C_mc * cos(phi_deg*pi/180.0)
aa   = sqrt(1.0 + kk*kk)
py   = (p_T + cc/aa) / (1.0 - kk/aa)
Ry   = py - p_T
pd_  = py - Ry * kk/aa
taud = kk * pd_ + cc

# Cap arc (upper half of the circle centered at (py, 0) with radius Ry),
# parametrised by pressure:
cap_tau(P_SI) = sqrt(Ry*Ry - (P_SI - py)*(P_SI - py))
# Drucker-Prager envelope:
DP_tau(P_SI)  = kk * P_SI + cc

# ----- plot -----------------------------------------------------------------

set multiplot layout 2,2 title "Popov et al. (2025) Fig. 5 — 0D stress integration\nMDOODZ 7.0 output, grid 21 x 15, {/Symbol D}t = 2 yr, paper Table 1 \"0D\" parameters (G=10^{10}, K=2x10^{11}, {/Symbol f}=30{/Symbol \260}, {/Symbol y}=10{/Symbol \260}, c_{MC}=10^6, p_T=-5x10^5, {/Symbol h}^{vp}=0)"

set grid
set xlabel "time [yr]"
set ylabel "stress [MPa]"
# Each plotted point is one MDOODZ timestep (centre-cell value at end of step).
# The lines connect successive timesteps.

set title "(a) Volumetric extension — 20 steps (40 yr),  {/Symbol e}_{xx}={/Symbol e}_{zz}=2.333e-15 s^{-1}"
set xrange [0:20]
set yrange [-0.6:0.1]
plot VOL using ($2/yr):(clean_sII($3, $1)/1e6)  with linespoints pt 7 ps 0.5 lw 2 lc rgb "black" title "{/Symbol t}_{II}",\
     VOL using ($2/yr):($4/1e6)                  with linespoints pt 7 ps 0.5 lw 2 lc rgb "red"   title "P"

set title "(b) Deviatoric shear — 15 steps (30 yr),  {/Symbol e}_{xy}=7e-14 s^{-1}"
set xrange [0:30]
set yrange [-0.2:1.0]
plot DEV using ($2/yr):(clean_sII($3, $1)/1e6)  with linespoints pt 7 ps 0.5 lw 2 lc rgb "black" title "{/Symbol t}_{II}",\
     DEV using ($2/yr):($4/1e6)                  with linespoints pt 7 ps 0.5 lw 2 lc rgb "red"   title "P"

set title "(c) Mixed strain — 10 steps (20 yr),  volumetric + shear superposition"
set xrange [0:20]
set yrange [-0.5:1.0]
plot MIX using ($2/yr):(clean_sII($3, $1)/1e6)  with linespoints pt 7 ps 0.5 lw 2 lc rgb "black" title "{/Symbol t}_{II}",\
     MIX using ($2/yr):($4/1e6)                  with linespoints pt 7 ps 0.5 lw 2 lc rgb "red"   title "P"

# Panel (d): meridional P–τII trajectories with yield-surface envelope overlaid.
# The envelope is the cap arc for p_T ≤ P ≤ p_d joined to the DP line for P ≥ p_d.
# Explicitly set ranges BEFORE the plot command (auto-expansion is disabled so
# gnuplot doesn't shrink the x-axis to the trajectory extents).
unset autoscale
set xrange [-0.7:2.0]
set yrange [0:1.8]
set samples 400
set title "(d) Meridional P–{/Symbol t}_{II} trajectories + yield surface"
set xlabel "P [MPa]"
set ylabel "{/Symbol t}_{II} [MPa]"

# Combined envelope: cap arc for P ≤ p_d, DP line for P ≥ p_d.
envelope(P_MPa) = (P_MPa*1e6 <= pd_ ? cap_tau(P_MPa*1e6)/1e6 : DP_tau(P_MPa*1e6)/1e6)

plot envelope(x) with lines lw 2 dt 2 lc rgb "dark-grey" title "yield envelope",\
     VOL using ($7/1e6):(clean_sII($3, $1)/1e6) with linespoints pt 7 ps 0.6 lw 2 lc rgb "gold"  title "Extension (p_{local})",\
     DEV using ($4/1e6):(clean_sII($3, $1)/1e6) with linespoints pt 7 ps 0.6 lw 2 lc rgb "blue"  title "Shear",\
     MIX using ($4/1e6):(clean_sII($3, $1)/1e6) with linespoints pt 7 ps 0.6 lw 2 lc rgb "red"   title "Mixed"

unset multiplot
print "Wrote out/fig5_0d_stress.png"
