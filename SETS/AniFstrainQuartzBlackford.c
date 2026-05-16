// AniFstrainQuartzBlackford — homogeneous simple-shear calibration set for the
// Blackford+24 (Tectonics 43, e2023TC008166) Northern Snake Range quartzite
// dataset at T ~ 500-650°C.
//
// Reproduces the 38-sample γ_eq ∈ [1.28, 5.08] corpus under aniso_db = 6
// (case 6, anisoDelta_Quartz_Blackford, M_∞ = 0.20, γ_e = 4.00, slope = 11.0,
// δ_∞ = 3.20). Slope is rescaled relative to case 3 (3.0 → 11.0) to compensate
// for the smaller M values from pole-figure J vs Pennacchioni's ODF J — both
// modes target the same physical δ_∞ ≈ 3.2 for natural quartz mylonite.
//
// γ_e = 4.0 is centered in the LSQ envelope (saturation is poorly constrained
// by the Blackford γ-range [1.3, 5.1]) — see calibrate/calibrate_quartz.py
// and the FlowLaws.c case-6 docstring. The MDOODZ-vs-datapoint plot shows
// pfJ-converted M values bracketing the saturating-exp curve over the
// available γ range; the curve continues to saturate beyond Blackford's
// upper-γ coverage.
//
// One phase, ani_fstrain = 2, no inclusion, periodic-in-x simple shear.

#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    return 0;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
}

int main() {
    MdoodzSetup setup = {
            .SetParticles = &(SetParticles_ff){
                    .SetPhase   = SetPhase,
                    .SetDensity = SetDensity,
            },
            .SetBCs = &(SetBCs_ff){
                    .SetBCVx = SetSimpleShearBCVx,
                    .SetBCVz = SetSimpleShearBCVz,
            },
    };
    RunMDOODZ("AniFstrainQuartzBlackford.txt", &setup);
}
