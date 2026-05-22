// AniFstrainQuartzPennacchioni — homogeneous simple-shear calibration set for the
// Pennacchioni+10 (JGR 115, B12405) natural-shear-zone quartz dataset at ~500°C.
//
// Reproduces the audit-suspect 17 (γ, J) corpus γ ∈ [0.3, 15] under aniso_db = 3
// (case 3, anisoDelta_Quartz, M_∞ = 0.75, γ_e = 3.50, slope = 3.0, δ_∞ = 3.25).
// The 17 J values are NOT real per-sample MTEX data (Pennacchioni+10 tabulates
// no J-index); residual smoothness σ/M_∞ = 0.045 is 3–4× below the real-data
// floor (Hansen 0.142, Blackford 0.191), the signature of a hand-fitted
// saturating curve. The 17 datapoints are still useful as a calibration target,
// but their defensibility is "physical-prior", not "data-driven" — see
// misc/aniso_fstrain/audit/Pennacchioni_etal_2010_audit.md.
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
    RunMDOODZ("AniFstrainQuartzPennacchioni.txt", &setup);
}
