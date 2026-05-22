// AniFstrainCalcitePieri — homogeneous simple-shear calibration set for the
// Pieri+01 (Tectonophysics 330, 119) Carrara marble torsion dataset at 727°C.
//
// Reproduces the lab-data corpus γ ∈ {0, 1, 2, 5, 11} at 727°C and produces
// HDF5 outputs that mdoodz_vs_pieri.py overlays on Pieri's J-converted-to-δ
// datapoints. Uses aniso_db = 5 (case 5, anisoDelta_Calcite_HighT) which is
// the LSQ-fitted high-T saturating-exp form: M(γ) = 0.73·(1 − exp(−γ/5.28)),
// δ = 6·M + 1, asymptote δ_∞ = 5.40.
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
    RunMDOODZ("AniFstrainCalcitePieri.txt", &setup);
}
