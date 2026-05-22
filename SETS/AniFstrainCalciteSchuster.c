// AniFstrainCalciteSchuster — homogeneous simple-shear calibration set for the
// Schuster, Habler, Schafler & Abart 2018 (J. Struct. Geol.;
// doi:10.1016/j.jsg.2018.09.003) Carrara marble high-pressure-torsion (HPT)
// experiments at 450°C, 1-4 GPa.
//
// Annotation-only data — Schuster+18 reports γ-coverage in Table 2 (16 EBSD
// maps spanning γ ∈ [3, 63]) but does NOT tabulate J / M / peak-m.r.d. for
// any sample. The MDOODZ run + analytical curve are plotted alongside
// Schuster's γ-rug as a data-coverage annotation.
//
// Uses aniso_db = 4 (case 4, anisoDelta_Calcite_LowT, mid-low-T regime
// fitted to Barnhoorn+04 ≤650°C: M_∞ = 0.22, γ_e = 1.37, slope = 6.0,
// asymptote δ_∞ = 2.32). Schuster's 450°C is closer to the LowT regime than
// to Pieri's 727°C, but the regime caveat (HPT 1-4 GPa is off-spec for
// MDOODZ lithospheric calibration) applies — see annotation memo §"Regime
// caveat" for the full discussion.
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
    RunMDOODZ("AniFstrainCalciteSchuster.txt", &setup);
}
