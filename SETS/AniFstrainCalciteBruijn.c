// AniFstrainCalciteBruijn — homogeneous simple-shear calibration set for the
// Bruijn, Burlini & Drury 2011 (Tectonophysics 503, 75; doi:10.1016/j.tecto.2010.09.029)
// Carrara marble torsion experiments at 727°C.
//
// ⚠ Audit caveat — the 3 (γ, J) values stored at
// misc/aniso_fstrain/data/calcite/raw_Bruijn_etal_2011_Tectonophysics.csv
// CANNOT BE SOURCED FROM THE PUBLISHED PAPER (audit memo embedded in the CSV
// header; Bruijn's Fig. 8 left column reproduces Pieri+01 / Barnhoorn+04 pole
// figures with no numerical J labels). The values appear to be visual greyscale
// estimates. This SET therefore doubles as an **audit-control** mirroring
// AniFstrainQuartzPennacchioni: the MDOODZ-vs-datapoint plot visualises that
// the Bruijn-attributed points sit far enough off case 2 to flag the audit
// finding. Useful for the methodology paper §S1.4 audit-trail story.
//
// Uses aniso_db = 2 (case 2, anisoDelta_Calcite, T-averaged Carrara Paterson
// rig: M_∞ = 0.41, γ_e = 3.18, slope = 6.0, asymptote δ_∞ = 3.46).
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
    RunMDOODZ("AniFstrainCalciteBruijn.txt", &setup);
}
