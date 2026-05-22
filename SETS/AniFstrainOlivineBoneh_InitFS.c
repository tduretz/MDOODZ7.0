// AniFstrainOlivineBoneh_InitFS — Hansen-saturation regression SET for the
// F-init helper (aniso-init-from-finite-strain). Single-phase, homogeneous
// simple shear with γ̇ = 2.0, dt = 0.5, Nt = 14 → γ ≈ 14, well past Hansen
// γ_e = 3.96 saturation. aniso_factor = 4 with aniso_angle = 80° configures
// inherited cratonic fabric; ani_fstrain = 3 with aniso_db = 1 uses the
// Hansen olivine δ-law and the Boneh DiSRX relaxation; cold T (300 K) keeps
// the strain-rate gate active (τ_eff → ∞) so the cold-limit telescope applies.
// With ani_relax_eps_max = -1 (gate OFF), the operator split collapses to
// ani_fstrain == 2 (Hansen-strict).
//
// Acceptance: at step 14, δ_max ≤ 14.132 + 1e-6 globally. Without the F-init
// helper the operator split telescopes to δ_max ≈ 17.13 (Hansen ceiling
// + (aniso_factor − 1) offset), violating the bound by a factor-3 margin.
//
// One phase, ani_fstrain = 3, no inclusion, periodic-in-x simple shear.
//
// (Note: this scaffold mirrors AniFstrainQuartzBlackford.c — identical setup
//  callback, only the .txt configuration differs.)

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
    RunMDOODZ("AniFstrainOlivineBoneh_InitFS.txt", &setup);
}
