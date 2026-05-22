// AnisoFstrainBox — unified single-phase homogeneous-box callback for the
// `ani_fstrain` visual-tests suite (formerly 29 sister SETs). Drives Tests
// 1-6 of `VISUAL_TESTS/aniso_fstrain/` via Python runners that generate per-
// cell `.txt` files from `VISUAL_TESTS/aniso_fstrain/templates/base.txt` and
// invoke the binary with `argv[1] = <path-to-.txt>`. The C code itself is
// trivial: phase 0 everywhere, density from materials, T from
// `bkg_temperature`, simple-shear BCs.
//
// Replaces (in SETS/CMakeLists.txt):
//   AniFstrainOlivineCompare_F{1,2,3}            (Test 1, 3 sister SETs)
//   AniFstrainOlivineFInit_a{000,030,060,080,090,120}  (Test 2, 6 sets)
//   AniFstrainOlivineFactorSweep_f{01..14}      (Test 3, 9 sets)
//   AniFstrainOlivineInheritSweep_f{1,2,4,8}    (Test 4, 4 sets)
//   AniFstrainOlivineTSweep_T{0500..1500}       (Test 5, 6 sets)
//   AniFstrainOlivineRegimeCell                 (Test 6 single binary)
//
// `AniFstrainOlivineBoneh_InitFS` is the Phase-1 regression SET and is
// deliberately NOT consolidated here — separate scope.

#include "mdoodz.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    return 0;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
}

double SetT(MdoodzInput *input, Coordinates coordinates) {
    return input->model.bkg_temperature;
}

int main(int nargs, char *args[]) {
    char *input_file;
    if (nargs < 2) {
        asprintf(&input_file, DefaultTextFilename(__FILE__)); // default = AnisoFstrainBox.txt
    } else {
        asprintf(&input_file, "%s", args[1]);                 // explicit .txt from Python runner
    }
    printf("Running MDoodz7.0 (AnisoFstrainBox) using %s\n", input_file);

    MdoodzSetup setup = {
            .SetParticles = &(SetParticles_ff){
                    .SetPhase       = SetPhase,
                    .SetDensity     = SetDensity,
                    .SetTemperature = SetT,
            },
            .SetBCs = &(SetBCs_ff){
                    .SetBCVx = SetSimpleShearBCVx,
                    .SetBCVz = SetSimpleShearBCVz,
            },
    };
    RunMDOODZ(input_file, &setup);
    free(input_file);
}
