#include "mdoodz.h"


int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  if (coordinates.z < 2e3 / input->scaling.L && coordinates.z > -2e3 / input->scaling.L) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureShearBCVx,
                  .SetBCVz = SetPureShearBCVz,
          },
  };
  RunMDOODZ("AnisotropyLayer.txt", &setup);
}
