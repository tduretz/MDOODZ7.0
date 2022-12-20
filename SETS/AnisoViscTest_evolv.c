#include "mdoodz.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (zeroC) / input->scaling.T;
  if (input->model.eqn_state > 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

void SetDefGrad(MdoodzInput *input, Coordinates coordinates, int phase, Tensor2D* F) {
  const double radius = input->model.user1 / input->scaling.L;
  F->xx=1.; F->xz=0.; F->zx=0.; F->zz=1.;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    // nothing special in the inclusion
  }
  else {
    // a bit more pure shear strain in the matrix
    F->xx=1.1; F->xz=0.; F->zx=0.; F->zz=0.9;
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
                   .SetDefGrad            = SetDefGrad,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("AnisoViscTest_evolv.txt", &setup);
}
