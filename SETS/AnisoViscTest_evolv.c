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

Tensor2D SetDefGrad(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    // nothing special in the inclusion
    return (Tensor2D) {.xx = 1., .xz = 0., .zx = 0., .zz = 1.,}; // initializes Fxx, Fxz, Fzx, Fzz 
  }
  else {
    // a bit more pure shear strain in the matrix
    return (Tensor2D) {.xx = 1.1, .xz = 0., .zx = 0., .zz = 0.9,}; // initializes Fxx, Fxz, Fzx, Fzz 
  }
}

double SetAnisoAngle(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    // nothing special in the inclusion
    return 90.0;
  }
  else {
    // a bit more pure shear strain in the matrix
    return 135.0;
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
                   .SetDefGrad            = SetDefGrad,
                   .SetAnisoAngle         = SetAnisoAngle,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("AnisoViscTest_evolv.txt", &setup);
}
