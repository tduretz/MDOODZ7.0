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
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lx  = input->model.xmax - input->model.xmin;
  const double VxW = -input->model.user2 / input->scaling.L * input->scaling.t;
  const double VxE =  input->model.user2 / input->scaling.L * input->scaling.t;
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = constant_shear_stress;
  } else if (position == W) {
    bc.value = VxW;
    bc.type  = constant_velocity;
  } else if (position == E)
  {
    bc.value = VxE;
    bc.type  = constant_velocity;
  } else {
    bc.value = 0.0;
    bc.type  = inside;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lz = input->model.zmax - input->model.zmin;
  const double VzS = -input->model.user3 / input->scaling.L * input->scaling.t;
  const double VzN =  input->model.user3 / input->scaling.L * input->scaling.t;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0;
    bc.type  = constant_shear_stress;
  } else if (position == S) {
    bc.value = VzS;
    bc.type  = constant_velocity;
  } else if (position == N)
  {
    bc.value = VzN;
    bc.type  = constant_velocity;
  } else {
    bc.value = 0;
    bc.type  = inside;
  }
  return bc;
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  // .SetBCVx = SetPureOrSimpleShearBCVx,
                  // .SetBCVz = SetPureOrSimpleShearBCVz,
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ("Popov2025_Pureshear_VEVP.txt", &setup);
}
