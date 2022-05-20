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
  if (input->model.eqn_state > 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (input->model.shear_style == 0) {
    if (position == S || position == N || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == W || position == E) {
      bc.value = -coordinates.x * input->model.EpsBG;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  } else {
    const double Lz = (double) (input->model.zmax - input->model.zmin);
    if (position == S || position == SE || position == SW) {
      bc.value = -input->model.EpsBG * Lz;
      bc.type  = 11;
    } else if (position == N || position == NE || position == NW) {
      bc.value = input->model.EpsBG * Lz;
      bc.type  = 11;
    } else if (position == E) {
      bc.value = 0.0;
      bc.type  = -12;
    } else if (position == W) {
      bc.value = 0.0;
      bc.type  = -2;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (input->model.shear_style == 0) {
    if (position == W || position == E || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == S || position == N) {
      bc.value = coordinates.z * input->model.EpsBG;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  } else {
    if (position == E || position == W || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type = -12;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type = 0;
    } else {
      bc.value = 0.0;
      bc.type = -1;
    }
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
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ("ShearTemplate.txt", &setup);
}
