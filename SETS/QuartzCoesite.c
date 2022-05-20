#include "mdoodz.h"


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

SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == S || position == SE || position == SW) {
    bc.value = particleTemperature;
    bc.type  = 1;
  } else if (position == N || position == NE || position == NW) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else if (position == W || position == E) {
    bc.value = mantleTemperature;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = 0;
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
