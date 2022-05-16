#include "visual-tests.h"


int SetPhase2(MdoodzInstance *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity2(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

SetBC SetBCVx2(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (instance->model.shear_style == 0) {
    if (position == S || position == N || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == W || position == E) {
      bc.value = -coordinates.x * instance->model.EpsBG;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == S || position == SE || position == SW) {
      bc.value = -instance->model.EpsBG * Lz;
      bc.type  = 11;
    } else if (position == N || position == NE || position == NW) {
      bc.value = instance->model.EpsBG * Lz;
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

SetBC SetBCVz2(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (instance->model.shear_style == 0) {
    if (position == W || position == E || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == S || position == N) {
      bc.value = coordinates.z * instance->model.EpsBG;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  } else {
    if (position == E || position == W || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type  = -12;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  }
  return bc;
}

MdoodzInstance CreateShearTemplateInstance(int shear_style) {
  char *inputFileName;
  if (shear_style) {
    inputFileName = "ShearTemplate1.txt";
  } else {
    inputFileName = "ShearTemplate.txt";
  }
  return (MdoodzInstance){
          .inputFileName = inputFileName,
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase   = SetPhase2,
                   .SetDensity = SetDensity2,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetBCVx2,
                  .SetBCVz = SetBCVz2,
          },
  };
}
