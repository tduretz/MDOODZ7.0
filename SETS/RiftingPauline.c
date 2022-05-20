#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"


double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel = -0.0e3 / instance->scaling.L;
  const double h_pert    = instance->model.user3 / instance->scaling.L;
  return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness  = instance->model.user1 / instance->scaling.L;
  const double crustThickness        = instance->model.user2 / instance->scaling.L;
  const double perturbationAmplitude = instance->model.user3 / instance->scaling.L;
  const double mohoLevel             = -crustThickness - perturbationAmplitude * cos(2 * M_PI * coordinates.x / (instance->model.xmax - instance->model.xmin));
  const bool   isBelowLithosphere    = coordinates.z < -lithosphereThickness;
  const bool   isAboveMoho           = coordinates.z > mohoLevel;

  if (instance->model.user4 && isAboveMoho) {
    const bool is2500MAboveMoho = coordinates.z > mohoLevel + 2500 / instance->scaling.L;
    const bool is4500MAboveMoho = coordinates.z > mohoLevel + 4500 / instance->scaling.L;
    const bool is7000MAboveMoho = coordinates.z > mohoLevel + 7000 / instance->scaling.L;
    const bool is9000MAboveMoho = coordinates.z > mohoLevel + 9000 / instance->scaling.L;
    if (is2500MAboveMoho && !is4500MAboveMoho) {
      return 7;
    } else if (is7000MAboveMoho && !is9000MAboveMoho) {
      return 8;
    } else {
      return 1;
    }
  } else if (isAboveMoho) {
    return 1;
  } else if (isBelowLithosphere) {
    return 3;
  } else {
    return 2;
  }
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (1330.0 + 273.15) / instance->scaling.T;
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature;
  } else {
    return particleTemperature;
  }
}

double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const int asthenospherePhase = 3;
  return instance->materials.gs_ref[asthenospherePhase];
}

double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.EpsBG;
}

double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.EpsBG;
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W || position == E) {
    bc.value = -coordinates.x * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}


SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = coordinates.z * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}

char SetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else {
    bc.value = 0.0;
    bc.type  = 0;
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

void MutateInput(MdoodzInput *instance) {
  int asthenospherePhases[1]   = {3};
  instance->crazyConductivity = &(CrazyConductivity){
          .multiplier = 1000,
          .nPhases    = 1,
          .phases     = asthenospherePhases,
  };
}

int main() {
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetGrainSize          = SetGrainSize,
                  .SetHorizontalVelocity = SetHorizontalVelocity,
                  .SetVerticalVelocity   = SetVerticalVelocity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx,
                  .SetBCVz    = SetBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
          .MutateInput = MutateInput,

  };
  RunMDOODZ("RiftingPauline.txt", &setup);
}
