#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"


double SetSurfaceZCoord(MdoodzInstance *instance, double x_coord) {
  const double TopoLevel = -0.0e3 / instance->scaling.L;
  const double h_pert    = instance->model.user3 / instance->scaling.L;
  return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
}

int SetPhase(MdoodzInstance *instance, Coordinates coordinates) {
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

double SetTemperature(MdoodzInstance *instance, Coordinates coordinates) {
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

double SetGrainSize(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  const int astenospherePhase = 3;
  return instance->materials.gs_ref[astenospherePhase];
}

double SetHorizontalVelocity(MdoodzInstance *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.EpsBG;
}

double SetVerticalVelocity(MdoodzInstance *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.EpsBG;
}

char SetBCVxType(MdoodzInstance *instance, POSITION position) {
  if (position == NORTH || position == SOUTH || position == NORTHWEST || position == SOUTHWEST || position == NORTHEAST || position == SOUTHEAST) {
    return 13;
  } else if (position == WEST || position == EAST) {
    return 0;
  } else {
    return -1;
  }
}

double SetBCVxValue(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  if (position == WEST || position == EAST) {
    return -coordinates.x * instance->model.EpsBG;
  } else {
    return 0;
  }
}

char SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (position == WEST || position == EAST || position == SOUTHWEST || position == SOUTHEAST || position == NORTHWEST || position == NORTHEAST) {
    return 13;
  } else if (position == SOUTH || position == NORTH) {
    return 0;
  } else {
    return -1;
  }
}

double SetBCVzValue(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  if (position == NORTH || position == SOUTH || position == NORTHWEST || position == SOUTHWEST || position == NORTHEAST || position == SOUTHEAST) {
    return coordinates.z * instance->model.EpsBG;
  } else {
    return 0;
  }
}

char SetBCPType(MdoodzInstance *instance, POSITION position) {
  if (position == NORTHEAST || position == NORTHWEST) {
    return 0;
  } else {
    return -1;
  }
}

double SetBCTValue(MdoodzInstance *instance, POSITION position, double gridTemperature) {
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    return surfaceTemperature;
  } else {
    return 0;
  }
}

double SetBCTValueNew(MdoodzInstance *instance, POSITION position, double gridTemperature) {
  double surfaceTemperature = zeroC / instance->scaling.T;
  double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {
    return gridTemperature;
  } else if (position == NORTH || position == NORTHEAST || position == NORTHWEST) {
    return surfaceTemperature;
  } else if (position == WEST || position == EAST) {
    return mantleTemperature;
  } else {
    return 0;
  }
}

char SetBCTType(MdoodzInstance *instance, POSITION position) {
  if (position == FREE_SURFACE) {
    return 1;
  } else {
    return 0;
  }
}

char SetBCTTypeNew(MdoodzInstance *instance, POSITION position) {
  if (position == NORTH || position == SOUTH || position == NORTHWEST || position == SOUTHWEST || position == NORTHEAST || position == SOUTHEAST) {
    return 1;
  } else {
    return 0;
  }
}

int main(int nargs, char *args[]) {
  int            astenospherePhases[1] = {3};
  MdoodzInstance instance              = {
                       .inputFileName          = GetSetupFileName(nargs, args),
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
                               .SetBCVxType    = SetBCVxType,
                               .SetBCVzType    = SetBCVzType,
                               .SetBCVxValue   = SetBCVxValue,
                               .SetBCVzValue   = SetBCVzValue,
                               .SetBCPType     = SetBCPType,
                               .SetBCTType     = SetBCTType,
                               .SetBCTTypeNew  = SetBCTTypeNew,
                               .SetBCTValue    = SetBCTValue,
                               .SetBCTValueNew = SetBCTValueNew,
          },
                       .crazyConductivity = &(CrazyConductivity){
                               .multiplier = 1000,
                               .nPhases    = 1,
                               .phases     = astenospherePhases,
          },

  };
  RunMDOODZ(&instance);
}
