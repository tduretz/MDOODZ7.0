#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInstance *instance, Coordinates coordinates) {
    const double A          = 2e-3/instance->scaling.L;
    const double layer_bot0 =-5e-2/instance->scaling.L;
    const double layer_top0 = 5e-2/instance->scaling.L;
    const double Lx         = (instance->model.xmax - instance->model.xmin);
    const double layer_top  = layer_top0 - A*cos(coordinates.x*2.0*M_PI/Lx);
    const double layer_bot  = layer_bot0 + A*cos(coordinates.x*2.0*M_PI/Lx);
    if (coordinates.z>layer_bot && coordinates.z<layer_top) {
        return 1;
    } 
    else {
        return 0;
  }
}
double SetGrainSize(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  return instance->materials.gs_ref[phase];
}

double SetDensity(MdoodzInstance *instance, Coordinates coordinates, int phase) {  // phase
    return instance->materials.rho[phase];
}

double SetTemperature(MdoodzInstance *instance, Coordinates coordinates) {
  const double T = (instance->model.user0 + zeroC) / instance->scaling.T;
  return T;
}

char SetBCVxType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST) {
      return 0;
    } else if (position == SOUTH || position == NORTH) {
      return 13;
    } else {
      return -1;
    }
  }
}

double SetBCVxValue(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST) {
      return -coordinates.x * instance->model.EpsBG;
    } else {
      return 0;
    }
  }
}

char SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST) {
      return 13;
    } else if (position == SOUTH || position == NORTH) {
      return 0;
    } else {
      return -1;
    }
  } 
}

double SetBCVzValue(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  if (instance->model.shear_style == 0) {
    if (position == NORTH || position == NORTHEAST || position == NORTHWEST || position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {
      return coordinates.z * instance->model.EpsBG;
    } else {
      return 0;
    }
  }
}

//----------------------------- THERMAL BC -----------------------------//

char SetBCTType(MdoodzInstance *instance, POSITION position) {
    return 0;
}

char SetBCTTypeNew(MdoodzInstance *instance, POSITION position) {
    return 0;
}

double SetBCTValue(MdoodzInstance *instance, POSITION position, double gridTemperature) {
    return gridTemperature;
}

double SetBCTValueNew(MdoodzInstance *instance, POSITION position, double gridTemperature) {
    return gridTemperature;
}

//----------------------------- MAIN -----------------------------//

int main(int nargs, char *args[]) {
  MdoodzInstance instance = {
          .inputFileName = GetSetupFileName(nargs, args),
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
                   .SetGrainSize          = SetGrainSize,
                   .SetTemperature        = SetTemperature,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVxType  = SetBCVxType,
                  .SetBCVxValue = SetBCVxValue,
                  .SetBCVzType  = SetBCVzType,
                  .SetBCVzValue = SetBCVzValue,
                  .SetBCTType   = SetBCTType,
                  .SetBCTValue  = SetBCTValue,
                  .SetBCTTypeNew   = SetBCTTypeNew,
                  .SetBCTValueNew  = SetBCTValueNew,
          },
  };
  RunMDOODZ(&instance);
}
