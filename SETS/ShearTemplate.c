#include "mdoodz.h"


int SetPhase(MdoodzInstance *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
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
  } else {
    if (position == SOUTH || position == NORTH || position == NORTHWEST || position == SOUTHWEST || position == NORTHEAST || position == SOUTHEAST) {
      return 11;
    } else if (position == EAST ) {
      return -12;
    } else if (position == WEST) {
      return -2;
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
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {
      return -instance->model.EpsBG * Lz;
    } else if (position == NORTH || position == NORTHEAST || position == NORTHWEST) {
      return instance->model.EpsBG * Lz;
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
  } else {
    if (position == EAST || position == WEST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST ) {
      return -12;
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
  } else {
    return 0;
  }
}

int main(int nargs, char *args[]) {
  MdoodzInstance instance = {
          .inputFileName = GetSetupFileName(nargs, args),
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVxType  = SetBCVxType,
                  .SetBCVxValue = SetBCVxValue,
                  .SetBCVzType  = SetBCVzType,
                  .SetBCVzValue = SetBCVzValue,
          },
  };
  RunMDOODZ(&instance);
}
