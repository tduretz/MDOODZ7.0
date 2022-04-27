#include "assert.h"
#include "hdf5.h"
#include "math.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"
#define FILENAME "Output00001.gzip.h5"


double SetHorizontalVelocity(MdoodzInstance *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.EpsBG;
}

double SetVerticalVelocity(MdoodzInstance *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.EpsBG;
}

int SetPhase(MdoodzInstance *instance, Coordinates coordinates) {
  double       xc     = 0.0;
  double       zc     = 0.0;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double X      = coordinates.x - xc;
  const double Z      = coordinates.z - zc;
  if (X * X + Z * Z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInstance *instance, Coordinates coordinates) {
  const double T_init = (instance->model.user0 + zeroC) / instance->scaling.T;
  double       xc     = 0.0;
  double       zc     = 0.0;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double X      = coordinates.x - xc;
  const double Z      = coordinates.z - zc;
  int          phase;
  if (X * X + Z * Z < radius * radius) {
    phase = 1;
  } else {
    phase = 0;
  }
  if (instance->model.eqn_state > 0) {
    return instance->materials.rho[phase] * (1 - instance->materials.alp[phase] * (T_init - instance->materials.T0[phase]));
  } else {
    return instance->materials.rho[phase];
  }
}

int SetBCVxType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST) {
      return 0;
    } else if (position == SOUTH || position == NORTH) {
      return 13;
    } else {
      return -1;
    }
  } else {
    if (position == WEST || position == NORTHWEST || position == SOUTHWEST) {
      return -2;
    } else if (position == EAST || position == NORTHEAST || position == SOUTHEAST) {
      return -12;
    } else if (position == SOUTH || position == NORTH) {
      return 11;
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
    if (position == SOUTH) {
      return -instance->model.EpsBG * Lz;
    } else if (position == NORTH) {
      return instance->model.EpsBG * Lz;
    } else {
      return 0;
    }
  }
}

int SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST) {
      return 13;
    } else if (position == SOUTH || position == NORTH) {
      return 0;
    } else {
      return -1;
    }
  } else {
    if (position == WEST || position == EAST) {
      return 0;
    } else if (position == SOUTH || position == NORTH) {
      return -12;
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
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == WEST || position == EAST || position == NORTHWEST || position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {
      return 0.0 * instance->model.EpsBG * Lz;
    } else {
      return 0;
    }
  }
}

int main() {
  MdoodzInstance instance = NewMdoodzInstance();
  instance.inputFileName  = "ShearTemplate.txt";
  instance.SetParticles   = &(SetParticles_ff){
            .SetPhase              = SetPhase,
            .SetVerticalVelocity   = SetVerticalVelocity,
            .SetHorizontalVelocity = SetHorizontalVelocity,
            .SetDensity            = SetDensity,
  };
  instance.SetBCs = &(SetBCs_ff){
          .SetBCVxType  = SetBCVxType,
          .SetBCVxValue = SetBCVxValue,
          .SetBCVzType  = SetBCVzType,
          .SetBCVzValue = SetBCVzValue,
  };
  instance.RunMDOODZ(&instance);

  hid_t File            = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t IterationsGroup = H5Gopen(File, "Iterations", H5P_DEFAULT);
  hid_t NumberStepsDataset =
          H5Dopen(IterationsGroup, "NumberSteps", H5P_DEFAULT);
  int NumberStepsArray[1] = {0};
  H5Dread(NumberStepsDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          NumberStepsArray);
  int StepsCount = NumberStepsArray[0];

  assert(StepsCount == 1);
}
