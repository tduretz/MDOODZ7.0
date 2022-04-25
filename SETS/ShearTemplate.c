#include "mdoodz.h"


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

double SetGrainSize(MdoodzInstance *instance, Coordinates coordinates) {
  return 0.0;
}

double SetPorosity(MdoodzInstance *instance, Coordinates coordinates) {
  return 0.0;
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

double SetTemperature(MdoodzInstance *instance, Coordinates coordinates) {
  return (instance->model.user0 + zeroC) / instance->scaling.T;
}

double SetXComponent(MdoodzInstance *instance, Coordinates coordinates) {
  return 0.0;
}

int SetBCVxType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == LEFT || position == RIGHT) {
      return 0;
    } else if (position == BOTTOM || position == TOP) {
      return 13;
    } else {
      return -1;
    }
  } else {
    if (position == LEFT) {
      return -2;
    } else if (position == RIGHT) {
      return -12;
    } else if (position == BOTTOM || position == TOP) {
      return 11;
    } else {
      return -1;
    }
  }
}

double SetBCVxValue(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  if (instance->model.shear_style == 0) {
    if (position == LEFT || position == RIGHT) {
      return -coordinates.x * instance->model.EpsBG;
    } else {
      return 0;
    }
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == BOTTOM) {
      return -instance->model.EpsBG * Lz;
    } else if (position == TOP) {
      return instance->model.EpsBG * Lz;
    } else {
      return 0;
    }
  }
}

int SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == LEFT || position == RIGHT) {
      return 13;
    } else if (position == BOTTOM || position == TOP) {
      return 0;
    } else {
      return -1;
    }
  } else {
    if (position == LEFT || position == RIGHT) {
      return 0;
    } else if (position == BOTTOM || position == TOP) {
      return -12;
    } else {
      return -1;
    }
  }
}

double SetBCVzValue(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  if (instance->model.shear_style == 0) {
    if (position == TOP || position == BOTTOM) {
      return coordinates.z * instance->model.EpsBG;
    } else {
      return 0;
    }
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == LEFT || position == RIGHT) {
      return 0.0 * instance->model.EpsBG * Lz;
    } else {
      return 0;
    }
  }
}

int SetBCPType(MdoodzInstance *instance, POSITION position) {
  return -1;
}

int SetBCTType(MdoodzInstance *instance, POSITION position) {
  return 0;
}

double SetBCTValue(MdoodzInstance *instance, POSITION position, double particleTemperature) {
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    return surfaceTemperature;
  } else {
    return 0;
  }
}

int main(int nargs, char *args[]) {
  MdoodzInstance instance = NewMdoodzInstance();
  instance.inputFileName  = GetSetupFileName(nargs, args);
  instance.SetParticles   = &(SetParticles_ff){
            .SetPhase              = SetPhase,
            .SetPorosity           = SetPorosity,
            .SetGrainSize          = SetGrainSize,
            .SetVerticalVelocity   = SetVerticalVelocity,
            .SetHorizontalVelocity = SetHorizontalVelocity,
            .SetDensity            = SetDensity,
            .SetXComponent         = SetXComponent,
            .SetTemperature        = SetTemperature,
  };
  instance.SetBCs = &(SetBCs_ff){
          .SetBCVxType  = SetBCVxType,
          .SetBCVxValue = SetBCVxValue,
          .SetBCVzType  = SetBCVzType,
          .SetBCVzValue = SetBCVzValue,
          .SetBCPType   = SetBCPType,
          .SetBCTType   = SetBCTType,
          .SetBCTValue  = SetBCTValue,
  };
  instance.RunMDOODZ(&instance);
}
