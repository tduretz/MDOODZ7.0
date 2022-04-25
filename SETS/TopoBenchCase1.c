#include "math.h"
#include "mdoodz.h"


double SetSurfaceZCoord(MdoodzInstance *instance, double x_coord) {
  double Amplitude  = 7e3 / instance->scaling.L;
  double Wavelength = 2800e3 / instance->scaling.L;
  return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
}

int SetSurfacePhase(MdoodzInstance *instance, double x_coord) {
  return 0;
}

int SetPhase(MdoodzInstance *instance, Coordinates coordinates) {
  double lithosphereBottomDepth = -100.0e3 / instance->scaling.L;
  if (coordinates.z > lithosphereBottomDepth) {
    return 1;
  } else {
    return 0;
  }
}

int SetBCVxType(MdoodzInstance *instance, POSITION position) {
  if (position == LEFT || position == RIGHT) {
    return 0;
  } else if (position == BOTTOM) {
    return 11;
  } else if (position == TOP) {
    return 13;
  } else {
    return -1;
  }
}

int SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (position == LEFT || position == RIGHT) {
    return 13;
  } else if (position == BOTTOM || position == TOP) {
    return 0;
  } else {
    return -1;
  }
}

int SetBCPType(MdoodzInstance *instance, POSITION position) {
  if (position == TOPRIGHT || position == TOPLEFT) {
    return 0;
  } else {
    return -1;
  }
}

double SetBCTValue(MdoodzInstance *instance, POSITION position, double particleTemperature) {
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    return surfaceTemperature;
  } else {
    return 0;
  }
}

double SetBCTValueNew(MdoodzInstance *instance, POSITION position, double particleTemperature) {
  double surfaceTemperature = zeroC / instance->scaling.T;
  double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == BOTTOM || position == BOTTOMRIGHT || position == BOTTOMLEFT) {
    return particleTemperature;
  } else if (position == TOP || position == TOPRIGHT || position == TOPLEFT) {
    return surfaceTemperature;
  } else if (position == LEFT || position == RIGHT) {
    return mantleTemperature;
  } else {
    return 0;
  }
}

int SetBCTType(MdoodzInstance *instance, POSITION position) {
  if (position == FREE_SURFACE) {
    return 1;
  } else {
    return 0;
  }
}

int SetBCTTypeNew(MdoodzInstance *instance, POSITION position) {
  if (position == TOP || position == BOTTOM) {
    return 1;
  } else {
    return 0;
  }
}

int main(int nargs, char *args[]) {
  MdoodzInstance instance         = NewMdoodzInstance();
  instance.inputFileName          = GetSetupFileName(nargs, args);
  instance.BuildInitialTopography = &(BuildInitialTopography_ff){
          .SetSurfacePhase  = SetSurfacePhase,
          .SetSurfaceZCoord = SetSurfaceZCoord,
  };
  instance.SetParticles = &(SetParticles_ff){
          .SetPhase = SetPhase,
  };
  instance.SetBCs = &(SetBCs_ff){
          .SetBCVxType    = SetBCVxType,
          .SetBCVzType    = SetBCVzType,
          .SetBCPType     = SetBCPType,
          .SetBCTType     = SetBCTType,
          .SetBCTTypeNew  = SetBCTTypeNew,
          .SetBCTValue    = SetBCTValue,
          .SetBCTValueNew = SetBCTValueNew,
  };
  int               astenospherePhases[2] = {2, 4};
  CrazyConductivity crazyConductivity     = {
              .multiplier = 1000,
              .nPhases    = 2,
              .phases     = astenospherePhases,
  };

  instance.crazyConductivity = &crazyConductivity;
  instance.RunMDOODZ(&instance);
}
