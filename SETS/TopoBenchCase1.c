#include "math.h"
#include "mdoodz.h"


double SetSurfaceZCoord(MdoodzInstance *instance, double x_coord) {
  double Amplitude  = 7e3 / instance->scaling.L;
  double Wavelength = 2800e3 / instance->scaling.L;
  return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
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
   if (position == BOTTOM || position == BOTTOMLEFT || position == BOTTOMRIGHT) {
    return 11;
  } else if (position == TOP || position == TOPLEFT || position == TOPRIGHT) {
    return 13;
  } else if (position == LEFT || position == RIGHT) {
    return 0;
  } else {
    return -1;
  }
}

int SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (position == LEFT || position == RIGHT || position == BOTTOMLEFT || position == BOTTOMRIGHT || position == TOPLEFT || position == TOPRIGHT) {
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

double SetBCVxValue() {
  return 0.0;
}

double SetBCVzValue() {
  return 0.0;
}

int main(int nargs, char *args[]) {
  MdoodzInstance instance         = NewMdoodzInstance();
  instance.inputFileName          = GetSetupFileName(nargs, args);
  instance.BuildInitialTopography = &(BuildInitialTopography_ff){
          .SetSurfaceZCoord = SetSurfaceZCoord,
  };
  instance.SetParticles = &(SetParticles_ff){
          .SetPhase = SetPhase,
  };
  instance.SetBCs = &(SetBCs_ff){
          .SetBCVxType    = SetBCVxType,
          .SetBCVxValue   = SetBCVxValue,
          .SetBCVzValue   = SetBCVzValue,
          .SetBCVzType    = SetBCVzType,
          .SetBCPType     = SetBCPType,
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
