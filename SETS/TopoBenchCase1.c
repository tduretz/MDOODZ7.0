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

char SetBCVxType(MdoodzInstance *instance, POSITION position) {
   if (position == SOUTH || position == SOUTHWEST || position == SOUTHEAST) {
    return 11;
  } else if (position == NORTH || position == NORTHWEST || position == NORTHEAST) {
    return 13;
  } else if (position == WEST || position == EAST) {
    return 0;
  } else {
    return -1;
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

char SetBCPType(MdoodzInstance *instance, POSITION position) {
  if (position == NORTHEAST || position == NORTHWEST) {
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
  int            astenospherePhases[2] = {2, 4};
  MdoodzInstance instance              = {
                       .inputFileName          = GetSetupFileName(nargs, args),
                       .BuildInitialTopography = &(BuildInitialTopography_ff){
                               .SetSurfaceZCoord = SetSurfaceZCoord,
          },
                       .SetParticles = &(SetParticles_ff){
                               .SetPhase = SetPhase,
          },
                       .SetBCs = &(SetBCs_ff){
                               .SetBCVxType  = SetBCVxType,
                               .SetBCVxValue = SetBCVxValue,
                               .SetBCVzValue = SetBCVzValue,
                               .SetBCVzType  = SetBCVzType,
                               .SetBCPType   = SetBCPType,
          },
                       .crazyConductivity = &(CrazyConductivity){
                               .multiplier = 1000,
                               .nPhases    = 2,
                               .phases     = astenospherePhases,
          }};
  RunMDOODZ(&instance);
}
