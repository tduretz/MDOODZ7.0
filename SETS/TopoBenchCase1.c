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

BC SetBCVx(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  BC bc = {.value = 0.0};
  if (position == S || position == SW || position == SE) {
    bc.type = 11;
  } else if (position == N || position == NW || position == NE) {
    bc.type = 13;
  } else if (position == W || position == E) {
    bc.type = 0;
  } else {
    bc.type = -1;
  }
  return bc;
}

BC SetBCVz(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  BC bc = {.value = 0.0};
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.type = 13;
  } else if (position == S || position == N) {
    bc.type = 0;
  } else {
    bc.type = -1;
  }
  return bc;
}

BC SetBCP(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  BC bc = {.value = 0.0};
  if (position == NE || position == NW) {
    bc.type = 0;
  } else {
    bc.type = -1;
  }
  return bc;
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
                               .SetBCVx = SetBCVx,
                               .SetBCVz = SetBCVz,
                               .SetBCP  = SetBCP,
          },
                       .crazyConductivity = &(CrazyConductivity){
                               .multiplier = 1000,
                               .nPhases    = 2,
                               .phases     = astenospherePhases,
          }};
  RunMDOODZ(&instance);
}
