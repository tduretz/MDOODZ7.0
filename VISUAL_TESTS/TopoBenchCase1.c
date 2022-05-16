#include "math.h"
#include "visual-tests.h"


double SetSurfaceZCoord1(MdoodzInstance *instance, double x_coord) {
  double Amplitude  = 7e3 / instance->scaling.L;
  double Wavelength = 2800e3 / instance->scaling.L;
  return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
}

int SetPhase1(MdoodzInstance *instance, Coordinates coordinates) {
  double lithosphereBottomDepth = -100.0e3 / instance->scaling.L;
  if (coordinates.z > lithosphereBottomDepth) {
    return 1;
  } else {
    return 0;
  }
}

SetBC SetBCVx1(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == S || position == SW || position == SE) {
    bc.value = 0.0;
    bc.type  = 11;
  } else if (position == N || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if (position == W || position == E) {
    bc.value = 0.0;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz1(MdoodzInstance *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

char SetBCPType1(MdoodzInstance *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

MdoodzInstance CreateTopoBenchCase1Instance() {
  int astenospherePhases[2] = {2, 4};
  return (MdoodzInstance){
          .inputFileName          = "TopoBenchCase1.txt",
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord1,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase = SetPhase1,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx1,
                  .SetBCVz    = SetBCVz1,
                  .SetBCPType = SetBCPType1,
          },
          .crazyConductivity = &(CrazyConductivity){
                  .multiplier = 1000,
                  .nPhases    = 2,
                  .phases     = astenospherePhases,
          }};
}
