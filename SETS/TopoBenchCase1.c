#include "math.h"
#include "mdoodz.h"


double SetSurfaceZCoord(MdoodzInput *input, double x_coord) {
  double Amplitude  = 7e3 / input->scaling.L;
  double Wavelength = 2800e3 / input->scaling.L;
  return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  double lithosphereBottomDepth = -100.0e3 / input->scaling.L;
  if (coordinates.z > lithosphereBottomDepth) {
    return 1;
  } else {
    return 0;
  }
}

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == S || position == SW || position == SE) {
    bc.value = 0.0;
    bc.type = 11;
  } else if (position == N || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type = 13;
  } else if (position == W || position == E) {
    bc.value = 0.0;
    bc.type = 0;
  } else {
    bc.value = 0.0;
    bc.type = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type = 13;
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type = 0;
  } else {
    bc.value = 0.0;
    bc.type = -1;
  }
  return bc;
}

char SetBCPType(MdoodzInput *input, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

int main() {
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase = SetPhase,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx,
                  .SetBCVz    = SetBCVz,
                  .SetBCPType = SetBCPType,
          },
  };
  RunMDOODZ("TopoBenchCase1.txt", &setup);
}
