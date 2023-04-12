#include "math.h"
#include "mdoodz.h"


double SetSurfaceZCoord(MdoodzInput *input, double x_coord) {
  double Amplitude  = 1.0 / input->scaling.L;
  double Wavelength = 0.5 / input->scaling.L;
  return Amplitude * exp(-x_coord*x_coord/Wavelength/Wavelength);
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  return 0;
}

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == S || position == SW || position == SE) {
    bc.value = 0.0;
    bc.type = 11;
  } else if (position == N || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type = 11;
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
    bc.type = 11;
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type = 0;
  } else {
    bc.value = 0.0;
    bc.type = -1;
  }
  return bc;
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
          },
  };
  RunMDOODZ("FreeSurfaceBenchmarkFEM.txt", &setup);
}
