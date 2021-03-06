#include "math.h"
#include "mdoodz.h"

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double A          = 2e-3 / instance->scaling.L;
  const double layer_bot0 = -5e-2 / instance->scaling.L;
  const double layer_top0 = 5e-2 / instance->scaling.L;
  const double Lx         = (instance->model.xmax - instance->model.xmin);
  const double layer_top  = layer_top0 - A * cos(coordinates.x * 2.0 * M_PI / Lx);
  const double layer_bot  = layer_bot0 + A * cos(coordinates.x * 2.0 * M_PI / Lx);
  if (coordinates.z > layer_bot && coordinates.z < layer_top) {
    return 1;
  } else {
    return 0;
  }
}
double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  return instance->materials.gs_ref[phase];
}

double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {// phase
  return instance->materials.rho[phase];
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double T = (instance->model.user0 + zeroC) / instance->scaling.T;
  return T;
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == NE || position == NW || position == SE || position == SW) {
    bc.value = -coordinates.x * instance->model.EpsBG;
    bc.type  = 0;
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type  = 13;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == NE || position == NW || position == SE || position == SW) {
    bc.value = coordinates.z * instance->model.EpsBG;
    bc.type  = 13;
  } else if (position == W || position == E) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = coordinates.z * instance->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

//----------------------------- THERMAL SetBC -----------------------------//


SetBC SetBCT(MdoodzInput *instance, POSITION position, double gridTemperature) {
  return (SetBC){.value = gridTemperature, .type = 0};
}

SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double gridTemperature) {
  return (SetBC){.value = gridTemperature, .type = 0};
}

//----------------------------- MAIN -----------------------------//

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase       = SetPhase,
                   .SetDensity     = SetDensity,
                   .SetGrainSize   = SetGrainSize,
                   .SetTemperature = SetTemperature,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVz   = SetBCVz,
                  .SetBCVx   = SetBCVx,
                  .SetBCT    = SetBCT,
                  .SetBCTNew = SetBCTNew,
          },
  };
  RunMDOODZ( "PinchSwellGSE.txt", &setup);
}
