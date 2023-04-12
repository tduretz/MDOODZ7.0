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

//----------------------------- THERMAL SetBC -----------------------------//


SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  if (position == W || position == E || position == S || position == N) {
    bc.type  = constant_heatflux;
    bc.value = 0.;
  }
  return bc;
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
                  .SetBCVz   = SetPureShearBCVz,
                  .SetBCVx   = SetPureShearBCVx,
                  .SetBCT    = SetBCT,
          },
  };
  RunMDOODZ( "PinchSwellGSE.txt", &setup);
}
