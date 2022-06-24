#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"


int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double HLit  = 70e3 / input->scaling.L;
  const double X0    = (75e3 * (double) rand() / RAND_MAX - 75e3 / 2) / input->scaling.L;
  const double Z0    = -HLit * (double) rand() / RAND_MAX;
  const double X     = coordinates.x;
  const double Z     = coordinates.z;
  const double a     = 1 / 0.005 / input->scaling.L;
  const double b     = 1 / 0.000016 / input->scaling.L;
  const double theta = 0 * M_PI / 180;
  if ((a * (X - X0) * (X - X0) + b * (Z - Z0) * (Z - Z0)) * pow(cos(theta), 2) + (b * (X - X0) * (X - X0) + a * (Z - Z0) * (Z - Z0)) * pow(sin(theta), 2) + (a - b) * (X - X0) * (Z - Z0) * sin(2 * theta) - 1 < 0) {
    if (Z0 > -HLit / 2) {
      return 1;
    }
    if (Z0 <= -HLit / 2) {
      return 3;
    }
  }
  if (coordinates.z < -HLit / 2) {
    return 2;
  }
  return 0;
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  const double T0    = 273.15 / input->scaling.T;
  const double TmaxW = (1250 + 273.15) / input->scaling.T + input->model.Tinit;
  const double HLit  = 70e3 / input->scaling.L;
  const double Htot  = 2 * HLit;
  return ((TmaxW - T0) / Htot) * (-coordinates.z - 0e3 / input->scaling.L) + T0;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double TPart = SetTemperature(input, coordinates);
  if ( input->model.eqn_state ==  1 ) {
    return input->materials.rho[phase] * (1 -  input->materials.alp[phase] * (TPart - input->materials.T0[phase]) );
  } else {
    return input->materials.rho[phase];
  }
}

char SetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC  bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else {
    bc.value = 0.0;
    bc.type  = 0;
  }
  return bc;
}


SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC  bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == S || position == SE || position == SW) {
    bc.value = particleTemperature;
    bc.type  = 1;
  } else if (position == N || position == NE || position == NW) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else if (position == W || position == E) {
    bc.value = mantleTemperature;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = 0;
  }
  return bc;
}

void MutateInput(MdoodzInput *instance) {
  int asthenospherePhases[1]  = {3};
  instance->crazyConductivity = &(CrazyConductivity){
          .multiplier = 1000,
          .nPhases    = 1,
          .phases     = asthenospherePhases,
  };
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetTemperature = SetTemperature,
                  .SetDensity     = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetPureShearBCVx,
                  .SetBCVz    = SetPureShearBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
          .MutateInput = MutateInput,

  };
  RunMDOODZ("RiftingPauline.txt", &setup);
}
