#include "mdoodz.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double HCrust  = 30e3 / input->scaling.L;
  const double HMantle = 40e3 / input->scaling.L;
  if (coordinates.z > -HCrust) {
    return 0;
  } else if (coordinates.z > -HCrust - HMantle) {
    return 1;
  } else {
    return 2;
  }
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  const double T0    = 273.15 / input->scaling.T;
  const double TmaxW = (1250 + 273.15) / input->scaling.T;
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


int main() {
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetTemperature = SetTemperature,
                  .SetDensity     = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetPureShearBCVx,
                  .SetBCVz    = SetPureShearBCVz,
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
  };
  RunMDOODZ("RiftingAnisotropy.txt", &setup);
}
