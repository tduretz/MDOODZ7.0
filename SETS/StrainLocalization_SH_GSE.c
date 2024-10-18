#include "mdoodz.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  return (input->model.user0 + zeroC)/input->scaling.T;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates,  double particleTemperature) {
  SetBC     bc;
  if (position == W || position == E || position == S || position == N) {
    bc.type  = constant_heatflux;
    bc.value = 0.;
  }
  return bc;
}

int main() {
  MdoodzSetup input = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetDensity     = SetDensity,
                  .SetTemperature = SetTemperature,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx   = SetPureOrSimpleShearBCVx,
                  .SetBCVz   = SetPureOrSimpleShearBCVz,
                  .SetBCT    = SetBCT,
          },
  };
  RunMDOODZ("StrainLocalization_SH_GSE.txt", &input);
}
