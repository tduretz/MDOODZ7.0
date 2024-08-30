#include "mdoodz.h"
#include <stdlib.h>
#include "stdio.h"
#include <math.h>

int SetPhase(MdoodzInput *input, Coordinates coord) {
  const double r = input->model.user1 / input->scaling.L;
  int phase = 0;
  if ( fabs(coord.x)>0.5 ) phase = 1;
  if (coord.x * coord.x + coord.z * coord.z < r*r ) phase = 2;
  return phase;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

double SetPressure(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.bkg_pressure;
}

double SetTxx(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.preload_sxxd;
}

double SetTzz(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.preload_szzd;
}

double SetTxz(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.preload_sxz;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "ShearBoxLaeti.txt"); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
                   .SetPressure           = SetPressure,
                   .SetTxx                = SetTxx,
                   .SetTzz                = SetTzz,
                   .SetTxz                = SetTxz,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
