#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

double FixTemperature(MdoodzInput *input, double p) {
  int thermal_evolution = (int) input->model.user0;
  if (thermal_evolution) {
    double T = 273.15/input->scaling.T + 0.1*input->model.time;
 // Tfix = (log( mesh->p_in[0]*scaling.S/1e9 /b) / a + zeroC)/scaling.T;
    return T; //(log(p * input->scaling.S / 1e9 / b) / a + zeroC) / input->scaling.T;
  } else {
    return input->model.bkg_temperature;
  }
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, DefaultTextFilename(__FILE__)); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
                  .FixTemperature = FixTemperature,
          },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
