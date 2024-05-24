#include "mdoodz.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double inner_radius = 15. / input->scaling.L;
  const double outer_radius = 50. / input->scaling.L;
  const double x = coordinates.x, z =  coordinates.z;
  int phase = 0;
  if (x * x + z * z < outer_radius * outer_radius) {
    phase = 1;
  }
  if (x * x + z * z < inner_radius * inner_radius) {
    phase = 2;
  }
  return phase;
}

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // Set checkerboard for phase 0
    Ax = cos( 6.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( 6.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Ax = cos( 24.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( 24.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 ) {
        dual_phase += input->model.Nb_phases;
    }

  return dual_phase;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  const double P_init = (input->model.bkg_pressure         ) / input->scaling.S;
  if (1 == 0) {
    return input->materials.rho[phase] * exp(input->materials.bet[phase]*P_init -  input->materials.alp[phase] * T_init);
  } else {
    return input->materials.rho[phase];
  }
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  return instance->model.user0 / instance->scaling.T;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "MeltingOverpressure.txt"); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDualPhase          = SetDualPhase,
                   .SetDensity            = SetDensity,
                   .SetTemperature        = SetTemperature,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);

}
