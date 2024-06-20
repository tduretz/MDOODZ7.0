#include "mdoodz.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  int phase = 0;
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
  return input->materials.rho[phase];
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double x = coordinates.x, z =  coordinates.z;
  const double Tcountry   = instance->model.user0 / instance->scaling.T;
  const double Tintrusion = instance->model.user1 / instance->scaling.T;
  const double radius     = instance->model.user2 / instance->scaling.L;
  double T = Tcountry;
  // if (x * x + z * z < radius * radius) {
  //   T = Tintrusion;
  // }
  if (fabs(x)<radius) T = Tintrusion;
  return T;
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  const double Tcountry   = instance->model.user0 / instance->scaling.T;
  if (position == W || position == E) {
    bc.type  = constant_temperature;
    bc.value = Tcountry;
  }
  if (position == S || position == N) {
    bc.type  = constant_heatflux;
    bc.value = 0.0;
  }
  return bc;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "CrystallisationStuewe1995.txt"); // Default
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
                  .SetBCT    = SetBCT,
          },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);

}
