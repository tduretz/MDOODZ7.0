#include "mdoodz.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"



// for the ellipse, we have several arguments: effective radius, aspect ratio, orientation, X0,Z0
int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    double radius = $RADIUS
    double AR = $AR
    double a,b, cos_theta = cos($THETA*PI/180.0), sin_theta=sin($THETA*PI/180.0)
    double $X0,$Z0

    if (AR>1){
      // compute half axis lengths so that we have the same "volume" as a circle with radius
      b   = radius / sqrt(AR);
      a   = radius * sqrt(AR);
    
      Xeff = (coordinates.x-X0)*cos_theta + (coordinates.z-Z0)*sin_theta;
      Zeff = (coordinates.x-X0)*sin_theta - (coordinates.z-Z0)*cos_theta;
    } 
    else {
      b = radius;
      a = radius;
      Xeff = coordinates.x-X0;
      Zeff = coordinates.z-Z0;
    }

    if ( Xeff*Xeff/(a*a) +Zeff*Zeff/(b*b) <1) {
        return 1;
    } else {
        return 0;
    }
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

int main(int argc, char *args[]) {

// Input file name
  char *input_file;
// effective radius, aspect ratio, orientation, X0,Z0 
  double radius, AR, theta, X0, Z0
  if ( nargs < 2 ) { // no input given
    asprintf(&input_file, "Shrinking.txt"); // Default
    radius = 0.1;
    AR     = 1;
    theta  = 0;
    X0     = 0;
    Z0     = 0;
  }
  else if ( nargs == 2) { // only input file given
    asprintf(&input_file, "%s", args[1]);
    radius = 0.1;
    AR     = 1;
    theta  = 0;
    X0     = 0;
    Z0     = 0;
  }
  else if ( nargs == 3) { // only input file and inclusion radius given --> circular inclusion
    asprintf(&input_file, "%s", args[1]);
    radius = atof(args[2]);
    AR     = 1;
    theta  = 0;
    X0     = 0;
    Z0     = 0;
  }
  else if ( nargs == 5) { // input file, inclusion radius, aspect ratio and orientation are given
    asprintf(&input_file, "%s", args[1]);
    radius = atof(args[2]);
    AR     = atof(args[3]);
    theta  = atof(args[4]);
    X0     = 0;
    Z0     = 0;
  }
  else if ( nargs == 7) { // input file, inclusion radius, aspect ratio, orientation, X0, Z0 are given
    asprintf(&input_file, "%s", args[1]);
    radius = atof(args[2]);
    AR     = atof(args[3]);
    theta  = atof(args[4]);
    X0     = atof(args[5]);
    Z0     = atof(args[6]);
  }
  else {
    asprintf(&input_file, "Shrinking.txt"); // Default
    radius = 0.1;
    AR     = 1;
    theta  = 0;
    X0     = 0;
    Z0     = 0;
    printf("Wrong number of arguments, running MDoodz7.0 using default %s\n", input_file);
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  
  if (AR==1.0) { // circular inclusion
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhaseCircle,
                   .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };

  }
  else { // elliptic inclusion
    MdoodzSetup setup = {
        .SetParticles  = &(SetParticles_ff){
                .SetPhase              = SetPhaseEllipse,
                .SetDensity            = SetDensity,
        },
        .SetBCs = &(SetBCs_ff){
              .SetBCVx = SetPureOrSimpleShearBCVx,
              .SetBCVz = SetPureOrSimpleShearBCVz,
        },
  };
  }
  RunMDOODZ("Shrinking.txt", &setup);
}
