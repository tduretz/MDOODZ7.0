#include "mdoodz.h"
#include "math.h"

/*
int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = 0.01 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}
*/

// for the ellipse, we have several arguments: effective radius, aspect ratio, orientation, X0,Z0
int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    double radius = 0.01 / input->scaling.L;
    double AR = 2;
    double a,b, Xeff,Zeff;
    double cos_theta = cos(45*M_PI/180.0); 
    double sin_theta=sin(45*M_PI/180.0);
    double X0=0;
    double Z0=0;

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

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDualPhase          = SetDualPhase,
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
  RunMDOODZ("Shrinking_preload.txt", &setup);
}
