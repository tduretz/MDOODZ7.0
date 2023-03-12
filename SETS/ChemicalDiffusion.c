#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = 0.25 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetXComponent(MdoodzInput *input, Coordinates coordinates, int phase) {
  double X = 0.0;
  if (phase == 1) X = 0.1;
  return X;
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

SetBC SetBCC(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  // if (position == W  || position == N || position == E || position == S) {
  //   bc.type  = constant_X;
  //   bc.value = 0.1;
  // }
  if (position == S || position == N) {
    bc.type  = constant_X;
    bc.value = 0.0;
  }
  if (position == E || position == W) {
    bc.type  = periodic_X_x;
    bc.value = 0.0;
  }
  return bc;
}

double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return 2.0e1;
}

int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetHorizontalVelocity = SetHorizontalVelocity,
                   .SetPhase              = SetPhase,
                   .SetDualPhase          = SetDualPhase,
                   .SetDensity            = SetDensity,
                   .SetXComponent         = SetXComponent,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
                  .SetBCC  = SetBCC,
          },
  };
  RunMDOODZ("ChemicalDiffusion.txt", &setup);
}
