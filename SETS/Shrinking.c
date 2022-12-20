#include "mdoodz.h"
<<<<<<< HEAD
#include "math.h"
=======
>>>>>>> main

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = 0.25 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

<<<<<<< HEAD
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

=======
>>>>>>> main
double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  const double P_init = (input->model.PrBG         ) / input->scaling.S;
  if (input->model.eqn_state > 0) {
    return input->materials.rho[phase] * exp(input->materials.bet[phase]*P_init -  input->materials.alp[phase] * T_init);
  } else {
    return input->materials.rho[phase];
  }
}

<<<<<<< HEAD
=======
void SetDefGrad(MdoodzInput *input, Coordinates coordinates, int phase, Tensor2D* F) {
  const double radius = input->model.user1 / input->scaling.L;
  F->xx=1.; F->xz=0.; F->zx=0.; F->zz=1.;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    // nothing special in the inclusion
  }
  else {
    // a bit more pure shear strain in the matrix
    F->xx=1.1; F->xz=0.; F->zx=0.; F->zz=0.9;
  }
}

>>>>>>> main
int main() {
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
<<<<<<< HEAD
                   .SetDualPhase          = SetDualPhase,
                   .SetDensity            = SetDensity,
=======
                   .SetDensity            = SetDensity,
                   .SetDefGrad            = SetDefGrad,
>>>>>>> main
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("Shrinking.txt", &setup);
}
