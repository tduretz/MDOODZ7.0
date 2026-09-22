#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"


double Vx_S = -0.5;
double Vx_N = 0.5;

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

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC           bc;  
  if (position == W || position == E) {
    bc.type  = 0;      
    bc.value = coordinates.z;
  } else if (position == S || position == SE || position == SW) {
    bc.type  = 11;
    bc.value = Vx_S; // vitesse non-nulle en x au bord S
  } else if (position == N || position == NE || position == NW) {
    bc.type  = 11;
    bc.value = Vx_N; // vitesse non-nulle en x au bord N
  } else {
    bc.type  = -1;
    bc.value = 0.0;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  if (position == N || position == S || position == NE || position == NW || position == SE || position == SW) {
    bc.type  = 0;
    bc.value = 0.;
  }
  else if (position == W) {
    bc.type  = 11;
    bc.value = 0.; // vitesse nulle en z sur le bord W
  }
  else if (position == E) {
    bc.type  = 11;
    bc.value = 0.; // vitesse nulle en z sur le bord E
  } else {
    bc.type  = -1;
    bc.value = 0.0;
  }
  return bc;
}


int main() {
  MdoodzSetup instance = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase   = SetPhase,
                  .SetDensity = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ("noPerBC_ShearTemplateAnisotropy.txt", &instance);
}