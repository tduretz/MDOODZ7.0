#include "complex.h"
#include "mdoodz.h"
#include "math.h"

int SetPhase(MdoodzInput *instance, Coordinates coords) {

  int phase = 0.;
  const double x0     = -500e3/instance->scaling.L;
  const double x1     = -400e3/instance->scaling.L;
  const double z_LAB  = -100e3/instance->scaling.L;
  const double z_slab = -200e3/instance->scaling.L;

  if (coords.x>x0 && coords.z>z_LAB) phase = 1;

  if (coords.x>x0 && coords.x<x1 && coords.z>z_slab) phase = 1;
  
  if (coords.z > 0.) phase = 2;
  
  // Return
  return phase;
    
}

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    // Passive tracer function. Useful for visualisation
    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;
    double f = 4.;

    // Set checkerboard for phase 0
    Ax = cos( f*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( f*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 && phase==0 ) {
        dual_phase = input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Az = sin( 3*f*2.0*M_PI*coordinate.z / Lz  );
    Ax = cos( 3*f*2.0*M_PI*coordinate.x / Lx  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 && phase==1 ) {
        dual_phase = input->model.Nb_phases;
    }

    return dual_phase;
}

double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  return instance->materials.rho[phase];
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;  
  if (position == W || position == E) {
    bc.type  = 0;
    bc.value = 0.;
  } else if (position == S || position == SE || position == SW) {
    bc.type  = 13;
    bc.value = 0.0;
  } else if (position == N || position == NE || position == NW) {
    bc.type  = 13;
    bc.value = 0.0;
  } else {
    bc.type  = -1;
    bc.value = 0.0;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  if (position == N || position == NE || position == NW || position == S || position == SE || position == SW) {
    bc.type  = 0;
    bc.value = 0.;
  }
  else if (position == W) {
    bc.type  = 13;
    bc.value = 0.;
  }
  else if (position == E) {
    bc.type  = 13;
    bc.value = 0.;
  } else {
    bc.type  = -1;
    bc.value = 0.0;
  }
  return bc;
}

int main() {
  MdoodzSetup instance = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase     = SetPhase,
                  .SetDualPhase = SetDualPhase,
                  .SetDensity   = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ("SubductionBenchmarkCase2.txt", &instance);
}
