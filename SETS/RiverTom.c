#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"

double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  double TopoLevel = 0.0;
  if (fabs(x_coord) < instance->model.surf_Winc/2.0) TopoLevel = -500.0/instance->scaling.L;
  return TopoLevel;
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double H_crust  = instance->model.user1 / instance->scaling.L;
  if (coordinates.z> -H_crust) return 0;
  else return 1;
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

  return dual_phase;
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (instance->model.user0 + 273.15) / instance->scaling.T;
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature;
  } else {
    return particleTemperature;
  }
}

double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.bkg_strain_rate;
}

double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.bkg_strain_rate;
}

char SetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == free_surface) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else {
    bc.value = 0.0;
    bc.type  = 0;
  }
  return bc;
}


SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC     bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  bc.value = 0.;
  bc.type  = 0;
  if (position == S || position == SE || position == SW) {
    bc.value = particleTemperature;
    bc.type  = 1;
  } 
  if (position == N || position == NE || position == NW) {
    bc.value = particleTemperature;
    bc.type  = 1;
  } 
  if (position == W || position == E ) {
    bc.value = mantleTemperature;
    bc.type  = 0;
  }
  return bc;
}

void AddCrazyConductivity(MdoodzInput *input) {
  int               *asthenospherePhases = (int *) malloc(sizeof(int));
  CrazyConductivity *crazyConductivity   = (CrazyConductivity *) malloc(sizeof(CrazyConductivity));
  asthenospherePhases[0]                 = 3;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000;
  input->crazyConductivity               = crazyConductivity;
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lx    = instance->model.xmax - instance->model.xmin;
  const double V_tot =  Lx * instance->model.bkg_strain_rate; 
  const double x = coordinates.x, z = coordinates.z;

  // Evaluate velocity of W and E boundaries
  const double VxW =  0.5*V_tot;
  const double VxE = -0.5*V_tot;
  
  // Assign BC values
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W) {
    bc.value = VxW;
    bc.type  = 0;
  } else if (position == E) {
    bc.value = VxE;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lx = instance->model.xmax - instance->model.xmin;
  double V_tot    =  Lx * instance->model.bkg_strain_rate; // |VxW| + |VxE|
 
  // Compute bottom boundary velocity and account for erosion
  const double surf_Winc = instance->model.surf_Winc;
  const double surf_Vinc = instance->model.surf_Vinc;
  const double Vz_corr   = surf_Winc*surf_Vinc / Lx;
  const double VzS = -V_tot * (0.0 - instance->model.zmin) / Lx + Vz_corr;

  // Set boundary nodes types and values
  if (position == W || position == SW || position == NW ) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if ( position == E || position == SE || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = VzS;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}

int main() {
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetHorizontalVelocity = SetHorizontalVelocity,
                  .SetVerticalVelocity   = SetVerticalVelocity,
                  .SetDualPhase          = SetDualPhase,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx, // manuel
                  .SetBCVz    = SetBCVz,
                  // .SetBCVx    = SetPureShearBCVx, // automatic
                  // .SetBCVz    = SetPureShearBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ("RiverTom.txt", &setup);
}
