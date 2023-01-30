#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // Set checkerboard for phase 0
    Az = sin( 48.0*2.0*M_PI*coordinate.z / Lz  );
    if ( Az>0.0 && dual_phase==0 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Az = sin( 48.0*2.0*M_PI*coordinate.z / Lz  );
    if ( Az>0.0 && dual_phase==1 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 2
    Az = sin( 48.0*2.0*M_PI*coordinate.z / Lz  );
    if ( Az>0.0 && dual_phase==2 ) {
        dual_phase += input->model.Nb_phases;
    }

  return dual_phase;
}


double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel   = -0.0e3 / instance->scaling.L;
  const double basin_width = 30.0e3 / instance->scaling.L;
  const double h_pert      = instance->model.user3 / instance->scaling.L;
  const double x           = x_coord - instance->model.user4 / instance->scaling.L;
  return TopoLevel + 0*h_pert * exp( - (x*x) / (2.0*basin_width*basin_width)   );
}

double SetNoise(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const double x        = coordinates.x - instance->model.user4 / instance->scaling.L;
  const double z        = coordinates.z;
  const double basin_width = 30.0e3 / instance->scaling.L;
  const double noise = ((double) rand() / (double) RAND_MAX) - 0.5;
  const double filter_x = exp( - (x*x)/ (2.0*basin_width*basin_width) );
  const double filter_z = exp( - (z*z)/ (2.0*basin_width*basin_width*4.0) );
  return  noise * filter_x * filter_z;
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double basin_width           = 30.0e3 / instance->scaling.L;
  const double h_pert                = instance->model.user3 / instance->scaling.L;
  const double lithosphereThickness  = instance->model.user1 / instance->scaling.L;
  const double crustThickness        = instance->model.user2 / instance->scaling.L;
  const double perturbationAmplitude = instance->model.user3 / instance->scaling.L;
  const double x                     = coordinates.x - instance->model.user4 / instance->scaling.L;
  const double mohoLevel             = -crustThickness + 0*h_pert * exp( - (x*x) / (2.0*basin_width*basin_width)   );
  const bool   isBelowLithosphere    = coordinates.z < -lithosphereThickness;
  const bool   isAboveMoho           = coordinates.z > mohoLevel;
  const bool   isLowerCrust          = coordinates.z < (mohoLevel + 0.5*crustThickness);

  if (isAboveMoho) {
    if (isLowerCrust) {
    return 1;
    }
    else {
    return 0;
    }
  } else if (isBelowLithosphere) {
    return 3;
  } else {
    return 2;
  }
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (1330.0 + 273.15) / instance->scaling.T;
  
  const double Tamp               =  50.0 / instance->scaling.T;
  const double x                  = coordinates.x - instance->model.user4 / instance->scaling.L;
  const double z                  = coordinates.z;
  const double basin_width        = 30.0e3 / instance->scaling.L;
  const double filter_x           = Tamp * exp( - (x*x)/ (2.0*basin_width*basin_width) );
  
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature + filter_x;
  } else {
    return particleTemperature;
  }
}

double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const int asthenospherePhase = 3;
  return instance->materials.gs_ref[asthenospherePhase];
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
  if (position == FREE_SURFACE) {
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
  const double surfaceTemperature = zeroC / instance->scaling.T;
  const double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == S || position == SE || position == SW) {
    bc.value = particleTemperature; // + filter_x;
    bc.type  = 1;
  } else if (position == N || position == NE || position == NW) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else if (position == W || position == E) {
    bc.value = 0.0; // mantleTemperature
    bc.type  = 0;
  } else {
    bc.value = 0.0;
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

int main() {
  srand(69); // Force random generator seed for reproducibility 
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetGrainSize          = SetGrainSize,
                  .SetNoise              = SetNoise,
                  .SetDualPhase          = SetDualPhase,

          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetPureShearBCVx,
                  .SetBCVz    = SetPureShearBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ("NeckingReview.txt", &setup);
}
