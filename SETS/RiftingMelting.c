#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"
#include "stdio.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double basin_width           = 30.0e3 / input->scaling.L;
  const double h_pert                = input->model.user3 / input->scaling.L;
  const double lithosphereThickness  = input->model.user1 / input->scaling.L;
  const double crustThickness        = input->model.user2 / input->scaling.L;
  const double perturbationAmplitude = input->model.user3 / input->scaling.L;
  const double x                     = coordinates.x;// - input->model.user4 / input->scaling.L;
  const double mohoLevel             = -crustThickness + 0*h_pert * exp( - (x*x) / (2.0*basin_width*basin_width)   );
  const double w                     = input->model.user5 / input->scaling.L;
  const double h_LAB                 = 0.0e3 / input->scaling.L;
  const double LABLevel              = -lithosphereThickness + h_LAB * exp( -(x*x) / (2.0*w*w)   );
  const bool   isBelowLithosphere    = coordinates.z < LABLevel;
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

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // Set checkerboard for phase 0
    Ax = cos(         36.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( (Lz/Lx)*36.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Ax = cos(         36.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( (Lz/Lx)*36.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 2
    Ax = cos(         36.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( (Lz/Lx)*36.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==2 ) {
        dual_phase += input->model.Nb_phases;
    }

  return dual_phase;
}


double SetSurfaceZCoord(MdoodzInput *input, double x_coord) {
  const double TopoLevel   = -0.0e3 / input->scaling.L;
  const double basin_width = 30.0e3 / input->scaling.L;
  const double h_pert      = input->model.user3 / input->scaling.L;
  const double x           = x_coord - input->model.user4 / input->scaling.L;
  return TopoLevel + h_pert * exp( - (x*x) / (2.0*basin_width*basin_width)   );
}

double SetNoise(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double x        = coordinates.x - input->model.user4 / input->scaling.L;
  const double z        = coordinates.z;
  const double basin_width = 30.0e3 / input->scaling.L;
  const double noise = ((double) rand() / (double) RAND_MAX) - 0.5;
  const double filter_x = exp( - (x*x)/ (2.0*basin_width*basin_width) );
  const double filter_z = exp( - (z*z)/ (2.0*basin_width*basin_width*4.0) );
  return  noise * filter_x * filter_z;
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  const double lithosphereThickness = input->model.user1 / input->scaling.L;
  const double surfaceTemperature   = 273.15 / input->scaling.T;
  const double mantleTemperature    = (1330.0 + 273.15) / input->scaling.T;
  const double Tamp               = 50.0 / input->scaling.T;
  const double x                  = coordinates.x - input->model.user4 / input->scaling.L;
  const double z                  = coordinates.z;
  const double basin_width        = 30.0e3 / input->scaling.L;
  const double filter_x           = Tamp * exp( - (x*x)/ (2.0*basin_width*basin_width) );
  
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature + filter_x;
  } else {
    return particleTemperature;
  }
}

double SetGrainSize(MdoodzInput *input, Coordinates coordinates, int phase) {
  const int asthenospherePhase = 3;
  return input->materials.gs_ref[asthenospherePhase];
}

char SetBCPType(MdoodzInput *input, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC SetBCT(MdoodzInput *input, POSITION position, Coordinates coordinates, double particleTemperature) {
  SetBC     bc;
  const double surface_temperature =          zeroC  / input->scaling.T;
  const double mantle_temperature  = (1330. + zeroC) / input->scaling.T;
  const double qz0 =  -18e-3 / (input->scaling.W/pow(input->scaling.L,2));
  const double qz1 =   input->model.user6  / (input->scaling.W/pow(input->scaling.L,2));
  const double w   =   input->model.user5 /input->scaling.L; 
  const double qz  = qz0 + qz1 * exp( -(coordinates.x*coordinates.x)/(2*w*w) );
  // printf("%2.2e %2.2e\n", coordinates.x, qz*(input->scaling.W/pow(input->scaling.L,2)));
  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = particleTemperature;
    // bc.type  = constant_heatflux;
    // bc.value = qz;
  }
  if (position == free_surface || position == N) {
    bc.type  = constant_temperature;
    bc.value = surface_temperature;
  } 
  if (position == W || position == E) {
    bc.type  = constant_heatflux;
    bc.value = 0.;
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

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "RiftingMelting.txt"); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  srand(69); // Force random generator seed for reproducibility 
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetGrainSize          = SetGrainSize,
                  // .SetNoise              = SetNoise,
                  .SetDualPhase          = SetDualPhase,

          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetPureShearBCVx,
                  .SetBCVz    = SetPureShearBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}