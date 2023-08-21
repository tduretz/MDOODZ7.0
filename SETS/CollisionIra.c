#include "math.h"
#include "mdoodz.h"
#include "stdlib.h"
#include "stdio.h"


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
  Rectangle westernContinent = {
          .sizeZ   = 180e3,
          .sizeX   = 600e3,
          .centreZ = -90e3,
          .centreX = -500e3,
          .angle   = 0,
  };
  Rectangle easternContinent = {
          .sizeZ   = 180e3,
          .sizeX   = 600e3,
          .centreZ = -90e3,
          .centreX = 500e3,
          .angle   = 0,
  };
  Rectangle HOceanicPlate = {
          .sizeZ   = 60e3,
          .sizeX   = 400e3,
          .centreZ = -30e3,
          .centreX = 0e3,
          .angle   = 0,
  };

  if (IsRectangleCoordinates(coordinates, westernContinent, instance->scaling.L)) {
    if (coordinates.z > -40e3 / instance->scaling.L) {
      return 1;
    } else {
      return 2;
    }
  } else if (IsRectangleCoordinates(coordinates, easternContinent, instance->scaling.L)) {
    if (coordinates.z > -30e3 / instance->scaling.L) {
      return 1;
    } else {
      return 2;
    }
  } else if (IsRectangleCoordinates(coordinates, HOceanicPlate, instance->scaling.L)) {
    if (coordinates.z > -20e3 / instance->scaling.L) {
      return 4;
    } else {
      return 2;
    }
  } else {
    return 3;
  }
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (1330.0 + 273.15) / instance->scaling.T;

  const double Tamp               =  50.0 / instance->scaling.T;
  const double x                  = coordinates.x - instance->model.user4 / instance->scaling.L;
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
  double surface_temperature =          zeroC  / instance->scaling.T;
  double mantle_temperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = mantle_temperature;
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

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  const double plateThickness = 180e3;
  SetBC bc;
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W) {
    if (coordinates.z > -plateThickness / input->scaling.L) {
      bc.value = -coordinates.x * input->model.bkg_strain_rate;
    } else {
      bc.value = 0.0;
    }
    bc.type  = 0;
  } else if (position == E) {
    if (coordinates.z > -200e3 / input->scaling.L) {
      bc.value = (-coordinates.x / 1) * input->model.bkg_strain_rate;
    } else if (coordinates.z < -250e3 / input->scaling.L && coordinates.z > -400e3 / input->scaling.L) {
      bc.value = -coordinates.x * input->model.bkg_strain_rate;
    } else {
      bc.value = 0.0;
    }
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == N) {
    bc.value = coordinates.z * input->model.bkg_strain_rate;
    bc.type  = 0;
  } else if (position == S) {
    if (coordinates.x > -100e3 / input->scaling.L && coordinates.x < 100e3 / input->scaling.L) {
      bc.value = -coordinates.z * input->model.bkg_strain_rate;
    } else {
      bc.value = 0.0;
    }
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
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
    asprintf(&input_file, "CollisionIra.txt"); // Default
  }
  else {
    printf("dodo %s\n", args[1]);
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  srand(69); // Force random generator seed for reproducibility
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff) {
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetGrainSize          = SetGrainSize,
                  .SetNoise              = SetNoise,

          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx,
                  .SetBCVz    = SetBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
