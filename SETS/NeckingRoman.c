#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"

double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel = -0.0e3 / instance->scaling.L;
  return TopoLevel;
}

double SetNoise(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const double x           = coordinates.x - instance->model.user4 / instance->scaling.L;
  const double z           = coordinates.z;
  const double basin_width = 30.0e3 / instance->scaling.L;
  const double noise       = ((double) rand() / (double) RAND_MAX) - 0.5;
  const double filter_x    = exp(-(x * x) / (2.0 * basin_width * basin_width));
  const double filter_z    = exp(-(z * z) / (2.0 * basin_width * basin_width * 4.0));
  return noise * filter_x * filter_z;
}

static Ellipse GetMicaEllipse(double centreX, double centreZ) {
  return (Ellipse) {
          .radiusZ = 10e3,
          .radiusX = 50e3,
          .centreX = centreX,
          .centreZ = centreZ,
          .angle   = 30,
  };
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  // remove perturbation
  // moho temperature should be around 500 degrees
  // add weak inclusion
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double crustThickness       = instance->model.user2 / instance->scaling.L;
  const double mohoLevel            = -crustThickness;
  const bool   isBelowLithosphere   = coordinates.z < -lithosphereThickness;
  const bool   isAboveMoho          = coordinates.z > mohoLevel;
  // const bool   isMicaEllipse         = IsEllipseCoordinates(coordinates, GetMicaEllipse(0e3, -120e3), instance->scaling.L);

  bool   isSerpentineEllipse = false;

// radius
//     if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) isSerpentineEllipse = true;

  const double a_ell = 10e3/ instance->scaling.L;
  const double b_ell = 20e3/ instance->scaling.L;
  const double x0 = 0.0, z0 = -120e3/instance->scaling.L;
  const double angle = 30*M_PI/180.0;
  const double x_ell = (coordinates.x - x0)*cos(angle) + (coordinates.z - z0)*sin(angle);
  const double z_ell =-(coordinates.x - x0)*sin(angle) + (coordinates.z - z0)*cos(angle);
  if (pow(x_ell/a_ell,2.0) + pow(z_ell/b_ell,2.0) < 1.0) isSerpentineEllipse = true;


  if (isSerpentineEllipse) {
    return 1;
  } else if (isAboveMoho) {
    return 0;
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

  const double Tamp                 = 50.0 / instance->scaling.T;
  const double x                    = coordinates.x - instance->model.user4 / instance->scaling.L;
  const double z                    = coordinates.z;
  const double basin_width          = 30.0e3 / instance->scaling.L;
  const double filter_x             = Tamp * exp(-(x * x) / (2.0 * basin_width * basin_width));

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
  SetBC  bc;
  double surface_temperature = zeroC / instance->scaling.T;
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
  if (nargs < 2) {
    asprintf(&input_file, "NeckingRoman.txt");// Default
  } else {
    printf("dodo %s\n", args[1]);
    asprintf(&input_file, "%s", args[1]);// Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  srand(69);// Force random generator seed for reproducibility
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetTemperature = SetTemperature,
                  .SetGrainSize   = SetGrainSize,
                  .SetNoise       = SetNoise,

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
