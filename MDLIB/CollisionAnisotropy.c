#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

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
          .radiusX = 10e3,
          .centreX = centreX,
          .centreZ = centreZ,
          .angle   = 0,
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
  const bool   isMicaEllipse         = IsEllipseCoordinates(coordinates, GetMicaEllipse(0e3, -55e3), instance->scaling.L);

  if (isMicaEllipse) {
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

double SetAnisoAngle(MdoodzInput *input, Coordinates coordinates, int phase, double predefined_angle) {
  static unsigned int seedIncrement = 0;
  srand(time(NULL) + seedIncrement++);

  // Use coordinates to determine the base angle, ensuring smooth variation
  double baseAngle = fmod((coordinates.x + coordinates.z), 180.0);

  // Generate a small random offset, for example within [-10, 10] degrees
  double offsetRange = 10.0; // Adjust this value to control the degree of randomness
  double randomOffset = ((rand() / (double)RAND_MAX) * 2 * offsetRange) - offsetRange;

  // Combine base angle with random offset
  double angle = baseAngle + randomOffset;

  // Ensure angle is within desired range, e.g., [0, 180]
  if (angle < 0.0) angle += 180.0;
  else if (angle > 180.0) angle -= 180.0;

  return angle;
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
  srand(time(NULL));
  char *input_file;
  if (nargs < 2) {
    asprintf(&input_file, "RiftingAnisotropy.txt");// Default
  } else {
    printf("dodo %s\n", args[1]);
    asprintf(&input_file, "%s", args[1]);// Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetTemperature = SetTemperature,
                  .SetGrainSize   = SetGrainSize,
                  .SetNoise       = SetNoise,
                  .SetAnisoAngle  = SetAnisoAngle,

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