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
          .radiusZ = 20e3,
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
  const double decol = instance->model.user1 / instance->scaling.L;
  const double xmax  = instance->model.xmax;
  const double zmin  = instance->model.zmin;
  const double basement_step = instance->model.user3 / instance->scaling.L;
  const double angle = (instance->model.user0 * M_PI/180.);
  const double basement_top_ordinate = zmin -  xmax*tan(angle);
  const double decol_top_ordinate    = zmin -  xmax*tan(angle) + decol;
  bool   is_decol    =  coordinates.z < coordinates.x*tan(angle) + decol_top_ordinate;
  bool   is_basement =  coordinates.z < coordinates.x*tan(angle) + basement_top_ordinate;
  
  // basement steps
  const double step1 =  25e3 / instance->scaling.L;
  const double step2 = -25e3 / instance->scaling.L;
  const double step3 = -75e3 / instance->scaling.L;

  if (coordinates.x<step1 && coordinates.x>step2) is_basement =  coordinates.z < coordinates.x*tan(angle) + basement_top_ordinate + 1*basement_step;
  if (coordinates.x<step2 && coordinates.x>step3) is_basement =  coordinates.z < coordinates.x*tan(angle) + basement_top_ordinate + 2*basement_step;
  if (coordinates.x<step3                       ) is_basement =  coordinates.z < coordinates.x*tan(angle) + basement_top_ordinate + 3*basement_step;

  if (coordinates.x<step1 && coordinates.x>step2) is_decol    =  coordinates.z < coordinates.x*tan(angle) + decol_top_ordinate    + 1*basement_step;
  if (coordinates.x<step2 && coordinates.x>step3) is_decol    =  coordinates.z < coordinates.x*tan(angle) + decol_top_ordinate    + 2*basement_step;
  if (coordinates.x<step3                       ) is_decol    =  coordinates.z < coordinates.x*tan(angle) + decol_top_ordinate    + 3*basement_step;

  int phase = 0;               // sediments
  if (is_decol)    phase = 1; // décollement
  if (is_basement) phase = 2; // basement
  return phase;
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

SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates,  double particleTemperature) {
  SetBC  bc;
  double gradient = instance->model.user2 / (instance->scaling.T/instance->scaling.L);
  const double zmin  = instance->model.zmin;
  double surface_temperature = zeroC / instance->scaling.T;
  double base_temperature  = zmin*gradient + zeroC/instance->scaling.T;
  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = base_temperature;
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


SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;

  const double L     = instance->model.xmax - instance->model.xmin;
  const double Edot  = instance->model.bkg_strain_rate;
  const double angle = (instance->model.user0 * M_PI/180.);
  const double VxE = -L*Edot*cos(angle);

  // Assign BC values
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W) {
    bc.value = 0.0;
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

  const double L     = instance->model.xmax - instance->model.xmin;
  const double Edot  = instance->model.bkg_strain_rate;
  const double angle = (instance->model.user0 * M_PI/180.);
  const double VzE = -L*Edot*sin(angle);

  if (position == W || position == SW || position == NW ) {
    bc.value = 0.0;
    bc.type  = 13;
  } else if ( position == E || position == SE || position == NE) {
    bc.value = VzE;
    bc.type  = 11;
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}


int main(int nargs, char *args[]) {
  srand(time(NULL));
  char *input_file;
  if (nargs < 2) {
    asprintf(&input_file, "WedgeAnisotropy.txt");// Default
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
                  .SetDualPhase   = SetDualPhase,
                  .SetTemperature = SetTemperature,
                  .SetGrainSize   = SetGrainSize,
                  .SetNoise       = SetNoise,
                  .SetAnisoAngle  = SetAnisoAngle,

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