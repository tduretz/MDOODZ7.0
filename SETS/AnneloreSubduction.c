#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"
#include "stdio.h"

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
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


double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel   = -0.0e3 / instance->scaling.L;
  // const double basin_width = 30.0e3 / instance->scaling.L;
  // const double h_pert      = instance->model.user3 / instance->scaling.L;
  // const double x           = x_coord - instance->model.user4 / instance->scaling.L;
  return TopoLevel;// + h_pert * exp( - (x*x) / (2.0*basin_width*basin_width)   );
}

int SetPhase(MdoodzInput *instance, Coordinates X) {
  // Define parameters
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double weakZoneWidth        = 10e3/instance->scaling.L;
  const double mohoLevel            = -5e3 / instance->scaling.L;
  const bool   isBelowLithosphere   = X.z < -lithosphereThickness;
  const bool   isAboveMoho          = X.z > mohoLevel;
  int phase = 0;
  if (X.z > - lithosphereThickness) phase = 1;
  return phase;
  
}

double SetTemperature(MdoodzInput *instance, Coordinates X) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (instance->model.user0 + 273.15) / instance->scaling.T;
  
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-X.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature;
  } else {
    return particleTemperature;
  }
}

double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const int asthenospherePhase = 0;
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

void AddCrazyConductivity(MdoodzInput *input) {
  int               *asthenospherePhases = (int *) malloc(sizeof(int));
  CrazyConductivity *crazyConductivity   = (CrazyConductivity *) malloc(sizeof(CrazyConductivity));
  asthenospherePhases[0]                 = 0;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000;
  input->crazyConductivity               = crazyConductivity;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "AnneloreSubduction.txt"); // Default
  }
  else {
    printf("dodo %s\n", args[1]);
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetGrainSize          = SetGrainSize,
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
