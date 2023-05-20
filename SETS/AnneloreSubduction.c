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

  return TopoLevel;
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {

  // Define parameters
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double weakZoneWidth        = 10e3/instance->scaling.L;
  const double mohoLevel            = -5e3 / instance->scaling.L;
  const bool   isBelowLithosphere   = coordinates.z < -lithosphereThickness;
  const bool   isAboveMoho          = coordinates.z > mohoLevel;
  int phase = 0;
  
  // Set all lithosphere
  if (coordinates.z>-lithosphereThickness) 
  {
    phase = 1;
  }
  
  // Set weak zone
   if (coordinates.z>-lithosphereThickness && coordinates.x>-weakZoneWidth + coordinates.z*1.3 && coordinates.x<weakZoneWidth + coordinates.z*1.3)
   {
     phase = 2;
   }

  // Set crustal layer
  if (isAboveMoho)
  {
    phase = 2;
  }
  
  // Return
  return phase;
  
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (instance->model.user0 + 273.15) / instance->scaling.T;
  
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature;
  } 
  else {
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
  double surface_temperature = (0.0 + 273.15) / instance->scaling.T ;
  double mantle_temperature  = (instance->model.user0 + 273.15) / instance->scaling.T;
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
    bc.value = 0.0;
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
