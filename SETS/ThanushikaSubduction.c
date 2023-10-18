#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"
#include "stdio.h"

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


double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel   = -0.0e3 / instance->scaling.L;

  return TopoLevel;
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {

  // Define parameters and initialise phases according to coordinates
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

double SetTemperature(MdoodzInput *instance, Coordinates X) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (instance->model.user0 + 273.15) / instance->scaling.T;
  
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-X.z) + surfaceTemperature;
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

// Boundary conditions
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

// Mimicking heat loss by mantle convection via artificially high thermal conductivity during thermal equilibration steps
void AddCrazyConductivity(MdoodzInput *input) {
  int               *asthenospherePhases = (int *) malloc(sizeof(int));
  CrazyConductivity *crazyConductivity   = (CrazyConductivity *) malloc(sizeof(CrazyConductivity));
  asthenospherePhases[0]                 = 0;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000;
  input->crazyConductivity               = crazyConductivity;
}

// Smooth transition of vertical boundary velocity profile
double BoundaryVelocityProfile(double V, double z, double z_LAB, double dz_smooth) {
  return -0.5*V*erfc( (z-z_LAB)/dz_smooth ) + V;
}

double BoundaryVelocityProfilePrimitive(double V, double z, double z_LAB, double dz_smooth) {
  return V*(z - 0.5*(z-z_LAB) * erfc( (z-z_LAB)/dz_smooth ) );
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double dz_smooth = 10e3/instance->scaling.L;
  double V_tot    =  instance->model.user2/instance->scaling.V; // |VxW| + |VxE|
  double VxW, VxE;
  double z_LAB = -instance->model.user1/instance->scaling.L, z_top = instance->model.zmax;
  const double x = coordinates.x, z = coordinates.z;

  // Evaluate velocity of W and E boundaries
    VxW =  -0.5*V_tot;
    VxE =  0.5*V_tot;

  // Apply smooth transition with depth
  VxW = BoundaryVelocityProfile(VxW, z, z_LAB, dz_smooth);
  VxE = BoundaryVelocityProfile(VxE, z, z_LAB, dz_smooth);

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
  const double Lx = instance->model.xmax - instance->model.xmin, dz_smooth = 10e3/instance->scaling.L;
  double V_tot    =  instance->model.user2/instance->scaling.V; // |VxW| + |VxE|
  double tet_W, alp_W, tet_E, alp_E;
  double VxW, VxE, VzW, VzE, VzS=0.0;
  double z_LAB = -instance->model.user1/instance->scaling.L, z_top = instance->model.zmax;
  const double x = coordinates.x, z = coordinates.z;  

  // Evaluate velocity of W and E boundaries
    VxW = -0.5*V_tot;
    VxE =  0.5*V_tot;
    VzW =  0.0;
    VzE =  0.0;

  // Compute compensating VzS value
  const double z_min       = instance->model.zmin; 
  const double prim_W_zmax = BoundaryVelocityProfilePrimitive(VxW, z_top, z_LAB, dz_smooth);
  const double prim_W_zmin = BoundaryVelocityProfilePrimitive(VxW, z_min, z_LAB, dz_smooth);
  const double prim_E_zmax = BoundaryVelocityProfilePrimitive(VxE, z_top, z_LAB, dz_smooth);
  const double prim_E_zmin = BoundaryVelocityProfilePrimitive(VxE, z_min, z_LAB, dz_smooth);
  const double intW = prim_W_zmax - prim_W_zmin;
  const double intE = prim_E_zmax - prim_E_zmin;
  VzS = - ( fabs(intW) + fabs(intE) ) / Lx;

  // Apply smooth transition with depth
  VzW = -BoundaryVelocityProfile(VzW, z, z_LAB, dz_smooth);
  VzE = -BoundaryVelocityProfile(VzE, z, z_LAB, dz_smooth);

  if (position == W || position == SW || position == NW ) {
    bc.value = VzW;
    bc.type  = 11;
  } else if ( position == E || position == SE || position == NE) {
    bc.value = VzE;
    bc.type  = 11;
  } else if (position == S || position == N) {
    bc.value = VzS;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}

// Main function applies all of the above defined
int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "ThanushikaSubduction.txt"); // Default
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
                  .SetBCVx    = SetBCVx, // SetPureShearBCVx
                  .SetBCVz    = SetBCVz, // SetPureShearBCVz
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}