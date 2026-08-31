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
    double f = 8.;

    // Set checkerboard for phase 0
    Ax = cos( 4*f*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( f*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 && phase==0 ) {
        dual_phase = input->model.Nb_phases;
    }

    // Set checkerboard for phase 2
    Az = sin( 3*f*2.0*M_PI*coordinate.z / Lz  );
    Ax = cos( 4*3*f*2.0*M_PI*coordinate.x / Lx  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 && phase==1 ) {
        dual_phase = input->model.Nb_phases+1;
    }

  return dual_phase;
}


// Initial free-surface level. Also used by SetBCVz to bound the boundary inflow integral,
// so that the velocity boundary conditions conserve volume.
double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel   = 0.0e3 / instance->scaling.L;

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

SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates, double particleTemperature) {
  SetBC     bc;
  double surface_temperature = (0.0 + 273.15) / instance->scaling.T ;
  double mantle_temperature  = (instance->model.user0 + 273.15) / instance->scaling.T;
  (void)coordinates;
  (void)particleTemperature;
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

// Smooth transition of horizontal boundary velocity with depth: V in the plate (z > z_LAB), 0 below
double BoundaryVelocityProfile(double V, double z, double z_LAB, double dz_smooth) {
  return -0.5*V*erfc( (z-z_LAB)/dz_smooth ) + V;
}

// Exact antiderivative of BoundaryVelocityProfile with respect to z
double BoundaryVelocityProfilePrimitive(double V, double z, double z_LAB, double dz_smooth) {
  const double u = (z-z_LAB)/dz_smooth;
  return V*( z - 0.5*dz_smooth*( u*erfc(u) - exp(-u*u)/sqrt(M_PI) ) );
}

// Plate velocities imposed on the W and E boundaries.
// user2 = V_tot = VxE - VxW: total convergence (<0) or divergence (>0) rate [m/s]
// user3 = fraction of V_tot applied on the W boundary (0.5: symmetric push)
void BoundaryPlateVelocities(MdoodzInput *instance, double *VxW, double *VxE) {
  const double V_tot = instance->model.user2/instance->scaling.V;
  const double fW    = instance->model.user3;
  *VxW = -fW*V_tot;
  *VxE = (1.0-fW)*V_tot;
}

// Vertical velocity on the S boundary that balances the net volume flux through the W and E boundaries.
// The inflow profiles are integrated from the base of the model up to the free surface (not up to zmax):
// the nodes above the free surface are "air" and do not transport material.
double BalancingBottomVelocity(MdoodzInput *instance) {
  const double Lx        = instance->model.xmax - instance->model.xmin;
  const double z_min     = instance->model.zmin;
  const double z_LAB     = -instance->model.user1/instance->scaling.L;
  const double dz_smooth = 10e3/instance->scaling.L;
  const double z_surf_W  = SetSurfaceZCoord(instance, instance->model.xmin);
  const double z_surf_E  = SetSurfaceZCoord(instance, instance->model.xmax);
  double VxW, VxE;
  BoundaryPlateVelocities(instance, &VxW, &VxE);
  // Volume fluxes (per unit length in y), positive in +x direction
  const double fluxW = BoundaryVelocityProfilePrimitive(VxW, z_surf_W, z_LAB, dz_smooth) - BoundaryVelocityProfilePrimitive(VxW, z_min, z_LAB, dz_smooth);
  const double fluxE = BoundaryVelocityProfilePrimitive(VxE, z_surf_E, z_LAB, dz_smooth) - BoundaryVelocityProfilePrimitive(VxE, z_min, z_LAB, dz_smooth);
  const double net_inflow = fluxW - fluxE;
  // Mass conservation: what enters through the sides must leave through the base (VzS < 0 for net inflow)
  return -net_inflow / Lx;
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double dz_smooth = 10e3/instance->scaling.L;
  const double z_LAB     = -instance->model.user1/instance->scaling.L;
  const double z         = coordinates.z;
  double VxW, VxE;

  // Evaluate velocity of W and E boundaries
  BoundaryPlateVelocities(instance, &VxW, &VxE);

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
  const double VzW = 0.0, VzE = 0.0;
  const double VzS = BalancingBottomVelocity(instance);

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
    asprintf(&input_file, "AnneloreSubduction.txt"); // Default
  }
  else {
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
