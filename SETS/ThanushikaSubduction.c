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
  double TopoLevel           = -0.0e3 / instance->scaling.L;
  const double seaLevel      = 0.0e3 / instance->scaling.L;
  const double ht_iso        =  4.8e3 / instance->scaling.L;
  const double weakZoneWidth = 10e3/instance->scaling.L;
  const double x_ext_start   =  0.0e3 / instance->scaling.L;
  const double x_ext_end     =  100.0e3 / instance->scaling.L;
  const double marginWidth   =  x_ext_end - x_ext_start;

  if (x_coord >= x_ext_start && x_coord <= x_ext_end)
  {
    // TopoLevel = ht_iso / 2.0 * ( 1.0 + tanh( ((x_coord - marginWidth / 2.0) / (2.0 * marginWidth)) ) );
    TopoLevel = seaLevel + (ht_iso - seaLevel) / marginWidth * x_coord;
  }
  if (x_coord >= x_ext_end)
  {
    TopoLevel = ht_iso;
  }

  return TopoLevel;
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {

  // Define parameters and initialise phases according to coordinates
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double weakZoneWidth        =  10.0e3/instance->scaling.L;
  // const double weakZoneDepth        = -120e3 / instance->scaling.L;
        double mohoLevel            = -instance->model.user0 / instance->scaling.L;
  const double seaLevel             = 0.0e3 / instance->scaling.L;
  const double ht_iso               = 4.8e3 / instance->scaling.L;
  const bool   isBelowLithosphere   = coordinates.z < -lithosphereThickness;
  const bool   isAboveMoho          = coordinates.z > mohoLevel;
  const bool   isContinentalPlate   = coordinates.x > 0.0e3 / instance->scaling.L;
  const double x_ext_start          =  0.0e3 / instance->scaling.L;
  const double x_ext_end            =  100.0e3 / instance->scaling.L;
  const double marginWidth          =  x_ext_end - x_ext_start;
  // const double subductionAngle      = tan((90.0 - instance->model.user4)/180*PI);
  const double subductionAngle      = instance->model.user4/180*PI;
  const double ocCrustThickness     = instance->model.user7 / instance->scaling.L;
  const double ocLithThickness      = instance->model.user5 / instance->scaling.L;
  const double sedimentThickness    = instance->model.user8 / instance->scaling.L;
  const double weakZoneDepth        = -fmin(lithosphereThickness, ocLithThickness) - 20.0e3 / instance->scaling.L;
  int phase = 0; // Upper mantle
  
  // Calculate distance of current point from center line of weak layer
  double d = (coordinates.x - x_ext_start) * sin(subductionAngle) - (-coordinates.z - 0.0) * cos(subductionAngle);

  // Set continental lithosphere
  if (coordinates.z >= -lithosphereThickness && d > weakZoneWidth / 2.0) 
  {
    phase = 1;
  }

  // Set continental crust
  if (coordinates.x >= x_ext_start && coordinates.x <= x_ext_end)
  {
    // mohoLevel = mohoLevel * ( 1.0 + tanh( ((coordinates.x - marginWidth / 2.0)) / (2.0 * marginWidth) ) );
    mohoLevel = seaLevel + (mohoLevel + ht_iso - seaLevel) / marginWidth * coordinates.x;
  }
  else
  {
    mohoLevel = mohoLevel + ht_iso;
  }
  if (coordinates.z >= mohoLevel && coordinates.x >= x_ext_start)
  {
    phase = 4;
  }
  
  // Set weak zone
  if (coordinates.z > weakZoneDepth && fabs(d) <= weakZoneWidth / 2.0)
  {
    phase = 3;
  }

  // Set oceanic lithosphere
  if (coordinates.z > -ocLithThickness && d < -weakZoneWidth / 2.0)
  {
    phase = 2;
  }

  // Set oceanic crust
  if (coordinates.z > -(ocCrustThickness + sedimentThickness) && d < -weakZoneWidth / 2.0)
  {
    phase = 5;
  }

  // Set sediments
  if (coordinates.z > -sedimentThickness && d < -weakZoneWidth / 2.0)
  {
    phase = 6;
  }
   
  // Return
  return phase;
  
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (1350.0 + 273.15) / instance->scaling.T;
  const double weakZoneWidth        =  10.0e3/instance->scaling.L;
  const double subductionAngle      = tan((90.0 - instance->model.user4)/180*PI);
  const double ocLithThickness      = instance->model.user5 / instance->scaling.L;
        double particleTemperature  = mantleTemperature;
        double thermalDiffusivity   = 1.0e-6 / (instance->scaling.L*instance->scaling.L / instance->scaling.t);
  const double adiabaticGradient    = instance->model.user6 / (instance->scaling.T / instance->scaling.L);
  const double crustThickness       = instance->model.user0 / instance->scaling.L;
  const double mohoTemperature      = (400.0 + 273.15) / instance->scaling.T;

  // Set conductive gradient in oceanic lithosphere (lithospheric thickness is precalculated based on thermal age)
  particleTemperature  = ((mantleTemperature - surfaceTemperature) / ocLithThickness) * (-coordinates.z) + surfaceTemperature;
  
  // Set temperature in the continental crust
  if (coordinates.x > weakZoneWidth - coordinates.z / tan(subductionAngle) && coordinates.z >= -crustThickness )
  {
    particleTemperature  = ((mohoTemperature - surfaceTemperature) / crustThickness) * (-coordinates.z) + surfaceTemperature;
  }
  
  // Set conductive thermal gradient in the continental lithosphere
  if (coordinates.x > weakZoneWidth - coordinates.z / tan(subductionAngle) && coordinates.z >= -lithosphereThickness && coordinates.z < -crustThickness)
  {
      particleTemperature  = ((mantleTemperature - mohoTemperature) / (lithosphereThickness - crustThickness)) * (-coordinates.z) + mohoTemperature;
  }
  
  // Set mantle adiabat
  if (particleTemperature >= mantleTemperature)
  {
    particleTemperature = adiabaticGradient * (-coordinates.z) + mantleTemperature;
  }

  // Return
  return particleTemperature;
}

int AdjustPhaseToTemperature(MdoodzInput *instance, Coordinates coordinates, double particleTemperature, int phase) {

  // Define parameters and initialise phases according to coordinates
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double weakZoneWidth        =  10.0e3 / instance->scaling.L;
  const double seaLevel             =  0.0e3 / instance->scaling.L;
  const double mohoLevel            = -instance->model.user0 / instance->scaling.L;
  const double subductionAngle      = tan(instance->model.user4/180*PI);
  const double mantleTemperature    = (1350.0 + 273.15) / instance->scaling.T;
  
  // Correct oceanic lithosphere for temperature
  if (particleTemperature < mantleTemperature - 0.01*mantleTemperature && phase == 0 && coordinates.x <= -weakZoneWidth - coordinates.z * subductionAngle)
  {
    phase = 2;
  }
  
  // Return
  return phase;
  
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

SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates,  double particleTemperature) {
  SetBC     bc;
  double surface_temperature = (0.0 + 273.15) / instance->scaling.T ;
  double mantle_temperature  = (1350.0 + 273.15) / instance->scaling.T;
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double ocLithThickness      = instance->model.user5 / instance->scaling.L;
  const double minLithThickness     = fmin(ocLithThickness, lithosphereThickness);
  const double upperMantleThickness = fabs(instance->model.zmin) - minLithThickness;
  const double adiabaticGradient    = instance->model.user6 / (instance->scaling.T / instance->scaling.L);

  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = mantle_temperature + upperMantleThickness * adiabaticGradient;
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
  const double ocLithThickness           = input->model.user5 / input->scaling.L;
  const double ocRefLithThick            = 120.0e3 / input->scaling.L;
  asthenospherePhases[0]                 = 0;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000.0;//ocRefLithThick * 25.0 / ocLithThickness;
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
  const double pushFactor = instance->model.user3;

  // Evaluate velocity of W and E boundaries
    VxW =  -pushFactor * V_tot;
    VxE =   (1 - pushFactor) * V_tot;

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
    VxW = -V_tot;
    VxE =  0.0;
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
                  .SetPhase                 = SetPhase,
                  .SetTemperature           = SetTemperature,
                  // .AdjustPhaseToTemperature = AdjustPhaseToTemperature,
                  .SetGrainSize             = SetGrainSize,
                  // .SetDualPhase          = SetDualPhase,

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