#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"
#include "stdio.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  double TopoLevel = 0.0;
  // if (fabs(x_coord) < instance->model.surf_Winc/2.0) TopoLevel = -500.0/instance->scaling.L;
  return TopoLevel;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double H_crust  = instance->model.user1 / instance->scaling.L;
  const double H_uppercrust  = instance->model.user2 / instance->scaling.L;
    int phase = 2;
    
    if (coordinates.z> -H_crust) phase = 1;
    if (coordinates.z> -H_uppercrust) phase = 0;
  
  return phase;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    int    il, dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // Set checkerboard for phase 0
    Ax = cos( 50.0*M_PI*coordinate.x / Lx  );
    Az = sin( 80.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
        dual_phase += input->model.Nb_phases;
    }

    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 ) {
        dual_phase += input->model.Nb_phases;
    }
    
  return dual_phase;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (instance->model.user0 + 273.15) / instance->scaling.T;
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature;
  } else {
    return particleTemperature;
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.bkg_strain_rate;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.bkg_strain_rate;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lx    = instance->model.xmax - instance->model.xmin;
  const double V_tot =  Lx * instance->model.bkg_strain_rate; 
  const double x = coordinates.x, z = coordinates.z;

  // Evaluate velocity of W and E boundaries
  const double VxW =  0.5*V_tot;
  const double VxE = -0.5*V_tot;
  
  // Assign BC values
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.type  = constant_shear_stress;
    bc.value = 0.0;
  } else if (position == W) {
    bc.type  = constant_velocity;
    bc.value = VxW;
  } else if (position == E) {
    bc.type  = constant_velocity;
    bc.value = VxE;    
  } else {
    bc.type  = inside;
    bc.value = 0.0;  
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lx = instance->model.xmax - instance->model.xmin;
  double V_tot    =  Lx * instance->model.bkg_strain_rate; // |VxW| + |VxE|
 
  // Compute bottom boundary velocity and account for erosion
  const double surf_Winc = instance->model.surf_Winc;
  const double surf_Vinc = instance->model.surf_Vinc;
  const double Vz_corr   = surf_Winc*surf_Vinc / Lx;
  // const double hW        = instance->flux->west - instance->model.zmin;
  // const double hE        = instance->flux->east - instance->model.zmin;
  // const double VzS       = -0.5 * V_tot * (hW + hE) / Lx + Vz_corr;
  const double VxW_H     = instance->flux->west;
  const double VxE_H     = instance->flux->east;
  const double VzS       = -(fabs(VxW_H) + fabs(VxE_H))/Lx + Vz_corr;

  // Here we need the total inflow flux: ncell_W * dz * VxW + ncell_E * dz * VxE
  // ncell_W is the number of Vx points for which we apply a BC value on the west 
  // ncell_E is the number of Vx points for which we apply a BC value on the east
  // Then:  VzS = -(ncell_W * dz * VxW + ncell_E * dz * VxE)/Lx + Vz_corr

  // Set boundary nodes types and values
  if (position == W || position == SW || position == NW ) {
    bc.type  = constant_shear_stress;
    bc.value = 0.0;
  } else if ( position == E || position == SE || position == NE) {
    bc.type  = constant_shear_stress;
    bc.value = 0.0;
  } else if (position == S || position == N) {
    bc.type  = constant_velocity;
    bc.value = VzS;
  } else {
    bc.type  = inside;
    bc.value = 0.0;
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddCrazyConductivity(MdoodzInput *input) {
  int               *asthenospherePhases = (int *) malloc(sizeof(int));
  CrazyConductivity *crazyConductivity   = (CrazyConductivity *) malloc(sizeof(CrazyConductivity));
  asthenospherePhases[0]                 = 3;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000;
  input->crazyConductivity               = crazyConductivity;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int main(int nargs, char *args[]) {
// Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "RiverTomLowerCrust.txt"); // Default
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
                  .SetHorizontalVelocity = SetHorizontalVelocity,
                  .SetVerticalVelocity   = SetVerticalVelocity,
                  .SetDualPhase          = SetDualPhase,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx, // manuel
                  .SetBCVz    = SetBCVz,
                  .SetBCT     = SetBCT,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
