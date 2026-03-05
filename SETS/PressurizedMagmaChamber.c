#include "mdoodz.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel   = -0.0e3 / instance->scaling.L;
  const double basin_width = 30.0e3 / instance->scaling.L;
  const double h_pert      = instance->model.user3 / instance->scaling.L;
  const double x           = x_coord - instance->model.user4 / instance->scaling.L;
  return TopoLevel;
//   return TopoLevel + h_pert * exp( - (x*x) / (2.0*basin_width*basin_width)   );
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  // Intrusion
  const double a_intr  = input->model.therm_perturb_rad_x;
  const double b_intr  = input->model.therm_perturb_rad_z;
  const double xc_intr = input->model.therm_perturb_x0;
  const double zc_intr = input->model.therm_perturb_z0;
  const double x = coordinates.x, z =  coordinates.z;
  const double position_intr = (x - xc_intr) * (x - xc_intr) / (a_intr * a_intr) + (z - zc_intr) * (z - zc_intr) / (b_intr * b_intr);
  
  // Initialize phase for current marker
  int phase = 0;

  // Set second phase if needed
  // if (position_intr <= 1.0)
  // {
  //   phase = 1;
  // }
  
  return phase;
}

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
    
    int    dual_phase = phase;
    double Lx = input->model.xmax - input->model.xmin;
    double Lz = input->model.zmax - input->model.zmin;
    double Ax, Az;

    // Set checkerboard for phase 0
    Ax = cos( 6.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( 6.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
        dual_phase += input->model.Nb_phases;
    }

    // Set checkerboard for phase 1
    Ax = cos( 24.0*2.0*M_PI*coordinate.x / Lx  );
    Az = sin( 24.0*2.0*M_PI*coordinate.z / Lz  );
    if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==1 ) {
        dual_phase += input->model.Nb_phases;
    }

  return dual_phase;
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  const double surfaceTemperature = ( 20.0 + zeroC) / input->scaling.T;
  const double mohoTemperature    = (450.0 + zeroC) / input->scaling.T;
  const double gradT              = (mohoTemperature - surfaceTemperature) / fabs(input->model.zmin);
  const double Tamp               = (750.0 + zeroC) / input->scaling.T;

  // Intrusion
  const double a_intr  = input->model.therm_perturb_rad_x;
  const double b_intr  = input->model.therm_perturb_rad_z;
  const double xc_intr = input->model.therm_perturb_x0;
  const double zc_intr = input->model.therm_perturb_z0;
  const double x = coordinates.x, z =  coordinates.z;
  const double position_intr = (x - xc_intr) * (x - xc_intr) / (a_intr * a_intr) + (z - zc_intr) * (z - zc_intr) / (b_intr * b_intr);

  double particleTemperature = gradT * (-z) + surfaceTemperature;

  // if (position_intr <= 1.0) {
  //   particleTemperature = Tamp;
  // }
  
  return particleTemperature;
}

double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.bkg_strain_rate;
}

double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.bkg_strain_rate;
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates,  double particleTemperature) {
  SetBC     bc;
  double surface_temperature =          zeroC  / instance->scaling.T;
  double moho_temperature    = (450.0 + zeroC) / instance ->scaling.T;
  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = moho_temperature;
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

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lx  = input->model.xmax - input->model.xmin;
  const double VxW = -input->model.user2 / input->scaling.L * input->scaling.t;
  const double VxE =  input->model.user2 / input->scaling.L * input->scaling.t;
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0.0;
    bc.type  = constant_shear_stress;
  } else if (position == W) {
    bc.value = VxW;
    bc.type  = constant_velocity;
  } else if (position == E)
  {
    bc.value = VxE;
    bc.type  = constant_velocity;
  } else {
    bc.value = 0.0;
    bc.type  = inside;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Lz = input->model.zmax - input->model.zmin;
  const double VzS = -input->model.user3 / input->scaling.L * input->scaling.t;
  const double VzN =  input->model.user3 / input->scaling.L * input->scaling.t;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = constant_shear_stress;
  } else if (position == S) {
    bc.value = VzS;
    bc.type  = constant_velocity;
  } else if (position == N)
  {
    bc.value = VzN;
    bc.type  = constant_velocity;
  } else {
    bc.value = 0;
    bc.type  = inside;
  }
  return bc;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "PressurizedMagmaChamber.txt"); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
        .BuildInitialTopography = &(BuildInitialTopography_ff){
            .SetSurfaceZCoord = SetSurfaceZCoord,
        },
        .SetParticles  = &(SetParticles_ff){
            .SetPhase                 = SetPhase,
            .SetTemperature           = SetTemperature,
        },
        .SetBCs = &(SetBCs_ff){
            .SetBCVx = SetPureOrSimpleShearBCVx,
            .SetBCVz = SetPureOrSimpleShearBCVz,
            // .SetBCVx = SetBCVx,
            // .SetBCVz = SetBCVz,
            .SetBCT  = SetBCT,
        },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);

}
