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
  const double a_out   =   4.0e3 / input->scaling.L;
  const double b_out   =   2.0e3 / input->scaling.L;
  const double a_inn   =   2.0e3 / input->scaling.L;
  const double b_inn   =   1.0e3 / input->scaling.L;
  const double xc_e    =   0.0e3 / input->scaling.L;  // X-coord of ellipse center
  const double zc_e    = -10.0e3 / input->scaling.L;  // Z-coord of ellipse center
  const double d_f     =   1.0e3 / input->scaling.L;  // Fault zone width (true)
  const double xc_fl   =   0.0e3 / input->scaling.L;  // X-coord of left fault boundary
  const double zc_fl   =  -5.0e3 / input->scaling.L;  // Z-coord of left fault boundary
  const double angle_f = 60.0 / 180.0 * M_PI;           // Fault zone angle
  const double xc_fr   = cos(angle_f) * d_f + xc_fl;  // X-coord of right fault boundary
  const double zc_fr   = sin(angle_f) * d_f + zc_fl;  // Z-coord of right fault boundary
  const double x = coordinates.x, z =  coordinates.z; // Coords of current marker

  // Determine position of marker w.r.t. magma chamber (ellipse)
  const double position_out = (x - xc_e) * (x - xc_e) / (a_out * a_out) + (z - zc_e) * (z - zc_e) / (b_out * b_out);
  const double position_inn = (x - xc_e) * (x - xc_e) / (a_inn * a_inn) + (z - zc_e) * (z - zc_e) / (b_inn * b_inn);
  
  // Determine position of marker w.r.t. fault zone
  const double s_l = -(sin(angle_f) * (x - xc_fl) + cos(angle_f) * (z - zc_fl)); // Signed value: if s >= 0 -> point is above fault boundary, s <= 0 point is below fault boundary
  const double s_r = -(sin(angle_f) * (x - xc_fr) + cos(angle_f) * (z - zc_fr));
  
  // Initialize phase for current marker
  int phase = 0;

  // Set faulted zone
  if (s_l <= 0.0 && s_r >= 0.0 && z >= zc_e)
  {
    phase = 0;
  }
  
  // Set chamber and mush
  // if (position_out <= 1.0) {
  //   phase = 1;
  // }
  // if (position_inn <= 1.0) {
  //   phase = 2;
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

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  const double P_init = (input->model.bkg_pressure         ) / input->scaling.S;
  if (0 == 0) {
    return input->materials.rho[phase] * exp(input->materials.bet[phase]*P_init -  input->materials.alp[phase] * T_init);
  } else {
    return input->materials.rho[phase];
  }
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double surfaceTemperature = (15.0+273.15) / instance->scaling.T;
  const double mohoTemperature    = (600.0 + zeroC) / instance->scaling.T;
  const double gradT              = (mohoTemperature - surfaceTemperature) / fabs(instance->model.zmin);
  const double Tamp               = (900.0 + zeroC) / instance->scaling.T;

  // Ellipse
  const double a_out = 4.0e3 / instance->scaling.L;
  const double b_out = 2.0e3 / instance->scaling.L;
  const double a_inn = 2.0e3 / instance->scaling.L;
  const double b_inn = 1.0e3 / instance->scaling.L;
  const double xc_e  = 0.0 / instance->scaling.L;
  const double zc_e  = -10.0e3 / instance->scaling.L;
  const double x = coordinates.x, z =  coordinates.z;
  const double position_out = (x - xc_e) * (x - xc_e) / (a_out * a_out) + (z - zc_e) * (z - zc_e) / (b_out * b_out);
  const double position_inn = (x - xc_e) * (x - xc_e) / (a_inn * a_inn) + (z - zc_e) * (z - zc_e) / (b_inn * b_inn);
  
  double particleTemperature = gradT * (-z) + surfaceTemperature;

  // if (position_out <= 1.0) {
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
  double moho_temperature    = (600.0 + zeroC) / instance ->scaling.T;
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
    bc.value = 0.;
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
            .SetDensity               = SetDensity,
            .SetTemperature           = SetTemperature,
            .SetHorizontalVelocity    = SetHorizontalVelocity,
            .SetVerticalVelocity      = SetVerticalVelocity
        },
        .SetBCs = &(SetBCs_ff){
            .SetBCVx = SetPureShearBCVx,
            .SetBCVz = SetPureShearBCVz,
            .SetBCT  = SetBCT,
        },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);

}
