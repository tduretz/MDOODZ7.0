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
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double a_out = 4.0e3 / input->scaling.L;
  const double b_out = 2.0e3 / input->scaling.L;
  const double xc_e  = 0.0 / input->scaling.L;
  const double zc_e  = -10.0e3 / input->scaling.L;
  const double x = coordinates.x, z =  coordinates.z;
  const double position_out = (x - xc_e) * (x - xc_e) / (a_out * a_out) + (z - zc_e) * (z - zc_e) / (b_out * b_out);
  int phase = 0;

  if (position_out <= 1.0) {
    phase = 1;
  }

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

double SetTemperature(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const double surfaceTemperature = (15.0+273.15) / instance->scaling.T;
  const double mantleTemperature  = (1000.0) / instance->scaling.T;
  const double gradT              = (mantleTemperature - surfaceTemperature) / fabs(instance->model.zmin);
  const double Tamp               = (900.0 + zeroC) / instance->scaling.T;
  const double d                  =  coordinates.z;

  // Initial temperature: set the non linear function here
  //double particleTemperature      = (500.0 + zeroC) / instance->scaling.T;

  // Linear gradient example:
  //double particleTemperature = gradT * (-d) + surfaceTemperature;  

  // O'reilly and Davies Heat pipe geotherm
  // Load and set values
  const double T_s = surfaceTemperature;
  const double T_D = mantleTemperature;
  const double rho = instance->materials.rho[1];
  const double c_p = instance->materials.Cp[1];
  const double v_z = (instance->model.user4)/(instance->scaling.V);
  const double H = instance->materials.Qr[1];
  const double t_cond = instance->materials.k[1];

  // Spatial vals
  const double z = (-d);
  const double D = fabs(instance->model.zmin);

  // Dependent equations
  const double kappa = ((t_cond)/(rho*c_p)); // thermal diffusivity (m^2/s)
  // u = downward advection velocity / thermal diffusivity (m-1)
  const double u = (v_z/kappa);
  // Q (volumetric heat flux) = Volumetric heat production rate/ (density * heat capacity * Thermal diffusivity) (K m2)
  const double Q = ((H)/(rho*c_p*kappa));

  double particleTemperature = (T_s + Q * z/u + ((T_D-T_s-Q*D/u)/(1-exp(-u*D)))*(exp(u*(z-D))-exp(-u*D)));

  if (particleTemperature>1000)
  {
    particleTemperature == 1000/(instance->scaling.T);
  }
  
  //printf("checkvals: %e %e %e %e %e %e %e %e %e %e\n", particleTemperature, u, Q, v_z, rho, c_p, H, z, D, t_cond);

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
  const double surface_temperature =          zeroC  / instance->scaling.T;
  const double basalTemperature  = (1000.0 + zeroC) / instance->scaling.T;
  
  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = basalTemperature;
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
    asprintf(&input_file, "HeatPipeRelaxation.txt"); // Default
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