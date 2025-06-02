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
  const double a_out = 4.0e3 / input->scaling.L;
  const double b_out = 2.0e3 / input->scaling.L;
  const double a_inn = 2.0e3 / input->scaling.L;
  const double b_inn = 1.0e3 / input->scaling.L;
  const double xc_e  = 0.0 / input->scaling.L;
  const double zc_e  = -10.0e3 / input->scaling.L;
  const double x = coordinates.x, z =  coordinates.z;
  const double position_out = (x - xc_e) * (x - xc_e) / (a_out * a_out) + (z - zc_e) * (z - zc_e) / (b_out * b_out);
  const double position_inn = (x - xc_e) * (x - xc_e) / (a_inn * a_inn) + (z - zc_e) * (z - zc_e) / (b_inn * b_inn);
  int phase = 0;
  if (position_out <= 1.0) {
    phase = 1;
  }
  if (position_inn <= 1.0) {
    phase = 2;
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

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double surfaceTemperature = (15.0+273.15) / instance->scaling.T;
  const double mohoTemperature    = (600.0 + zeroC) / instance->scaling.T;
  const double gradT              = (mohoTemperature - surfaceTemperature) / fabs(instance->model.zmin);
  const double Tamp               = 50.0 / instance->scaling.T;
  const double x                  = coordinates.x - instance->model.user4 / instance->scaling.L;
  const double z                  = coordinates.z;
  const double basin_width        = 30.0e3 / instance->scaling.L;
  
  double particleTemperature = gradT * (-z) + surfaceTemperature;

    return particleTemperature;
}

double AdjustTemperatureToPhase(MdoodzInput *instance, Coordinates coordinates, double particleTemperature, int phase){
    const double chamberTemperature = (900.0 + zeroC) / instance->scaling.T;
    if (phase > 0)
    {
        return chamberTemperature;
    }
    else
    {
        return particleTemperature;
    }
    
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
            .AdjustTemperatureToPhase = AdjustTemperatureToPhase,
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
