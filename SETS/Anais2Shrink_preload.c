#include "mdoodz.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  int phase = 0;
  // Layers
  const int nb_layers = 50;
  const double dL    = (input->model.user7)/input->scaling.L;   // distance entre 2 droites [m]
  const double space_factor = 2.0;
  const double alpha = (input->model.user9)*M_PI/180.0;   // angle => pente de la droite [°=>rad]
  const double shift = dL/cos(alpha);   // epaisseur en y entre 2 droite [m]
  const double space = (input->model.user8)*shift;       //espace entre 2 bandes [m]
  double yref  = 0.0/input->scaling.L;    // coordonnée y de reference pour la droite inférieure [m]
  const double yref0 = yref-nb_layers/2.0*shift-nb_layers/2.0*space+space/2.0;
  double yp_ls, yp_us;
    
  yref = yref0;
  for (int il=0;il<nb_layers;il++) {
    yp_ls = tan(alpha)*coordinates.x + yref;
    yp_us = tan(alpha)*coordinates.x + yref + shift;
    if (coordinates.z>yp_ls && coordinates.z<yp_us) {
        if (input->model.layering==1) {
            phase = 3;
        }
    }           
    yref += shift + space;
  }
  // Inclusion
  const double radius = (input->model.user3) / input->scaling.L;
  const double b_ell = (input->model.user5)*radius;
  const double a_ell = (input->model.user4)*radius;
  const double x0 = 0.0, z0 = 0.0/input->scaling.L;
  const double angle = -(input->model.user6)*M_PI/180.0;
  const double x_ell = (coordinates.x - x0)*cos(angle) + (coordinates.z - z0)*sin(angle);
  const double z_ell =-(coordinates.x - x0)*sin(angle) + (coordinates.z - z0)*cos(angle);
  if (pow(x_ell/b_ell,2.0) + pow(z_ell/a_ell,2.0) < 1.0) {
    phase = 1;
  } 
  if (phase >= input->model.Nb_phases) {printf("Non existing phase number, update .txt file\n"); exit(11);} 
  else { return phase; }
}

int SetDualPhase(MdoodzInput *input, Coordinates coordinate, int phase) {
  int    dual_phase = phase;
  double Lx = input->model.xmax - input->model.xmin;
  double Lz = input->model.zmax - input->model.zmin;
  double Ax, Az;
  // Set checkerboard for phase 0
  Ax = cos( 10.0*2.0*M_PI*coordinate.x / Lx  );
  Az = cos( 10.0*2.0*M_PI*coordinate.z / Lz  );
  if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==0 ) {
      dual_phase += input->model.Nb_phases;
  }
  // Set checkerboard for phase 2
  Ax = sin( 10.0*2.0*M_PI*coordinate.x / Lx  );
  Az = sin( 10.0*2.0*M_PI*coordinate.z / Lz  );
  if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==2 ) {
      dual_phase += input->model.Nb_phases;
  }
  return dual_phase;
}

// double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
//   const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
//   const double P_init = (input->model.bkg_pressure         ) / input->scaling.S;
//   if (1 == 0) {
//     return input->materials.rho[phase] * exp(input->materials.bet[phase]*P_init -  input->materials.alp[phase] * T_init);
//   } else {
//     return input->materials.rho[phase];
//   }
// }

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  return (input->model.user0 + zeroC) / input->scaling.T;
}

double SetPressure(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.bkg_pressure;
}

double SetTxx(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.preload_sxxd;
}

double SetTzz(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.preload_szzd;
}

double SetTxz(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.preload_sxz;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if (nargs < 2) {
    asprintf(&input_file, DefaultTextFilename(__FILE__)); // Default
  } else {
    asprintf(&input_file, "%s", args[1]);// Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase              = SetPhase,
                   .SetDualPhase          = SetDualPhase,
                   .SetTemperature        = SetTemperature,
                   .SetPressure           = SetPressure,
                   .SetTxx                = SetTxx,
                   .SetTzz                = SetTzz,
                   .SetTxz                = SetTxz,
                  //  .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
