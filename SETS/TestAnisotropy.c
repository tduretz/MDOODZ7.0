#include "complex.h"
#include "mdoodz.h"
#include "math.h"
#include "stdio.h"
#include "time.h"
#include "stdlib.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double rad = instance->model.user1 / instance->scaling.L;
  const double x = coordinates.x;
  const double z = coordinates.z;
  const bool inclusion = instance->model.user0;
  if (inclusion) {
    if ( instance->model.shear_style==0 ) {
      if ( pow(x,2) + pow(z,2) < pow(rad,2) ) {
        return 1;
      } else {
        return 0;
      }
    }
    else {
      if ( (pow(x-instance->model.xmin,2) + pow(z,2) < pow(rad,2)) || (pow(x-instance->model.xmax,2) + pow(z,2) < pow(rad,2)) ) {
        return 1;
      } else {
        return 0;
      }
    }
  }
  else return 0;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  return instance->materials.rho[phase];
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm     = 1.0;
  const double mc     = 1e3;
  double       Vx, Vz, P, eta, sxx, szz, x, z;
  if (instance->model.shear_style == 0) {
    printf("SETTING PURE SHEAR BC\n");
    if (position == W || position == E) {
      x = coord.x;
      z = coord.z;
      bc.type  = 0;
      bc.value = coord.x*instance->model.bkg_strain_rate;
    } else if (position == S || position == SE || position == SW) {
      x = coord.x;
      z = coord.z + instance->model.dx / 2.0;// make sure it lives on the boundary
      bc.type  = 13;
      bc.value = 0.0;
    } else if (position == N || position == NE || position == NW) {
      x = coord.x;
      z = coord.z - instance->model.dx / 2.0;// make sure it lives on the boundary
      bc.type  = 13;
      bc.value = 0.0;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  } else {
    printf("SETTING SIMPLE SHEAR BC\n");
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == S || position == SE || position == SW) {
      x = coord.x;
      z = coord.z + instance->model.dx / 2.0;// make sure it lives on the boundary
      bc.type  = 13;
      bc.value = 0.0;
    } else if (position == N || position == NE || position == NW) {
      x = coord.x;
      z = coord.z - instance->model.dx / 2.0;// make sure it lives on the boundary
      bc.type  = 13;
      bc.value = 0.0;
    } else if (position == E) {
      bc.value = 0.0;
      bc.type  = -12;
    } else if (position == W) {
      bc.value = 0.0;
      bc.type  = -2;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm     = 1.0;
  const double mc     = 1e3;
  double       Vx, Vz, P, eta, sxx, szz, x, z;
  if (instance->model.shear_style == 0) {
    if (position == N || position == NE || position == NW || position == S || position == SE || position == SW) {
      x = coord.x;
      z = coord.z;
      bc.type  = 0;
      bc.value = -coord.z*instance->model.bkg_strain_rate;
    }
    else if (position == W) {
      x = coord.x + instance->model.dx / 2.0;
      z = coord.z;
      bc.type  = 11;
      bc.value = 0.0;
    }
    else if (position == E) {
      x = coord.x + instance->model.dx / 2.0;
      z = coord.z;
      bc.type  = 13;
      bc.value = 0.0;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  } 
  else {
    if (position == E || position == W || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type = -12;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type = 0;
    } else {
      bc.value = 0.0;
      bc.type = -1;
    }
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int main(int nargs, char *args[]) {
  srand(time(NULL));
  char *input_file;
  if (nargs < 2) {
    asprintf(&input_file, "TestAnisotropy.txt");// Default
  } else {
    asprintf(&input_file, "%s", args[1]);// Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup instance = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase   = SetPhase,
                  .SetDensity = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ(input_file, &instance);
  free(input_file);
}