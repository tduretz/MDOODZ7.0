#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel = 10 / instance->scaling.L;
  const double h_pert    = instance->model.user3 / instance->scaling.L;
  return TopoLevel + 0.0*x_coord;// + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double radius = input->model.user1 / input->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double T_init = (input->model.user0 + zeroC) / input->scaling.T;
  if (1 == 0) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (T_init - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}


SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  double       x, z;
  // ----------- Pure shear ----------- //
  if (instance->model.shear_style==0) {
    if (position == W) {
      x = coord.x;
      z = coord.z;
      bc.type  = 2;
      bc.value = -1.;//instance->stress->west*0.0;//1*z*4;
    } else if (position == E) {
      x = coord.x;
      z = coord.z;
      bc.type  = 2;
      bc.value = -1.;//instance->stress->east*0.0;//1*z*4;
    } else if (position == S || position == SE || position == SW) {
        x = coord.x;
        z = coord.z + instance->model.dx / 2.0;// make sure it lives on the boundary
        bc.type  = 13;
        bc.value = 0*3e-3;
        // Pick a point in the middle
        if ( (fabs(x)-0.0) < instance->model.dx/2) {
          printf("COUCOU x S\n");
          bc.type  = 11;
          bc.value = 0*3e-3;
        }
    } else if (position == N || position == NE || position == NW) {
      x = coord.x;
      z = coord.z - instance->model.dx / 2.0;// make sure it lives on the boundary
      bc.type  = 13;
      bc.value = 0.0;
        // Pick a point in the middle
        if ( (fabs(x)-0.0) < instance->model.dx/2) {
          printf("COUCOU x N\n");
          bc.type  = 11;
          bc.value = 0*3e-3;
        }
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  }
  // ----------- Simple shear ----------- //
  else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == S || position == SE || position == SW) {
      bc.type  = 11;
      bc.value = -instance->model.bkg_strain_rate * Lz;
    } else if (position == N || position == NE || position == NW) {
      bc.type  = 11;
      bc.value =  instance->model.bkg_strain_rate * Lz;
    } else if (position == E) {
      // bc.type  = 0;
      // bc.value = 0.0;
      bc.type  = -12;
      bc.value = 0.0;
    } else if (position == W) {
      // bc.type  = 0;
      // bc.value = 0.0;
      bc.type  = -2;
      bc.value = 0.0;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  double       x, z;
  // ----------- Pure shear ----------- //
  if (instance->model.shear_style==0) {
    if (position == N) {
      x = coord.x;
      z = coord.z;
      bc.type  = 2;
      bc.value = 1.;
    }
    else if (position == S) {
      x = coord.x;
      z = coord.z;
      bc.type  = 2;
      bc.value = 1.;
    }
    else if (position == W || position == SE || position == SW) {
      x = coord.x + instance->model.dx / 2.0;
      z = coord.z;
      bc.type  = 13;
      bc.value = 0.0;
      // Pick a point in the middle
      if ( (fabs(z)-0.0) < instance->model.dz/2) {
        printf("COUCOU z W\n");
        bc.type  = 11;
        bc.value = 0;
      }
    }
    else if (position == E || position == NE || position == NW) {
      x = coord.x - instance->model.dx / 2.0;
      z = coord.z;
      bc.type  = 13;
      bc.value = 0.0;
      // Pick a point in the middle
      if ( (fabs(z)-0.0) < instance->model.dz/2) {
        printf("COUCOU z E\n");
        bc.type  = 11;
        bc.value = 0;
      }
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  }
  // ----------- Simple shear ----------- //
  else {
    if ( position == W || position == NW || position == SW) {
      bc.type  = -12;
      bc.value = 0.0;
      // bc.type  = 13;
      // bc.value = 0.0;
      // // Pick a point in the middle
      // if ( (fabs(coord.z)-0.0) < instance->model.dz/2) {
      //   printf("COUCOU z W\n");
      //   bc.type  = 11;
      //   bc.value = 0;
      // }
    }
    else if (position == E || position == NE || position == SE) {
      bc.type  = -12;
      bc.value = 0.0;
      // bc.type  = 13;
      // bc.value = 0.0;
    } else if (position == S) {
      bc.type  = 2;
      bc.value = -1.0;
      //  if ( (fabs(coord.x)-0.0) < instance->model.dx/2) {
      //   bc.type  = 2;
      //   bc.value = 1.0;
      // }
    }
      else if (position == N) {
      bc.type  = 0;
      bc.value = 0.0;
      // if ( (fabs(coord.x)-0.0) < instance->model.dx/2) {
      //   bc.type  = 2;
      //   bc.value = 1.0;
      // }

    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  }
  return bc;
}

int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;

  if ( nargs < 2 ) {
    asprintf(&input_file, DefaultTextFilename(__FILE__)); // Default
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
                   .SetPhase              = SetPhase,
                   .SetDensity            = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
