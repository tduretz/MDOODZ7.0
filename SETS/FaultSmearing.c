#include "math.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"

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

    // // Set checkerboard for phase 2
    // Ax = cos( 10.0*2.0*M_PI*coordinate.x / Lx  );
    // Az = cos( 10.0*2.0*M_PI*coordinate.z / Lz  );
    // if ( ( (Az<0.0 && Ax<0.0) || (Az>0.0 && Ax>0.0) ) && dual_phase==2 ) {
    //     dual_phase += input->model.Nb_phases;
    // }

  return dual_phase;
}

char *ReadGeometryFile(char *inputFile, int nb_elems) {
  char *ph_hr = malloc((nb_elems) * sizeof(char));
  FILE *fid   = fopen(inputFile, "rb");
  if (!fid) {
    fprintf(stderr, "\nUnable to open file %s. I will exit here ... \n", inputFile);
    exit(2);
  }
  fread(ph_hr, sizeof(char), nb_elems, fid);
  fclose(fid);
  return ph_hr;
}

void MutateInput(MdoodzInput *input) {
  if (input->model.user0 == 0) {
    printf("Phase map to be drawn in a SetParticles\n");
    return;
  } else if (input->model.user0 == 1) {
    char *fileName = input->model.import_file;
    char  inputFilePath[255];
    sprintf(inputFilePath, "%s/%s", input->model.import_files_dir, input->model.import_file);
    printf("Phase map will be built based on %s\n", fileName);
    const int nx       = (int) (input->model.user7);
    const int nz       = nx;
    const int nb_elems = nx * nz;
    input->geometry    = malloc(sizeof(Geometry));
    input->geometry[0] = (Geometry){
            .nx       = nx,
            .nz       = nz,
            .nb_elems = nb_elems,
            .ph_hr    = ReadGeometryFile(inputFilePath, nb_elems),
    };
  }
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  if (input->geometry) {
    const int nx    = input->geometry->nx+1;
    const int nz    = input->geometry->nz+1;
    // ------------------------- //
    // Locate markers in the image files
    // Find index of minimum/west temperature node
    double    dx_hr = (input->model.xmax - input->model.xmin) / (nx-1);
    double    dz_hr = (input->model.zmax - input->model.zmin) / (nz-1);
    double    dstx  = (coordinates.x - input->model.xmin);
    int       ix    = (int) ceil(dstx / dx_hr) - 1;
    // Find index of minimum/west pressure node
    double    dstz  = (coordinates.z - input->model.zmin);
    int       iz    = (int) ceil(dstz / dz_hr) - 1;
    // Attribute phase
    if (ix < 0) {
      printf("sauceisse!!!\n");
      exit(1);
    }
    if (iz < 0) {
      printf("puréee!!!\n");
      exit(1);
    }
    if (ix + iz * (nx-1) > input->geometry->nb_elems) {
      printf("puréee!!!\n");
      exit(1);
    }
    // Default: phase from input file
    int phase = input->geometry->ph_hr[ix + iz * (nx-1)];
    // Add seed
    const double radius = input->model.user1 / input->scaling.L;
    const double theta  = -0.0 * M_PI / 180.0;
    const double Xn     = coordinates.x * cos(theta) - coordinates.z * sin(theta);
    const double Zn     = coordinates.x * sin(theta) + coordinates.z * cos(theta);
    if ( (pow(Xn / radius / 1.0, 2) + pow(Zn / radius / 1.0, 2) - 1.0) < 0) {
      phase = 1;
    } 
    return phase;
  } 
  else {
    // Simple setup with only one seed
    const double radius = input->model.user1 / input->scaling.L;
    const double theta  = -0.0 * M_PI / 180.0;
    const double Xn     = coordinates.x * cos(theta) - coordinates.z * sin(theta);
    const double Zn     = coordinates.x * sin(theta) + coordinates.z * cos(theta);
    if ( (pow(Xn / radius / 1.0, 2) + pow(Zn / radius / 1.0, 2) - 1.0) < 0) {
      return 1;
    } else {
      return 0;
    }
  }
}

double SetNoise(MdoodzInput *input, Coordinates coordinates, int phase) {
  // srand(69);
  return ((double) rand() / (double) RAND_MAX) - 0.5;
}

double SetPressure(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->model.bkg_pressure;
}

// SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
//   SetBC bc;
//   const double x = coordinates.x;
  
//   // Assign BC values
//   if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
//     bc.value = 0;
//     bc.type  = 13;
//   } else if (position == W || position == E) {
//     bc.value =  -x * (instance->model.bkg_strain_rate - instance->model.bkg_div_rate/3.0);
//     bc.type  = 0;
//   } else {
//     bc.value = 0.0;
//     bc.type  = -1;
//   }
//   return bc;
// }

// SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
//   SetBC bc;
//   const double z = coordinates.z;

//   // Set boundary nodes types and values
//   if (position == W || position == SW || position == NW || position == E || position == SE || position == NE ) {
//     bc.value = 0.0;
//     bc.type  = 13;
//   } else if (position == S || position == N) {
//     bc.value = z * (instance->model.bkg_strain_rate + instance->model.bkg_div_rate/3.0);
//     bc.type  = 0;
//   } else {
//     bc.value = 0;
//     bc.type  = -1;
//   }
//   return bc;
// }

int main(int nargs, char *args[]) {
// Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "FaultSmearing.txt"); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetNoise       = SetNoise,
                  .SetPressure    = SetPressure,
                  .SetDualPhase   = SetDualPhase,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
          .MutateInput = MutateInput,
  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
