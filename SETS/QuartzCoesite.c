#include "math.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"


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
  if (input->model.user1 == 1) {
    printf("Phase map be drawn in a SetParticles\n");
    return;
  } else if (input->model.user1 == 2) {
    char *fileName = input->model.input_file;
    printf("Phase map will be built based on %s\n", fileName);
    const int nx       = 1921;
    const int nz       = 1921;
    const int nb_elems = nx * nz;
    input->geometry    = malloc(sizeof(Geometry));
    input->geometry[0] = (Geometry){
            .nx       = nx,
            .nz       = nz,
            .nb_elems = nb_elems,
            .ph_hr    = ReadGeometryFile(fileName, nb_elems),
    };
  }
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  if (input->geometry) {
    srand(69);
    const int nx    = input->geometry->nx;
    const int nz    = input->geometry->nz;
    // ------------------------- //
    // Locate markers in the image files
    // Find index of minimum/west temperature node
    double    dx_hr = (input->model.xmax - input->model.xmin) / nx;
    double    dz_hr = (input->model.zmax - input->model.zmin) / nz;
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
    if (ix + iz * (input->geometry->nx) > input->geometry->nb_elems) {
      printf("puréee!!!\n");
      exit(1);
    }
    return (int) input->geometry->ph_hr[ix + iz * nx];
  } else {
    const double radius = input->model.user2 / input->scaling.L;
    const double theta  = -30.0 * M_PI / 180.0;
    const double Xn     = coordinates.x * cos(theta) - coordinates.z * sin(theta);
    const double Zn     = coordinates.x * sin(theta) + coordinates.z * cos(theta);
    if (pow(Xn / radius * 2.0, 2) + pow(Zn / radius / 2.0, 2) - 1 < 0) {
      return 1;
    } else {
      return 0;
    }
  }
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase = SetPhase,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
          .MutateInput = MutateInput,
  };
  RunMDOODZ("QuartzCoesite.txt", &setup);
}
