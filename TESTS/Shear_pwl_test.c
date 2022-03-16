#include "hdf5.h"
#include "header_MDOODZ.h"
#include "assert.h"
#define FILENAME "Output00001.gzip.h5"

int main() {
  printf("Test: Shear_pwl model with OPT=ON and OMP=ON should finish and produce output\n");
  //printf("Number of threads used: %d\n", omp_get_num_threads());

  const char *args[] = {
      "Some string",
      "Setup03_strongCircleAnisoConst.txt",
  };
  RunMDOODZ(2, (char **)args);

  hid_t File = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

  assert(!!File);
}