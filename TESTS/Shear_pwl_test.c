#include "hdf5.h"
#include "header_MDOODZ.h"
#include "assert.h"
#define FILENAME "Output00001.gzip.h5"

int main() {
  const char *args[] = {
      "Some string",
      "Setup03_strongCircleAnisoConst.txt",
  };
  RunMDOODZ(2, (char **)args);

  hid_t File = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

  assert(!!File);
}