#include "hdf5.h"
#include "header_MDOODZ.h"
#include "assert.h"
#define FILENAME "Output00001.gzip.h5"

int main() {
  const char *args[] = {
      "Some string",
      "ShearTemplate.txt",
  };
  RunMDOODZ(2, (char **)args);

  hid_t File = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t IterationsGroup = H5Gopen(File, "Iterations", H5P_DEFAULT);
  hid_t NumberStepsDataset =
      H5Dopen(IterationsGroup, "NumberSteps", H5P_DEFAULT);
  int NumberStepsArray[1] = {0};
  H5Dread(NumberStepsDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          NumberStepsArray);
  int StepsCount = NumberStepsArray[0];

  assert(StepsCount == 1);
}