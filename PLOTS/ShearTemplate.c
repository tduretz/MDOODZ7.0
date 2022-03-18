#include "hdf5.h"
#include "stdio.h"
#define FILENAME "ShearTemplateReference.gzip.h5"

static void ReshapeArray(const int x, const float *flat_array) {
  const float(*reshaped)[x] = (void *)flat_array;
  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      printf("%f\t", reshaped[i][j]);
    }
    printf("\n");
  }
}

int main() {
  hid_t File = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t CentersGroup = H5Gopen(File, "Centers", H5P_DEFAULT);
  hid_t PressureDataset = H5Dopen(CentersGroup, "P", H5P_DEFAULT);
  float NumberStepsArray[10000];
  H5Dread(PressureDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          NumberStepsArray);
  ReshapeArray(100, NumberStepsArray);
}