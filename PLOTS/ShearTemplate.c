#include "hdf5read.h"
#include "matrices.h"
#include "stdio.h"
#define FILENAME "ShearTemplateReference.gzip.h5"

int main() {
  Hdf5File hdf5File = ReadHdf5("ShearTemplateReference.gzip.h5");
  float *pressureArray = GetPressureArray(hdf5File);
  Reshape reshapePressure = InitReshape(hdf5File.nx, hdf5File.nz);
  float(*pressureMatrix)[hdf5File.nx - 1] = ReshapeFloatArray(pressureArray, &reshapePressure);

  for (int i = 0; i < hdf5File.nx - 1; i++) {
    for (int j = 0; j < hdf5File.nz - 1; j++) {
      printf("%d\t%d\t%f\n", i, j, pressureMatrix[j][i]);
    }
  }
}