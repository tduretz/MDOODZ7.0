#include "hdf5read.h"
#include "assert.h"
#include "stdlib.h"

bool Equals(float a, float b) {
  float epsilon = 0.0000000001;
  return a - b < epsilon;
}

void ShouldReadPressure() {
  Hdf5File hdf5File = ReadHdf5("ShearTemplateReference.gzip.h5");
  float *pressureArray = GetPressureArray(hdf5File);
  assert(Equals(pressureArray[0], 0));
  assert(Equals(pressureArray[1], -0.000087));
  assert(Equals(pressureArray[3493], 0.046443));
  assert(Equals(pressureArray[9998], -0.000087));
  assert(Equals(pressureArray[9999], 0));
  free(pressureArray);
}

int main() {
  ShouldReadPressure();
}