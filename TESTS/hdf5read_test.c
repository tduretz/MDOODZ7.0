#include "hdf5read.h"

#include "assert.h"
#include "stdio.h"
#include "math.h"
#include "stdlib.h"

bool EqualsF(float a, float b) {
  float epsilon = 0.000001;
  return fabs(a - b) < epsilon;
}

void ShouldReadPressure() {
  Hdf5File hdf5File = ReadHdf5("ShearTemplateReference.gzip.h5");
  float *pressureArray = GetPressureArray(hdf5File);
  assert(EqualsF(pressureArray[0], 0));
  assert(EqualsF(pressureArray[1], -0.000087));
  assert(EqualsF(pressureArray[3493], 0.046443));
  assert(EqualsF(pressureArray[9998], -0.000087));
  assert(EqualsF(pressureArray[9999], 0));
  free(pressureArray);
}

void ShouldReadCoordinates() {
  Hdf5File hdf5File = ReadHdf5("ShearTemplateReference.gzip.h5");
  Coordinates coordinates = GetCoordinates(hdf5File);
  assert(EqualsF(coordinates.zcCoordArray[0], -0.495000));
  assert(EqualsF(coordinates.zcCoordArray[99], 0.495000));
  assert(EqualsF(coordinates.zcCoordArray[50], 0.005000));
  assert(EqualsF(coordinates.xcCoordArray[0], -0.495000));
  assert(EqualsF(coordinates.xcCoordArray[99], 0.495000));
  assert(EqualsF(coordinates.xcCoordArray[50], 0.005000));
  free(coordinates.zcCoordArray);
  free(coordinates.xcCoordArray);
}

void ShouldReadNodes() {
  Hdf5File hdf5File = ReadHdf5("ShearTemplateReference.gzip.h5");
  Nodes nodes = GetNodes(hdf5File);
  printf("%f\n", nodes.VzArray[0]);
  printf("%f\n", nodes.VzArray[10200]);
  printf("%f\n", nodes.VxArray[312]);
  printf("%f\n", nodes.VzArray[8099]);
  printf("%f\n", nodes.VxArray[1]);
  assert(EqualsF(nodes.VzArray[0], -0.500000));
  assert(EqualsF(nodes.VzArray[10200], 0.500000));
  assert(EqualsF(nodes.VxArray[312], 0.412373));
  assert(EqualsF(nodes.VzArray[8099], 0.283529));
  assert(EqualsF(nodes.VxArray[1], 0.490271));
  assert(EqualsF(nodes.VxArray[0], 0.500000));
}

int main() {
  ShouldReadCoordinates();
  ShouldReadPressure();
  ShouldReadNodes();
}