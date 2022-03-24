#include "hdf5read.h"
#include "matrices.h"
#include "stdio.h"
#include "stdlib.h"
#define FILENAME "ShearTemplateReference.gzip.h5"

int main() {
  Hdf5File hdf5File = ReadHdf5("ShearTemplateReference.gzip.h5");
  float *pressureArray = GetPressureArray(hdf5File);
  Reshape reshapePressure = InitReshape(hdf5File.nx, hdf5File.nz);
  float(*pressureMatrix)[hdf5File.nx - 1] = ReshapeFloatArray(pressureArray, &reshapePressure);
  Coordinates coordinates = GetCoordinates(hdf5File);

  Nodes nodes = GetNodes(hdf5File);
  Reshape reshapeVx = InitReshape(hdf5File.nx, hdf5File.nz + 1);
  reshapeVx.xMaxLimit = -1;
  reshapeVx.yMinLimit = 1;
  reshapeVx.yMaxLimit = -1;
  float(*VxMatrix)[hdf5File.nx - 1] = ReshapeFloatArray(nodes.VxArray, &reshapeVx);
  reshapeVx.xMinLimit = 1;
  reshapeVx.yMinLimit = 1;
  reshapeVx.yMaxLimit = -1;
  float(*VxMatrix2)[hdf5File.nx - 1] = ReshapeFloatArray(nodes.VxArray, &reshapeVx);
  Sum sumVx = InitSum(hdf5File.nx - 1, hdf5File.nz - 1);
  sumVx.multiplier = 0.5;
  float(*VxSumMatrix)[hdf5File.nx - 1] = SumFloatMatrices(VxMatrix, VxMatrix2, &sumVx);

  Reshape reshapeVz = InitReshape(hdf5File.nx + 1, hdf5File.nz);
  reshapeVz.xMinLimit = -1;
  reshapeVz.xMaxLimit = -1;
  reshapeVz.yMaxLimit = -1;
  float(*VzMatrix)[hdf5File.nx] = ReshapeFloatArray(nodes.VxArray, &reshapeVz);
  reshapeVz.xMinLimit = 1;
  reshapeVz.xMaxLimit = -1;
  reshapeVz.yMinLimit = 1;
  float(*VzMatrix2)[hdf5File.nx] = ReshapeFloatArray(nodes.VxArray, &reshapeVz);
  float(*VzSumMatrix)[hdf5File.nx - 1] = SumFloatMatrices(VzMatrix, VzMatrix2, &sumVx);

  for (int i = 0; i < hdf5File.nx - 1; i++) {
    for (int j = 0; j < hdf5File.nz - 1; j++) {
      //printf("%f\t%f\t%f\t%f\t%f\n", coordinates.xcCoordArray[i], coordinates.zcCoordArray[j], pressureMatrix[j][i], VxSumMatrix[j][i], VzSumMatrix[j][i]);
      printf("%f\t", VxMatrix[j][i]);
    }
    printf("\n");
  }

  free(pressureArray);
  free(pressureMatrix);
  free(hdf5File.filePtr);
}