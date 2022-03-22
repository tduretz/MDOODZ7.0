#include "hdf5.h"
#include "stdio.h"
#define FILENAME "ShearTemplateReference.gzip.h5"

int main() {
  hid_t File = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

  hid_t ModelGroup = H5Gopen(File, "Model", H5P_DEFAULT);
  hid_t paramsDataset = H5Dopen(ModelGroup, "Params", H5P_DEFAULT);
  double params[8];
  H5Dread(paramsDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          params);
  int nx = (int)params[3];
  int nz = (int)params[4];

  hid_t CentersGroup = H5Gopen(File, "Centers", H5P_DEFAULT);
  hid_t PressureDataset = H5Dopen(CentersGroup, "P", H5P_DEFAULT);
  float PressureArray[(nx - 1) * (nz - 1)];
  H5Dread(PressureDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          PressureArray);

  hid_t VxNodesGroup = H5Gopen(File, "VxNodes", H5P_DEFAULT);
  hid_t VxDataset = H5Dopen(VxNodesGroup, "Vx", H5P_DEFAULT);
  float xComponentArray[nx * (nz + 1)];
  H5Dread(VxDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          xComponentArray);

  hid_t VzNodesGroup = H5Gopen(File, "VzNodes", H5P_DEFAULT);
  hid_t VzDataset = H5Dopen(VzNodesGroup, "Vz", H5P_DEFAULT);
  float zComponentArray[(nx + 1) * nz];
  H5Dread(VzDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          zComponentArray);

  hid_t xcCoord = H5Dopen(ModelGroup, "xc_coord", H5P_DEFAULT);
  float xcCoordArray[nx - 1];
  H5Dread(xcCoord, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          xcCoordArray);

  hid_t zcCoord = H5Dopen(ModelGroup, "zc_coord", H5P_DEFAULT);
  float zcCoordArray[nz - 1];
  H5Dread(zcCoord, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          zcCoordArray);

  int flatIndex = 0;
  for (int i = 0; i < nx - 1; i++) {
    for (int j = 0; j < nz - 1; j++) {
      // printf("%f\t%f\t%f\t\n", xcCoordArray[j], zcCoordArray[i],
      // PressureArray[flatIndex]);
      flatIndex++;
    }
  }

  float arr1[(nx - 1) * (nz - 1)];
  flatIndex = 0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nz; j++) {
      if (i < nx && j > 0 && j < nx) {
        arr1[flatIndex] = xComponentArray[j + i * j];
        flatIndex++;
      }
    }
  }

  float arr2[(nx - 1) * (nz - 1)];
  flatIndex = 0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < nz; j++) {
      if (i > 0 && j > 0 && j < nx) {
        arr2[flatIndex] = xComponentArray[j + i * j];
        flatIndex++;
      }
    }
  }

  float averagedXcCoordArray[(nx - 1) * (nz - 1)];
  flatIndex = 0;
  for (int i = 0; i < nx - 1; i++) {
    for (int j = 0; j < nz - 1; j++) {
      averagedXcCoordArray[flatIndex] = (arr1[flatIndex] + arr2[flatIndex]) / 2;
      flatIndex++;
    }
  }

  float averagedZcCoordArray[(nx - 1) * (nz - 1)];

  flatIndex = 0;
  for (int i = 0; i < nx - 1; i++) {
    for (int j = 0; j < nz - 1; j++) {
      printf("%f\t%f\t%f\t\n", xcCoordArray[j], zcCoordArray[i],
             PressureArray[flatIndex]);
      flatIndex++;
    }
  }
}