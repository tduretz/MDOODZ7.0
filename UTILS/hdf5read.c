#include "hdf5read.h"

#include "hdf5.h"
#include "stdlib.h"
#include "string.h"

Hdf5File ReadHdf5(char *filename) {
  hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t modelGroup = H5Gopen(file, "Model", H5P_DEFAULT);
  hid_t paramsDataset = H5Dopen(modelGroup, "Params", H5P_DEFAULT);
  double params[8];
  H5Dread(paramsDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          params);
  int nx = (int) params[3];
  int nz = (int) params[4];
  hid_t *file_p = &file;
  Hdf5File hdf5File = {
          .fileName = filename,
          .nz = nz,
          .nx = nx,
          .file_p = file_p};
  return hdf5File;
}

float *GetPressureArray(Hdf5File hdf5File) {
  int arraySize = (hdf5File.nx - 1) * (hdf5File.nz - 1);
  hid_t centersGroup = H5Gopen(*hdf5File.file_p, "Centers", H5P_DEFAULT);
  hid_t pressureDataset = H5Dopen(centersGroup, "P", H5P_DEFAULT);
  float pressureArray[arraySize];
  H5Dread(pressureDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          pressureArray);
  printf("%f\n", pressureArray[10]);
  float(*pressureArrayPtr) = malloc(sizeof(float) * arraySize);
  memcpy(pressureArrayPtr, pressureArray, sizeof(float) * arraySize);
  return pressureArrayPtr;
}