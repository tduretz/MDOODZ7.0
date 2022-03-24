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

  hid_t *filePtr = malloc(sizeof(file));
  memcpy(filePtr, &file, sizeof(file));
  Hdf5File hdf5File = {
          .fileName = filename,
          .nz = nz,
          .nx = nx,
          .filePtr = filePtr};
  return hdf5File;
}

float *GetPressureArray(Hdf5File hdf5File) {
  int arraySize = (hdf5File.nx - 1) * (hdf5File.nz - 1);
  hid_t centersGroup = H5Gopen(*hdf5File.filePtr, "Centers", H5P_DEFAULT);
  hid_t pressureDataset = H5Dopen(centersGroup, "P", H5P_DEFAULT);
  float pressureArray[arraySize];
  H5Dread(pressureDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          pressureArray);
  float *pressureArrayPtr = malloc(sizeof(float) * arraySize);
  memcpy(pressureArrayPtr, pressureArray, sizeof(float) * arraySize);
  return pressureArrayPtr;
}

Coordinates GetCoordinates(Hdf5File hdf5File) {
  hid_t modelGroup = H5Gopen(*hdf5File.filePtr, "Model", H5P_DEFAULT);

  hid_t xcCoord = H5Dopen(modelGroup, "xc_coord", H5P_DEFAULT);
  float xcCoordArray[hdf5File.nx - 1];
  H5Dread(xcCoord, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          xcCoordArray);

  hid_t zcCoord = H5Dopen(modelGroup, "zc_coord", H5P_DEFAULT);
  float zcCoordArray[hdf5File.nz - 1];
  H5Dread(zcCoord, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          zcCoordArray);

  float *xcCoordArrayPtr = malloc(sizeof(float) * hdf5File.nx - 1);
  memcpy(xcCoordArrayPtr, xcCoordArray, sizeof(float) * hdf5File.nx - 1);
  float *zcCoordArrayPtr = malloc(sizeof(float) * hdf5File.nz - 1);
  memcpy(zcCoordArrayPtr, zcCoordArray, sizeof(float) * hdf5File.nz - 1);

  Coordinates coordinates = {.xcCoordArray = xcCoordArrayPtr, .zcCoordArray = zcCoordArrayPtr};
  return coordinates;
}

Nodes GetNodes(Hdf5File hdf5File) {
  int nxVx = hdf5File.nx;
  int nzVx = hdf5File.nz + 1;
  hid_t VxNodesGroup = H5Gopen(*hdf5File.filePtr, "VxNodes", H5P_DEFAULT);
  hid_t VxDataset = H5Dopen(VxNodesGroup, "Vx", H5P_DEFAULT);
  float VxArray[nxVx * nzVx];
  H5Dread(VxDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          VxArray);
  float *VxArrayPtr = malloc(sizeof(float) * nxVx * nzVx);
  memcpy(VxArrayPtr, VxArray, sizeof(float) * nxVx * nzVx);

  int nxVz = hdf5File.nx + 1;
  int nzVz = hdf5File.nz;
  hid_t VzNodesGroup = H5Gopen(*hdf5File.filePtr, "VzNodes", H5P_DEFAULT);
  hid_t VzDataset = H5Dopen(VzNodesGroup, "Vz", H5P_DEFAULT);
  float VzArray[nxVz * nzVz];
  H5Dread(VzDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          VzArray);
  float *VzArrayPtr = malloc(sizeof(float) * nxVz * nzVz);
  memcpy(VzArrayPtr, VzArray, sizeof(float) * nxVz * nzVz);

  Nodes nodes = {.VxArray = VxArrayPtr, .VzArray = VzArrayPtr};
  return nodes;
}