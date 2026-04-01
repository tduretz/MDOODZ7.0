#ifndef MDOODZ_TEST_HELPERS_H
#define MDOODZ_TEST_HELPERS_H

#include <hdf5.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>

// Read the Newton iteration count from the HDF5 output
static int getStepsCount(const char *hdf5FileName) {
  hid_t file            = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t iterationsGroup = H5Gopen(file, "Iterations", H5P_DEFAULT);
  hid_t dataset         = H5Dopen(iterationsGroup, "NumberSteps", H5P_DEFAULT);
  int numberSteps[1]    = {-1};
  H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, numberSteps);
  H5Dclose(dataset);
  H5Gclose(iterationsGroup);
  H5Fclose(file);
  printf("Number of iterations: %d\n", numberSteps[0]);
  return numberSteps[0];
}

// Read the final absolute residual from Iterations/rx_abs or rp_abs
static double getFinalResidual(const char *hdf5FileName, const char *datasetName) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, "Iterations", H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, datasetName, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  double *buf = (double *)malloc(dims[0] * sizeof(double));
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  // Find the last non-zero residual
  int nSteps = getStepsCount(hdf5FileName);
  double res = buf[nSteps > 0 ? nSteps - 1 : 0];
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return res;
}

// Read a float field from HDF5 and return the maximum value
static double getMaxFieldValue(const char *hdf5FileName, const char *groupName, const char *datasetName) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, groupName, H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, datasetName, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  float *buf = (float *)malloc(dims[0] * sizeof(float));
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  double maxVal = -DBL_MAX;
  for (hsize_t i = 0; i < dims[0]; i++) {
    if (buf[i] > maxVal) maxVal = buf[i];
  }
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return maxVal;
}

// Read a float field from HDF5 and return the minimum value
static double getMinFieldValue(const char *hdf5FileName, const char *groupName, const char *datasetName) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, groupName, H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, datasetName, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  float *buf = (float *)malloc(dims[0] * sizeof(float));
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  double minVal = DBL_MAX;
  for (hsize_t i = 0; i < dims[0]; i++) {
    if (buf[i] < minVal) minVal = buf[i];
  }
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return minVal;
}

// Read a float field from HDF5 and return the mean value
static double getMeanFieldValue(const char *hdf5FileName, const char *groupName, const char *datasetName) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, groupName, H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, datasetName, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  float *buf = (float *)malloc(dims[0] * sizeof(float));
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  double sum = 0.0;
  for (hsize_t i = 0; i < dims[0]; i++) {
    sum += buf[i];
  }
  double mean = sum / (double)dims[0];
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return mean;
}

// Read a double param from /Model/Params at given index
static double getModelParam(const char *hdf5FileName, int index) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, "Model", H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, "Params", H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  double *buf = (double *)malloc(dims[0] * sizeof(double));
  H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  double val = buf[index];
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return val;
}

// Read a single float value at a specific index from a dataset
static double getFieldValueAt(const char *hdf5FileName, const char *groupName, const char *datasetName, hsize_t index) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, groupName, H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, datasetName, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  float *buf = (float *)malloc(dims[0] * sizeof(float));
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  double val = (index < dims[0]) ? buf[index] : 0.0;
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return val;
}

// Read max of absolute values from a float field
static double getMaxAbsFieldValue(const char *hdf5FileName, const char *groupName, const char *datasetName) {
  hid_t file  = H5Fopen(hdf5FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t group = H5Gopen(file, groupName, H5P_DEFAULT);
  hid_t dset  = H5Dopen(group, datasetName, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  float *buf = (float *)malloc(dims[0] * sizeof(float));
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
  double maxAbsVal = 0.0;
  for (hsize_t i = 0; i < dims[0]; i++) {
    double absVal = fabs((double)buf[i]);
    if (absVal > maxAbsVal) maxAbsVal = absVal;
  }
  free(buf);
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(group);
  H5Fclose(file);
  return maxAbsVal;
}

#endif
