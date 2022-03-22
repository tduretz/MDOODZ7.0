#include "stdlib.h"
#include "matrices.h"
#include "stdbool.h"

bool IsCutOff(int i, int j, Reshape *reshape) {
  if (j + 1 <= reshape->xMinLimit) {
    return true;
  }
  if (j >= reshape->cols + reshape->xMaxLimit) {
    return true;
  }
  if (i + 1 <= reshape->yMinLimit) {
    return true;
  }
  if (i >= reshape->cols + reshape->yMaxLimit) {
    return true;
  }
  return false;
}

int *ReshapeIntArray(int *array, Reshape *reshape) {
  int(*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0 ; i < reshape->rows ; i++) {
    for (int j = 0 ; j < reshape->cols ; j++) {
      int flatIndex = j + i * reshape->cols;
      if (IsCutOff(i, j, reshape)) {
        continue;
      }
      arr[updatedIndex] = array[flatIndex];
      updatedIndex++;
    }
  }
  return arr;
}

float *ReshapeFloatArray(float *array, Reshape *reshape) {
  float(*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0 ; i < reshape->rows ; i++) {
    for (int j = 0 ; j < reshape->cols ; j++) {
      int flatIndex = j + i * reshape->cols;
      if (IsCutOff(i, j, reshape)) {
        continue;
      }
      arr[updatedIndex] = array[flatIndex];
      updatedIndex++;
    }
  }
  return arr;
}

double *ReshapeDoubleArray(double *array, Reshape *reshape) {
  double (*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0 ; i < reshape->rows ; i++) {
    for (int j = 0 ; j < reshape->cols ; j++) {
      int flatIndex = j + i * reshape->cols;
      if (IsCutOff(i, j, reshape)) {
        continue;
      }
      arr[updatedIndex] = array[flatIndex];
      updatedIndex++;
    }
  }
  return arr;
}

int *SumIntMatrices(int *matrix, int *matrix2, Sum *sum) {
  int(*arr) = malloc(sizeof *arr);
  for (int i = 0 ; i < sum->rows ; i++) {
    for (int j = 0 ; j < sum->cols ; j++) {
      int flatIndex = j + i * sum->cols;
      arr[flatIndex] = matrix[flatIndex] + matrix2[flatIndex];
    }
  }
  return arr;
}

double *SumDoubleMatrices(double *matrix, double *matrix2, Sum *sum) {
  double(*arr) = malloc(sizeof *arr);
  for (int i = 0 ; i < sum->rows ; i++) {
    for (int j = 0 ; j < sum->cols ; j++) {
      int flatIndex = j + i * sum->cols;
      double doubleSum = matrix[flatIndex] + matrix2[flatIndex];
      arr[flatIndex] = doubleSum;
    }
  }
  return arr;
}