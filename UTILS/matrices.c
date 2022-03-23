#include "stdlib.h"
#include "matrices.h"
#include "stdbool.h"

static bool IsCutOff(int i, int j, Reshape *reshape) {
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

Reshape InitReshape(int rows, int cols) {
  Reshape const e = {.rows = rows,
                     .cols = cols,
                     .xMinLimit = 0,
                     .xMaxLimit = 0,
                     .yMinLimit = 0,
                     .yMaxLimit = 0};
  return e;
}

Sum InitSum(int rows, int cols) {
  Sum const e = {
          .rows = rows,
          .cols = cols,
          .multiplier = 1,
  };
  return e;
}

int *ReshapeIntArray(const int *array, Reshape *reshape) {
  int(*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0; i < reshape->rows; i++) {
    for (int j = 0; j < reshape->cols; j++) {
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

float *ReshapeFloatArray(const float *array, Reshape *reshape) {
  float(*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0; i < reshape->rows; i++) {
    for (int j = 0; j < reshape->cols; j++) {
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

double *ReshapeDoubleArray(const double *array, Reshape *reshape) {
  double *arr =
          (double *) malloc(sizeof(double) * reshape->cols * reshape->rows);
  int updatedIndex = 0;
  for (int i = 0; i < reshape->rows; i++) {
    for (int j = 0; j < reshape->cols; j++) {
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

int *SumIntMatrices(const int *matrix, const int *matrix2, Sum *sum) {
  int(*arr) = malloc(sizeof *arr);
  for (int i = 0; i < sum->rows; i++) {
    for (int j = 0; j < sum->cols; j++) {
      int flatIndex = j + i * sum->cols;
      arr[flatIndex] = matrix[flatIndex] + matrix2[flatIndex];
    }
  }
  return arr;
}

double *SumDoubleMatrices(const double *matrix, const double *matrix2,
                          Sum *sum) {
  double(*arr) = malloc(sizeof *arr);
  for (int i = 0; i < sum->rows; i++) {
    for (int j = 0; j < sum->cols; j++) {
      int flatIndex = j + i * sum->cols;
      double doubleSum =
              sum->multiplier * (matrix[flatIndex] + matrix2[flatIndex]);
      arr[flatIndex] = doubleSum;
    }
  }
  return arr;
}