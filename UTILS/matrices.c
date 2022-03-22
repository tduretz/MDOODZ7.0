#include "stdlib.h"
#include "matrices.h"
#include "stdbool.h"

bool IsCutOff(int i, int j, CutCommand *cutCommand) {
  if (j + 1 <= cutCommand->xMinLimit) {
    return true;
  }
  if (j >= cutCommand->cols + cutCommand->xMaxLimit) {
    return true;
  }
  if (i + 1 <= cutCommand->yMinLimit) {
    return true;
  }
  if (i >= cutCommand->cols + cutCommand->yMaxLimit) {
    return true;
  }
  return false;
}

int *ReshapeIntArray(int *array, CutCommand *cutCommand) {
  int(*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0 ; i < cutCommand->rows ; i++) {
    for (int j = 0 ; j < cutCommand->cols ; j++) {
      int flatIndex = j + i * cutCommand->cols;
      if (IsCutOff(i, j, cutCommand)) {
        continue;
      }
      arr[updatedIndex] = array[flatIndex];
      updatedIndex++;
    }
  }
  return arr;
}

int *ReshapeFloatArray(float *array, CutCommand *cutCommand) {
  float(*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0 ; i < cutCommand->rows ; i++) {
    for (int j = 0 ; j < cutCommand->cols ; j++) {
      int flatIndex = j + i * cutCommand->cols;
      if (IsCutOff(i, j, cutCommand)) {
        continue;
      }
      arr[updatedIndex] = array[flatIndex];
      updatedIndex++;
    }
  }
  return arr;
}

int *ReshapeDoubleArray(double *array, CutCommand *cutCommand) {
  double (*arr) = malloc(sizeof *arr);
  int updatedIndex = 0;
  for (int i = 0 ; i < cutCommand->rows ; i++) {
    for (int j = 0 ; j < cutCommand->cols ; j++) {
      int flatIndex = j + i * cutCommand->cols;
      if (IsCutOff(i, j, cutCommand)) {
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