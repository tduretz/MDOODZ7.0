#ifndef MDOODZ_MATRICES_H
#define MDOODZ_MATRICES_H

typedef struct Reshape {
  int rows;
  int cols;
  int xMinLimit;
  int xMaxLimit;
  int yMinLimit;
  int yMaxLimit;
} Reshape;

int *ReshapeIntArray(int *array, Reshape *reshape);

int *ReshapeFloatArray(float *array, Reshape *reshape);

int *ReshapeDoubleArray(double *array, Reshape *reshape);

typedef struct Sum {
  int rows;
  int cols;
} Sum;

int *SumIntMatrices(int *matrix, int *matrix2, Sum *sum);

#endif
