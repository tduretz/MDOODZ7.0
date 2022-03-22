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

float *ReshapeFloatArray(float *array, Reshape *reshape);

double *ReshapeDoubleArray(double *array, Reshape *reshape);

typedef struct Sum {
  int rows;
  int cols;
} Sum;

int *SumIntMatrices(int *matrix, int *matrix2, Sum *sum);

double *SumDoubleMatrices(double *matrix, double *matrix2, Sum *sum);

#endif
